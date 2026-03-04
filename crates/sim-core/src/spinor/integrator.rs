use num_complex::Complex64;
use rustfft::FftPlanner;
use std::f64::consts::PI;
use std::sync::Arc;

use crate::potential::Potential;
use super::wavefunction::SpinorWavefunction;

/// Magnetic field configuration.
#[derive(Clone)]
pub struct MagneticField {
    /// Uniform longitudinal field B_z (Zeeman splitting, no spin mixing).
    pub bz: f64,
    /// Transverse field B_x (spin precession, mixes ↑ and ↓).
    pub bx: f64,
    /// Longitudinal field gradient dB_z/dx (for Stern-Gerlach effect).
    /// Total B_z at position x = bz + gradient_bz * x.
    pub gradient_bz: f64,
}

impl MagneticField {
    pub fn zero() -> Self {
        Self { bz: 0.0, bx: 0.0, gradient_bz: 0.0 }
    }

    pub fn longitudinal(bz: f64) -> Self {
        Self { bz, bx: 0.0, gradient_bz: 0.0 }
    }

    pub fn transverse(bx: f64) -> Self {
        Self { bz: 0.0, bx, gradient_bz: 0.0 }
    }

    /// B_z(x) = bz + gradient_bz * x
    pub fn bz_at(&self, x: f64) -> f64 {
        self.bz + self.gradient_bz * x
    }
}

/// Split-operator integrator for a spin-1/2 particle.
///
/// Handles the spatial Hamiltonian (V + T via FFT) and magnetic coupling
/// (Zeeman splitting via spin-dependent potential, Larmor precession via
/// spin rotation).
pub struct SpinorIntegrator {
    /// e^{-iV_↑(x) dt/2} for spin-up: V_↑ = V(x) + B_z/2
    potential_half_up: Vec<Complex64>,
    /// e^{-iV_↓(x) dt/2} for spin-down: V_↓ = V(x) - B_z/2
    potential_half_down: Vec<Complex64>,
    /// e^{-iT(k) dt} for each frequency (same for both spins).
    kinetic_full: Vec<Complex64>,
    /// Spin rotation half-step matrix element: cos(B_x dt/4)
    spin_cos: f64,
    /// Spin rotation half-step matrix element: sin(B_x dt/4)
    spin_sin: f64,
    /// Forward FFT plan.
    fft: Arc<dyn rustfft::Fft<f64>>,
    /// Inverse FFT plan.
    ifft: Arc<dyn rustfft::Fft<f64>>,
    /// Scratch buffer.
    scratch: Vec<Complex64>,
    n: usize,
}

impl SpinorIntegrator {
    pub fn new(
        wf: &SpinorWavefunction,
        potential: &Potential,
        field: &MagneticField,
        dt: f64,
    ) -> Self {
        let n = wf.n;
        let length = wf.x_max - wf.x_min;

        // Potential half-step for each spin component
        // V_↑(x) = V(x) + B_z(x)/2,  V_↓(x) = V(x) - B_z(x)/2
        let potential_half_up: Vec<Complex64> = (0..n)
            .map(|i| {
                let x = wf.x(i);
                let v = potential.value_at(x) + field.bz_at(x) / 2.0;
                Complex64::new(0.0, -v * dt / 2.0).exp()
            })
            .collect();

        let potential_half_down: Vec<Complex64> = (0..n)
            .map(|i| {
                let x = wf.x(i);
                let v = potential.value_at(x) - field.bz_at(x) / 2.0;
                Complex64::new(0.0, -v * dt / 2.0).exp()
            })
            .collect();

        // Kinetic full step (spin-independent)
        let kinetic_full: Vec<Complex64> = (0..n)
            .map(|i| {
                let j = if i <= n / 2 {
                    i as f64
                } else {
                    i as f64 - n as f64
                };
                let k = 2.0 * PI * j / length;
                let t_k = k * k / 2.0;
                Complex64::new(0.0, -t_k * dt).exp()
            })
            .collect();

        // Spin rotation half-step: exp(-i B_x σ_x dt/4)
        // = cos(B_x dt/4) I - i sin(B_x dt/4) σ_x
        // Matrix: (cos   -i sin)
        //         (-i sin  cos )
        let angle = field.bx * dt / 4.0;
        let spin_cos = angle.cos();
        let spin_sin = angle.sin();

        let mut planner = FftPlanner::new();
        let fft = planner.plan_fft_forward(n);
        let ifft = planner.plan_fft_inverse(n);
        let scratch_len = fft.get_inplace_scratch_len().max(ifft.get_inplace_scratch_len());

        Self {
            potential_half_up,
            potential_half_down,
            kinetic_full,
            spin_cos,
            spin_sin,
            fft,
            ifft,
            scratch: vec![Complex64::new(0.0, 0.0); scratch_len],
            n,
        }
    }

    /// Advance the spinor wavefunction by one time step.
    pub fn step(&mut self, wf: &mut SpinorWavefunction) {
        let n = self.n;

        // 1. Potential half-step (spin-dependent)
        for i in 0..n {
            wf.up[i] *= self.potential_half_up[i];
            wf.down[i] *= self.potential_half_down[i];
        }

        // 2. Spin rotation half-step (if B_x ≠ 0)
        if self.spin_sin.abs() > 0.0 {
            self.apply_spin_rotation(wf);
        }

        // 3. FFT both components to momentum space
        self.fft
            .process_with_scratch(&mut wf.up, &mut self.scratch);
        self.fft
            .process_with_scratch(&mut wf.down, &mut self.scratch);

        // 4. Kinetic full step (same for both)
        for i in 0..n {
            wf.up[i] *= self.kinetic_full[i];
            wf.down[i] *= self.kinetic_full[i];
        }

        // 5. IFFT back to position space
        self.ifft
            .process_with_scratch(&mut wf.up, &mut self.scratch);
        self.ifft
            .process_with_scratch(&mut wf.down, &mut self.scratch);

        // 6. Normalize (rustfft convention)
        let inv_n = 1.0 / n as f64;
        for c in wf.up.iter_mut() {
            *c *= inv_n;
        }
        for c in wf.down.iter_mut() {
            *c *= inv_n;
        }

        // 7. Spin rotation half-step
        if self.spin_sin.abs() > 0.0 {
            self.apply_spin_rotation(wf);
        }

        // 8. Potential half-step
        for i in 0..n {
            wf.up[i] *= self.potential_half_up[i];
            wf.down[i] *= self.potential_half_down[i];
        }
    }

    /// Apply the spin rotation exp(-i B_x σ_x dt/4) at each grid point.
    ///
    /// (ψ_↑')   (cos     -i sin) (ψ_↑)
    /// (ψ_↓') = (-i sin   cos  ) (ψ_↓)
    fn apply_spin_rotation(&self, wf: &mut SpinorWavefunction) {
        let cos = Complex64::new(self.spin_cos, 0.0);
        let neg_i_sin = Complex64::new(0.0, -self.spin_sin);

        for i in 0..self.n {
            let u = wf.up[i];
            let d = wf.down[i];
            wf.up[i] = cos * u + neg_i_sin * d;
            wf.down[i] = neg_i_sin * u + cos * d;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::SpinorSimulation;

    #[test]
    fn free_spinor_preserves_norm() {
        let mut sim = SpinorSimulation::new(
            512, -20.0, 20.0, 0.01,
            Potential::free(), MagneticField::zero(),
        );
        sim.wf.set_gaussian(
            0.0, 2.0, 3.0,
            Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0),
        );

        let n0 = sim.wf.norm();
        sim.step_n(500);
        let n1 = sim.wf.norm();

        assert!(
            (n1 - n0).abs() < 1e-10,
            "norm drift: {n0} → {n1}"
        );
    }

    #[test]
    fn zeeman_preserves_norm() {
        let mut sim = SpinorSimulation::new(
            512, -20.0, 20.0, 0.01,
            Potential::harmonic(0.0, 1.0),
            MagneticField::longitudinal(2.0),
        );
        // Start in superposition of spin states
        sim.wf.set_gaussian(
            0.0, 1.0, 0.0,
            Complex64::new(1.0, 0.0), Complex64::new(1.0, 0.0),
        );

        let n0 = sim.wf.norm();
        sim.step_n(500);
        let n1 = sim.wf.norm();

        assert!(
            (n1 - n0).abs() < 1e-10,
            "norm drift with B_z: {n0} → {n1}"
        );
    }

    #[test]
    fn larmor_precession() {
        // Spin-up in B_x field should precess: ⟨σ_z⟩ oscillates
        let bx = 2.0;
        let dt = 0.005;
        let mut sim = SpinorSimulation::new(
            512, -20.0, 20.0, dt,
            Potential::harmonic(0.0, 1.0), // trap to prevent spreading
            MagneticField::transverse(bx),
        );
        sim.wf.set_gaussian(
            0.0, 1.0, 0.0,
            Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0),
        );

        // Initially ⟨σ_z⟩ = 1
        assert!((sim.wf.expected_sz() - 1.0).abs() < 1e-6);

        // Evolve to quarter period: t = π/(2 B_x)
        let quarter_period = PI / (2.0 * bx);
        let steps = (quarter_period / dt).round() as usize;
        sim.step_n(steps);

        // ⟨σ_z⟩ should be ≈ 0
        let sz = sim.wf.expected_sz();
        assert!(
            sz.abs() < 0.1,
            "at quarter period ⟨σ_z⟩ = {sz}, expected ≈ 0"
        );

        // Evolve to half period: ⟨σ_z⟩ ≈ -1
        sim.step_n(steps);
        let sz = sim.wf.expected_sz();
        assert!(
            (sz + 1.0).abs() < 0.1,
            "at half period ⟨σ_z⟩ = {sz}, expected ≈ -1"
        );
    }

    #[test]
    fn pure_spin_up_in_bz_stays_up() {
        // Spin-up is an eigenstate of σ_z, so B_z shouldn't change ⟨σ_z⟩
        let mut sim = SpinorSimulation::new(
            512, -20.0, 20.0, 0.01,
            Potential::harmonic(0.0, 1.0),
            MagneticField::longitudinal(3.0),
        );
        sim.wf.set_gaussian(
            0.0, 1.0, 0.0,
            Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0),
        );

        sim.step_n(500);
        let sz = sim.wf.expected_sz();
        assert!(
            (sz - 1.0).abs() < 1e-6,
            "spin-up in B_z: ⟨σ_z⟩ = {sz}, expected 1.0"
        );
    }
}
