use num_complex::Complex64;
use rustfft::FftPlanner;
use std::f64::consts::PI;
use std::sync::Arc;

use crate::potential::Potential;
use crate::spinor::wavefunction::SpinorWavefunction;

/// Split-operator integrator for the 1D Dirac equation.
///
/// H = -ic σ_x ∂/∂x + mc² σ_z + V(x)
///
/// The kinetic + mass step is a 2×2 unitary rotation at each momentum mode:
///   exp(-i H_k dt) = cos(E_k dt) I - i sin(E_k dt) (H_k / E_k)
/// where E_k = √(c²k² + m²c⁴).
pub struct DiracIntegrator {
    /// exp(-i V(x) dt/2) for each grid point (same for both components).
    potential_half: Vec<Complex64>,
    /// Precomputed kinetic rotation parameters at each frequency.
    /// Each entry: (cos(E_k dt), sin(E_k dt), ck/E_k, mc²/E_k)
    kinetic_params: Vec<KineticRotation>,
    /// Forward FFT plan.
    fft: Arc<dyn rustfft::Fft<f64>>,
    /// Inverse FFT plan.
    ifft: Arc<dyn rustfft::Fft<f64>>,
    /// Scratch buffer.
    scratch: Vec<Complex64>,
    n: usize,
}

/// Precomputed rotation parameters for one momentum mode.
#[derive(Clone, Copy)]
struct KineticRotation {
    /// cos(E_k dt)
    cos_e: f64,
    /// sin(E_k dt)
    sin_e: f64,
    /// ck / E_k  (off-diagonal mixing strength)
    ck_over_e: f64,
    /// mc² / E_k  (diagonal splitting)
    mc2_over_e: f64,
}

impl DiracIntegrator {
    /// Create a new Dirac integrator.
    ///
    /// - `c`: speed of light (in atomic units, physical value ≈ 137)
    /// - `mass`: particle mass (1.0 for electron in atomic units)
    pub fn new(
        wf: &SpinorWavefunction,
        potential: &Potential,
        c: f64,
        mass: f64,
        dt: f64,
    ) -> Self {
        let n = wf.n;
        let length = wf.x_max - wf.x_min;
        let mc2 = mass * c * c;

        // Potential half-step (electrostatic, same for both components)
        let potential_half: Vec<Complex64> = (0..n)
            .map(|i| {
                let v = potential.value_at(wf.x(i));
                Complex64::new(0.0, -v * dt / 2.0).exp()
            })
            .collect();

        // Kinetic rotation at each momentum mode
        let kinetic_params: Vec<KineticRotation> = (0..n)
            .map(|i| {
                let j = if i <= n / 2 {
                    i as f64
                } else {
                    i as f64 - n as f64
                };
                let k = 2.0 * PI * j / length;
                let ck = c * k;
                let e_k = (ck * ck + mc2 * mc2).sqrt();

                if e_k < 1e-30 {
                    // k = 0 and m = 0: identity
                    KineticRotation {
                        cos_e: 1.0,
                        sin_e: 0.0,
                        ck_over_e: 0.0,
                        mc2_over_e: 0.0,
                    }
                } else {
                    KineticRotation {
                        cos_e: (e_k * dt).cos(),
                        sin_e: (e_k * dt).sin(),
                        ck_over_e: ck / e_k,
                        mc2_over_e: mc2 / e_k,
                    }
                }
            })
            .collect();

        let mut planner = FftPlanner::new();
        let fft = planner.plan_fft_forward(n);
        let ifft = planner.plan_fft_inverse(n);
        let scratch_len = fft.get_inplace_scratch_len().max(ifft.get_inplace_scratch_len());

        Self {
            potential_half,
            kinetic_params,
            fft,
            ifft,
            scratch: vec![Complex64::new(0.0, 0.0); scratch_len],
            n,
        }
    }

    /// Advance the wavefunction by one time step.
    pub fn step(&mut self, wf: &mut SpinorWavefunction) {
        let n = self.n;

        // 1. Potential half-step (same for both components)
        for i in 0..n {
            wf.up[i] *= self.potential_half[i];
            wf.down[i] *= self.potential_half[i];
        }

        // 2. FFT both components
        self.fft
            .process_with_scratch(&mut wf.up, &mut self.scratch);
        self.fft
            .process_with_scratch(&mut wf.down, &mut self.scratch);

        // 3. Dirac kinetic + mass step: 2×2 rotation at each k
        //
        // exp(-i H_k dt) = cos(E dt) I - i sin(E dt) H_k/E
        //
        // H_k/E = ( mc²/E    ck/E )
        //         ( ck/E   -mc²/E )
        //
        // So: ψ_up'  = [cos(E dt) - i sin(E dt) mc²/E] ψ_up  - i sin(E dt) ck/E  ψ_down
        //     ψ_down' = -i sin(E dt) ck/E  ψ_up + [cos(E dt) + i sin(E dt) mc²/E] ψ_down
        for i in 0..n {
            let r = &self.kinetic_params[i];
            let a = Complex64::new(r.cos_e, -r.sin_e * r.mc2_over_e); // diagonal: cos - i sin mc²/E
            let b = Complex64::new(0.0, -r.sin_e * r.ck_over_e); // off-diag: -i sin ck/E

            let u = wf.up[i];
            let d = wf.down[i];
            wf.up[i] = a * u + b * d;
            wf.down[i] = b * u + a.conj() * d; // lower diagonal: cos + i sin mc²/E = conj(a)
        }

        // 4. IFFT both components
        self.ifft
            .process_with_scratch(&mut wf.up, &mut self.scratch);
        self.ifft
            .process_with_scratch(&mut wf.down, &mut self.scratch);

        // 5. Normalize (rustfft convention)
        let inv_n = 1.0 / n as f64;
        for c in wf.up.iter_mut() {
            *c *= inv_n;
        }
        for c in wf.down.iter_mut() {
            *c *= inv_n;
        }

        // 6. Potential half-step
        for i in 0..n {
            wf.up[i] *= self.potential_half[i];
            wf.down[i] *= self.potential_half[i];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::DiracSimulation;
    use crate::potential::Potential;
    use num_complex::Complex64;

    #[test]
    fn free_dirac_preserves_norm() {
        let mut sim = DiracSimulation::new(512, -30.0, 30.0, 0.005, Potential::free(), 10.0, 1.0);
        sim.wf.set_gaussian(
            0.0, 2.0, 3.0,
            Complex64::new(1.0, 0.0),
            Complex64::new(0.0, 0.0),
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
    fn positive_energy_packet_moves_at_group_velocity() {
        // For a relativistic particle, group velocity v_g = dE/dk = c²k/E_k
        let c: f64 = 10.0;
        let mass = 1.0;
        let k0 = 2.0;
        let dt = 0.005;
        let steps = 400;

        let mc2 = mass * c * c;
        let e_k = (c * c * k0 * k0 + mc2 * mc2).sqrt();
        let v_group = c * c * k0 / e_k;

        let mut sim =
            DiracSimulation::new(2048, -50.0, 50.0, dt, Potential::free(), c, mass);
        sim.init_positive_energy_packet(0.0, 2.0, k0);

        let x0 = sim.wf.expected_x();
        sim.step_n(steps);
        let x1 = sim.wf.expected_x();

        let expected_disp = v_group * dt * steps as f64;
        let actual_disp = x1 - x0;

        assert!(
            (actual_disp - expected_disp).abs() < 0.5,
            "displacement: expected {expected_disp:.3}, got {actual_disp:.3}"
        );
    }

    #[test]
    fn zitterbewegung_with_naive_init() {
        // Upper-component-only init with nonzero momentum creates a superposition
        // of positive and negative energy states, producing Zitterbewegung.
        let c = 5.0;
        let mut sim = DiracSimulation::new(1024, -30.0, 30.0, 0.002, Potential::free(), c, 1.0);
        sim.wf.set_gaussian(
            0.0, 2.0, 3.0, // nonzero momentum for stronger mixing
            Complex64::new(1.0, 0.0),
            Complex64::new(0.0, 0.0),
        );

        // The lower component should acquire amplitude from Zitterbewegung
        sim.step_n(200);
        let prob_down = sim.wf.prob_down();
        assert!(
            prob_down > 0.01,
            "expected lower component excitation, got P(↓) = {prob_down}"
        );
    }

    #[test]
    fn positive_energy_init_suppresses_zitterbewegung() {
        let c = 10.0;
        let mut sim = DiracSimulation::new(1024, -30.0, 30.0, 0.005, Potential::free(), c, 1.0);
        sim.init_positive_energy_packet(0.0, 2.0, 3.0);

        // Lower component should be small but present (it's the small component)
        let p_down_0 = sim.wf.prob_down();

        sim.step_n(200);
        let p_down_1 = sim.wf.prob_down();

        // The lower component probability should stay roughly constant
        // (no Zitterbewegung oscillation)
        assert!(
            (p_down_1 - p_down_0).abs() < 0.02,
            "lower component changed: {p_down_0:.4} → {p_down_1:.4}"
        );
    }
}
