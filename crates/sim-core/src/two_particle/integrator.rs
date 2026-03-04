use num_complex::Complex64;
use rustfft::FftPlanner;
use std::f64::consts::PI;
use std::sync::Arc;

use super::potential::Potential2D;
use super::wavefunction::Wavefunction2D;

/// 2D split-operator (split-step Fourier) time integrator.
///
/// Advances ψ(x₁, x₂) by dt using:
///   ψ(t+dt) ≈ e^{-iV dt/2} · IFFT₂D[ e^{-iT(k₁,k₂) dt} · FFT₂D[ e^{-iV dt/2} · ψ ] ]
///
/// The 2D FFT is performed as row-wise 1D FFTs followed by column-wise 1D FFTs.
pub struct SplitOperator2D {
    /// e^{-iV(x₁,x₂) dt/2} for each grid point (n² entries, row-major).
    potential_half: Vec<Complex64>,
    /// e^{-iT(k₁,k₂) dt} for each frequency pair (n² entries, row-major).
    kinetic_full: Vec<Complex64>,
    /// Forward FFT plan (length n).
    fft: Arc<dyn rustfft::Fft<f64>>,
    /// Inverse FFT plan (length n).
    ifft: Arc<dyn rustfft::Fft<f64>>,
    /// Scratch buffer for FFT.
    scratch: Vec<Complex64>,
    /// Contiguous buffer for column FFTs.
    col_buf: Vec<Complex64>,
    /// Grid points per axis.
    n: usize,
}

impl SplitOperator2D {
    pub fn new(wf: &Wavefunction2D, potential: &Potential2D, dt: f64) -> Self {
        let n = wf.n;
        let length = wf.x_max - wf.x_min;
        let xs = wf.xs();

        // Precompute potential half-step: exp(-i V(x₁, x₂) dt/2)
        let potential_half: Vec<Complex64> = {
            let mut v = Vec::with_capacity(n * n);
            for i1 in 0..n {
                for i2 in 0..n {
                    let val = potential.value_at(xs[i1], xs[i2]);
                    v.push(Complex64::new(0.0, -val * dt / 2.0).exp());
                }
            }
            v
        };

        // Precompute kinetic full step: exp(-i (k₁² + k₂²)/2 dt)
        // Frequency indices in FFT standard order
        let freq = |i: usize| -> f64 {
            let j = if i <= n / 2 { i as f64 } else { i as f64 - n as f64 };
            2.0 * PI * j / length
        };

        let kinetic_full: Vec<Complex64> = {
            let mut k = Vec::with_capacity(n * n);
            for i1 in 0..n {
                let k1 = freq(i1);
                for i2 in 0..n {
                    let k2 = freq(i2);
                    let t_k = (k1 * k1 + k2 * k2) / 2.0;
                    k.push(Complex64::new(0.0, -t_k * dt).exp());
                }
            }
            k
        };

        let mut planner = FftPlanner::new();
        let fft = planner.plan_fft_forward(n);
        let ifft = planner.plan_fft_inverse(n);
        let scratch_len = fft.get_inplace_scratch_len().max(ifft.get_inplace_scratch_len());

        Self {
            potential_half,
            kinetic_full,
            fft,
            ifft,
            scratch: vec![Complex64::new(0.0, 0.0); scratch_len],
            col_buf: vec![Complex64::new(0.0, 0.0); n],
            n,
        }
    }

    /// Advance the wavefunction by one time step.
    pub fn step(&mut self, wf: &mut Wavefunction2D) {
        let n = self.n;
        let psi = &mut wf.psi;

        // 1. Potential half-step in position space
        for i in 0..n * n {
            psi[i] *= self.potential_half[i];
        }

        // 2. Forward 2D FFT (row-wise then column-wise)
        self.fft_2d(psi, true);

        // 3. Kinetic full step in momentum space
        for i in 0..n * n {
            psi[i] *= self.kinetic_full[i];
        }

        // 4. Inverse 2D FFT
        self.fft_2d(psi, false);

        // 5. Normalize (1/n per axis = 1/n² total, from rustfft convention)
        let inv_n2 = 1.0 / (n * n) as f64;
        for c in psi.iter_mut() {
            *c *= inv_n2;
        }

        // 6. Potential half-step again
        for i in 0..n * n {
            psi[i] *= self.potential_half[i];
        }
    }

    /// Perform a 2D FFT (forward or inverse) via row-wise + column-wise 1D FFTs.
    fn fft_2d(&mut self, psi: &mut [Complex64], forward: bool) {
        let n = self.n;
        let plan = if forward { &self.fft } else { &self.ifft };

        // Row-wise FFTs (each row is contiguous in memory)
        for i1 in 0..n {
            let row = &mut psi[i1 * n..(i1 + 1) * n];
            plan.process_with_scratch(row, &mut self.scratch);
        }

        // Column-wise FFTs (extract to contiguous buffer, FFT, copy back)
        for i2 in 0..n {
            // Extract column
            for i1 in 0..n {
                self.col_buf[i1] = psi[i1 * n + i2];
            }
            // FFT
            plan.process_with_scratch(&mut self.col_buf, &mut self.scratch);
            // Copy back
            for i1 in 0..n {
                psi[i1 * n + i2] = self.col_buf[i1];
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::two_particle::potential::{Interaction, Potential2D};
    use crate::two_particle::TwoParticleSimulation;
    use crate::potential::Potential;

    #[test]
    fn free_particles_preserve_norm() {
        let pot = Potential2D::new(Potential::free(), Potential::free(), Interaction::None);
        let mut sim = TwoParticleSimulation::new(128, -15.0, 15.0, 0.005, pot);
        sim.wf.set_product_gaussian(-3.0, 1.5, 2.0, 3.0, 1.5, -2.0);

        let n0 = sim.norm();
        sim.step_n(500);
        let n1 = sim.norm();

        assert!(
            (n1 - n0).abs() < 1e-10,
            "norm drift: {n0} → {n1}"
        );
    }

    #[test]
    fn interacting_particles_preserve_norm() {
        let pot = Potential2D::new(
            Potential::free(),
            Potential::free(),
            Interaction::SoftCoulomb { g: 1.0, epsilon: 0.5 },
        );
        let mut sim = TwoParticleSimulation::new(128, -15.0, 15.0, 0.005, pot);
        sim.wf.set_product_gaussian(-5.0, 1.5, 2.0, 5.0, 1.5, -2.0);

        let n0 = sim.norm();
        sim.step_n(500);
        let n1 = sim.norm();

        assert!(
            (n1 - n0).abs() < 1e-10,
            "norm drift with interaction: {n0} → {n1}"
        );
    }

    #[test]
    fn non_interacting_marginals_match_1d() {
        // Two non-interacting particles: marginals should match independent 1D evolution
        use crate::Simulation;

        let x0_1 = -3.0;
        let sigma_1 = 1.5;
        let k0_1 = 2.0;
        let x0_2 = 3.0;
        let sigma_2 = 1.0;
        let k0_2 = -1.0;

        let n = 128;
        let x_min = -15.0;
        let x_max = 15.0;
        let dt = 0.005;
        let steps = 200;

        // 2D simulation
        let pot2d = Potential2D::new(Potential::free(), Potential::free(), Interaction::None);
        let mut sim2d = TwoParticleSimulation::new(n, x_min, x_max, dt, pot2d);
        sim2d.wf.set_product_gaussian(x0_1, sigma_1, k0_1, x0_2, sigma_2, k0_2);

        // Independent 1D simulations
        let mut sim1 = Simulation::new(n, x_min, x_max, dt, Potential::free());
        sim1.wf.set_gaussian(x0_1, sigma_1, k0_1);
        let mut sim2 = Simulation::new(n, x_min, x_max, dt, Potential::free());
        sim2.wf.set_gaussian(x0_2, sigma_2, k0_2);

        sim2d.step_n(steps);
        sim1.step_n(steps);
        sim2.step_n(steps);

        // Compare marginals to 1D probability densities
        let m1 = sim2d.wf.marginal_1();
        let d1 = sim1.wf.probability_density();
        let max_diff_1: f64 = m1.iter().zip(d1.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);

        let m2 = sim2d.wf.marginal_2();
        let d2 = sim2.wf.probability_density();
        let max_diff_2: f64 = m2.iter().zip(d2.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);

        assert!(
            max_diff_1 < 1e-6,
            "marginal_1 differs from 1D: max |Δ| = {max_diff_1}"
        );
        assert!(
            max_diff_2 < 1e-6,
            "marginal_2 differs from 1D: max |Δ| = {max_diff_2}"
        );
    }

    #[test]
    fn energy_conservation_with_interaction() {
        let pot = Potential2D::new(
            Potential::harmonic(0.0, 1.0),
            Potential::harmonic(0.0, 1.0),
            Interaction::SoftCoulomb { g: 0.5, epsilon: 0.5 },
        );
        let mut sim = TwoParticleSimulation::new(128, -15.0, 15.0, 0.005, pot);
        sim.wf.set_product_gaussian(-2.0, 1.0, 0.0, 2.0, 1.0, 0.0);

        let e0 = sim.expected_energy();
        sim.step_n(500);
        let e1 = sim.expected_energy();

        assert!(
            (e1 - e0).abs() < 0.05,
            "energy drift: {e0} → {e1}"
        );
    }
}
