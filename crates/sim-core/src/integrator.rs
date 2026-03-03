use num_complex::Complex64;
use rustfft::FftPlanner;
use std::f64::consts::PI;
use std::sync::Arc;

use crate::potential::Potential;
use crate::wavefunction::Wavefunction;

/// Split-operator (split-step Fourier) time integrator.
///
/// Advances the wavefunction by dt using:
///   ψ(t+dt) ≈ e^{-iV dt/2} · IFFT[ e^{-iT(k) dt} · FFT[ e^{-iV dt/2} · ψ(t) ] ]
///
/// This is second-order accurate in dt and exactly preserves unitarity
/// (total probability is conserved).
pub struct SplitOperator {
    /// e^{-iV(x) dt/2} for each grid point (potential half-step in position space).
    potential_half: Vec<Complex64>,
    /// e^{-iT(k) dt} for each frequency (kinetic full step in momentum space).
    kinetic_full: Vec<Complex64>,
    /// Forward FFT plan.
    fft: Arc<dyn rustfft::Fft<f64>>,
    /// Inverse FFT plan.
    ifft: Arc<dyn rustfft::Fft<f64>>,
    /// Scratch buffer for FFT.
    scratch: Vec<Complex64>,
    /// Grid size (for FFT normalization).
    n: usize,
}

impl SplitOperator {
    /// Create a new integrator for the given wavefunction grid, potential, and time step.
    pub fn new(wf: &Wavefunction, potential: &Potential, dt: f64) -> Self {
        let n = wf.n;
        let length = wf.x_max - wf.x_min;

        // Pre-compute potential half-step: exp(-i V(x_j) dt/2)
        let potential_half: Vec<Complex64> = (0..n)
            .map(|i| {
                let v = potential.value_at(wf.x(i));
                Complex64::new(0.0, -v * dt / 2.0).exp()
            })
            .collect();

        // Pre-compute kinetic full step: exp(-i k²/2 dt)
        // Frequencies in FFT standard order: 0, 1, ..., N/2-1, -N/2, -N/2+1, ..., -1
        // k_j = 2π j / L
        let kinetic_full: Vec<Complex64> = (0..n)
            .map(|i| {
                let j = if i <= n / 2 {
                    i as f64
                } else {
                    i as f64 - n as f64
                };
                let k = 2.0 * PI * j / length;
                let t_k = k * k / 2.0; // kinetic energy = k²/2 (atomic units, m=1)
                Complex64::new(0.0, -t_k * dt).exp()
            })
            .collect();

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
            n,
        }
    }

    /// Advance the wavefunction by one time step.
    pub fn step(&mut self, wf: &mut Wavefunction) {
        let psi = &mut wf.psi;

        // 1. Apply potential half-step in position space
        for i in 0..self.n {
            psi[i] *= self.potential_half[i];
        }

        // 2. FFT to momentum space
        self.fft.process_with_scratch(psi, &mut self.scratch);

        // 3. Apply kinetic full step in momentum space
        for i in 0..self.n {
            psi[i] *= self.kinetic_full[i];
        }

        // 4. Inverse FFT back to position space
        self.ifft.process_with_scratch(psi, &mut self.scratch);

        // 5. Normalize (rustfft doesn't normalize the inverse FFT)
        let inv_n = 1.0 / self.n as f64;
        for c in psi.iter_mut() {
            *c *= inv_n;
        }

        // 6. Apply potential half-step again
        for i in 0..self.n {
            psi[i] *= self.potential_half[i];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Simulation;

    #[test]
    fn free_particle_preserves_norm() {
        let mut sim = Simulation::new(512, -20.0, 20.0, 0.01, Potential::free());
        sim.wf.set_gaussian(0.0, 2.0, 3.0);

        let initial_norm = sim.norm();
        sim.step_n(1000);
        let final_norm = sim.norm();

        assert!(
            (final_norm - initial_norm).abs() < 1e-10,
            "norm drift: {initial_norm} → {final_norm}"
        );
    }

    #[test]
    fn free_particle_moves_with_group_velocity() {
        // A wave packet with momentum k0 should move at velocity v = k0/m = k0 (atomic units).
        let k0 = 3.0;
        let dt = 0.01;
        let steps = 500;

        let mut sim = Simulation::new(2048, -40.0, 40.0, dt, Potential::free());
        sim.wf.set_gaussian(0.0, 2.0, k0);

        let x0 = sim.expected_x();
        sim.step_n(steps);
        let x1 = sim.expected_x();

        let expected_displacement = k0 * dt * steps as f64;
        let actual_displacement = x1 - x0;

        assert!(
            (actual_displacement - expected_displacement).abs() < 0.1,
            "displacement: expected {expected_displacement}, got {actual_displacement}"
        );
    }

    #[test]
    fn harmonic_oscillator_ground_state_is_stationary() {
        // The ground state of a harmonic oscillator (ω=1) is a Gaussian with σ = 1/√(2ω).
        // ψ₀(x) = (ω/π)^{1/4} exp(-ωx²/2), which matches our Gaussian with σ² = 1/(2ω).
        // It should not change shape over time (it's an eigenstate).
        let omega: f64 = 1.0;
        let sigma = 1.0 / (2.0 * omega).sqrt();

        let mut sim = Simulation::new(1024, -15.0, 15.0, 0.005, Potential::harmonic(0.0, omega));
        sim.wf.set_gaussian(0.0, sigma, 0.0);

        let initial_density: Vec<f64> = sim.wf.probability_density();
        sim.step_n(2000); // 10 time units
        let final_density: Vec<f64> = sim.wf.probability_density();

        // The shape should be essentially unchanged.
        let max_diff: f64 = initial_density
            .iter()
            .zip(final_density.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);

        assert!(
            max_diff < 1e-4,
            "ground state changed: max |Δρ| = {max_diff}"
        );
    }

    #[test]
    fn energy_conservation() {
        let mut sim = Simulation::new(1024, -20.0, 20.0, 0.01, Potential::harmonic(0.0, 1.0));
        sim.wf.set_gaussian(-3.0, 1.0, 0.0); // displaced from equilibrium

        let e0 = sim.expected_energy();
        sim.step_n(1000);
        let e1 = sim.expected_energy();

        assert!(
            (e1 - e0).abs() < 0.01,
            "energy drift: {e0} → {e1}"
        );
    }

    #[test]
    fn tunneling_through_barrier() {
        // A particle with momentum hitting a finite barrier should partially tunnel.
        let k0 = 4.0;
        let barrier = Potential::barrier(-0.5, 0.5, 10.0); // barrier height > kinetic energy k²/2 = 8

        let mut sim = Simulation::new(2048, -30.0, 30.0, 0.005, barrier);
        sim.wf.set_gaussian(-8.0, 1.5, k0);

        sim.step_n(3000);

        // Some probability should have tunneled to x > 0.5
        let prob_right: f64 = sim
            .wf
            .probability_density()
            .iter()
            .enumerate()
            .filter(|&(i, _)| sim.wf.x(i) > 0.5)
            .map(|(_, p)| p * sim.wf.dx)
            .sum();

        assert!(
            prob_right > 0.001,
            "expected some tunneling, got P(x>0.5) = {prob_right}"
        );
        assert!(
            prob_right < 0.99,
            "expected partial reflection, got P(x>0.5) = {prob_right}"
        );
    }
}
