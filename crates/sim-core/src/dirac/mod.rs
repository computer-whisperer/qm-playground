pub mod integrator;

use integrator::DiracIntegrator;
use num_complex::Complex64;
use std::f64::consts::PI;

use crate::potential::Potential;
use crate::spinor::wavefunction::SpinorWavefunction;

/// Complete 1D Dirac equation simulation state.
///
/// Uses the same `SpinorWavefunction` as the spin module — the Dirac spinor
/// has two components (upper/lower) on the same spatial grid.
pub struct DiracSimulation {
    pub wf: SpinorWavefunction,
    pub potential: Potential,
    pub c: f64,
    pub mass: f64,
    pub integrator: DiracIntegrator,
    pub time: f64,
    pub dt: f64,
}

impl DiracSimulation {
    pub fn new(
        n: usize,
        x_min: f64,
        x_max: f64,
        dt: f64,
        potential: Potential,
        c: f64,
        mass: f64,
    ) -> Self {
        let wf = SpinorWavefunction::new(n, x_min, x_max);
        let integrator = DiracIntegrator::new(&wf, &potential, c, mass, dt);
        Self {
            wf,
            potential,
            c,
            mass,
            integrator,
            time: 0.0,
            dt,
        }
    }

    /// Advance by one time step.
    pub fn step(&mut self) {
        self.integrator.step(&mut self.wf);
        self.time += self.dt;
    }

    /// Advance by n time steps.
    pub fn step_n(&mut self, n: usize) {
        for _ in 0..n {
            self.step();
        }
    }

    /// Recompute integrator (call after changing potential, c, mass, or dt).
    pub fn rebuild_integrator(&mut self) {
        self.integrator =
            DiracIntegrator::new(&self.wf, &self.potential, self.c, self.mass, self.dt);
    }

    pub fn norm(&self) -> f64 {
        self.wf.norm()
    }

    pub fn expected_x(&self) -> f64 {
        self.wf.expected_x()
    }

    /// Expected momentum via finite differences (both components).
    pub fn expected_p(&self) -> f64 {
        let n = self.wf.n;
        let dx = self.wf.dx;
        let mut sum = Complex64::new(0.0, 0.0);
        for i in 0..n {
            let ip = (i + 1) % n;
            let im = (i + n - 1) % n;
            let du = (self.wf.up[ip] - self.wf.up[im]) / (2.0 * dx);
            sum += self.wf.up[i].conj() * du;
            let dd = (self.wf.down[ip] - self.wf.down[im]) / (2.0 * dx);
            sum += self.wf.down[i].conj() * dd;
        }
        (Complex64::new(0.0, -1.0) * sum * dx).re
    }

    /// Initialize a positive-energy Dirac wave packet.
    ///
    /// For each momentum mode k, the positive-energy Dirac spinor is:
    ///   u(k) = N ( 1, ck/(E_k + mc²) )
    ///
    /// This method creates a Gaussian in the upper component, FFTs it, projects
    /// each mode onto the positive-energy spinor, then IFFTs both components back.
    /// The result is a wave packet with no negative-energy content, so it
    /// propagates cleanly without Zitterbewegung.
    pub fn init_positive_energy_packet(&mut self, x0: f64, sigma: f64, k0: f64) {
        let n = self.wf.n;
        let length = self.wf.x_max - self.wf.x_min;
        let mc2 = self.mass * self.c * self.c;

        // Start with a spatial Gaussian in the upper component
        let spatial_norm = (2.0 * PI * sigma * sigma).powf(-0.25);
        for i in 0..n {
            let x = self.wf.x(i);
            let dx = x - x0;
            let envelope = (-dx * dx / (4.0 * sigma * sigma)).exp();
            let phase = Complex64::new(0.0, k0 * x).exp();
            self.wf.up[i] = spatial_norm * envelope * phase;
            self.wf.down[i] = Complex64::new(0.0, 0.0);
        }

        // FFT the upper component
        use rustfft::FftPlanner;
        let mut planner = FftPlanner::new();
        let fft = planner.plan_fft_forward(n);
        let ifft = planner.plan_fft_inverse(n);
        let mut scratch = vec![Complex64::new(0.0, 0.0); fft.get_inplace_scratch_len().max(ifft.get_inplace_scratch_len())];

        fft.process_with_scratch(&mut self.wf.up, &mut scratch);

        // At each k, project onto the positive-energy spinor:
        // u(k) ∝ (1, ck/(E_k + mc²))
        // The upper component amplitude f(k) gets multiplied by the spinor:
        //   up[k]   = f(k) · 1 / norm_factor
        //   down[k] = f(k) · ck/(E_k + mc²) / norm_factor
        // where norm_factor = √(1 + (ck/(E_k+mc²))²) to preserve the total norm.
        for i in 0..n {
            let j = if i <= n / 2 {
                i as f64
            } else {
                i as f64 - n as f64
            };
            let k = 2.0 * PI * j / length;
            let ck = self.c * k;
            let e_k = (ck * ck + mc2 * mc2).sqrt();

            let ratio = if e_k + mc2 > 1e-30 {
                ck / (e_k + mc2)
            } else {
                0.0
            };
            let norm_factor = (1.0 + ratio * ratio).sqrt();

            let f = self.wf.up[i];
            self.wf.up[i] = f / norm_factor;
            self.wf.down[i] = f * ratio / norm_factor;
        }

        // IFFT both components
        ifft.process_with_scratch(&mut self.wf.up, &mut scratch);
        ifft.process_with_scratch(&mut self.wf.down, &mut scratch);

        // Normalize (rustfft convention)
        let inv_n = 1.0 / n as f64;
        for c in self.wf.up.iter_mut() {
            *c *= inv_n;
        }
        for c in self.wf.down.iter_mut() {
            *c *= inv_n;
        }

        self.wf.normalize();
    }
}
