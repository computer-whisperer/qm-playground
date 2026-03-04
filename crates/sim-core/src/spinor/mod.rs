pub mod integrator;
pub mod wavefunction;

use integrator::{MagneticField, SpinorIntegrator};
use num_complex::Complex64;
use wavefunction::SpinorWavefunction;

use crate::potential::Potential;

/// Complete spin-1/2 single-particle simulation state.
pub struct SpinorSimulation {
    pub wf: SpinorWavefunction,
    pub potential: Potential,
    pub field: MagneticField,
    pub integrator: SpinorIntegrator,
    pub time: f64,
    pub dt: f64,
}

impl SpinorSimulation {
    pub fn new(
        n: usize,
        x_min: f64,
        x_max: f64,
        dt: f64,
        potential: Potential,
        field: MagneticField,
    ) -> Self {
        let wf = SpinorWavefunction::new(n, x_min, x_max);
        let integrator = SpinorIntegrator::new(&wf, &potential, &field, dt);
        Self {
            wf,
            potential,
            field,
            integrator,
            time: 0.0,
            dt,
        }
    }

    pub fn step(&mut self) {
        self.integrator.step(&mut self.wf);
        self.time += self.dt;
    }

    pub fn step_n(&mut self, n: usize) {
        for _ in 0..n {
            self.step();
        }
    }

    /// Recompute integrator arrays (call after changing potential, field, or dt).
    pub fn rebuild_integrator(&mut self) {
        self.integrator =
            SpinorIntegrator::new(&self.wf, &self.potential, &self.field, self.dt);
    }

    pub fn norm(&self) -> f64 {
        self.wf.norm()
    }

    pub fn expected_x(&self) -> f64 {
        self.wf.expected_x()
    }

    /// Expected momentum ⟨p⟩ via finite differences (sum of both components).
    pub fn expected_p(&self) -> f64 {
        let n = self.wf.n;
        let dx = self.wf.dx;
        let mut sum = Complex64::new(0.0, 0.0);
        for i in 0..n {
            let ip = (i + 1) % n;
            let im = (i + n - 1) % n;
            // spin-up contribution
            let du = (self.wf.up[ip] - self.wf.up[im]) / (2.0 * dx);
            sum += self.wf.up[i].conj() * du;
            // spin-down contribution
            let dd = (self.wf.down[ip] - self.wf.down[im]) / (2.0 * dx);
            sum += self.wf.down[i].conj() * dd;
        }
        (Complex64::new(0.0, -1.0) * sum * dx).re
    }

    /// Expected energy ⟨H⟩ = ⟨T⟩ + ⟨V⟩ + ⟨H_mag⟩.
    pub fn expected_energy(&self) -> f64 {
        let n = self.wf.n;
        let dx = self.wf.dx;

        // ⟨T⟩: kinetic energy via finite differences (both components)
        let mut t_sum = Complex64::new(0.0, 0.0);
        for i in 0..n {
            let ip = (i + 1) % n;
            let im = (i + n - 1) % n;
            let d2u = (self.wf.up[ip] - 2.0 * self.wf.up[i] + self.wf.up[im]) / (dx * dx);
            let d2d =
                (self.wf.down[ip] - 2.0 * self.wf.down[i] + self.wf.down[im]) / (dx * dx);
            t_sum += self.wf.up[i].conj() * d2u + self.wf.down[i].conj() * d2d;
        }
        let expected_t = (-0.5 * t_sum * dx).re;

        // ⟨V⟩: potential energy (same for both components)
        let mut v_sum = 0.0;
        for i in 0..n {
            let v = self.potential.value_at(self.wf.x(i));
            v_sum += v * (self.wf.up[i].norm_sqr() + self.wf.down[i].norm_sqr());
        }
        let expected_v = v_sum * dx;

        // ⟨H_mag⟩: magnetic energy
        // B_z term: (B_z/2)(|ψ_↑|² - |ψ_↓|²) = (B_z/2)⟨σ_z⟩
        let mag_z = self.field.bz / 2.0 * self.wf.expected_sz();
        // B_x term: (B_x/2)⟨σ_x⟩
        let mag_x = self.field.bx / 2.0 * self.wf.expected_sx();

        expected_t + expected_v + mag_z + mag_x
    }
}
