pub mod entanglement;
pub mod integrator;
pub mod potential;
pub mod wavefunction;

use integrator::SplitOperator2D;
use num_complex::Complex64;
use potential::Potential2D;
use wavefunction::Wavefunction2D;

pub use wavefunction::ParticleSymmetry;

/// Complete two-particle simulation state.
pub struct TwoParticleSimulation {
    pub wf: Wavefunction2D,
    pub potential: Potential2D,
    pub integrator: SplitOperator2D,
    pub time: f64,
    pub dt: f64,
}

impl TwoParticleSimulation {
    pub fn new(n: usize, x_min: f64, x_max: f64, dt: f64, potential: Potential2D) -> Self {
        let wf = Wavefunction2D::new(n, x_min, x_max);
        let integrator = SplitOperator2D::new(&wf, &potential, dt);
        Self {
            wf,
            potential,
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

    /// Recompute integrator arrays (call after changing potential or dt).
    pub fn rebuild_integrator(&mut self) {
        self.integrator = SplitOperator2D::new(&self.wf, &self.potential, self.dt);
    }

    pub fn norm(&self) -> f64 {
        self.wf.norm()
    }

    pub fn purity(&self) -> f64 {
        entanglement::purity(&self.wf)
    }

    pub fn expected_x1(&self) -> f64 {
        self.wf.expected_x1()
    }

    pub fn expected_x2(&self) -> f64 {
        self.wf.expected_x2()
    }

    /// Total energy ⟨H⟩ = ⟨T⟩ + ⟨V⟩ via direct summation.
    ///
    /// ⟨T⟩ uses the kinetic energy operator in momentum space (computed via finite differences).
    /// ⟨V⟩ = ∫∫ V(x₁,x₂) |ψ|² dx₁ dx₂.
    pub fn expected_energy(&self) -> f64 {
        let n = self.wf.n;
        let dx = self.wf.dx;
        let dx2 = dx * dx;
        let psi = &self.wf.psi;

        // ⟨V⟩
        let xs = self.wf.xs();
        let mut v_sum = 0.0;
        for i1 in 0..n {
            for i2 in 0..n {
                let v = self.potential.value_at(xs[i1], xs[i2]);
                v_sum += v * psi[i1 * n + i2].norm_sqr();
            }
        }
        let expected_v = v_sum * dx2;

        // ⟨T⟩ via finite differences: T = -½(∂²/∂x₁² + ∂²/∂x₂²)
        // ⟨T⟩ = -½ ∫∫ ψ* (∂²ψ/∂x₁² + ∂²ψ/∂x₂²) dx₁ dx₂
        let mut t_sum = Complex64::new(0.0, 0.0);
        for i1 in 0..n {
            let i1p = (i1 + 1) % n;
            let i1m = (i1 + n - 1) % n;
            for i2 in 0..n {
                let i2p = (i2 + 1) % n;
                let i2m = (i2 + n - 1) % n;

                let idx = i1 * n + i2;
                let psi_c = psi[idx];

                // ∂²ψ/∂x₁²
                let d2_x1 = (psi[i1p * n + i2] - 2.0 * psi_c + psi[i1m * n + i2]) / (dx * dx);
                // ∂²ψ/∂x₂²
                let d2_x2 = (psi[i1 * n + i2p] - 2.0 * psi_c + psi[i1 * n + i2m]) / (dx * dx);

                t_sum += psi_c.conj() * (d2_x1 + d2_x2);
            }
        }
        let expected_t = (-0.5 * t_sum * dx2).re;

        expected_t + expected_v
    }
}
