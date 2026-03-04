pub mod integrator;
pub mod potential;
pub mod spinor;
pub mod two_particle;
pub mod units;
pub mod wavefunction;

use integrator::SplitOperator;
use potential::Potential;
use wavefunction::Wavefunction;

/// Complete simulation state: wavefunction + potential + integrator.
pub struct Simulation {
    pub wf: Wavefunction,
    pub potential: Potential,
    pub integrator: SplitOperator,
    pub time: f64,
    pub dt: f64,
}

impl Simulation {
    pub fn new(n: usize, x_min: f64, x_max: f64, dt: f64, potential: Potential) -> Self {
        let wf = Wavefunction::new(n, x_min, x_max);
        let integrator = SplitOperator::new(&wf, &potential, dt);
        Self {
            wf,
            potential,
            integrator,
            time: 0.0,
            dt,
        }
    }

    /// Advance the simulation by one time step.
    pub fn step(&mut self) {
        self.integrator.step(&mut self.wf);
        self.time += self.dt;
    }

    /// Advance the simulation by n time steps.
    pub fn step_n(&mut self, n: usize) {
        for _ in 0..n {
            self.step();
        }
    }

    /// Recompute integrator arrays (call after changing potential or dt).
    pub fn rebuild_integrator(&mut self) {
        self.integrator = SplitOperator::new(&self.wf, &self.potential, self.dt);
    }

    /// Total probability (should stay ≈ 1.0 if things are working).
    pub fn norm(&self) -> f64 {
        self.wf.norm()
    }

    /// Expected position ⟨x⟩.
    pub fn expected_x(&self) -> f64 {
        self.wf.expected_x()
    }

    /// Expected momentum ⟨p⟩ (computed via finite differences).
    pub fn expected_p(&self) -> f64 {
        self.wf.expected_p()
    }

    /// Total energy ⟨H⟩ = ⟨T⟩ + ⟨V⟩.
    pub fn expected_energy(&self) -> f64 {
        self.wf.expected_kinetic_energy() + self.wf.expected_potential_energy(&self.potential)
    }
}
