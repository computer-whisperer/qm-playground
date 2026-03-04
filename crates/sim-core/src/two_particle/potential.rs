use std::f64::consts::PI;

use crate::potential::Potential;

/// Interaction potential between two particles.
#[derive(Clone)]
pub enum Interaction {
    /// No interaction (particles are independent).
    None,
    /// Soft-Coulomb repulsion: V_int = g / sqrt((x₁ - x₂)² + ε²).
    /// The softening parameter ε prevents the 1/r singularity.
    SoftCoulomb { g: f64, epsilon: f64 },
    /// Contact (delta-like) interaction: V_int = g · δ_w(x₁ - x₂)
    /// approximated as a narrow Gaussian of width `width`.
    Contact { g: f64, width: f64 },
}

impl Interaction {
    pub fn value_at(&self, x1: f64, x2: f64) -> f64 {
        match self {
            Interaction::None => 0.0,
            Interaction::SoftCoulomb { g, epsilon } => {
                let dx = x1 - x2;
                g / (dx * dx + epsilon * epsilon).sqrt()
            }
            Interaction::Contact { g, width } => {
                let dx = x1 - x2;
                let w = *width;
                g * (-dx * dx / (2.0 * w * w)).exp() / (w * (2.0 * PI).sqrt())
            }
        }
    }
}

/// Two-particle potential: V(x₁, x₂) = V₁(x₁) + V₂(x₂) + V_int(x₁, x₂).
#[derive(Clone)]
pub struct Potential2D {
    /// Single-particle potential for particle 1.
    pub v1: Potential,
    /// Single-particle potential for particle 2.
    pub v2: Potential,
    /// Interaction between the two particles.
    pub interaction: Interaction,
}

impl Potential2D {
    pub fn new(v1: Potential, v2: Potential, interaction: Interaction) -> Self {
        Self { v1, v2, interaction }
    }

    /// Evaluate V(x₁, x₂) = V₁(x₁) + V₂(x₂) + V_int(x₁, x₂).
    pub fn value_at(&self, x1: f64, x2: f64) -> f64 {
        self.v1.value_at(x1) + self.v2.value_at(x2) + self.interaction.value_at(x1, x2)
    }

    /// Precompute V on the full n×n grid (row-major).
    pub fn values_on_grid(&self, xs: &[f64]) -> Vec<f64> {
        let n = xs.len();
        let mut vals = vec![0.0; n * n];
        for (i1, &x1) in xs.iter().enumerate() {
            for (i2, &x2) in xs.iter().enumerate() {
                vals[i1 * n + i2] = self.value_at(x1, x2);
            }
        }
        vals
    }
}
