/// Potential energy function V(x).
///
/// Multiple potential types can be combined (summed) to build complex scenarios.
#[derive(Clone)]
pub struct Potential {
    terms: Vec<PotentialTerm>,
}

#[derive(Clone)]
pub enum PotentialTerm {
    /// No potential (free particle).
    Free,
    /// Infinite square well approximated by steep walls.
    /// V = 0 for x in [left, right], V = wall_height outside.
    Well {
        left: f64,
        right: f64,
        wall_height: f64,
    },
    /// Rectangular barrier: V = height for x in [left, right], 0 elsewhere.
    Barrier {
        left: f64,
        right: f64,
        height: f64,
    },
    /// Harmonic oscillator: V = ½ω²(x - x0)².
    Harmonic { x0: f64, omega: f64 },
    /// Double well: two harmonic wells separated by a barrier.
    /// V = min(½ω²(x - x1)², ½ω²(x - x2)²) capped at barrier_height.
    DoubleWell {
        x1: f64,
        x2: f64,
        omega: f64,
        barrier_height: f64,
    },
    /// Arbitrary potential defined by a closure.
    Custom(CustomPotential),
}

/// Wrapper to make closures Clone-able by boxing and using Arc.
#[derive(Clone)]
pub struct CustomPotential {
    func: std::sync::Arc<dyn Fn(f64) -> f64 + Send + Sync>,
}

impl CustomPotential {
    pub fn new(f: impl Fn(f64) -> f64 + Send + Sync + 'static) -> Self {
        Self {
            func: std::sync::Arc::new(f),
        }
    }

    pub fn eval(&self, x: f64) -> f64 {
        (self.func)(x)
    }
}

impl PotentialTerm {
    fn value_at(&self, x: f64) -> f64 {
        match self {
            PotentialTerm::Free => 0.0,
            PotentialTerm::Well {
                left,
                right,
                wall_height,
            } => {
                if x >= *left && x <= *right {
                    0.0
                } else {
                    *wall_height
                }
            }
            PotentialTerm::Barrier {
                left,
                right,
                height,
            } => {
                if x >= *left && x <= *right {
                    *height
                } else {
                    0.0
                }
            }
            PotentialTerm::Harmonic { x0, omega } => {
                let dx = x - x0;
                0.5 * omega * omega * dx * dx
            }
            PotentialTerm::DoubleWell {
                x1,
                x2,
                omega,
                barrier_height,
            } => {
                let v1 = 0.5 * omega * omega * (x - x1).powi(2);
                let v2 = 0.5 * omega * omega * (x - x2).powi(2);
                v1.min(v2).min(*barrier_height)
            }
            PotentialTerm::Custom(cp) => cp.eval(x),
        }
    }
}

impl Potential {
    pub fn free() -> Self {
        Self {
            terms: vec![PotentialTerm::Free],
        }
    }

    pub fn well(left: f64, right: f64, wall_height: f64) -> Self {
        Self {
            terms: vec![PotentialTerm::Well {
                left,
                right,
                wall_height,
            }],
        }
    }

    pub fn barrier(left: f64, right: f64, height: f64) -> Self {
        Self {
            terms: vec![PotentialTerm::Barrier {
                left,
                right,
                height,
            }],
        }
    }

    pub fn harmonic(x0: f64, omega: f64) -> Self {
        Self {
            terms: vec![PotentialTerm::Harmonic { x0, omega }],
        }
    }

    pub fn double_well(x1: f64, x2: f64, omega: f64, barrier_height: f64) -> Self {
        Self {
            terms: vec![PotentialTerm::DoubleWell {
                x1,
                x2,
                omega,
                barrier_height,
            }],
        }
    }

    /// Evaluate the total potential at position x.
    pub fn value_at(&self, x: f64) -> f64 {
        self.terms.iter().map(|t| t.value_at(x)).sum()
    }

    /// Evaluate the potential at each grid point.
    pub fn values_on_grid(&self, xs: &[f64]) -> Vec<f64> {
        xs.iter().map(|&x| self.value_at(x)).collect()
    }

    /// Add an additional potential term.
    pub fn add_term(&mut self, term: PotentialTerm) {
        self.terms.push(term);
    }

    /// The current list of terms (for UI display).
    pub fn terms(&self) -> &[PotentialTerm] {
        &self.terms
    }
}
