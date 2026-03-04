use num_complex::Complex64;
use std::f64::consts::PI;

/// A spin-1/2 wavefunction on a 1D uniform grid.
///
/// The state is a 2-component spinor: [ψ_↑(x), ψ_↓(x)].
/// Both components share the same spatial grid.
pub struct SpinorWavefunction {
    /// Spin-up component ψ_↑(x).
    pub up: Vec<Complex64>,
    /// Spin-down component ψ_↓(x).
    pub down: Vec<Complex64>,
    pub n: usize,
    pub x_min: f64,
    pub x_max: f64,
    pub dx: f64,
}

impl SpinorWavefunction {
    pub fn new(n: usize, x_min: f64, x_max: f64) -> Self {
        let dx = (x_max - x_min) / n as f64;
        Self {
            up: vec![Complex64::new(0.0, 0.0); n],
            down: vec![Complex64::new(0.0, 0.0); n],
            n,
            x_min,
            x_max,
            dx,
        }
    }

    /// Position of the i-th grid point.
    pub fn x(&self, i: usize) -> f64 {
        self.x_min + (i as f64 + 0.5) * self.dx
    }

    /// All grid positions.
    pub fn xs(&self) -> Vec<f64> {
        (0..self.n).map(|i| self.x(i)).collect()
    }

    /// Initialize as a Gaussian wave packet with a definite spin state.
    ///
    /// `spin_up_amp` and `spin_down_amp` set the spin superposition (will be normalized).
    /// For pure spin-up: (1, 0). For spin-down: (0, 1). For +x: (1, 1).
    pub fn set_gaussian(
        &mut self,
        x0: f64,
        sigma: f64,
        k0: f64,
        spin_up_amp: Complex64,
        spin_down_amp: Complex64,
    ) {
        let spatial_norm = (2.0 * PI * sigma * sigma).powf(-0.25);
        let spin_norm = (spin_up_amp.norm_sqr() + spin_down_amp.norm_sqr()).sqrt();
        let alpha = spin_up_amp / spin_norm;
        let beta = spin_down_amp / spin_norm;

        for i in 0..self.n {
            let x = self.x(i);
            let dx = x - x0;
            let envelope = (-dx * dx / (4.0 * sigma * sigma)).exp();
            let phase = Complex64::new(0.0, k0 * x).exp();
            let spatial = spatial_norm * envelope * phase;

            self.up[i] = alpha * spatial;
            self.down[i] = beta * spatial;
        }
        self.normalize();
    }

    /// Total probability ∫ (|ψ_↑|² + |ψ_↓|²) dx.
    pub fn norm(&self) -> f64 {
        let sum: f64 = self
            .up
            .iter()
            .zip(self.down.iter())
            .map(|(u, d)| u.norm_sqr() + d.norm_sqr())
            .sum();
        sum * self.dx
    }

    /// Normalize so total probability = 1.
    pub fn normalize(&mut self) {
        let n = self.norm().sqrt();
        if n > 0.0 {
            for c in &mut self.up {
                *c /= n;
            }
            for c in &mut self.down {
                *c /= n;
            }
        }
    }

    /// Total probability density |ψ_↑(x)|² + |ψ_↓(x)|² at each grid point.
    pub fn probability_density(&self) -> Vec<f64> {
        self.up
            .iter()
            .zip(self.down.iter())
            .map(|(u, d)| u.norm_sqr() + d.norm_sqr())
            .collect()
    }

    /// Spin-up probability density |ψ_↑(x)|².
    pub fn density_up(&self) -> Vec<f64> {
        self.up.iter().map(|c| c.norm_sqr()).collect()
    }

    /// Spin-down probability density |ψ_↓(x)|².
    pub fn density_down(&self) -> Vec<f64> {
        self.down.iter().map(|c| c.norm_sqr()).collect()
    }

    /// Expected position ⟨x⟩ = ∫ x (|ψ_↑|² + |ψ_↓|²) dx.
    pub fn expected_x(&self) -> f64 {
        let mut sum = 0.0;
        for i in 0..self.n {
            sum += self.x(i) * (self.up[i].norm_sqr() + self.down[i].norm_sqr());
        }
        sum * self.dx
    }

    /// Expected spin along z: ⟨σ_z⟩ = ∫ (|ψ_↑|² - |ψ_↓|²) dx.
    pub fn expected_sz(&self) -> f64 {
        let sum: f64 = self
            .up
            .iter()
            .zip(self.down.iter())
            .map(|(u, d)| u.norm_sqr() - d.norm_sqr())
            .sum();
        sum * self.dx
    }

    /// Expected spin along x: ⟨σ_x⟩ = 2 Re ∫ ψ_↑* ψ_↓ dx.
    pub fn expected_sx(&self) -> f64 {
        let sum: Complex64 = self
            .up
            .iter()
            .zip(self.down.iter())
            .map(|(u, d)| u.conj() * d)
            .sum();
        2.0 * (sum * self.dx).re
    }

    /// Expected spin along y: ⟨σ_y⟩ = 2 Im ∫ ψ_↑* ψ_↓ dx.
    pub fn expected_sy(&self) -> f64 {
        let sum: Complex64 = self
            .up
            .iter()
            .zip(self.down.iter())
            .map(|(u, d)| u.conj() * d)
            .sum();
        2.0 * (sum * self.dx).im
    }

    /// Spin-up probability: ∫ |ψ_↑|² dx.
    pub fn prob_up(&self) -> f64 {
        self.up.iter().map(|c| c.norm_sqr()).sum::<f64>() * self.dx
    }

    /// Spin-down probability: ∫ |ψ_↓|² dx.
    pub fn prob_down(&self) -> f64 {
        self.down.iter().map(|c| c.norm_sqr()).sum::<f64>() * self.dx
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spin_up_gaussian_is_normalized() {
        let mut wf = SpinorWavefunction::new(1024, -20.0, 20.0);
        wf.set_gaussian(0.0, 1.0, 0.0, Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0));
        assert!((wf.norm() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn pure_spin_up_expectations() {
        let mut wf = SpinorWavefunction::new(1024, -20.0, 20.0);
        wf.set_gaussian(0.0, 1.0, 0.0, Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0));

        assert!((wf.expected_sz() - 1.0).abs() < 1e-10, "⟨σ_z⟩ = {}", wf.expected_sz());
        assert!(wf.expected_sx().abs() < 1e-10);
        assert!(wf.expected_sy().abs() < 1e-10);
        assert!((wf.prob_up() - 1.0).abs() < 1e-10);
        assert!(wf.prob_down().abs() < 1e-10);
    }

    #[test]
    fn superposition_expectations() {
        let mut wf = SpinorWavefunction::new(1024, -20.0, 20.0);
        // |+x⟩ = (|↑⟩ + |↓⟩)/√2
        wf.set_gaussian(0.0, 1.0, 0.0, Complex64::new(1.0, 0.0), Complex64::new(1.0, 0.0));

        assert!(wf.expected_sz().abs() < 1e-10, "⟨σ_z⟩ = {}", wf.expected_sz());
        assert!((wf.expected_sx() - 1.0).abs() < 1e-10, "⟨σ_x⟩ = {}", wf.expected_sx());
        assert!(wf.expected_sy().abs() < 1e-10, "⟨σ_y⟩ = {}", wf.expected_sy());
        assert!((wf.prob_up() - 0.5).abs() < 1e-10);
        assert!((wf.prob_down() - 0.5).abs() < 1e-10);
    }
}
