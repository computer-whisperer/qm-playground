use num_complex::Complex64;
use std::f64::consts::PI;

use crate::potential::Potential;

/// A 1D wavefunction on a uniform grid.
///
/// ψ is stored as a vector of complex amplitudes at each grid point.
/// The grid runs from x_min to x_max with n points and spacing dx.
pub struct Wavefunction {
    pub psi: Vec<Complex64>,
    pub n: usize,
    pub x_min: f64,
    pub x_max: f64,
    pub dx: f64,
}

impl Wavefunction {
    /// Create a zero wavefunction on a grid of n points.
    pub fn new(n: usize, x_min: f64, x_max: f64) -> Self {
        let dx = (x_max - x_min) / n as f64;
        Self {
            psi: vec![Complex64::new(0.0, 0.0); n],
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

    /// Initialize as a Gaussian wave packet.
    ///
    /// ψ(x) = (2πσ²)^{-1/4} exp(-(x-x0)²/(4σ²)) exp(i k0 x)
    ///
    /// - x0: center position
    /// - sigma: width (spatial uncertainty)
    /// - k0: central wave number (determines average momentum p = ℏk0 = k0 in atomic units)
    pub fn set_gaussian(&mut self, x0: f64, sigma: f64, k0: f64) {
        let norm = (2.0 * PI * sigma * sigma).powf(-0.25);
        for i in 0..self.n {
            let x = self.x(i);
            let dx = x - x0;
            let envelope = (-dx * dx / (4.0 * sigma * sigma)).exp();
            let phase = Complex64::new(0.0, k0 * x).exp();
            self.psi[i] = norm * envelope * phase;
        }
        // Renormalize to ensure ∫|ψ|² dx = 1 on the discrete grid.
        self.normalize();
    }

    /// Probability density |ψ(x)|² at each grid point.
    pub fn probability_density(&self) -> Vec<f64> {
        self.psi.iter().map(|c| c.norm_sqr()).collect()
    }

    /// Total probability ∫|ψ|² dx (should be 1.0).
    pub fn norm(&self) -> f64 {
        self.psi.iter().map(|c| c.norm_sqr()).sum::<f64>() * self.dx
    }

    /// Normalize the wavefunction so ∫|ψ|² dx = 1.
    pub fn normalize(&mut self) {
        let n = self.norm().sqrt();
        if n > 0.0 {
            for c in &mut self.psi {
                *c /= n;
            }
        }
    }

    /// Expected position ⟨x⟩ = ∫ x |ψ|² dx.
    pub fn expected_x(&self) -> f64 {
        let mut sum = 0.0;
        for i in 0..self.n {
            sum += self.x(i) * self.psi[i].norm_sqr();
        }
        sum * self.dx
    }

    /// Expected momentum ⟨p⟩ = -i ∫ ψ* (∂ψ/∂x) dx.
    ///
    /// Computed via central finite differences.
    pub fn expected_p(&self) -> f64 {
        let mut sum = Complex64::new(0.0, 0.0);
        for i in 0..self.n {
            // Central difference with periodic boundary
            let ip = (i + 1) % self.n;
            let im = (i + self.n - 1) % self.n;
            let dpsi_dx = (self.psi[ip] - self.psi[im]) / (2.0 * self.dx);
            sum += self.psi[i].conj() * dpsi_dx;
        }
        // ⟨p⟩ = -i ∫ ψ* ∂ψ/∂x dx  (in atomic units, ℏ = 1)
        let result = Complex64::new(0.0, -1.0) * sum * self.dx;
        result.re
    }

    /// Expected kinetic energy ⟨T⟩ = -½ ∫ ψ* (∂²ψ/∂x²) dx.
    ///
    /// Computed via central finite differences.
    pub fn expected_kinetic_energy(&self) -> f64 {
        let mut sum = Complex64::new(0.0, 0.0);
        for i in 0..self.n {
            let ip = (i + 1) % self.n;
            let im = (i + self.n - 1) % self.n;
            let d2psi = (self.psi[ip] - 2.0 * self.psi[i] + self.psi[im]) / (self.dx * self.dx);
            sum += self.psi[i].conj() * d2psi;
        }
        // ⟨T⟩ = -½ ∫ ψ* ∂²ψ/∂x² dx
        let result = -0.5 * sum * self.dx;
        result.re
    }

    /// Expected potential energy ⟨V⟩ = ∫ V(x) |ψ|² dx.
    pub fn expected_potential_energy(&self, potential: &Potential) -> f64 {
        let mut sum = 0.0;
        for i in 0..self.n {
            sum += potential.value_at(self.x(i)) * self.psi[i].norm_sqr();
        }
        sum * self.dx
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gaussian_is_normalized() {
        let mut wf = Wavefunction::new(1024, -20.0, 20.0);
        wf.set_gaussian(0.0, 1.0, 0.0);
        let norm = wf.norm();
        assert!((norm - 1.0).abs() < 1e-10, "norm = {norm}");
    }

    #[test]
    fn gaussian_centered_at_x0() {
        let mut wf = Wavefunction::new(1024, -20.0, 20.0);
        wf.set_gaussian(3.0, 1.0, 0.0);
        let ex = wf.expected_x();
        assert!((ex - 3.0).abs() < 1e-6, "⟨x⟩ = {ex}");
    }

    #[test]
    fn gaussian_with_momentum() {
        let mut wf = Wavefunction::new(2048, -30.0, 30.0);
        let k0 = 5.0;
        wf.set_gaussian(0.0, 2.0, k0);
        let ep = wf.expected_p();
        assert!(
            (ep - k0).abs() < 0.05,
            "⟨p⟩ = {ep}, expected {k0}"
        );
    }
}
