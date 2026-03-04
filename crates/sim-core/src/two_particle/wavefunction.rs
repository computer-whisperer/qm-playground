use num_complex::Complex64;
use std::f64::consts::PI;

/// A 2D wavefunction ψ(x₁, x₂) on a uniform grid.
///
/// Both particles share the same spatial grid (x_min, x_max, n points).
/// Stored in row-major order: psi[i1 * n + i2] = ψ(x₁_i1, x₂_i2).
pub struct Wavefunction2D {
    pub psi: Vec<Complex64>,
    /// Grid points per axis (total storage = n²).
    pub n: usize,
    pub x_min: f64,
    pub x_max: f64,
    pub dx: f64,
}

impl Wavefunction2D {
    /// Create a zero wavefunction on an n×n grid.
    pub fn new(n: usize, x_min: f64, x_max: f64) -> Self {
        let dx = (x_max - x_min) / n as f64;
        Self {
            psi: vec![Complex64::new(0.0, 0.0); n * n],
            n,
            x_min,
            x_max,
            dx,
        }
    }

    /// Position of the i-th grid point along either axis.
    pub fn x(&self, i: usize) -> f64 {
        self.x_min + (i as f64 + 0.5) * self.dx
    }

    /// All grid positions along one axis.
    pub fn xs(&self) -> Vec<f64> {
        (0..self.n).map(|i| self.x(i)).collect()
    }

    /// Initialize as a product of two Gaussian wave packets (unentangled).
    ///
    /// ψ(x₁, x₂) = φ₁(x₁) · φ₂(x₂) where each factor is a normalized Gaussian:
    ///   φ(x) = (2πσ²)^{-1/4} exp(-(x-x0)²/(4σ²)) exp(i k0 x)
    pub fn set_product_gaussian(
        &mut self,
        x0_1: f64,
        sigma_1: f64,
        k0_1: f64,
        x0_2: f64,
        sigma_2: f64,
        k0_2: f64,
    ) {
        let norm_1 = (2.0 * PI * sigma_1 * sigma_1).powf(-0.25);
        let norm_2 = (2.0 * PI * sigma_2 * sigma_2).powf(-0.25);

        for i1 in 0..self.n {
            let x1 = self.x(i1);
            let dx1 = x1 - x0_1;
            let env1 = (-dx1 * dx1 / (4.0 * sigma_1 * sigma_1)).exp();
            let phase1 = Complex64::new(0.0, k0_1 * x1).exp();
            let phi1 = norm_1 * env1 * phase1;

            for i2 in 0..self.n {
                let x2 = self.x(i2);
                let dx2 = x2 - x0_2;
                let env2 = (-dx2 * dx2 / (4.0 * sigma_2 * sigma_2)).exp();
                let phase2 = Complex64::new(0.0, k0_2 * x2).exp();
                let phi2 = norm_2 * env2 * phase2;

                self.psi[i1 * self.n + i2] = phi1 * phi2;
            }
        }

        self.normalize();
    }

    /// Marginal probability density for particle 1: ρ₁(x₁) = ∫ |ψ(x₁, x₂)|² dx₂.
    pub fn marginal_1(&self) -> Vec<f64> {
        let n = self.n;
        let mut rho = vec![0.0; n];
        for i1 in 0..n {
            let mut sum = 0.0;
            for i2 in 0..n {
                sum += self.psi[i1 * n + i2].norm_sqr();
            }
            rho[i1] = sum * self.dx;
        }
        rho
    }

    /// Marginal probability density for particle 2: ρ₂(x₂) = ∫ |ψ(x₁, x₂)|² dx₁.
    pub fn marginal_2(&self) -> Vec<f64> {
        let n = self.n;
        let mut rho = vec![0.0; n];
        for i2 in 0..n {
            let mut sum = 0.0;
            for i1 in 0..n {
                sum += self.psi[i1 * n + i2].norm_sqr();
            }
            rho[i2] = sum * self.dx;
        }
        rho
    }

    /// 2D probability density |ψ(x₁, x₂)|² as a flat array (row-major, n×n).
    pub fn probability_density(&self) -> Vec<f64> {
        self.psi.iter().map(|c| c.norm_sqr()).collect()
    }

    /// Total probability ∫∫ |ψ|² dx₁ dx₂ (should be 1.0).
    pub fn norm(&self) -> f64 {
        self.psi.iter().map(|c| c.norm_sqr()).sum::<f64>() * self.dx * self.dx
    }

    /// Normalize the wavefunction so ∫∫ |ψ|² dx₁ dx₂ = 1.
    pub fn normalize(&mut self) {
        let n = self.norm().sqrt();
        if n > 0.0 {
            for c in &mut self.psi {
                *c /= n;
            }
        }
    }

    /// Expected position of particle 1: ⟨x₁⟩ = ∫∫ x₁ |ψ|² dx₁ dx₂.
    pub fn expected_x1(&self) -> f64 {
        let n = self.n;
        let mut sum = 0.0;
        for i1 in 0..n {
            let x1 = self.x(i1);
            for i2 in 0..n {
                sum += x1 * self.psi[i1 * n + i2].norm_sqr();
            }
        }
        sum * self.dx * self.dx
    }

    /// Expected position of particle 2: ⟨x₂⟩ = ∫∫ x₂ |ψ|² dx₁ dx₂.
    pub fn expected_x2(&self) -> f64 {
        let n = self.n;
        let mut sum = 0.0;
        for i1 in 0..n {
            for i2 in 0..n {
                let x2 = self.x(i2);
                sum += x2 * self.psi[i1 * n + i2].norm_sqr();
            }
        }
        sum * self.dx * self.dx
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn product_gaussian_is_normalized() {
        let mut wf = Wavefunction2D::new(128, -15.0, 15.0);
        wf.set_product_gaussian(0.0, 1.0, 0.0, 0.0, 1.0, 0.0);
        let norm = wf.norm();
        assert!((norm - 1.0).abs() < 1e-10, "norm = {norm}");
    }

    #[test]
    fn product_gaussian_expected_positions() {
        let mut wf = Wavefunction2D::new(128, -15.0, 15.0);
        wf.set_product_gaussian(2.0, 1.0, 0.0, -3.0, 1.0, 0.0);
        let ex1 = wf.expected_x1();
        let ex2 = wf.expected_x2();
        assert!((ex1 - 2.0).abs() < 0.01, "⟨x₁⟩ = {ex1}");
        assert!((ex2 + 3.0).abs() < 0.01, "⟨x₂⟩ = {ex2}");
    }

    #[test]
    fn marginals_are_normalized() {
        let mut wf = Wavefunction2D::new(128, -15.0, 15.0);
        wf.set_product_gaussian(0.0, 1.5, 2.0, 3.0, 1.0, -1.0);

        let m1 = wf.marginal_1();
        let m2 = wf.marginal_2();
        let sum1: f64 = m1.iter().sum::<f64>() * wf.dx;
        let sum2: f64 = m2.iter().sum::<f64>() * wf.dx;

        assert!((sum1 - 1.0).abs() < 1e-6, "∫ρ₁ dx₁ = {sum1}");
        assert!((sum2 - 1.0).abs() < 1e-6, "∫ρ₂ dx₂ = {sum2}");
    }
}
