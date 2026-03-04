use num_complex::Complex64;

use super::wavefunction::Wavefunction2D;

/// Compute the purity Tr(ρ₁²) of the reduced density matrix for particle 1.
///
/// For a product (unentangled) state, purity = 1.0.
/// For a maximally entangled state of Schmidt rank r, purity = 1/r.
///
/// Computed directly as:
///   Tr(ρ₁²) = Σ_{i,k} |Σ_j ψ(i,j) ψ*(k,j)|² dx⁴
///
/// This is O(n³) — call sparingly (e.g., every ~10 frames).
pub fn purity(wf: &Wavefunction2D) -> f64 {
    let n = wf.n;
    let dx = wf.dx;
    let dx4 = dx * dx * dx * dx;
    let psi = &wf.psi;

    let mut sum = 0.0;
    for i in 0..n {
        for k in i..n {
            // Compute ρ₁(i, k) = Σ_j ψ(i,j) ψ*(k,j) dx
            let mut dot = Complex64::new(0.0, 0.0);
            let row_i = i * n;
            let row_k = k * n;
            for j in 0..n {
                dot += psi[row_i + j] * psi[row_k + j].conj();
            }
            let ns = dot.norm_sqr();
            if i == k {
                sum += ns; // diagonal: counted once
            } else {
                sum += 2.0 * ns; // off-diagonal: ρ is Hermitian, so (i,k) and (k,i) contribute equally
            }
        }
    }
    sum * dx4
}

/// Schmidt number K = 1 / Purity.
///
/// K = 1 for product states, K > 1 for entangled states.
pub fn schmidt_number(wf: &Wavefunction2D) -> f64 {
    1.0 / purity(wf)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn product_state_purity_is_one() {
        let mut wf = Wavefunction2D::new(128, -15.0, 15.0);
        wf.set_product_gaussian(0.0, 1.5, 2.0, 0.0, 1.5, -2.0);

        let p = purity(&wf);
        assert!(
            (p - 1.0).abs() < 1e-4,
            "product state purity = {p}, expected ~1.0"
        );
    }

    #[test]
    fn purity_drops_after_interacting_evolution() {
        use crate::two_particle::potential::{Interaction, Potential2D};
        use crate::two_particle::TwoParticleSimulation;
        use crate::potential::Potential;

        let pot = Potential2D::new(
            Potential::free(),
            Potential::free(),
            Interaction::SoftCoulomb { g: 2.0, epsilon: 0.5 },
        );
        let mut sim = TwoParticleSimulation::new(128, -15.0, 15.0, 0.005, pot);
        sim.wf.set_product_gaussian(-4.0, 1.5, 3.0, 4.0, 1.5, -3.0);

        let p0 = purity(&sim.wf);
        assert!(
            (p0 - 1.0).abs() < 1e-4,
            "initial purity = {p0}, expected ~1.0"
        );

        // Evolve until packets overlap and interact
        sim.step_n(800);

        let p1 = purity(&sim.wf);
        assert!(
            p1 < 0.99,
            "purity after interaction = {p1}, expected < 1.0 (entanglement)"
        );
    }
}
