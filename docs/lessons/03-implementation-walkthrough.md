# How the Identical Particle Simulation Works

A code-level walkthrough of the Stage 3 additions. Read alongside the
[physics primer](03-identical-particles.md) for the theory.

The good news: most of the infrastructure from Chapter 2 is unchanged. The
integrator, the grid, the entanglement measurement, the heatmap rendering —
all identical. What's new is initialization and a symmetry projection utility.

## What Changed

| Component | Chapter 2 | Chapter 3 |
|-----------|-----------|-----------|
| Wavefunction type | `Wavefunction2D` | Same |
| Initialization | `set_product_gaussian` | `set_symmetrized_gaussian` |
| Integrator | `SplitOperator2D` | Same (unchanged) |
| Entanglement | `purity()`, `schmidt_number()` | Same |
| New | — | `ParticleSymmetry` enum, `symmetrize()` |

The minimal diff reflects the physics: exchange symmetry is a constraint on the
*state*, not on the dynamics. The Hamiltonian doesn't change. The time evolution
doesn't change. Only the initial condition changes — and the consequences are
dramatic.

## The ParticleSymmetry Enum

```rust
enum ParticleSymmetry {
    Distinguishable,  // no constraint (Chapter 2 behavior)
    Boson,            // ψ(x₂, x₁) = +ψ(x₁, x₂)
    Fermion,          // ψ(x₂, x₁) = -ψ(x₁, x₂)
}
```

This is passed to the initialization function and stored in the viewer's app
state. The simulation engine itself doesn't track symmetry — it just evolves
whatever wavefunction it's given. Symmetry preservation is a consequence of the
Hamiltonian's structure, not enforcement by the integrator.

## Symmetrized Initialization

The key new function. Given two Gaussian parameter sets (a, b) and a symmetry type:

```rust
fn set_symmetrized_gaussian(&mut self,
    x0_a, sigma_a, k0_a,   // Gaussian "a"
    x0_b, sigma_b, k0_b,   // Gaussian "b"
    symmetry: ParticleSymmetry,
) {
    // For Distinguishable, this is just set_product_gaussian
    // For Boson/Fermion:
    let sign = match symmetry { Boson => +1.0, Fermion => -1.0 };

    for i1 in 0..n {
        let phi_a1 = gaussian(x(i1), x0_a, sigma_a, k0_a);  // "a" at position x₁
        let phi_b1 = gaussian(x(i1), x0_b, sigma_b, k0_b);  // "b" at position x₁

        for i2 in 0..n {
            let phi_a2 = gaussian(x(i2), x0_a, sigma_a, k0_a);  // "a" at position x₂
            let phi_b2 = gaussian(x(i2), x0_b, sigma_b, k0_b);  // "b" at position x₂

            // ψ = φ_a(x₁)φ_b(x₂) ± φ_a(x₂)φ_b(x₁)
            psi[i1 * n + i2] = phi_a1 * phi_b2 + sign * phi_a2 * phi_b1;
        }
    }
    self.normalize();
}
```

For each grid point (x₁, x₂), we compute both the "direct" term φ_a(x₁)φ_b(x₂)
and the "exchange" term φ_a(x₂)φ_b(x₁), then add them (bosons) or subtract them
(fermions).

### Why the direct construction instead of symmetrizing afterward?

We could also build a product state and then symmetrize it. The direct approach
avoids creating the unsymmetrized state and discarding half of it. More
importantly, it makes the Pauli exclusion case obvious: if a == b, the two
terms are identical, and the fermion subtraction gives exactly zero.

### The normalization factor

The physics primer says 1/√2, but the code just calls `normalize()`. Why?

For orthogonal single-particle states, the analytic normalization is indeed 1/√2.
But our Gaussians aren't generally orthogonal — two Gaussians at different
positions overlap. The discrete grid adds its own discretization error. Rather
than computing the exact overlap integral, we just normalize numerically. It's
one pass over the grid, and it's always correct.

### What happens when a == b (Pauli exclusion)

If both Gaussians have identical parameters (same x0, sigma, k0), then:

```
φ_a(x₁)φ_b(x₂) - φ_a(x₂)φ_b(x₁) = φ(x₁)φ(x₂) - φ(x₂)φ(x₁) = 0
```

Every element of psi is exactly zero. The norm is zero, so `normalize()` is a
no-op (it checks for zero norm). The viewer shows a blank heatmap with ‖ψ‖² = 0.
This is the Pauli exclusion principle: no wavefunction exists for two fermions
in the same state.

## The Re-Symmetrization Projection

In exact arithmetic, a symmetric Hamiltonian preserves symmetry perfectly. But
floating-point arithmetic introduces roundoff. Over thousands of steps, tiny
asymmetries can accumulate.

The `symmetrize()` method projects ψ back onto the symmetric or antisymmetric
subspace:

```rust
fn symmetrize(&mut self, symmetry: ParticleSymmetry) {
    let sign = match symmetry { Boson => +1.0, Fermion => -1.0 };

    for i1 in 0..n {
        for i2 in (i1+1)..n {
            let a = psi[i1 * n + i2];    // ψ(x₁, x₂)
            let b = psi[i2 * n + i1];    // ψ(x₂, x₁)
            let sym = (a + sign * b) / 2;
            psi[i1 * n + i2] = sym;
            psi[i2 * n + i1] = sign * sym;
        }
        // Fermion diagonal: force ψ(x, x) = 0
        if symmetry == Fermion {
            psi[i1 * n + i1] = 0;
        }
    }
    self.normalize();
}
```

This iterates only the upper triangle (i2 > i1) and sets each pair to the
average of what it should be. For bosons: ψ(x₁,x₂) = ψ(x₂,x₁) = (a+b)/2.
For fermions: ψ(x₁,x₂) = -ψ(x₂,x₁) = (a-b)/2, and the diagonal is zeroed.

Cost: O(n²), one pass. Cheap enough to run every frame if needed.

In practice, the test suite confirms that symmetry violation after 500 steps
of interacting evolution stays below 1e-10 — so re-symmetrization is available
as a safety net but rarely needed.

## Why the Integrator Doesn't Change

This deserves emphasis because it's the most important structural insight.

The split-operator method applies three operations per step:
1. exp(-iV dt/2) — pointwise multiply in position space
2. exp(-iT dt) — pointwise multiply in momentum space
3. exp(-iV dt/2) — pointwise multiply again

Each of these is a **multiplication by a function that's symmetric under
particle exchange**, when the Hamiltonian is exchange-symmetric:

- V(x₁, x₂) = V(x₂, x₁) when V₁ = V₂ and V_int depends only on |x₁-x₂|
- T(k₁, k₂) = (k₁²+k₂²)/2 = T(k₂, k₁) always (both particles have the same mass)

The 2D FFT also preserves symmetry: FFT is a linear operator that commutes with
index permutation (swapping the row and column axes).

So every operation in the time step preserves the symmetry of ψ. If ψ starts
symmetric, it stays symmetric. If antisymmetric, it stays antisymmetric. No code
changes needed.

This is not an accident — it's a consequence of the quantum mechanical theorem
that the time evolution operator commutes with any symmetry of the Hamiltonian.
The split-operator method inherits this property because each factor individually
has the symmetry.

## Entanglement for Identical Particles

The purity calculation is unchanged, but its interpretation shifts.

For distinguishable particles, a product state has purity = 1.0. For identical
particles, the symmetrized state:

```
ψ = [φ_a(x₁)φ_b(x₂) + φ_a(x₂)φ_b(x₁)] / √2
```

is *not* a product state — it has two terms. The reduced density matrix is not
rank-1. So purity starts below 1.0 even without interaction.

This initial purity depends on how much φ_a and φ_b overlap. Well-separated
Gaussians (orthogonal) give purity ≈ 0.5 for both bosons and fermions. As
they're brought closer together, the bosonic purity increases (toward 1.0 for
identical states) while the fermionic purity becomes undefined (the state
vanishes).

The viewer computes purity on reset for identical particle scenarios so the
initial display is accurate, rather than defaulting to 1.0.

## The Viewer Scenarios

Three new scenarios, each targeting a specific phenomenon:

**Boson bunching**: Two bosons in a harmonic well (ω = 0.5), starting at x = ±3
with weak repulsive interaction. The symmetrized initial state has enhanced
probability on the x₁ = x₂ diagonal. Under evolution, the bosons oscillate but
tend to cluster. Compare with the distinguishable "Trapped + interaction"
scenario using the same parameters to see the exchange effect.

**Fermion antibunching**: Same setup, but antisymmetrized. The initial state has
a node (zero probability) along the entire x₁ = x₂ diagonal — the exchange
hole. Even with zero interaction, the particles can't be found at the same
position. Under evolution, the nodal line persists because the Hamiltonian
preserves antisymmetry.

**Pauli exclusion**: Two fermions initialized with identical parameters (both at
x = 0, same σ and k). The antisymmetrized state is identically zero. ‖ψ‖² = 0.
Nothing to evolve. This is the starkest demonstration: the state literally
doesn't exist.

The "Identical Particles" section appears as a third group in the scenario
sidebar, below "One Particle" and "Two Particles." The observables panel shows
the symmetry type ("Bosons (symmetric)" or "Fermions (antisymmetric)") for
identical particle scenarios.
