# Identical Particles

In Chapters 1 and 2 we treated particles as distinguishable — "particle 1" and
"particle 2" were separate entities that happened to interact. But in quantum
mechanics, identical particles are *genuinely indistinguishable*. Two electrons
aren't just "very similar" — there is no experiment, even in principle, that can
tell which is which. This has profound consequences.

## Classical vs Quantum Identity

Classically, identical particles are still distinguishable: you can track their
trajectories. Label one red and one blue, follow them through time, and you always
know which is which. Their identity is bookkeeping.

Quantum mechanically, this fails. Particles don't have trajectories. Two electrons
in an atom don't take paths you can follow — their wavefunctions overlap, and
there's no fact of the matter about "which electron is where." The labels "1" and
"2" are artifacts of our notation, not features of reality.

This means the physics must be unchanged if we swap the labels. That constraint
alone — the **exchange symmetry** — gives us the Pauli exclusion principle, the
periodic table, the stability of matter, and the difference between metals and
insulators.

## Exchange Symmetry

If two particles are truly identical, swapping them cannot change any observable.
Since observables depend on |ψ|², we need:

```
|ψ(x₂, x₁)|² = |ψ(x₁, x₂)|²
```

This allows exactly two possibilities:

```
ψ(x₂, x₁) = +ψ(x₁, x₂)     (symmetric — bosons)
ψ(x₂, x₁) = -ψ(x₁, x₂)     (antisymmetric — fermions)
```

No other option works. Swapping could multiply ψ by some complex factor, but
since swapping twice must return to the original (swapping 1↔2 then 1↔2 again
is a no-op), that factor squared must equal 1. The only solutions are +1 and -1.

This is not a preference or approximation. It's a **superselection rule**: a
particle is either always a boson or always a fermion. There are no particles
that are "sometimes symmetric."

## Bosons and Fermions

### Bosons (symmetric)

Bosons have wavefunctions that are symmetric under exchange. Photons, gluons,
the Higgs boson, helium-4 atoms, and any composite particle with integer spin
are bosons.

Key property: **any number of bosons can occupy the same quantum state.** This
is what makes lasers work (many photons in the same mode) and enables
Bose-Einstein condensation (many atoms in the same ground state).

### Fermions (antisymmetric)

Fermions have wavefunctions that are antisymmetric under exchange. Electrons,
protons, neutrons, quarks, and any composite particle with half-integer spin
are fermions.

Key property: **no two fermions can occupy the same quantum state.** This is
the Pauli exclusion principle — and it falls directly out of the antisymmetry.
If two fermions are in the same state φ, the antisymmetrized wavefunction is:

```
ψ(x₁, x₂) = φ(x₁)φ(x₂) - φ(x₂)φ(x₁) = 0
```

Zero. Not "small" — identically zero. The state doesn't exist.

### The Spin-Statistics Theorem

Why are fermions antisymmetric and bosons symmetric? This is the content of the
**spin-statistics theorem**: particles with half-integer spin (1/2, 3/2, ...)
must be fermions, and particles with integer spin (0, 1, 2, ...) must be bosons.

The proof requires relativistic quantum field theory — it doesn't follow from
non-relativistic quantum mechanics alone. In our simulation we simply impose the
symmetry as a constraint. Nature enforces it automatically through the structure
of spacetime.

Our 1D simulation doesn't model spin. We treat the particles as spinless, and
impose bosonic or fermionic symmetry directly on the spatial wavefunction. In
full 3D with spin, the total wavefunction (space × spin) must have the right
symmetry — which allows fermions to share a spatial state if their spins are
opposite. In 1D without spin, the exclusion is absolute.

## Building Symmetric and Antisymmetric States

In Chapter 2, we initialized product states: ψ = φ_a(x₁) · φ_b(x₂). For
identical particles, we must symmetrize or antisymmetrize:

**Bosonic:**
```
ψ_B(x₁, x₂) = [φ_a(x₁)φ_b(x₂) + φ_a(x₂)φ_b(x₁)] / √2
```

**Fermionic:**
```
ψ_F(x₁, x₂) = [φ_a(x₁)φ_b(x₂) - φ_a(x₂)φ_b(x₁)] / √2
```

The 1/√2 maintains normalization (assuming φ_a and φ_b are orthogonal or at
least normalized). The fermionic case is a **Slater determinant**:

```
ψ_F(x₁, x₂) = (1/√2) det | φ_a(x₁)  φ_a(x₂) |
                            | φ_b(x₁)  φ_b(x₂) |
```

For N fermions, this generalizes to an N×N determinant (the full Slater
determinant). For two particles, it's just the formula above.

### What these look like on the heatmap

A distinguishable product state ψ = φ_a(x₁)φ_b(x₂) appears as a rectangular
blob centered at (a, b) in configuration space.

The **bosonic** state adds the (b, a) blob: two rectangular blobs, one at
(a, b) and one at (b, a), superposed constructively. Along the x₁ = x₂
diagonal, the two blobs overlap and *add* — bosons are enhanced near the diagonal.

The **fermionic** state subtracts the (b, a) blob: the same two blobs, but
with opposite signs. Along x₁ = x₂, they cancel exactly — the wavefunction
is zero everywhere on the diagonal. Fermions have a **node** along x₁ = x₂.
They cannot be found at the same position.

## Exchange Effects: Bunching and Antibunching

Even without any interaction potential, identical particles behave differently
from distinguishable ones. These are **exchange effects** — purely quantum
phenomena with no classical counterpart.

### Boson Bunching

Two bosons in the same potential tend to cluster together. Their symmetrized
wavefunction has enhanced probability along the x₁ = x₂ diagonal — they
"prefer" to be near each other. This isn't caused by an attractive force.
There is no force. It's a consequence of the symmetry constraint on the
wavefunction.

In the simulation: put two bosons in a harmonic well with slightly different
initial states and no interaction. Watch the marginal distributions — they'll
be pushed closer together compared to distinguishable particles in the same
setup.

### Fermion Antibunching (the Exchange Hole)

Two fermions in the same potential avoid each other. The antisymmetric
wavefunction has zero probability along x₁ = x₂ — they *cannot* be at the
same position. Around the diagonal, probability is suppressed. This creates
an **exchange hole**: a region of depleted probability around each fermion
where the other fermion cannot be found.

The exchange hole exists even with zero interaction potential. It's not
Coulomb repulsion (though that adds on top). It's a geometric consequence
of antisymmetry.

In the simulation: put two fermions in a harmonic well with no interaction.
The heatmap shows a clear node along the diagonal. The marginal distributions
are pushed apart compared to distinguishable particles. This is the same
effect that fills electron shells in atoms — each electron excludes others
from its state, forcing them into higher energy levels.

### Comparing All Three

The same physical setup — same well, same initial states, same interaction —
produces three different outcomes depending on particle statistics:

- **Distinguishable**: simple product, no correlation
- **Bosons**: clustered, enhanced diagonal, bunching
- **Fermions**: separated, nodal diagonal, antibunching

The difference is entirely in the initial symmetry. The Hamiltonian is the same.
The integrator is the same. Only the initial condition changes — and the
consequences are dramatic.

## Why the Integrator Doesn't Change

A critical implementation detail: the split-operator method from Chapter 2
works without modification. If ψ starts symmetric, it stays symmetric under
time evolution, because the Hamiltonian is symmetric:

```
H(x₁, x₂) = H(x₂, x₁)
```

when both particles have the same mass and the same external potential (V₁ = V₂),
and the interaction depends only on |x₁ - x₂|. All of our potentials satisfy
this. So symmetry is preserved automatically — we don't need to modify the
time-stepping algorithm.

In practice, numerical roundoff can break exact symmetry over many steps. A
periodic re-symmetrization step (projecting ψ back onto the symmetric or
antisymmetric subspace) keeps things clean. This is a cheap O(n²) operation
and doesn't change the physics.

## Entanglement in Identical Particles

Entanglement for identical particles is subtle. A symmetrized product state:

```
ψ_B = [φ_a(x₁)φ_b(x₂) + φ_a(x₂)φ_b(x₁)] / √2
```

looks entangled — it's not a product state. But is this "real" entanglement, or
just a consequence of the symmetry bookkeeping?

This is an active debate in quantum foundations. For our purposes: the purity
of the reduced density matrix still works as a measure, but its baseline changes.
A symmetrized product state has purity < 1 even without any interaction, because
the symmetrization itself creates correlations. The interesting physics is in
how purity *changes* from this baseline when interactions are turned on.

The simulation tracks this naturally — you can compare purity for non-interacting
identical particles (the symmetry baseline) against interacting identical particles
(symmetry plus dynamical correlations).

## The Pauli Exclusion Principle

The exclusion principle is not an additional postulate. It's a theorem:

**If ψ is antisymmetric and φ_a = φ_b, then ψ = 0.**

Proof:
```
ψ(x₁, x₂) = φ(x₁)φ(x₂) - φ(x₂)φ(x₁) = 0
```

That's it. The antisymmetry does all the work. Two fermions cannot occupy the
same single-particle state because the resulting wavefunction is identically zero —
not "very small," not "suppressed," but exactly zero everywhere.

In the simulation, you can try this: initialize two fermions with the same
Gaussian parameters (same x0, sigma, k0). The antisymmetrization produces
ψ = 0. The norm is zero. There is literally nothing there.

This is why the periodic table has shells. Each electron state (defined by
quantum numbers n, l, m, s) can hold at most one electron. Once a state is
occupied, additional electrons must go into higher-energy states. This filling
pattern gives atoms their chemical properties.

## What to Watch For in the Visualization

### Boson scenarios

- The heatmap is symmetric about the x₁ = x₂ diagonal (by construction)
- Probability is *enhanced* on the diagonal — bosons cluster
- Marginal distributions shift toward each other compared to distinguishable particles
- With interaction, entanglement develops on top of the symmetry baseline

### Fermion scenarios

- The heatmap is antisymmetric: probability is zero on the x₁ = x₂ diagonal
- A clear **nodal line** runs along the diagonal — the exchange hole
- Marginal distributions shift apart compared to distinguishable particles
- The "same state" initialization produces ψ = 0 — exclusion principle in action
- With repulsive interaction, the exchange hole and Coulomb repulsion combine

### Comparing statistics

Run the same setup with distinguishable, bosonic, and fermionic symmetry. The
dramatic differences — from the same Hamiltonian, the same integrator, the same
grid — demonstrate that exchange symmetry is not a small correction. It's a
fundamental structural feature of quantum mechanics that reshapes everything.

## How the Code Works

For a walkthrough of the (anti)symmetrized initialization, the re-symmetrization
step, and how the existing integrator handles identical particles automatically —
see the [implementation walkthrough](03-implementation-walkthrough.md).

## Why This Matters for the Project

Exchange symmetry is the gateway to chemistry. Without it, all atoms would have
the same ground state (all electrons in the 1s orbital), and there would be no
periodic table, no chemical bonds, no molecules. The exclusion principle forces
diversity.

For the project's long-term goal of learned coarse-graining:
- A neural net that approximates many-electron systems *must* respect exchange
  symmetry. This is a hard constraint, not a soft preference.
- Architectures that build in symmetry (equivariant neural networks, antisymmetric
  layers) learn faster and generalize better than those that have to discover it.
- The jump from two particles to many particles is where second quantization
  becomes necessary — the Slater determinant formalism doesn't scale, and field
  operators provide a cleaner framework.

## What Comes Next

Continue to [Chapter 4: Spin](04-spin.md) — where particles gain an internal
degree of freedom, magnetic fields split energy levels, and the connection between
exchange symmetry and the periodic table becomes complete.

With spin and identical particles in hand, the further steps:

- **Second quantization**: Rewrite the theory in terms of creation and
  annihilation operators. This handles antisymmetry automatically and extends
  naturally to systems where particle number changes (photon emission, pair
  creation).
- **Many-body approximations**: Hartree-Fock, density functional theory — the
  workhorses of computational chemistry, all built on exchange symmetry.
- **Coarse-graining**: Can a neural network learn the effective physics of
  many identical particles from few-particle training data?
