# Spin

Chapters 1–3 treated particles as featureless points — characterized entirely by
position and momentum. But real particles carry an internal degree of freedom
called **spin**. It has no classical analog. It doesn't mean the particle is
literally spinning. It's an intrinsic angular momentum that exists even for
point particles, and it changes everything.

Spin is why atoms have magnetic moments, why the periodic table has its structure,
why MRI machines work, and why quantum computers use qubits. It's also the
bridge to the Dirac equation and relativistic quantum mechanics.

## What Is Spin?

Every fundamental particle has a fixed amount of intrinsic angular momentum.
Electrons, protons, and neutrons have spin-1/2. Photons have spin-1. The Higgs
boson has spin-0.

For spin-1/2 (the case we'll simulate), measuring the spin along any axis gives
one of exactly two results: +ℏ/2 ("spin up") or −ℏ/2 ("spin down"). There is no
middle ground. This two-valuedness is the defining feature of spin-1/2 and the
basis of the qubit.

In atomic units (ℏ = 1), the spin measurement outcomes are +1/2 and −1/2.

### Not a tiny spinning ball

If the electron were a classical spinning sphere, its surface would need to move
faster than light to produce the observed angular momentum given its size (which
is, as far as we know, zero). Spin isn't rotation. It's a fundamentally quantum
property with no classical counterpart — a new kind of thing that the math
demands and experiments confirm.

## The Spinor

A spin-1/2 particle's spin state lives in a 2D complex vector space. The state
is a **spinor** — a two-component complex vector:

```
|χ⟩ = α|↑⟩ + β|↓⟩ = (α)
                      (β)
```

where |α|² + |β|² = 1. Here:
- |α|² is the probability of measuring spin-up
- |β|² is the probability of measuring spin-down
- The relative phase between α and β determines the spin orientation in the
  transverse (x-y) plane

### The full wavefunction

For a particle with both spatial and spin degrees of freedom, the complete state
is:

```
|Ψ⟩ = ψ_↑(x)|↑⟩ + ψ_↓(x)|↓⟩
```

This is a pair of spatial wavefunctions — one for each spin component. The
probability of finding the particle at position x with spin up is |ψ_↑(x)|²,
and with spin down is |ψ_↓(x)|². The total probability density at x is
|ψ_↑(x)|² + |ψ_↓(x)|², and integrating this over all x gives 1.

In the simulation, this means storing **two** arrays of complex numbers instead
of one. The computational cost doubles, but the grid structure is unchanged.

## The Pauli Matrices

Spin operators are represented by the **Pauli matrices** — three 2×2 matrices
that form a basis for spin measurements:

```
σ_x = (0 1)    σ_y = (0 -i)    σ_z = (1  0)
      (1 0)          (i  0)          (0 -1)
```

The spin operator along axis a is S_a = σ_a/2 (in units of ℏ = 1).

Key properties:
- Each σ has eigenvalues ±1 (so S has eigenvalues ±1/2)
- σ_z|↑⟩ = +|↑⟩ and σ_z|↓⟩ = −|↓⟩ — our basis states are σ_z eigenstates
- σ_x flips the spin: σ_x|↑⟩ = |↓⟩ and σ_x|↓⟩ = |↑⟩
- σ_y also flips spin but with a phase: σ_y|↑⟩ = i|↓⟩
- They anticommute: σ_x σ_y = iσ_z (and cyclic permutations)

### Spin expectation values

For a spinor (α, β):

```
⟨σ_x⟩ = 2 Re(α* β)     — spin component in the x direction
⟨σ_y⟩ = 2 Im(α* β)     — spin component in the y direction
⟨σ_z⟩ = |α|² - |β|²    — spin component in the z direction
```

These three numbers define the **Bloch vector** — a point on the unit sphere
that represents the spin state geometrically. Pure states live on the surface
of the Bloch sphere; mixed states are inside it.

## Magnetic Fields and the Zeeman Effect

A magnetic field couples to spin. The interaction Hamiltonian is:

```
H_mag = -μ · B = (g_s μ_B / ℏ) S · B
```

In atomic units with g_s ≈ 2 and the Bohr magneton μ_B = 1/2, this simplifies.
For a field along z:

```
H_mag = B_z σ_z / 2
```

This adds +B_z/2 to the energy of spin-up and −B_z/2 to spin-down. An energy
level that was degenerate (same energy for both spins) **splits** into two levels
separated by B_z. This is the **Zeeman effect** — first observed in 1896 and
one of the early clues that spin exists.

In the simulation: a longitudinal field (B_z) doesn't mix the spin components.
It just adds a different potential to each:

```
V_↑(x) = V(x) + B_z/2
V_↓(x) = V(x) - B_z/2
```

Each spin component evolves independently under its own effective potential.

## Larmor Precession

A transverse magnetic field (along x or y) does something more interesting:
it **mixes the spin components**. A spin-up particle in a B_x field will
oscillate between up and down at the **Larmor frequency** ω_L = B_x.

The time evolution under H = B_x σ_x / 2 rotates the spinor:

```
(α(t))   ( cos(B_x t/2)    -i sin(B_x t/2) ) (α(0))
(β(t)) = ( -i sin(B_x t/2)  cos(B_x t/2)   ) (β(0))
```

Starting from pure spin-up (α=1, β=0):
- At t=0: ⟨σ_z⟩ = 1 (fully up)
- At t=π/(2B_x): ⟨σ_z⟩ = 0 (equal superposition)
- At t=π/B_x: ⟨σ_z⟩ = −1 (fully down)
- At t=2π/B_x: ⟨σ_z⟩ = 1 (back to up)

This is **Larmor precession** — the spin precesses around the field direction,
exactly analogous to a gyroscope precessing in gravity, but quantized.

In the simulation: a transverse field requires a spin rotation step that mixes
ψ_↑ and ψ_↓ at each grid point. This is a 2×2 unitary matrix applied pointwise.

## The Hamiltonian with Spin

For a single particle with spin in a potential V(x) and magnetic field B:

```
H = (-½ ∂²/∂x² + V(x)) ⊗ I₂ + B_z σ_z/2 + B_x σ_x/2
```

The first term is the spatial Hamiltonian acting identically on both spin
components. The second and third terms are the magnetic coupling, which acts
on the spin at each spatial point.

The split-operator method extends naturally:
1. Apply V half-step (potentially spin-dependent if B_z ≠ 0)
2. Apply spin rotation half-step (if B_x ≠ 0)
3. FFT both components
4. Apply T full step (same for both — momentum doesn't depend on spin)
5. IFFT both components
6. Apply spin rotation half-step
7. Apply V half-step

The kinetic energy doesn't distinguish spin states. Only the potential and
magnetic terms do.

## Two Particles with Spin: Singlet and Triplet

This is where spin connects back to Chapter 3 and explains something deep.

Two spin-1/2 particles have a 4D spin space: |↑↑⟩, |↑↓⟩, |↓↑⟩, |↓↓⟩. These
combine into states of definite total spin:

### Spin triplet (S=1, symmetric spin)

```
|↑↑⟩                              m_s = +1
(|↑↓⟩ + |↓↑⟩)/√2                  m_s =  0
|↓↓⟩                              m_s = -1
```

Three states, all symmetric under particle exchange.

### Spin singlet (S=0, antisymmetric spin)

```
(|↑↓⟩ - |↓↑⟩)/√2                  m_s =  0
```

One state, antisymmetric under particle exchange.

### The connection to spatial symmetry

For fermions, the **total** wavefunction (space × spin) must be antisymmetric.
Since antisymmetric = symmetric × antisymmetric, this means:

- **Singlet** (antisymmetric spin) → **symmetric spatial** wavefunction
- **Triplet** (symmetric spin) → **antisymmetric spatial** wavefunction

This is the stunning revelation: **Chapter 3's "boson" mode was actually two
fermions in a spin singlet, and "fermion" mode was two fermions in a spin
triplet.** The spatial symmetry we imposed wasn't arbitrary — it was the
consequence of the spin state.

Two electrons *can* share a spatial orbital — if their spins are opposite
(singlet). This is why each atomic orbital holds exactly two electrons:
one spin-up, one spin-down. The third electron must go to a different orbital.
That's the Pauli exclusion principle in its full spin-aware form.

## Exchange Energy

Because singlet and triplet have different spatial symmetries, they have
different energies when there's an interaction:

```
E_singlet = ⟨ψ_sym | H | ψ_sym⟩       (particles can overlap → higher V_int)
E_triplet = ⟨ψ_anti | H | ψ_anti⟩     (particles avoid each other → lower V_int)
```

The difference E_singlet − E_triplet is the **exchange energy** J. For repulsive
interactions (like Coulomb), J > 0: the triplet state has lower energy because
the particles stay apart, reducing their repulsion.

This energy difference is the origin of:

- **Hund's first rule**: In atoms, electrons prefer to have parallel spins (triplet)
  because the exchange energy is lower. This is why the ground state of carbon has
  two unpaired electrons, not a paired pair.

- **Ferromagnetism**: In iron, the exchange energy between neighboring atoms
  favors parallel spins. Below the Curie temperature, all the spins align,
  creating a permanent magnet. Magnetism is an exchange effect.

- **Chemical bonding**: Covalent bonds are spin singlets. The two electrons in
  a bond have opposite spins, which allows them to share the bonding orbital
  (symmetric spatial state), lowering the kinetic energy.

## What to Watch For in the Visualization

### Single-particle spinor scenarios

**Zeeman splitting**: Start with a Gaussian in an eigenstate (pure spin-up).
Apply a B_z field. The spin-up and spin-down components separate in energy
but don't mix. If started in a superposition, the two components oscillate
at different frequencies, creating a beating pattern.

**Larmor precession**: Start with spin-up in a B_x field. Watch ⟨σ_z⟩ oscillate
between +1 and −1 as the spin precesses. The spatial wavefunction barely
changes (if V(x) is spin-independent), but the spin oscillates rapidly.
Both ψ_↑ and ψ_↓ are visible on the same plot.

**Stern-Gerlach analog**: Start a moving wave packet in a gradient magnetic
field (B_z varying with x). The spin-up and spin-down components are deflected
in opposite directions, splitting into two beams — the quantum version of the
historic experiment that discovered spin.

### Two-particle spin scenarios

**Singlet vs triplet**: Run the same potential with symmetric (singlet) and
antisymmetric (triplet) spatial states. Observe different energies, different
spatial distributions. The exchange energy is visible as a measurable energy
difference between the two.

## How the Code Works

For a walkthrough of the 2-component spinor data structure, the modified
split-operator method, and the spin rotation step — see the
[implementation walkthrough](04-implementation-walkthrough.md).

## What Comes Next

With spin in hand, we're ready for the **Dirac equation** — the relativistic
wave equation that *predicts* spin rather than imposing it. The Dirac equation
naturally produces 2-component spinors (or rather, 4-component ones that include
antimatter), and our spinor infrastructure carries over directly.

Continue to [Chapter 5: The Dirac Equation](05-dirac-equation.md).
