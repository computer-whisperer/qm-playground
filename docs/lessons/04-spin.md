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

Spin operators are represented by the **Pauli matrices** — three 2×2 matrices,
one for each spatial axis. They act on the spinor (α, β) by matrix
multiplication:

```
σ_x = (0 1)    σ_y = (0 -i)    σ_z = (1  0)
      (1 0)          (i  0)          (0 -1)
```

Each matrix corresponds to measuring spin along one axis. For example, σ_z
applied to a spinor returns +1 times spin-up and -1 times spin-down — it
distinguishes the two states. The physical spin operator along axis a is
S_a = σ_a/2 (in our units where ℏ = 1), so the measurement outcomes are ±1/2.

What each matrix does to a spinor:
- **σ_z** leaves the spinor unchanged except for sign: σ_z(α, β) = (α, -β).
  It measures "how much up vs down."
- **σ_x** swaps the components: σ_x(α, β) = (β, α). It flips spin-up to
  spin-down and vice versa — the spin analog of a NOT gate.
- **σ_y** also swaps but with imaginary phases: σ_y(α, β) = (-iβ, iα). It
  combines a flip with a 90° rotation in the complex plane.

The notation σ_z|↑⟩ = +|↑⟩ means: applying the σ_z matrix to the spin-up
vector (1, 0) gives back (1, 0) with a factor of +1. Similarly,
σ_z|↓⟩ = −|↓⟩ means applying σ_z to (0, 1) gives (0, -1) = -1 × (0, 1).
In linear algebra terms, |↑⟩ and |↓⟩ are the eigenvectors of σ_z with
eigenvalues +1 and -1.

### Spin expectation values

For a spinor (α, β), the average spin measurement along each axis is:

```
⟨σ_z⟩ = |α|² - |β|²    — spin component in the z direction
⟨σ_x⟩ = 2 Re(α* β)     — spin component in the x direction
⟨σ_y⟩ = 2 Im(α* β)     — spin component in the y direction
```

Here α\* means the complex conjugate of α (flip the sign of the imaginary
part), and Re/Im extract the real/imaginary parts.

**⟨σ_z⟩** is the simplest: if the particle is all spin-up (α=1, β=0), then
⟨σ_z⟩ = 1. All spin-down: ⟨σ_z⟩ = -1. Equal superposition: ⟨σ_z⟩ = 0.

**⟨σ_x⟩ and ⟨σ_y⟩** depend on the relative phase between α and β. Two states
with the same |α| and |β| can have different x and y components depending on
whether α and β are in phase, out of phase, or 90° apart.

These three numbers (⟨σ_x⟩, ⟨σ_y⟩, ⟨σ_z⟩) form the **Bloch vector** — a
point on or inside a unit sphere that represents the spin state geometrically.
When the particle is in a definite spin state (a pure state), the Bloch vector
has length 1 and sits on the sphere's surface. When spin is entangled with
position (as in the Stern-Gerlach effect), the vector shrinks toward the
center — the spin is no longer in a definite state on its own.

## Magnetic Fields and the Zeeman Effect

A magnetic field interacts with spin because spin produces a tiny magnetic
moment (like a compass needle). The energy of a compass needle depends on its
orientation relative to the field — aligned is low energy, anti-aligned is high.

The full formula involves physical constants (the g-factor, the Bohr magneton),
but in our atomic units it simplifies to a clean result. For a magnetic field
pointing along the z axis:

```
H_mag = B_z σ_z / 2
```

What this means concretely: σ_z gives +1 for spin-up and -1 for spin-down,
so this term adds **+B_z/2 to the energy of spin-up** and **-B_z/2 to
spin-down**. An energy level that was the same for both spins now **splits**
into two levels separated by B_z. This is the **Zeeman effect** — first
observed in 1896 and one of the early clues that spin exists.

In the simulation, a longitudinal field (B_z) doesn't mix the spin components.
It just adds a different potential to each:

```
V_↑(x) = V(x) + B_z/2       spin-up sees a slightly higher potential
V_↓(x) = V(x) - B_z/2       spin-down sees a slightly lower potential
```

Each spin component evolves independently under its own effective potential,
like two separate particles in slightly different wells.

## Larmor Precession

A transverse magnetic field (along x or y) does something more interesting:
it **mixes the spin components**. While a z-field just shifts energies, an
x-field actively flips spin-up into spin-down and vice versa. A spin-up
particle in a B_x field will oscillate between up and down at the **Larmor
frequency** ω_L = B_x.

The time evolution is a rotation — the B_x field literally rotates the spin
state over time:

```
α(t) =  cos(B_x t/2) · α(0)  -  i·sin(B_x t/2) · β(0)
β(t) = -i·sin(B_x t/2) · α(0)  +  cos(B_x t/2) · β(0)
```

This is a 2×2 matrix multiplying the spinor at each moment. The cos and sin
control how much of the original up-component stays up vs rotates into down
(and vice versa). The factor of i (the imaginary unit) is needed to keep the
rotation unitary (probability-preserving).

Starting from pure spin-up (α=1, β=0), the z-component of spin oscillates:
- At t=0: ⟨σ_z⟩ = 1 (fully up)
- At t=π/(2B_x): ⟨σ_z⟩ = 0 (equal superposition of up and down)
- At t=π/B_x: ⟨σ_z⟩ = −1 (fully down)
- At t=2π/B_x: ⟨σ_z⟩ = 1 (back to up — full cycle)

This is **Larmor precession** — the spin precesses around the field direction,
analogous to a gyroscope precessing in gravity, but quantized. The precession
period T = 2π/B_x, so stronger fields make it spin faster.

In the simulation: the transverse field requires a rotation step that mixes
ψ_↑ and ψ_↓ at each grid point. This is the 2×2 matrix above, applied to the
pair (ψ_↑(x), ψ_↓(x)) at every position x.

## The Full Picture: Spatial Motion + Spin

A particle with spin has three things happening simultaneously:
1. **Spatial kinetic energy** — the usual -½ ∂²/∂x² from Chapter 1, which
   affects both spin components identically (momentum doesn't care about spin)
2. **Potential energy** V(x) — also the same for both components (for an
   electrostatic potential; magnetic fields make it spin-dependent)
3. **Magnetic coupling** — B_z splits the energies, B_x rotates the spin

Written as a formula:

```
H = [-½ ∂²/∂x² + V(x)] · I  +  B_z σ_z/2  +  B_x σ_x/2
     └── same for both ──┘      └── splits ──┘  └ rotates ┘
```

Here I is the 2×2 identity matrix — it means "do the same thing to both
components." The ⊗ (tensor product) symbol you'll see in textbooks is just a
formal way of saying "this spatial operator acts on both spin components
independently."

The split-operator method extends naturally. Each time step:
1. Apply V half-step (spin-dependent if B_z ≠ 0: different phase for ↑ and ↓)
2. Apply spin rotation half-step (if B_x ≠ 0: mix ↑ and ↓)
3. FFT both components to momentum space
4. Apply kinetic full step (same for both — momentum doesn't depend on spin)
5. IFFT both components back to position space
6. Apply spin rotation half-step
7. Apply V half-step

The cost is roughly 2× Chapter 1: two FFTs instead of one, plus some
pointwise 2×2 operations.

## Two Particles with Spin: Singlet and Triplet

This is where spin connects back to Chapter 3 and explains something deep.

Two spin-1/2 particles each have a spin that can be up or down, giving four
possible combinations: both up, first up + second down, first down + second
up, and both down. Written in the |particle₁ particle₂⟩ notation:

```
|↑↑⟩    |↑↓⟩    |↓↑⟩    |↓↓⟩
```

These four basis states combine into groups based on what happens when you
swap the two particles' labels:

### Spin triplet (total spin S=1, symmetric)

Three combinations that **stay the same** when you swap particles 1 and 2:

```
|↑↑⟩                         both up (obviously symmetric)
(|↑↓⟩ + |↓↑⟩)/√2             one up one down, symmetric combination
|↓↓⟩                         both down (obviously symmetric)
```

The middle state is symmetric because swapping ↑↓ and ↓↑ gives back the
same sum. The 1/√2 keeps the total probability equal to 1.

### Spin singlet (total spin S=0, antisymmetric)

One combination that **picks up a minus sign** when you swap:

```
(|↑↓⟩ - |↓↑⟩)/√2             one up one down, antisymmetric combination
```

Swapping gives (|↓↑⟩ - |↑↓⟩)/√2 = minus the original. This is the only
antisymmetric possibility.

### The connection to spatial symmetry

Here's the key insight from Chapter 3: fermions require their **total**
wavefunction to be antisymmetric under exchange. The total wavefunction is
the spatial part times the spin part. For the product to be antisymmetric,
one factor must be symmetric and the other antisymmetric:

- **Singlet** (antisymmetric spin) → spatial part must be **symmetric**
- **Triplet** (symmetric spin) → spatial part must be **antisymmetric**

This retroactively explains Chapter 3: **the "boson" mode (symmetric spatial
wavefunction) was really two fermions in a spin singlet, and the "fermion"
mode (antisymmetric spatial wavefunction) was two fermions in a spin triplet.**
The spatial symmetry we imposed wasn't arbitrary — it was forced by the spin.

Two electrons *can* share a spatial orbital — if their spins are opposite
(singlet). This is why each atomic orbital holds exactly two electrons:
one spin-up, one spin-down. The third electron must go to a different orbital.
That's the Pauli exclusion principle in its full spin-aware form.

## Exchange Energy

Because singlet and triplet have different spatial symmetries, they have
different energies when there's an interaction potential between the particles.

The symmetric spatial state (singlet spin) allows both particles to be at the
same position — so for a repulsive interaction, the particles feel more
repulsion. The antisymmetric spatial state (triplet spin) has a node where
the particles overlap — they naturally stay apart, feeling less repulsion.

The energy difference between singlet and triplet is the **exchange energy** J.
For repulsive interactions (like Coulomb), J > 0: the triplet has lower energy
because the particles avoid each other.

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
