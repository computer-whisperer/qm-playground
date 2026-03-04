# Two-Particle Quantum Mechanics

In Chapter 1 we simulated a single particle on a line. That was enough to see tunneling,
quantization, and the uncertainty principle. But the truly strange part of quantum mechanics —
the part Einstein called "spooky action at a distance" — only appears when you have more
than one particle.

This chapter covers the physics behind Stage 2 of the simulation: two particles in 1D,
interacting via a potential.

## From One Particle to Two: Configuration Space

A single particle in 1D has a wavefunction ψ(x). Two particles in 1D have a wavefunction
ψ(x₁, x₂) — a complex number assigned to every *pair* of positions.

This is not two separate wavefunctions. It's a single function of two variables. The
probability of finding particle 1 near position a *and* particle 2 near position b is:

```
P(x₁ ≈ a, x₂ ≈ b) = |ψ(a, b)|² dx₁ dx₂
```

The domain of ψ is **configuration space** — the set of all possible configurations
of the system. For two particles in 1D, that's a 2D plane. For two particles in 3D,
it's 6D. For N particles in 3D, it's 3N-dimensional.

This is why quantum mechanics is hard: the wavefunction doesn't live in physical space.
It lives in configuration space, whose dimensionality grows linearly with particle count.
The computational cost grows *exponentially*.

### The Exponential Wall

Our single-particle simulation used n=1024 grid points. That's 1024 complex numbers.

Two particles on the same grid: n² = 1,048,576 complex numbers. We've dropped to n=256
per axis to keep this manageable (256² = 65,536), and already the simulation runs noticeably
slower.

Three particles would need n³. Four particles, n⁴. By about 10 particles, no computer on
Earth has enough memory. This exponential scaling is the fundamental obstacle in quantum
simulation — and the reason quantum computers are interesting.

## Product States and Entanglement

### Product States

If two particles are completely independent, their joint wavefunction factors:

```
ψ(x₁, x₂) = φ(x₁) · χ(x₂)
```

This is called a **product state**. Knowing about particle 1 tells you nothing about
particle 2. The heatmap of |ψ|² looks like a rectangle — the product of two independent
distributions.

In the simulation, every initial state is a product of two Gaussians. Before you press Run,
you're looking at a rectangular blob.

### Entanglement

Now turn on the interaction and let the particles scatter. The wavefunction develops
correlations:

```
ψ(x₁, x₂) ≠ φ(x₁) · χ(x₂)     for any choice of φ, χ
```

The state is **entangled**. Measuring particle 1's position now tells you something about
where particle 2 will be found. Not because of any signal between them — because the
wavefunction *itself* encodes a correlation that can't be described as "particle 1 does
its thing, particle 2 does its thing."

Watch the heatmap as the particles interact: the rectangular blob deforms into a shape
that can't be written as a product. The diagonal/anti-diagonal structure shows the
particles becoming correlated.

### Building Intuition with the Scenarios

**Free (no interaction):** Two Gaussians pass through each other. The heatmap stays
rectangular. The marginal distributions match what each particle would do alone.
Purity stays at 1.0. This is your control experiment.

**Scattering:** Two particles approach with opposite momenta. The soft-Coulomb repulsion
pushes them apart. After scattering, the wavefunction has a complicated shape: there's
some probability that both bounced back, some that they passed through each other, and
these possibilities are *correlated*. Purity drops — entanglement has developed.

**Trapped + interaction:** Both particles in a harmonic well with repulsion. They settle
into correlated bound states. Purity drops and stays low — persistent entanglement.

## Measuring Entanglement: Purity and Schmidt Number

How entangled is a state? We need a number.

### The Reduced Density Matrix

If you have the two-particle state ψ(x₁, x₂) but can only measure particle 1, what do
you know? You trace out particle 2 to get the **reduced density matrix**:

```
ρ₁(x, x') = ∫ ψ(x, x₂) ψ*(x', x₂) dx₂
```

This is a matrix (really an integral kernel) that captures everything about particle 1's
statistics. For a product state, ρ₁ = |φ⟩⟨φ| is a pure state (rank 1). For an entangled
state, ρ₁ is a *mixed* state — particle 1 isn't in any definite quantum state on its own.

### Purity

The **purity** Tr(ρ₁²) measures how "pure" the reduced state is:

```
Purity = Tr(ρ₁²) = ∫∫ |ρ₁(x, x')|² dx dx'
```

- Purity = 1.0 → product state (no entanglement)
- Purity = 1/n → maximally entangled

The simulation computes this by direct summation over the grid — O(n³) per evaluation,
so it updates every 10 frames rather than every frame.

### Schmidt Number

The **Schmidt number** K = 1/Purity gives the effective number of "terms" in the
entanglement:

- K = 1 → product state
- K = 2 → two significant terms (like a Bell pair)
- K > 10 → highly entangled

Watch K climb from 1.0 as the particles scatter. It quantifies what you're seeing in
the heatmap.

### The Schmidt Decomposition (Optional)

Any bipartite state can be written as:

```
ψ(x₁, x₂) = Σᵢ σᵢ uᵢ(x₁) vᵢ(x₂)
```

where σᵢ are the **Schmidt coefficients** (singular values of the wavefunction viewed as a
matrix), and uᵢ, vᵢ are orthonormal functions. A product state has one nonzero σ.
Entanglement means multiple nonzero σs.

This is literally the SVD of the wavefunction matrix. Purity = Σ σᵢ⁴. The Schmidt
number K = 1/Σσᵢ⁴ counts the effective number of nonzero σs.

## The Two-Particle Schrödinger Equation

The Hamiltonian for two interacting particles in 1D:

```
Ĥ = -½ ∂²/∂x₁² - ½ ∂²/∂x₂² + V₁(x₁) + V₂(x₂) + V_int(x₁, x₂)
```

Breaking this down:
- First two terms: kinetic energy of each particle (both with mass 1 in atomic units)
- V₁, V₂: external potentials on each particle (harmonic well, free, etc.)
- V_int: interaction between the particles

### Interaction Potentials

**Soft Coulomb:** V_int = g / √((x₁ - x₂)² + ε²)

The 1/r Coulomb potential, softened by ε to avoid the singularity when particles overlap.
g > 0 is repulsive, g < 0 is attractive. The ε parameter controls how sharp the repulsion
is at short range.

**Contact interaction:** V_int = g · δ_w(x₁ - x₂)

A narrow Gaussian approximation to a delta function. Models particles that only interact
when they're on top of each other — common in cold atom physics.

## The Split-Operator Method in 2D

The same split-operator approach from Chapter 1 extends to two particles. The key insight:
the 2D FFT separates into row-wise and column-wise 1D FFTs.

### Algorithm

```
1. Apply V half-step:  ψ *= exp(-i V(x₁,x₂) dt/2)     [position space]
2. Row-wise FFT:       transform x₂ → k₂
3. Column-wise FFT:    transform x₁ → k₁               [now in momentum space]
4. Apply T full step:  ψ *= exp(-i (k₁²+k₂²)/2 dt)     [momentum space]
5. Row-wise IFFT:      transform k₂ → x₂
6. Column-wise IFFT:   transform k₁ → x₁               [back to position space]
7. Normalize:          ψ /= n²                          [FFT convention]
8. Apply V half-step:  ψ *= exp(-i V(x₁,x₂) dt/2)
```

The potential step now involves the full 2D potential V(x₁, x₂) = V₁(x₁) + V₂(x₂) +
V_int(x₁, x₂), precomputed on the grid. The kinetic step uses T(k₁, k₂) = (k₁² + k₂²)/2,
also precomputed.

### Cost

Each step requires 4n 1D FFTs of length n (row-wise forward, column-wise forward,
row-wise inverse, column-wise inverse) plus O(n²) pointwise multiplications. Total cost
per step: O(n² log n). For n=256, this is about 500K operations — fast enough for real-time
interaction.

### Column FFTs

Row-wise FFTs are free because rows are contiguous in memory. Column FFTs require
extracting each column to a contiguous buffer, FFTing, and copying back. This is a
transpose-like operation and adds some overhead, but for n=256 it's negligible.

## What to Watch For in the Visualization

### The Heatmap

The main display shows |ψ(x₁, x₂)|² as a heatmap (inferno colormap). The horizontal
axis is x₁ (particle 1), the vertical axis is x₂ (particle 2).

- **Rectangle** → product state (unentangled)
- **Tilted/diagonal structure** → correlated → entangled
- **Anti-diagonal streak** → repulsive scattering (particles can't be in the same place)
- **Diagonal streak** → attractive correlation (particles prefer the same position)

### The Marginal Plots

The plot above the heatmap shows ρ₁(x₁) — the probability density of particle 1,
with particle 2 integrated out. The plot to the right shows ρ₂(x₂).

For non-interacting particles, these marginals match exactly what each particle would
do in isolation (verified by our test suite). For interacting particles, the marginals
spread and shift differently than independent particles would — the interaction affects
each particle's statistics.

### Purity and Schmidt Number

Watch these numbers in the observables panel:
- Start at Purity = 1.00, K = 1.00 (product state)
- During interaction: Purity drops, K rises
- After scattering: values stabilize at their new entangled values
- The amount of entanglement depends on interaction strength (g) and how much the
  particles overlap during evolution

## Why This Matters

Entanglement is not just a curiosity. It's the fundamental reason that:

1. **Quantum systems are hard to simulate classically.** An entangled state of N particles
   can't be described by N separate descriptions — you need the full n^N wavefunction.

2. **Quantum computing is possible.** Entanglement lets qubits share information in ways
   that classical bits can't, enabling exponential speedups for certain problems.

3. **Quantum chemistry is hard.** Electrons in molecules are entangled. This is why we
   can't just solve each electron independently — and why quantum chemistry is the
   leading near-term application of quantum computers.

4. **The measurement problem is sharp.** When you measure one particle of an entangled
   pair, the other particle's state instantly updates. This is nonlocal but can't transmit
   information (no FTL communication). Understanding why requires careful thinking about
   what the wavefunction means.

## How the Code Works

For a walkthrough of the 2D data structures, the row+column FFT trick, and how
entanglement is actually computed — see the
[implementation walkthrough](02-implementation-walkthrough.md).

## What Comes Next

After two-particle QM, the next steps toward the project's goal:

- **Identical particles**: Fermions (antisymmetric ψ) and bosons (symmetric ψ). This is
  where the Pauli exclusion principle comes from, and it's the foundation of all chemistry.
- **Second quantization**: Rewrite the theory in terms of creation/annihilation operators,
  which handles particle number changes and leads to quantum field theory.
- **Coarse-graining**: Train neural networks to approximate the many-body wavefunction
  or its effective dynamics at lower resolution — the bridge toward emergent behavior.
