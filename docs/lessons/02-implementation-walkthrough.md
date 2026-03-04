# How the Two-Particle Simulation Works

A code-level walkthrough of the Stage 2 simulation. Read alongside the
[physics primer](02-two-particle-qm.md) for the theory — this document covers
what changed when we went from one particle to two, and why.

## The Data: A 2D Wavefunction

One particle has ψ(x) — an array of n complex numbers. Two particles have
ψ(x₁, x₂) — a complex number for every *pair* of positions. That's an n×n grid:

```rust
struct Wavefunction2D {
    psi: Vec<Complex64>,  // n² complex amplitudes, row-major
    n: usize,             // grid points per axis (typically 256)
    x_min: f64,
    x_max: f64,
    dx: f64,              // same spacing on both axes
}
```

The storage layout is row-major: `psi[i1 * n + i2]` is ψ at position (x₁ = x(i1),
x₂ = x(i2)). Particle 1's coordinate varies along the "row" axis, particle 2's
along the "column" axis.

Both particles share the same spatial grid. This is a simplification — in principle
they could have different domains — but it keeps the FFT machinery uniform.

### The cost of a second particle

The grid dropped from n = 1024 (single particle) to n = 256 (two particles).
That's 256² = 65,536 complex numbers instead of 1,024. In memory it's manageable.
But every operation that was O(n) is now O(n²), and FFT operations that were
O(n log n) become O(n² log n). Three particles would be n³ — this exponential
scaling is the entire reason quantum mechanics is hard to simulate classically.

## Initialization: Product Gaussians

The initial state is always a product of two independent Gaussians:

```rust
fn set_product_gaussian(&mut self,
    x0_1: f64, sigma_1: f64, k0_1: f64,  // particle 1
    x0_2: f64, sigma_2: f64, k0_2: f64,  // particle 2
) {
    for i1 in 0..n {
        let phi1 = gaussian(x(i1), x0_1, sigma_1, k0_1);

        for i2 in 0..n {
            let phi2 = gaussian(x(i2), x0_2, sigma_2, k0_2);
            psi[i1 * n + i2] = phi1 * phi2;
        }
    }
    self.normalize();
}
```

This is a **product state**: ψ(x₁, x₂) = φ₁(x₁) · φ₂(x₂). The two particles
are independent — knowing about one tells you nothing about the other.

On the heatmap this looks like a rectangle: the 2D probability density is the
outer product of two 1D bell curves. The rectangular shape is the visual signature
of "no entanglement."

After the interaction turns on, the state generally stops being factorizable.
That's when entanglement appears, and the rectangle deforms.

## The 2D Potential

The potential has three parts:

```rust
struct Potential2D {
    v1: Potential,          // V₁(x₁) — external force on particle 1
    v2: Potential,          // V₂(x₂) — external force on particle 2
    interaction: Interaction,  // V_int(x₁, x₂) — between the particles
}

fn value_at(&self, x1: f64, x2: f64) -> f64 {
    self.v1.value_at(x1) + self.v2.value_at(x2) + self.interaction.value_at(x1, x2)
}
```

v1 and v2 reuse the single-particle `Potential` type from Chapter 1. Each can be
free, harmonic, a well, etc. The interaction depends on *both* positions —
specifically on how far apart the particles are.

### Interaction types

```rust
enum Interaction {
    None,                           // independent particles
    SoftCoulomb { g, epsilon },     // g / sqrt((x₁-x₂)² + ε²)
    Contact { g, width },           // g · gauss(x₁-x₂, width)
}
```

**Soft Coulomb** is a 1/r potential with a softening parameter ε that prevents
the singularity when particles overlap. The real Coulomb potential V = 1/|r| blows
up at r = 0; the soft version V = 1/√(r² + ε²) caps it at 1/ε. The `g` parameter
controls sign and strength: g > 0 is repulsive, g < 0 is attractive.

**Contact** is a narrow Gaussian approximation to a delta function. It's zero
unless the particles are practically on top of each other (within a few `width`
of separation). This models very short-range forces — common in cold atom physics
where atoms only interact when their wavefunctions overlap.

Both interactions depend only on the *difference* x₁ - x₂, which is physically
natural: the force between two particles depends on how far apart they are, not
on where they are in absolute space.

## Time Evolution: 2D Split-Operator

The algorithm is the same as Chapter 1 — Strang splitting with FFT — but the
FFTs are now two-dimensional.

### The integrator

```rust
struct SplitOperator2D {
    potential_half: Vec<Complex64>,  // exp(-i V(x₁,x₂) dt/2), n² entries
    kinetic_full:  Vec<Complex64>,   // exp(-i (k₁²+k₂²)/2 dt), n² entries
    fft:  Arc<dyn Fft>,              // 1D forward FFT plan (length n)
    ifft: Arc<dyn Fft>,              // 1D inverse FFT plan (length n)
    scratch: Vec<Complex64>,         // FFT workspace
    col_buf: Vec<Complex64>,         // contiguous buffer for column FFTs
    n: usize,
}
```

Notice the FFT plans are still **1D** (length n). We don't use a 2D FFT library.
Instead, a 2D FFT is decomposed into n row-wise 1D FFTs followed by n column-wise
1D FFTs. This is mathematically equivalent to a full 2D FFT.

### Precomputing the phase factors

At construction time, we build n² potential phase factors and n² kinetic phase
factors:

```rust
// Potential half-step: one entry per configuration-space point
for i1 in 0..n {
    for i2 in 0..n {
        let v = potential.value_at(xs[i1], xs[i2]);
        potential_half[i1*n + i2] = exp(-i * v * dt/2);
    }
}

// Kinetic full step: one entry per momentum-space point
for i1 in 0..n {
    let k1 = freq(i1);  // momentum for frequency index i1
    for i2 in 0..n {
        let k2 = freq(i2);
        let T = (k1*k1 + k2*k2) / 2.0;  // kinetic energy of this mode
        kinetic_full[i1*n + i2] = exp(-i * T * dt);
    }
}
```

The `freq()` function maps FFT indices to physical momenta, using the same
wrapping convention as the 1D case: indices 0..n/2 are positive frequencies,
n/2+1..n-1 are negative (wrapped) frequencies.

The kinetic energy T(k₁, k₂) = (k₁² + k₂²)/2 is the sum of each particle's
kinetic energy — they move independently in momentum space (it's only the
interaction potential that couples them).

### One time step

```rust
fn step(&mut self, wf: &mut Wavefunction2D) {
    let psi = &mut wf.psi;

    // 1. Potential half-step (position space)
    for i in 0..n*n {
        psi[i] *= potential_half[i];
    }

    // 2. Forward 2D FFT
    self.fft_2d(psi, forward=true);

    // 3. Kinetic full step (momentum space)
    for i in 0..n*n {
        psi[i] *= kinetic_full[i];
    }

    // 4. Inverse 2D FFT
    self.fft_2d(psi, forward=false);

    // 5. Normalize: divide by n² (rustfft doesn't normalize)
    for c in psi {
        *c /= n²;
    }

    // 6. Potential half-step again
    for i in 0..n*n {
        psi[i] *= potential_half[i];
    }
}
```

Structurally identical to the 1D version — only the FFT is different and we
divide by n² instead of n.

### How the 2D FFT works

The 2D FFT is the most interesting piece. It decomposes into row-wise and
column-wise passes:

```rust
fn fft_2d(&mut self, psi: &mut [Complex64], forward: bool) {
    let plan = if forward { &self.fft } else { &self.ifft };

    // Pass 1: FFT each row (rows are contiguous in memory)
    for i1 in 0..n {
        let row = &mut psi[i1*n .. (i1+1)*n];
        plan.process(row);
    }

    // Pass 2: FFT each column (columns are NOT contiguous — need a buffer)
    for i2 in 0..n {
        // Extract column into contiguous buffer
        for i1 in 0..n {
            col_buf[i1] = psi[i1*n + i2];
        }

        plan.process(&mut col_buf);

        // Scatter back to column
        for i1 in 0..n {
            psi[i1*n + i2] = col_buf[i1];
        }
    }
}
```

**Rows are free**: each row occupies a contiguous slice of memory, so we can hand
it directly to the FFT.

**Columns require a copy**: column elements are strided — `psi[0*n + i2]`,
`psi[1*n + i2]`, `psi[2*n + i2]`, etc. are n elements apart. FFT libraries
want contiguous data. So we extract each column into `col_buf`, FFT it there,
then scatter the results back.

This gather-FFT-scatter pattern is equivalent to a matrix transpose + row FFTs +
transpose back. Our approach avoids the full transpose. For n = 256, each column
copy is 256 complex numbers (4 KB) — well within cache and negligible compared
to the FFT cost.

**Total per step:** 4n FFTs of length n (n rows forward, n columns forward, n rows
inverse, n columns inverse) plus 3 pointwise multiplies of n² elements. Cost:
O(n² log n). At n = 256, that's roughly 500K multiply-add operations — comfortably
real-time.

## Computing Observables

### Marginal distributions

The heatmap shows the full 2D probability density. But we also want to know about
each particle individually. The **marginal** for particle 1 is what you get by
integrating out particle 2:

```rust
fn marginal_1(&self) -> Vec<f64> {
    let mut rho = vec![0.0; n];
    for i1 in 0..n {
        let mut sum = 0.0;
        for i2 in 0..n {
            sum += |psi[i1*n + i2]|²;
        }
        rho[i1] = sum * dx;   // integrate over x₂
    }
    rho
}
```

For each position x₁, we sum |ψ|² over all x₂ values. The result is the
probability density of particle 1 alone, regardless of where particle 2 is.
Marginal 2 is the same with the loops swapped.

For a product state, marginal_1 exactly matches what particle 1 would do
in an independent 1D simulation — the test suite verifies this by running
both and comparing element-by-element. When the interaction turns on, the
marginals diverge from independent predictions; the interaction affects each
particle's statistics through entanglement.

### Expected positions ⟨x₁⟩ and ⟨x₂⟩

Same idea as 1D, but summing over the 2D grid:

```rust
fn expected_x1(&self) -> f64 {
    let mut sum = 0.0;
    for i1 in 0..n {
        for i2 in 0..n {
            sum += x(i1) * |psi[i1*n + i2]|²;
        }
    }
    sum * dx²  // 2D integral: sum * dx₁ * dx₂
}
```

Note the `dx²` factor — the 2D discrete integral needs both differentials.

### Energy ⟨H⟩ = ⟨T⟩ + ⟨V⟩

The expected potential energy is the straightforward 2D weighted average:

```rust
// ⟨V⟩ = ∫∫ V(x₁,x₂) |ψ|² dx₁ dx₂
for i1 in 0..n {
    for i2 in 0..n {
        v_sum += V(xs[i1], xs[i2]) * |psi[i1*n + i2]|²;
    }
}
expected_v = v_sum * dx²;
```

Kinetic energy uses second derivatives in *both* directions — one for each
particle's momentum:

```rust
// T = -½(∂²/∂x₁² + ∂²/∂x₂²)
for i1 in 0..n {
    for i2 in 0..n {
        // second derivative in x₁ (central difference, periodic)
        let d2_x1 = (psi[(i1+1)*n + i2] - 2*psi[i1*n + i2] + psi[(i1-1)*n + i2]) / dx²;

        // second derivative in x₂ (central difference, periodic)
        let d2_x2 = (psi[i1*n + (i2+1)] - 2*psi[i1*n + i2] + psi[i1*n + (i2-1)]) / dx²;

        t_sum += conj(psi[i1*n + i2]) * (d2_x1 + d2_x2);
    }
}
expected_t = (-0.5 * t_sum * dx²).re;
```

The total energy ⟨T⟩ + ⟨V⟩ should be conserved. The test suite verifies drift
stays below 0.05 Hartree over 500 steps with a harmonic potential plus
soft-Coulomb interaction.

## Measuring Entanglement

This is the genuinely new physics in Stage 2. Everything else was "the same as 1D
but squared." Entanglement doesn't exist for one particle.

### The reduced density matrix

If you can only measure particle 1, the relevant object is the reduced density
matrix:

```
ρ₁(x, x') = ∫ ψ(x, x₂) ψ*(x', x₂) dx₂
```

On the grid, this becomes an n×n matrix where element (i, k) sums over the
"traced out" particle 2:

```rust
// ρ₁(i, k) = Σⱼ ψ(i,j) · conj(ψ(k,j)) · dx
let mut dot = Complex64::zero();
for j in 0..n {
    dot += psi[i*n + j] * conj(psi[k*n + j]);
}
// dot * dx gives one entry of the reduced density matrix
```

This is an inner product: row i of the wavefunction matrix dotted with row k.
If ψ is a product state ψ = φ₁ ⊗ φ₂, then ρ₁(i, k) = φ₁(i) · conj(φ₁(k)) ·
‖φ₂‖² — a rank-1 matrix. Entanglement makes it higher rank.

### Purity: how entangled are they?

The purity Tr(ρ₁²) tells you how "mixed" the reduced state is:

```rust
fn purity(wf: &Wavefunction2D) -> f64 {
    let mut sum = 0.0;
    for i in 0..n {
        for k in i..n {   // exploit Hermitian symmetry
            // compute |ρ₁(i,k)|² = |Σⱼ ψ(i,j) · conj(ψ(k,j))|²
            let mut dot = Complex64::zero();
            for j in 0..n {
                dot += psi[i*n + j] * conj(psi[k*n + j]);
            }
            let ns = dot.norm_sqr();

            if i == k {
                sum += ns;           // diagonal: count once
            } else {
                sum += 2.0 * ns;     // off-diagonal: (i,k) and (k,i) are conjugates
            }
        }
    }
    sum * dx⁴   // four dx factors: two from ρ, two from Tr(ρ²)
}
```

- Purity = 1.0 → product state, no entanglement
- Purity < 1.0 → entangled
- Purity = 1/n → maximally entangled (everything correlated with everything)

**The Hermitian symmetry optimization**: ρ₁ is Hermitian, so |ρ₁(i,k)|² = |ρ₁(k,i)|².
We only iterate the upper triangle (k ≥ i) and double-count the off-diagonals.
This cuts the work roughly in half.

**Why O(n³)?** Three nested loops: i, k, and j. Each (i, k) pair requires an
inner product over j. No way around it without computing the full SVD, which is
also O(n³). This is why we compute purity every ~10 frames instead of every frame.

**Why dx⁴?** The reduced density matrix has dimensions of 1/length (one dx from
the trace integral, one from the definition). Squaring it and tracing again
gives four dx factors. If you forget any of them, purity won't come out to 1.0
for a product state — the test suite catches this.

### Schmidt number

```rust
fn schmidt_number(wf: &Wavefunction2D) -> f64 {
    1.0 / purity(wf)
}
```

K = 1 for product states, K > 1 for entangled states. Physically, K tells you
"how many independent terms" are needed to describe the correlations. A Bell pair
has K = 2. Strongly interacting particles can reach K > 10.

### Connection to SVD

The physics primer mentions the Schmidt decomposition. Computationally, if you
reshape ψ from an n²-vector into an n×n matrix M where M(i1, i2) = ψ(i1, i2),
then the Schmidt decomposition is literally the singular value decomposition of
M. The singular values σᵢ are the Schmidt coefficients. Purity equals Σ σᵢ⁴.

We don't actually compute the SVD — it would be O(n³) and give us more information
than we need. The purity sum is equivalent and avoids constructing the full
decomposition.

## What Each Scenario Demonstrates

**Free (no interaction)** — the control experiment. Two Gaussians fly through
each other. Because there's no interaction, the state stays a product: ψ = φ₁ ⊗ φ₂
for all time. The heatmap stays rectangular. Purity stays at 1.0. The marginals
match independent 1D simulations exactly. This scenario validates that the 2D
code doesn't introduce spurious entanglement.

**Scattering (soft-Coulomb repulsion)** — two packets approach with opposite
momenta. When they overlap, the repulsive potential pushes them apart. But quantum
mechanics allows superposition: there's some amplitude for "both bounced back"
and some for "both passed through." These two outcomes are correlated — if particle
1 bounced, particle 2 probably did too. That correlation can't be factored into
independent descriptions, so entanglement develops. Watch the heatmap deform from
a rectangle into an anti-diagonal streak (particles avoiding the x₁ = x₂ diagonal
where the repulsion is strongest).

**Trapped + interaction** — both particles in harmonic wells with soft-Coulomb
repulsion. They settle into correlated bound states. The entanglement is persistent:
purity drops and stays low, unlike the scattering case where the particles
eventually fly apart and the entanglement "freezes" at its post-scattering value.

## Differences from the 1D Code

A summary of what changed structurally:

| Aspect | 1D (Chapter 1) | 2D (Chapter 2) |
|--------|-----------------|-----------------|
| Grid | n points | n² points |
| Default n | 1024 | 256 |
| Storage | Vec of n Complex64 | Vec of n² Complex64 |
| Potential | V(x) | V₁(x₁) + V₂(x₂) + V_int(x₁,x₂) |
| FFT | 1 forward + 1 inverse | 2n forward + 2n inverse (row + col) |
| Normalization divisor | n | n² |
| New observable | — | Purity, Schmidt number |
| Observables cost | O(n) each | O(n²) for position/energy, O(n³) for purity |
| Time per step | ~μs | ~ms |

The code reuses the single-particle `Potential` type for the external forces on
each particle. The interaction is the genuinely new piece — it's the only term
that depends on both coordinates simultaneously, and it's what makes entanglement
possible.
