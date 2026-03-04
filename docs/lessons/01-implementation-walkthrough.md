# How the Single-Particle Simulation Works

A code-level walkthrough of the Stage 1 simulation. Read alongside the
[physics primer](01-single-particle-qm.md) for the theory — this document covers
what the code actually does and why it's structured the way it is.

## The Data: A Wavefunction on a Grid

Quantum mechanics says ψ(x) is a continuous function. Computers need discrete arrays.
We represent ψ on a uniform grid of n points between x_min and x_max:

```rust
struct Wavefunction {
    psi: Vec<Complex64>,  // n complex amplitudes
    n: usize,             // grid points (typically 1024)
    x_min: f64,           // left edge of the domain
    x_max: f64,           // right edge
    dx: f64,              // spacing: (x_max - x_min) / n
}
```

Each entry `psi[i]` holds the complex amplitude at position `x(i)`:

```rust
fn x(&self, i: usize) -> f64 {
    self.x_min + (i as f64 + 0.5) * self.dx
}
```

The `+ 0.5` is a cell-centered convention: grid points sit at the *midpoints*
of their cells rather than at the boundaries. This avoids putting a grid point
exactly at x_min or x_max, which matters when boundary conditions are in play.

### Why Complex?

ψ is complex-valued. The *magnitude* |ψ(x)|² gives probability density. The
*phase* — the angle of ψ in the complex plane — encodes momentum. A wavefunction
`exp(i k x)` oscillates in phase with wavelength 2π/k, and that k *is* the
momentum (in atomic units where ℏ = 1).

So `psi[i]` isn't just a number telling you "how much probability is here." It's
a number telling you "how much probability is here, and which direction is it going."

## Initialization: The Gaussian Wave Packet

Every scenario starts with a Gaussian wave packet — a bell curve in probability,
with a phase twist that gives it momentum:

```rust
fn set_gaussian(&mut self, x0: f64, sigma: f64, k0: f64) {
    let norm = (2π σ²)^(-1/4);  // analytic normalization
    for i in 0..n {
        let x = self.x(i);
        let envelope = exp(-(x - x0)² / (4σ²));     // bell curve centered at x0
        let phase    = exp(i k0 x);                   // momentum kick
        psi[i] = norm * envelope * phase;
    }
    self.normalize();  // fix up for discrete grid
}
```

Three parameters control the initial state:
- **x0**: where the packet is centered
- **sigma**: how wide it is (position uncertainty Δx ≈ σ)
- **k0**: its average momentum (p = k0 in atomic units)

The `normalize()` call at the end ensures ∫|ψ|² dx = 1 on the discrete grid,
compensating for the fact that the analytic normalization constant doesn't
perfectly match the discrete sum.

### The uncertainty principle is baked in

A Gaussian is the *minimum uncertainty* state: Δx · Δp = ℏ/2. If you make sigma
small (tight position), the momentum-space representation is wide (uncertain
momentum), and vice versa. You can see this directly in the sim: a narrow Gaussian
spreads faster because its high-momentum components fly apart.

## The Potential: What the Particle Feels

The potential V(x) defines the forces on the particle. It's built from composable terms:

```rust
struct Potential {
    terms: Vec<PotentialTerm>,
}

enum PotentialTerm {
    Free,                                    // V = 0 everywhere
    Well { left, right, wall_height },       // V = 0 inside, wall_height outside
    Barrier { left, right, height },         // V = height inside, 0 outside
    Harmonic { x0, omega },                  // V = ½ω²(x - x0)²
    DoubleWell { x1, x2, omega, barrier },   // min of two parabolas, capped
    Custom(closure),                         // arbitrary V(x)
}
```

`value_at(x)` sums all terms — so you can combine a harmonic well with a barrier,
or any other mix. The simulation only ever calls `value_at(x)` at specific grid
points, so the potential doesn't need to be differentiable or even continuous. It's
just a lookup: given a position, what's the potential energy?

The Well variant deserves a note: quantum mechanics doesn't have "infinite walls"
in the way a textbook draws them. We approximate an infinite square well with a
very large `wall_height` (typically 1000 Hartree). The wavefunction decays
exponentially into the wall and is effectively zero there — close enough to the
ideal case for visualization.

## Time Evolution: The Split-Operator Method

This is the heart of the simulation. Given ψ at time t, we need ψ at time t + dt.

The Schrödinger equation says:

```
ψ(t + dt) = exp(-i Ĥ dt) ψ(t)
```

where Ĥ = T̂ + V̂ (kinetic plus potential energy). The problem: T̂ and V̂ don't commute,
so exp(-i(T̂+V̂)dt) ≠ exp(-iT̂ dt) · exp(-iV̂ dt). But we can *approximately* split them
using the Strang splitting:

```
exp(-i(T̂+V̂)dt) ≈ exp(-iV̂ dt/2) · exp(-iT̂ dt) · exp(-iV̂ dt/2) + O(dt³)
```

The key insight: each of these three factors is easy to apply, just in *different*
representations:

- **V̂ is diagonal in position space.** Multiplying by exp(-iV(x) dt/2) is a pointwise
  multiply — just scale each grid point by a complex phase.

- **T̂ is diagonal in momentum space.** In Fourier space, the kinetic energy of each
  frequency component is simply k²/2. So exp(-iT̂ dt) is also a pointwise multiply,
  but on the Fourier coefficients.

The FFT lets us switch between these representations in O(n log n).

### The Integrator

```rust
struct SplitOperator {
    potential_half: Vec<Complex64>,  // exp(-i V(xⱼ) dt/2) for each grid point
    kinetic_full:   Vec<Complex64>,  // exp(-i k²/2 dt) for each frequency
    fft:  Arc<dyn Fft>,              // forward FFT plan
    ifft: Arc<dyn Fft>,              // inverse FFT plan
    scratch: Vec<Complex64>,         // FFT workspace
    n: usize,
}
```

Both `potential_half` and `kinetic_full` are **precomputed once** when the integrator
is created. They don't change between time steps — they depend only on the grid,
the potential, and dt. This is a significant optimization: computing exp() is
expensive, and we'd be doing n of them per step otherwise.

### How the frequencies work

The FFT gives us Fourier coefficients in a specific order. For a grid of n points
spanning a domain of length L, the frequency corresponding to index i is:

```rust
let j = if i <= n/2 { i } else { i as f64 - n as f64 };
let k = 2π * j / L;
```

Indices 0 through n/2 correspond to positive frequencies (positive momenta). Indices
n/2+1 through n-1 correspond to negative frequencies — the FFT wraps them around.
This is standard FFT convention, not a quirk of our code.

The kinetic energy at each frequency is T(k) = k²/2 (in atomic units with m = 1),
and the phase factor we apply is exp(-i · k²/2 · dt).

### One time step

```rust
fn step(&mut self, wf: &mut Wavefunction) {
    // 1. Half-step in position space: ψ[i] *= exp(-i V(xᵢ) dt/2)
    for i in 0..n {
        psi[i] *= self.potential_half[i];
    }

    // 2. FFT to momentum space
    self.fft.process(psi);

    // 3. Full step in momentum space: ψ̃[i] *= exp(-i k²/2 dt)
    for i in 0..n {
        psi[i] *= self.kinetic_full[i];
    }

    // 4. IFFT back to position space
    self.ifft.process(psi);

    // 5. Normalize (rustfft convention: inverse FFT doesn't divide by n)
    for c in psi {
        *c /= n;
    }

    // 6. Half-step in position space again
    for i in 0..n {
        psi[i] *= self.potential_half[i];
    }
}
```

That's it. Six operations, four of which are pointwise multiplies on arrays, two
of which are FFTs. This runs at O(n log n) per step, which for n = 1024 is about
10,000 operations — fast enough for thousands of steps per frame.

### Why half-V, full-T, half-V?

This is the Strang splitting. The half-half sandwich makes the method **second-order
accurate**: the error per step is O(dt³), not O(dt²). You could do V-then-T or
T-then-V, but those are first-order — you'd need dt 10x smaller for the same accuracy.

The sandwich also makes the method **exactly unitary**: total probability is preserved
to machine precision, regardless of dt. This is a big deal — it means ∫|ψ|² dx stays
exactly 1.0 forever, even if the physics is slightly wrong due to time discretization.

### The rustfft normalization quirk

`rustfft` does not normalize either the forward or inverse FFT. The round-trip
FFT → IFFT gives you your original data multiplied by n. We divide by n after the
inverse FFT to compensate. This is just a library convention — mathematically,
the factor of 1/n can go on either the forward or inverse transform, or be split
as 1/√n on each. We put it all on the inverse side.

## Computing Observables

The simulation computes physical observables at each frame for display.

### Probability density

The simplest observable. For each grid point, compute |ψ(x)|²:

```rust
fn probability_density(&self) -> Vec<f64> {
    self.psi.iter().map(|c| c.norm_sqr()).collect()
}
```

This is what gets plotted as the blue filled curve in the viewer.

### Expected position ⟨x⟩

The quantum average of position — where you'd find the particle "on average":

```rust
fn expected_x(&self) -> f64 {
    let mut sum = 0.0;
    for i in 0..n {
        sum += x(i) * |psi[i]|²;
    }
    sum * dx  // discrete integral: sum * spacing
}
```

This is a weighted average of position, weighted by probability. For a Gaussian
centered at x0, this returns x0. For a packet bouncing off a wall, it tracks the
center of mass of the probability distribution.

### Expected momentum ⟨p⟩

Momentum isn't as simple as position. You can't just "read it off" the wavefunction.
The momentum operator in position space is -i ∂/∂x, so:

```rust
fn expected_p(&self) -> f64 {
    let mut sum = Complex64::zero();
    for i in 0..n {
        let dpsi_dx = (psi[i+1] - psi[i-1]) / (2 * dx);  // central difference
        sum += conj(psi[i]) * dpsi_dx;
    }
    (-i * sum * dx).re  // take real part (imaginary should be ~0)
}
```

The derivative is approximated by a central finite difference, with periodic
wrapping at the boundaries (index n wraps to 0, index -1 wraps to n-1). The
result should be purely real for a normalized state — the imaginary part is a
numerical artifact.

### Energy ⟨H⟩ = ⟨T⟩ + ⟨V⟩

Kinetic energy uses the second derivative (curvature of ψ):

```rust
fn expected_kinetic_energy(&self) -> f64 {
    let mut sum = Complex64::zero();
    for i in 0..n {
        let d2psi = (psi[i+1] - 2*psi[i] + psi[i-1]) / dx²;
        sum += conj(psi[i]) * d2psi;
    }
    (-0.5 * sum * dx).re
}
```

Potential energy is simpler — just the weighted average of V(x):

```rust
fn expected_potential_energy(&self, V: &Potential) -> f64 {
    let mut sum = 0.0;
    for i in 0..n {
        sum += V.value_at(x(i)) * |psi[i]|²;
    }
    sum * dx
}
```

Total energy is their sum. It should be conserved (constant over time) if the
simulation is working correctly. The test suite verifies this stays within 0.01
Hartree over 1000 steps for a displaced harmonic oscillator.

## The Simulation Loop

Putting it all together, the top-level `Simulation` struct bundles state and
provides a simple interface:

```rust
struct Simulation {
    wf: Wavefunction,
    potential: Potential,
    integrator: SplitOperator,
    time: f64,
    dt: f64,
}

impl Simulation {
    fn step(&mut self) {
        self.integrator.step(&mut self.wf);
        self.time += self.dt;
    }
}
```

The viewer calls `step()` some number of times per frame (configurable via a
"steps per frame" slider), then reads `wf.probability_density()` and the
observables for display. The simulation itself knows nothing about rendering —
it's a pure physics engine.

### Rebuilding the integrator

If the user changes the potential or dt, the precomputed phase factors are invalid.
The `rebuild_integrator()` method recomputes them:

```rust
fn rebuild_integrator(&mut self) {
    self.integrator = SplitOperator::new(&self.wf, &self.potential, self.dt);
}
```

This is a relatively expensive operation (n evaluations of exp()), so it only
happens on parameter changes, not every frame.

## What Each Scenario Demonstrates

The five single-particle scenarios are chosen to isolate specific quantum phenomena:

**Free particle** — no potential. The Gaussian spreads because different momentum
components travel at different speeds. This is pure dispersion, visible as the
probability curve widening while its peak drops. The total probability stays 1;
the area under the curve is constant.

**Infinite well** — tall walls on both sides. The wavefunction bounces back and
forth, interfering with itself. After enough time, standing wave patterns form —
these are the energy eigenstates. The specific pattern depends on the initial
conditions, but the quantized structure is visible in any case.

**Tunneling barrier** — a rectangular potential bump in the middle. Send a wave
packet at it. Some probability reflects, some tunnels through. Both outcomes
are present simultaneously — this is a superposition of "bounced back" and
"made it through." The transmission probability depends on the barrier height
and width relative to the particle's kinetic energy (k0²/2).

**Harmonic oscillator** — V = ½ω²x². The ground state Gaussian (σ = 1/√(2ω))
is stationary — its probability density doesn't change over time. A displaced
Gaussian oscillates back and forth like a classical spring, but also gradually
develops interference features as different energy components dephase.

**Double well** — two harmonic wells separated by a barrier. The wavefunction
tunnels between the wells on a timescale that depends on the barrier height.
Lower barriers → faster tunneling.

## Boundary Conditions and Periodic Artifacts

The FFT naturally imposes **periodic boundary conditions**: the right edge of the
grid wraps to the left edge. This means a wavefunction that hits the right boundary
reappears on the left.

For most scenarios this is invisible — the grid is wide enough (±15 to ±30 Bohr)
that ψ is effectively zero at the boundaries. But if you crank up the momentum or
let a free particle run for a long time, you'll eventually see it wrap around. This
isn't a bug; it's the correct behavior for a periodic domain. The "infinite well"
scenario handles this by putting large-V walls well inside the grid boundaries.
