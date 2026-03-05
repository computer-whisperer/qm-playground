# How the Dirac Simulation Works

A code-level walkthrough of the Stage 5 additions. Read alongside the
[physics primer](05-dirac-equation.md) for the theory.

## Reusing the Spinor Infrastructure

The Dirac equation is a two-component equation, just like the spin system from
Chapter 4. The upper component is the "particle" and the lower is the
"antiparticle" content. We reuse `SpinorWavefunction` from the spinor module —
same data structure, different physics.

```rust
pub struct DiracSimulation {
    pub wf: SpinorWavefunction,  // reused from spinor module
    pub potential: Potential,
    pub c: f64,                   // speed of light (tunable)
    pub mass: f64,
    pub integrator: DiracIntegrator,
    pub time: f64,
    pub dt: f64,
}
```

The key difference from the spinor integrator: instead of a non-relativistic
kinetic step (k²/2m) with magnetic field mixing, the Dirac integrator applies
a relativistic 2×2 unitary rotation at each momentum mode.

## The DiracIntegrator

### Precomputed data

```rust
struct DiracIntegrator {
    potential_half: Vec<Complex64>,     // exp(-iV dt/2) at each x
    kinetic_params: Vec<KineticRotation>, // rotation params at each k
    fft: Arc<dyn Fft<f64>>,
    ifft: Arc<dyn Fft<f64>>,
    scratch: Vec<Complex64>,
    n: usize,
}

struct KineticRotation {
    cos_e: f64,       // cos(E_k dt)
    sin_e: f64,       // sin(E_k dt)
    ck_over_e: f64,   // ck / E_k
    mc2_over_e: f64,  // mc² / E_k
}
```

At construction, we compute for each momentum mode k:
- The relativistic energy E_k = √(c²k² + m²c⁴)
- The trigonometric rotation parameters cos(E_k dt) and sin(E_k dt)
- The mixing ratios ck/E_k and mc²/E_k

These four numbers fully specify the 2×2 unitary rotation at each k.

### Why four parameters?

The free Dirac Hamiltonian in momentum space at a given k is:

```
H_k = ( mc²    ck  )
      ( ck    -mc² )
```

Its time evolution is exp(-iH_k dt) = cos(E_k dt) I - i sin(E_k dt) H_k/E_k.
Written as a matrix:

```
( cos - i sin·mc²/E      -i sin·ck/E    )
(    -i sin·ck/E      cos + i sin·mc²/E  )
```

The diagonal elements are complex conjugates of each other. The off-diagonal
elements are identical. So the entire 2×2 matrix is determined by cos(E_k dt),
sin(E_k dt), ck/E_k, and mc²/E_k.

### The time step

```rust
fn step(&mut self, wf: &mut SpinorWavefunction) {
    // 1. Potential half-step (same for both components)
    for i in 0..n {
        wf.up[i]   *= potential_half[i];
        wf.down[i] *= potential_half[i];
    }

    // 2. FFT both components to momentum space
    fft.process(&mut wf.up);
    fft.process(&mut wf.down);

    // 3. Dirac rotation at each k
    for i in 0..n {
        let r = &kinetic_params[i];
        let a = Complex64::new(r.cos_e, -r.sin_e * r.mc2_over_e);
        let b = Complex64::new(0.0, -r.sin_e * r.ck_over_e);

        let u = wf.up[i];
        let d = wf.down[i];
        wf.up[i]   = a * u + b * d;
        wf.down[i] = b * u + a.conj() * d;
    }

    // 4. IFFT back to position space
    ifft.process(&mut wf.up);
    ifft.process(&mut wf.down);

    // 5. Normalize (rustfft convention: divide by n)
    // 6. Potential half-step
}
```

The structure is V/2 → FFT → Dirac rotation → IFFT → V/2 (Strang splitting),
same sandwich as every other integrator in the project. The physics difference
is entirely in step 3: instead of multiplying by a scalar phase exp(-ik²dt/2m),
we apply a 2×2 matrix that mixes the two components.

### Comparing with the spinor integrator

| Feature | Spinor (Ch. 4) | Dirac (Ch. 5) |
|---------|---------------|---------------|
| Kinetic step | exp(-ik²/2m dt) (scalar, same for both) | 2×2 rotation (mixes components) |
| Mass term | Not present | mc²σ_z (splits positive/negative energy) |
| Potential | Spin-dependent via B_z | Same for both components |
| Dispersion | E = k²/2m (parabolic) | E = ±√(c²k² + m²c⁴) (hyperbolic) |
| Spin mixing | B_x rotation (separate step) | Built into kinetic rotation |

## Positive-Energy Initialization

The most important new method. A naive initialization (all amplitude in the
upper component) creates a superposition of positive and negative energy states,
producing Zitterbewegung. To avoid this, we project onto the positive-energy
sector.

### The Dirac spinor structure

At each momentum k, the positive-energy eigenstate of H_k is:

```
u(k) ∝ (     1        )
        ( ck/(E_k+mc²) )
```

The ratio ck/(E_k + mc²) is the "small component" — it's small when k << mc
(non-relativistic) and approaches 1 when k >> mc (ultra-relativistic).

### The algorithm

```rust
fn init_positive_energy_packet(&mut self, x0: f64, sigma: f64, k0: f64) {
    // 1. Create a Gaussian in the upper component only
    //    (standard spatial Gaussian with momentum k0)

    // 2. FFT the upper component to momentum space

    // 3. At each k, project onto positive-energy spinor:
    //    f(k) → f(k)/N · (1, ck/(E_k+mc²))
    //    where N = √(1 + (ck/(E_k+mc²))²) preserves the norm

    // 4. IFFT both components back to position space
    // 5. Normalize
}
```

After this initialization, each momentum mode has the correct ratio of upper
to lower components for a positive-energy state. The wave packet propagates
cleanly at its group velocity v_g = c²k₀/E_{k₀} without Zitterbewegung.

### What the small component looks like

In position space after the IFFT, the lower component is generally a modified
version of the upper component — slightly shifted and with different amplitude.
For a packet centered at k₀:
- At low k₀ (non-relativistic): the lower component is tiny (proportional to k₀/mc)
- At high k₀ (ultra-relativistic): the lower component approaches the upper in magnitude

This is the correct physical behavior: the "small component" of the Dirac
spinor carries the antiparticle content, which is negligible at low energies.

## The Speed of Light as a Knob

In atomic units, c ≈ 137. At this value, relativistic effects are tiny for
everyday momenta. We make c a tunable parameter:

- **c = 5**: strongly relativistic. Zitterbewegung is fast, Klein tunneling is
  visible with moderate barriers, the lower component has significant amplitude.
- **c = 10-20**: moderate relativistic effects. Good for seeing group velocity
  deviate from the non-relativistic p/m.
- **c = 30-50**: nearly non-relativistic. The lower component is tiny, the
  dispersion is approximately parabolic, and behavior matches Schrödinger.

The gap energy 2mc² determines the Klein tunneling threshold: barriers taller
than this become partially transparent.

## What Each Scenario Demonstrates

### Zitterbewegung

- **Init**: Naive (upper component only), nonzero momentum, c = 5
- **What happens**: The lower (orange) component lights up immediately as
  negative-energy modes are excited. The position ⟨x⟩ oscillates at
  frequency ~2mc².
- **Key observable**: P(↓) grows from 0 and fluctuates. Compare with
  positive-energy init where P(↓) stays constant.

### Klein tunneling

- **Init**: Positive-energy packet heading toward a barrier with V₀ > 2mc²
- **What happens**: Significant transmission through the barrier, in contrast
  to exponential suppression in non-relativistic tunneling.
- **Key observable**: The lower component lights up inside the barrier
  (antiparticle propagation), and a transmitted packet appears on the far side.
- **Why it works**: Inside the barrier, the particle energy aligns with the
  negative-energy branch, allowing propagating solutions.

### Non-relativistic limit

- **Init**: Positive-energy packet, high c (c = 30)
- **What happens**: Behavior matches the Schrödinger equation from Chapter 1.
  The lower component is nearly invisible. Dispersion is parabolic.
- **Key observable**: P(↓) ≈ 0, the group velocity is approximately k₀/m,
  and the packet spreads like a non-relativistic Gaussian.

## Cost Analysis

Same asymptotic cost as the spinor integrator: two FFTs + two IFFTs of length n,
plus O(n) pointwise operations. The Dirac rotation is slightly more expensive
per point (two complex multiplies and adds instead of one scalar multiply) but
the FFTs dominate.

| Grid | Per step | Memory |
|------|----------|--------|
| n = 1024 | ~0.1 ms | 2 × 1024 × 16 B ≈ 32 KB wavefunction |
| n = 2048 | ~0.2 ms | ~64 KB |
