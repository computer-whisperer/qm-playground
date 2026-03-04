# How the Spin Simulation Works

A code-level walkthrough of the Stage 4 additions. Read alongside the
[physics primer](04-spin.md) for the theory.

## The Data: A Two-Component Spinor

Where the scalar wavefunction had one complex number per grid point, the spinor
has two:

```rust
struct SpinorWavefunction {
    up:   Vec<Complex64>,   // ψ_↑(x) — spin-up component
    down: Vec<Complex64>,   // ψ_↓(x) — spin-down component
    n: usize,
    x_min: f64,
    x_max: f64,
    dx: f64,
}
```

Both components share the same spatial grid. Memory usage doubles compared to
Chapter 1, but the grid resolution stays at n = 1024. The total probability is:

```rust
fn norm(&self) -> f64 {
    self.up.iter().zip(self.down.iter())
        .map(|(u, d)| u.norm_sqr() + d.norm_sqr())
        .sum::<f64>() * self.dx
}
```

Both components contribute — the particle is *somewhere* with *some* spin, and
the probabilities must sum to 1.

## Initialization: Spatial Gaussian with Spin State

The initializer takes the usual spatial parameters (x0, sigma, k0) plus a spin
state given as two complex amplitudes:

```rust
fn set_gaussian(&mut self,
    x0: f64, sigma: f64, k0: f64,
    spin_up_amp: Complex64,      // α
    spin_down_amp: Complex64,    // β
) {
    // normalize spin amplitudes: |α|² + |β|² = 1
    let (alpha, beta) = normalize(spin_up_amp, spin_down_amp);

    for i in 0..n {
        let spatial = gaussian(x(i), x0, sigma, k0);
        up[i]   = alpha * spatial;
        down[i] = beta * spatial;
    }
    self.normalize();  // fix discrete grid normalization
}
```

The factored structure ψ_↑(x) = α · φ(x), ψ_↓(x) = β · φ(x) means the spin
and spatial parts are initially independent. Under evolution with a magnetic field,
this factorization can break — different spin components can develop different
spatial profiles (as in the Stern-Gerlach effect).

Common spin states:
- Pure spin-up: α = 1, β = 0
- Pure spin-down: α = 0, β = 1
- Spin along +x: α = 1, β = 1 (normalized to 1/√2 each)
- Spin along +y: α = 1, β = i

## The Magnetic Field

```rust
struct MagneticField {
    bz: f64,           // uniform longitudinal field
    bx: f64,           // uniform transverse field
    gradient_bz: f64,  // dB_z/dx (for Stern-Gerlach)
}
```

The total longitudinal field at position x is B_z(x) = bz + gradient_bz · x.

- **bz alone**: Zeeman effect. Each spin component sees a different potential.
  No mixing between ↑ and ↓.
- **bx alone**: Larmor precession. The spin rotates, mixing ↑ and ↓ at each
  spatial point. The spatial profiles don't separate.
- **gradient_bz**: Stern-Gerlach. The force F = -dV/dx differs for ↑ and ↓
  because V_↑ = V + B_z(x)/2 and V_↓ = V - B_z(x)/2 have different slopes.
  The two components are deflected in opposite directions.

## The Modified Integrator

The split-operator method extends naturally. The key insight: the kinetic energy
doesn't depend on spin (momentum doesn't care which way the spin points), but
the potential and magnetic terms do.

```rust
struct SpinorIntegrator {
    potential_half_up:   Vec<Complex64>,  // exp(-i V_↑(x) dt/2)
    potential_half_down: Vec<Complex64>,  // exp(-i V_↓(x) dt/2)
    kinetic_full:        Vec<Complex64>,  // exp(-i k²/2 dt) — same for both
    spin_cos: f64,                        // cos(B_x dt/4)
    spin_sin: f64,                        // sin(B_x dt/4)
    fft: ...,
    ifft: ...,
}
```

### Precomputed phase factors

The potential phase factors are now spin-dependent:

```rust
// V_↑(x) = V(x) + B_z(x)/2
potential_half_up[i] = exp(-i * (V(x_i) + B_z(x_i)/2) * dt/2)

// V_↓(x) = V(x) - B_z(x)/2
potential_half_down[i] = exp(-i * (V(x_i) - B_z(x_i)/2) * dt/2)
```

When B_z = 0, these are identical and the spin doesn't affect the spatial
evolution at all.

The kinetic phase factors are spin-independent — identical to Chapter 1.

### The spin rotation matrix

A transverse field B_x couples the spin components. The time evolution under
H_spin = B_x σ_x / 2 for a duration τ is:

```
exp(-i B_x σ_x τ/2) = cos(B_x τ/2) I - i sin(B_x τ/2) σ_x
                     = ( cos    -i sin )
                       ( -i sin  cos   )
```

We apply this as a half-step (τ = dt/2), so the precomputed values are
cos(B_x dt/4) and sin(B_x dt/4). When B_x = 0, sin = 0 and the rotation
is skipped entirely.

### One time step

```rust
fn step(&mut self, wf: &mut SpinorWavefunction) {
    // 1. Potential half-step (spin-dependent)
    for i in 0..n {
        wf.up[i]   *= potential_half_up[i];
        wf.down[i] *= potential_half_down[i];
    }

    // 2. Spin rotation half-step (if B_x ≠ 0)
    //    Pointwise 2×2 unitary on (up[i], down[i])
    if spin_sin != 0.0 {
        for i in 0..n {
            let u = wf.up[i];
            let d = wf.down[i];
            wf.up[i]   = cos * u - i_sin * d;
            wf.down[i] = -i_sin * u + cos * d;
        }
    }

    // 3. FFT both components
    fft.process(&mut wf.up);
    fft.process(&mut wf.down);

    // 4. Kinetic full step (same for both)
    for i in 0..n {
        wf.up[i]   *= kinetic_full[i];
        wf.down[i] *= kinetic_full[i];
    }

    // 5. IFFT both
    ifft.process(&mut wf.up);
    ifft.process(&mut wf.down);

    // 6. Normalize (1/n from rustfft convention)
    //    ... divide both up and down by n ...

    // 7. Spin rotation half-step (if B_x ≠ 0)
    // 8. Potential half-step
}
```

The structure is: V-spin-FFT-T-IFFT-spin-V, sandwiching the kinetic step with
potential and spin rotations on both sides (Strang splitting). This is second-
order accurate in dt and preserves unitarity, just like the scalar version.

**Cost**: two FFTs + two IFFTs of length n, plus O(n) pointwise operations.
Same asymptotic cost as Chapter 1, just with a constant factor of ~2×.

## Computing Spin Observables

The spin expectation values come from integrating the spinor components:

```rust
// ⟨σ_z⟩ = ∫ (|ψ_↑|² - |ψ_↓|²) dx
fn expected_sz(&self) -> f64 {
    self.up.iter().zip(self.down.iter())
        .map(|(u, d)| u.norm_sqr() - d.norm_sqr())
        .sum::<f64>() * self.dx
}

// ⟨σ_x⟩ = 2 Re ∫ ψ_↑* ψ_↓ dx
fn expected_sx(&self) -> f64 {
    let sum: Complex64 = self.up.iter().zip(self.down.iter())
        .map(|(u, d)| u.conj() * d)
        .sum();
    2.0 * (sum * self.dx).re
}

// ⟨σ_y⟩ = 2 Im ∫ ψ_↑* ψ_↓ dx
fn expected_sy(&self) -> f64 {
    // same inner product, take imaginary part
    2.0 * (sum * self.dx).im
}
```

These three numbers define the **Bloch vector** (⟨σ_x⟩, ⟨σ_y⟩, ⟨σ_z⟩). For a
pure spin state, this vector has length 1 and lies on the Bloch sphere. For a
mixed state (when spin and spatial degrees of freedom are entangled, as in
Stern-Gerlach after splitting), the Bloch vector is shorter — the spin is
no longer in a definite state.

The probabilities P(↑) = ∫|ψ_↑|² dx and P(↓) = ∫|ψ_↓|² dx give the measurement
outcome probabilities. They sum to 1 and relate to ⟨σ_z⟩ by
⟨σ_z⟩ = P(↑) - P(↓).

## The Connection to Chapter 3

The physics primer explains this in detail, but it bears repeating in code terms:

For two fermions with a spin-independent Hamiltonian, the total wavefunction
factors as (spatial) × (spin). The antisymmetry requirement on the total state
means:

- **Spin singlet** (antisymmetric spin, S=0) → **symmetric** spatial → Chapter 3's `Boson` mode
- **Spin triplet** (symmetric spin, S=1) → **antisymmetric** spatial → Chapter 3's `Fermion` mode

No new code is needed for this. Chapter 3's `ParticleSymmetry::Boson` already
simulates two fermions in a spin singlet, and `ParticleSymmetry::Fermion`
simulates two fermions in a spin triplet. The identical-particle scenarios
*were* spin physics all along.

## What Each Scenario Demonstrates

**Zeeman splitting**: Start in an equal spin superposition (α = β = 1/√2) in a
harmonic well with B_z ≠ 0. The ↑ and ↓ components see slightly different
potentials, so they oscillate at slightly different frequencies. The total
probability density shows a beating pattern as the two components drift in and
out of phase.

**Larmor precession**: Start spin-up in a B_x field. The transverse field
rotates the spin: ⟨σ_z⟩ oscillates between +1 and -1 at frequency B_x.
The spatial wavefunction barely changes (the harmonic trap keeps it still), but
the spin-up and spin-down probability densities swap back and forth. The period
is T = 2π/B_x.

**Stern-Gerlach**: Start with a moving wave packet in an equal spin superposition,
in a B_z gradient (B_z = gradient · x). The spin-up component feels a force
F = -gradient/2 (pulled one way), spin-down feels F = +gradient/2 (pulled the
other way). The packet splits into two beams, each with a definite spin. This is
the quantum mechanical version of the 1922 experiment that first revealed spin.
After splitting, ⟨σ_z⟩ ≈ 0 (equal parts up and down), but each spatial region
has a definite spin.

## Visualization

The spinor plot shows three curves:
- **Blue**: |ψ_↑(x)|² — spin-up probability density
- **Orange**: |ψ_↓(x)|² — spin-down probability density
- **White (dim)**: total |ψ_↑|² + |ψ_↓|² — total probability

The observables panel shows ⟨σ_z⟩, ⟨σ_x⟩, P(↑), and P(↓) alongside the usual
spatial observables. For Larmor precession, watch ⟨σ_z⟩ oscillate. For
Stern-Gerlach, watch the blue and orange curves separate spatially.
