# The Dirac Equation

Everything we've done so far has been non-relativistic: the Schrödinger equation
treats time and space asymmetrically (first derivative in t, second in x), and
the energy-momentum relation is E = p²/2m. This breaks down when particles move
at appreciable fractions of the speed of light.

The Dirac equation is the relativistic wave equation for spin-1/2 particles. It
doesn't just add corrections to Schrödinger — it's a fundamentally different
equation that **predicts** spin, **predicts** antimatter, and produces phenomena
that have no non-relativistic analog.

## Why Not Just Relativize Schrödinger?

The relativistic energy-momentum relation is:

```
E² = p²c² + m²c⁴
```

The obvious approach: replace E → iℏ∂/∂t and p → -iℏ∂/∂x to get the
**Klein-Gordon equation**:

```
-ℏ² ∂²ψ/∂t² = -ℏ²c² ∂²ψ/∂x² + m²c⁴ ψ
```

This is second-order in time (unlike Schrödinger, which is first-order). The
problem: |ψ|² is no longer a conserved probability density. It can go negative.
The Klein-Gordon equation describes relativistic *scalar* particles (spin-0),
but it's really a field equation, not a probability-amplitude equation.

Dirac's insight: take the **square root** of the energy-momentum relation.
Write E as a *linear* function of p, at the cost of making the coefficients
matrices:

```
E = c α p + β mc²
```

where α and β are matrices satisfying α² = β² = I and αβ + βα = 0. In 1D,
these can be 2×2 matrices — the Pauli matrices work: α = σ_x, β = σ_z.

The resulting equation is first-order in both time and space:

```
iℏ ∂ψ/∂t = (-iℏc σ_x ∂/∂x + mc² σ_z) ψ
```

where ψ is now a 2-component spinor — not by assumption, but as a consequence
of the mathematics. **Spin emerges from requiring a first-order relativistic
wave equation.**

## The 1D Dirac Hamiltonian

In our units (ℏ = 1, but keeping c and m explicit):

```
H = -ic σ_x ∂/∂x + mc² σ_z + V(x)

  = ( mc² + V(x)      -ic ∂/∂x   )
    ( -ic ∂/∂x      -mc² + V(x)  )
```

The two components of the spinor are coupled by the off-diagonal derivative
terms. The diagonal mc²σ_z term splits the spectrum: one component is
associated with positive energy (+mc²), the other with negative energy (−mc²).

### The dispersion relation

For a free particle (V = 0), plane wave solutions ψ ~ exp(ikx - iEt) satisfy:

```
E = ±√(c²k² + m²c⁴) = ±E_k
```

Two branches:
- **Positive energy** (+E_k ≥ +mc²): ordinary particles
- **Negative energy** (−E_k ≤ −mc²): antiparticle solutions

There's a **gap** of 2mc² between the branches. No states exist with |E| < mc².
This gap is where the mass lives — the rest energy.

Compare with Schrödinger: E = k²/2m has only positive energies and no gap. The
Dirac spectrum is qualitatively different.

### The speed of light as a parameter

In atomic units, c ≈ 137. At everyday momenta (k ~ 1), the kinetic energy
k²/2m ~ 0.5 Hartree while mc² ~ 137² / 2 ~ 9400 Hartree. Relativistic effects
are tiny.

In the simulation, we make c a tunable parameter. Setting c = 5 or c = 10
produces visible relativistic effects at moderate momenta. Setting c = 137
recovers essentially non-relativistic behavior. This lets you see exactly how
the Schrödinger equation emerges as a limit of the Dirac equation.

## The Split-Operator Method for Dirac

The split-operator approach works, but the kinetic step is different.

### Momentum space

In momentum space, the free Dirac Hamiltonian at each k is a 2×2 matrix:

```
H_k = ( mc²    ck  )
      ( ck    -mc² )
```

The time evolution for one step is exp(-iH_k dt), which is a 2×2 unitary
rotation:

```
exp(-iH_k dt) = cos(E_k dt) I - i sin(E_k dt) (H_k / E_k)
```

where E_k = √(c²k² + m²c⁴). This is computed analytically — no matrix
exponentiation needed.

### The full step

```
1. V half-step in position space:  ψ *= exp(-i V(x) dt/2)
2. FFT to momentum space
3. Dirac kinetic full step:       apply 2×2 rotation at each k
4. IFFT to position space
5. V half-step:                   ψ *= exp(-i V(x) dt/2)
```

The potential step is scalar (V(x) multiplies both components equally for an
electrostatic potential). The kinetic step is the 2×2 rotation that mixes the
components — this is where the relativistic physics lives.

## Key Phenomena

### Zitterbewegung (Trembling Motion)

If you initialize a Dirac particle with all amplitude in the upper component
(the "particle" component), you've implicitly created a superposition of
positive and negative energy states. These interfere, producing a rapid
oscillation in ⟨x⟩ at frequency 2mc²/ℏ — the **Zitterbewegung**.

This is not a physical oscillation of a real electron (its frequency is ~10²¹ Hz,
far too fast to observe directly). It's an artifact of projecting onto a single
component. But in the simulation with reduced c, it's clearly visible: the
expected position jitters even for a free particle.

Properly initializing in the positive-energy sector (both components with the
correct Dirac spinor structure) eliminates Zitterbewegung. The particle then
propagates cleanly at its group velocity.

### Klein Tunneling

The most dramatic relativistic effect. In non-relativistic quantum mechanics,
tunneling probability decreases exponentially with barrier height. In the Dirac
equation, the opposite can happen.

When a potential barrier V₀ exceeds the gap energy 2mc², something remarkable
occurs: the barrier allows **propagating solutions** because the particle's
energy inside the barrier matches the negative-energy (antiparticle) branch.
The particle effectively converts to an antiparticle inside the barrier and
back to a particle on the other side.

For V₀ >> 2mc², the transmission probability approaches **1.0** — the barrier
becomes transparent. This is the Klein paradox, and it's why you can't confine
relativistic particles with potential barriers alone.

In the simulation: send a wave packet at a barrier taller than 2mc². With the
Schrödinger equation, you'd see exponentially suppressed tunneling. With Dirac,
you see significant or even complete transmission.

### The Non-Relativistic Limit

As c → ∞ (or equivalently, for momenta k << mc), the Dirac equation reduces
to the Schrödinger equation. The lower component becomes small (proportional to
v/c), and the upper component satisfies:

```
i ∂ψ_↑/∂t ≈ (-1/(2m) ∂²/∂x² + V) ψ_↑
```

which is exactly the Schrödinger equation from Chapter 1.

In the simulation: increase c and watch the Dirac wavefunction converge to
Schrödinger behavior. The lower component shrinks, Zitterbewegung disappears,
and tunneling through barriers recovers the exponential suppression.

## Antiparticles

The negative-energy solutions of the Dirac equation were initially a puzzle.
Dirac proposed that the "vacuum" is a sea of filled negative-energy states
(the Dirac sea), and a hole in this sea appears as a particle with positive
charge — an antiparticle. This predicted the positron, discovered in 1932.

Modern quantum field theory replaces the Dirac sea with the concept of
antiparticle creation operators, but the mathematical structure is the same.
The negative-energy branch is real and physical — it's why pair creation
(e⁺e⁻ from a photon) is possible and why Klein tunneling works.

In our 1D simulation, the lower spinor component carries the antiparticle
content. When it has significant amplitude (as in Zitterbewegung or Klein
tunneling), antiparticle physics is active.

## What to Watch For in the Visualization

### Zitterbewegung scenario

- Initialize with amplitude only in the upper component
- Watch ⟨x⟩ oscillate rapidly even for a free particle
- The upper (blue) and lower (orange) components exchange amplitude
- Increase c to slow the oscillation frequency (it scales as mc²)
- Use positive-energy initialization to eliminate it

### Klein tunneling scenario

- A wave packet hits a potential barrier taller than 2mc²
- Non-relativistic expectation: exponentially suppressed tunneling
- Dirac reality: significant transmission through the barrier
- Increase barrier height → transmission *increases* (Klein paradox)
- The lower component lights up inside the barrier (antiparticle content)

### Non-relativistic limit scenario

- Compare Dirac and Schrödinger for the same setup
- At low c (c ~ 5): dramatic differences in dispersion and tunneling
- At high c (c ~ 50): behavior converges to Schrödinger
- The lower component amplitude decreases as c increases

## How the Code Works

For a walkthrough of the 2×2 momentum-space rotation, the relativistic
dispersion relation, and positive-energy initialization — see the
[implementation walkthrough](05-implementation-walkthrough.md).

## What Comes Next

The Dirac equation is the end of single-particle relativistic quantum mechanics.
To go further — many particles, particle creation and annihilation, the full
structure of matter — requires **quantum field theory**. That's where the
project's coarse-graining ambitions start to become relevant: can a neural
network learn to approximate the effective dynamics of a quantum field?
