# Quantum Mechanics Primer

This document covers the physics behind the Stage 1 simulation — enough to understand
what the code is doing and *why*.

## The Big Picture

Classical mechanics says a particle has a definite position and momentum at all times.
Quantum mechanics says: no. A particle is described by a **wavefunction** ψ(x,t) — a
complex-valued function spread over all of space.

The wavefunction is not a fuzzy cloud around a "real" position. There *is no* definite
position until you measure. The wavefunction is the complete description of the particle's state.

## The Wavefunction

ψ(x,t) is a complex number at every point x and time t.

**What it means:**
- |ψ(x)|² gives the **probability density** of finding the particle at position x
- The integral of |ψ|² over all space must equal 1 (the particle is *somewhere*)
- The phase of ψ (its angle in the complex plane) determines interference and momentum

**Key intuition:** ψ behaves like a wave — it can interfere constructively and
destructively. Two-slit experiments work because ψ passes through *both* slits, and the
two paths interfere.

## The Schrödinger Equation

The time evolution of ψ is governed by the Schrödinger equation:

```
iℏ ∂ψ/∂t = Ĥψ
```

The left side is i (the imaginary unit) times ℏ times the rate of change of ψ.
The right side is the **Hamiltonian operator** Ĥ applied to ψ. The Hamiltonian
encodes the total energy of the system — it determines how the wavefunction
evolves. The hat (^) is just notation to remind you it's an operator (something
that acts on functions), not a number.

For a single particle in a potential V(x), the Hamiltonian has two parts:

```
iℏ ∂ψ/∂t = -ℏ²/(2m) ∂²ψ/∂x² + V(x)ψ
```

Breaking this down:
- Left side: how ψ changes over time
- Right side, first term: **kinetic energy** — the second spatial derivative measures
  curvature. A highly curved wavefunction has high momentum and thus high kinetic energy.
- Right side, second term: **potential energy** — just V(x) multiplied by ψ at each point

In **atomic units** (ℏ = m = 1), this simplifies to:

```
i ∂ψ/∂t = -½ ∂²ψ/∂x² + V(x)ψ
```

This is what our code solves.

## Atomic Units

To avoid carrying around tiny numbers (ℏ ≈ 1.055 × 10⁻³⁴ J·s), we work in units where
the fundamental constants are 1:

| Quantity | Atomic unit | SI value |
|----------|-------------|----------|
| ℏ (action) | 1 | 1.055 × 10⁻³⁴ J·s |
| mₑ (electron mass) | 1 | 9.109 × 10⁻³¹ kg |
| e (charge) | 1 | 1.602 × 10⁻¹⁹ C |
| Length (Bohr radius) | 1 | 0.529 × 10⁻¹⁰ m |
| Energy (Hartree) | 1 | 4.360 × 10⁻¹⁸ J = 27.2 eV |
| Time (ℏ/Eₕ) | 1 | 2.419 × 10⁻¹⁷ s |

So "x = 5" in our simulation means 5 Bohr radii ≈ 2.6 Ångströms.

## Key Phenomena to Build Intuition For

### Quantization
In a confined space (particle in a box), only certain wavelengths fit. Since wavelength
determines energy, energy is **quantized** — only discrete values are allowed. This is
why atoms have specific energy levels.

### Tunneling
A classical ball can't roll over a hill if it doesn't have enough energy. A quantum
wavefunction can *leak through* a potential barrier, even when the particle's energy is
less than the barrier height. The probability of tunneling decreases exponentially with
barrier width and height.

### Superposition
A wavefunction can be a sum of multiple states simultaneously. A particle can be "in
two places at once" — not metaphorically, but as a direct consequence of the linearity
of the Schrödinger equation.

### Uncertainty Principle
A narrow wavefunction (well-defined position) must contain many wavelengths (spread
in momentum). You can't make both narrow simultaneously:

```
Δx · Δp ≥ ℏ/2
```

This isn't a measurement limitation — it's a fundamental property of waves.

## The Split-Operator Method

Our code solves the Schrödinger equation using the **split-operator** (split-step Fourier)
method. Here's the idea:

The Schrödinger equation says how ψ evolves in time. Over a short time step
dt, the new wavefunction is related to the old by a "time evolution operator":

```
ψ(t + dt) = e^{-iĤ dt} ψ(t)
```

The e^{...} here is a matrix exponential (or operator exponential) — a
generalization of the ordinary exponential function to operators. The hat on Ĥ
just marks it as an operator (something that acts on wavefunctions, not a
plain number). The key fact: e^{-iĤ dt} is **unitary**, which means it
preserves total probability. No probability leaks or appears — the particle
stays normalized.

The Hamiltonian splits into kinetic and potential parts: Ĥ = T̂ + V̂. These
don't commute (applying T then V gives a different result than V then T), but
for small dt we can approximate by alternating:

```
e^{-i(T̂+V̂)dt} ≈ e^{-iV̂ dt/2} · e^{-iT̂ dt} · e^{-iV̂ dt/2}
```

This is called **Strang splitting** and is accurate to order dt² (the error
per step is proportional to dt³).

The key trick:
- **V̂ is diagonal in position space**: e^{-iV dt/2} just multiplies each grid
  point by a complex phase factor e^{-iV(x)dt/2}. This is just one multiply
  per grid point.
- **T̂ is diagonal in momentum space**: e^{-iT dt} just multiplies each Fourier
  coefficient by e^{-ik²dt/2}. Also one multiply per grid point, but in
  Fourier space.
- So we: apply V half-step → FFT → apply T full step → inverse FFT → apply V half-step

Each half is simple pointwise multiplication. The FFT and inverse FFT (O(N log N))
are the expensive parts. The result preserves unitarity — total probability
stays exactly 1.

### Why FFT?

The kinetic energy operator involves ∂²/∂x². In Fourier space, derivatives become
multiplication: ∂²ψ/∂x² → -k²ψ̃(k). So the kinetic energy in momentum space is
simply T(k) = k²/2.

The FFT lets us switch between position space (where V is simple) and momentum space
(where T is simple) in O(N log N) time.

## What to Watch For in the Visualization

When running the simulation:

1. **Gaussian wave packet in free space**: Watch it spread. This *is* the uncertainty
   principle — the initial position uncertainty grows because different momentum
   components travel at different speeds.

2. **Particle hitting a barrier**: Part of ψ reflects, part tunnels through. The
   transmitted and reflected parts are *simultaneously real* — this is superposition.

3. **Particle in a box**: After initial transients, you'll see standing wave patterns —
   these are the energy eigenstates. Their frequencies are the quantized energy levels.

4. **Harmonic oscillator**: The ground state is a Gaussian that doesn't spread (it's
   an eigenstate). Displaced from center, it oscillates back and forth — the quantum
   analog of a classical spring.

## How the Code Works

For a walkthrough of how these ideas map to actual Rust code — the data structures,
the FFT loop, how observables are computed — see the
[implementation walkthrough](01-implementation-walkthrough.md).

## What Comes Next

Continue to [Chapter 2: Two-Particle QM](02-two-particle-qm.md) — where entanglement
appears and the exponential scaling of quantum mechanics becomes visceral.
