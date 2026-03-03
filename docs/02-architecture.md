# Architecture

## Crate Structure

```
crates/
├── sim-core/     Physics engine — no rendering dependencies
│   ├── lib.rs           Re-exports, top-level types
│   ├── wavefunction.rs  Complex wavefunction on a 1D grid
│   ├── potential.rs      Potential energy functions V(x)
│   ├── integrator.rs     Split-operator time evolution
│   └── units.rs          Physical constants, unit conversions
│
└── viewer/       Visualization app (egui + wgpu via eframe)
    ├── main.rs          Entry point, eframe setup
    └── app.rs           Simulation state, UI layout, rendering
```

## Design Principles

**sim-core is headless.** It knows nothing about rendering. This keeps it testable,
portable, and reusable when we later add ML training loops that need to run simulations
without a GUI.

**viewer drives the simulation.** It owns a sim-core `Simulation` and calls `step()`
each frame (or on demand). It reads the wavefunction state and renders it.

**Atomic units everywhere in sim-core.** No SI conversion inside the engine. The viewer
handles any unit display conversion for the UI.

## Time Evolution: Split-Operator Method

The integrator uses the split-step Fourier method:

```
ψ(t + dt) ≈ e^{-iV dt/2} · FFT⁻¹[ e^{-iT(k) dt} · FFT[ e^{-iV dt/2} · ψ(t) ] ]
```

Implementation:
1. Pre-compute `exp(-i V(x) dt/2)` for each grid point (potential half-step)
2. Pre-compute `exp(-i k²/2 dt)` for each frequency (kinetic full step)
3. Each time step: multiply by (1), FFT, multiply by (2), inverse FFT, multiply by (1)
4. Normalization: divide by N after inverse FFT (rustfft convention)

These pre-computed arrays are stored in the `SplitOperator` struct and only recomputed
when dt or the potential changes.

## Grid and Boundary Conditions

The 1D grid has N points from x_min to x_max. We use **periodic boundary conditions**
(implicit in the FFT). To simulate hard walls, we place a high potential barrier at the
edges. This is standard practice and simpler than implementing absorbing boundaries.

## Planned Extensions

As we progress through the pedagogical stages, the architecture will evolve:

- **Stage 2 (multi-particle)**: Wavefunction becomes multi-dimensional. Grid-based
  methods hit exponential scaling. This is where we'll explore tensor network or
  neural wavefunction representations.

- **Stage 3 (second quantization)**: State representation changes fundamentally —
  from particle positions to field occupation numbers. New crate likely needed.

- **Stage 4 (lattice gauge theory)**: Gauge fields on lattice links, matter fields
  on sites. Significant architectural rework expected.

Each stage may warrant its own crate within the workspace.

## WASM Support

eframe supports both native (wgpu) and web (WebGL/WebGPU) backends via feature flags.
The viewer's Cargo.toml has `native` and `web` features ready. To build for web:

```bash
# (future) requires trunk or wasm-pack setup
cargo build -p viewer --target wasm32-unknown-unknown --no-default-features --features web
```

sim-core has no platform dependencies and works on WASM as-is (rustfft is pure Rust).
