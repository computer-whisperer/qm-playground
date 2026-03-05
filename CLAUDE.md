# Emergence Sim

## Project Vision

Multi-scale physics simulation engine. Start from a quantum mechanical ground truth,
then train successive neural models to approximate emergent behavior at each higher scale.
The goal is to explore whether learned coarse-graining can bridge the gap from first-principles
QM to macro-scale phenomena.

## Current Stage

**Stage 5: Dirac Equation** — 1D relativistic Dirac equation with tunable speed of light.
Split-operator with 2×2 unitary rotation at each momentum mode. Zitterbewegung,
Klein tunneling, non-relativistic limit. Reuses SpinorWavefunction from Stage 4.
Prior stages (1-4) remain alongside.

## Repository Layout

Bare repo + worktree pattern:
- `.bare/` — bare git repository
- `main/` — main worktree (development happens here)
- Feature worktrees: `<branch-name>/` as siblings to `main/`

## Crate Structure

- `crates/sim-core/` — Physics simulation engine (no rendering dependencies)
- `crates/viewer/` — egui + wgpu visualization app (eframe)

## Building

```bash
cd main
cargo run -p viewer          # run the viewer
cargo test -p sim-core       # run physics tests
```

## Technical Notes

- **Units:** Atomic units (ℏ = m_e = e = 1). Energy in Hartree, length in Bohr radii, time in ℏ/E_h.
- **Time evolution:** Split-operator method via FFT (rustfft crate). 1D for single particle, 2D (row+column FFTs) for two particles.
- **Two-particle grid:** n=256 per axis (n²=65536 complex amplitudes). Domain [-30, 30] Bohr.
- **Entanglement:** Purity Tr(ρ₁²) via O(n³) direct reduced density matrix computation. Schmidt number K = 1/Purity.
- **Identical particles:** `ParticleSymmetry` enum (Distinguishable/Boson/Fermion). Symmetrized init via `set_symmetrized_gaussian()`. Optional `symmetrize()` projection for numerical drift.
- **Spinor:** `spinor` module for spin-1/2 particles. `SpinorWavefunction` (2-component), `SpinorIntegrator` (split-operator with Zeeman + precession). `MagneticField` supports uniform B_z, B_x, and position-dependent gradient.
- **Dirac:** `dirac` module for 1D Dirac equation. `DiracIntegrator` applies analytical 2×2 rotation exp(-iH_k dt) at each momentum mode. Tunable c (speed of light) and mass. `init_positive_energy_packet()` projects onto positive-energy Dirac spinors to suppress Zitterbewegung.
- **WASM pivot:** eframe feature flags (`native` / `web`) allow future browser deployment.

## Docs

See `docs/` for physics background and architecture notes. Start with `00-project-vision.md`.
