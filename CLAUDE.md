# Emergence Sim

## Project Vision

Multi-scale physics simulation engine. Start from a quantum mechanical ground truth,
then train successive neural models to approximate emergent behavior at each higher scale.
The goal is to explore whether learned coarse-graining can bridge the gap from first-principles
QM to macro-scale phenomena.

## Current Stage

**Stage 1: Single-particle quantum mechanics** — 1D Schrödinger equation with split-operator
time evolution. Building intuition for wavefunctions, tunneling, quantization before moving
to multi-particle systems and eventually QFT.

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
- **Time evolution:** Split-operator method via FFT (rustfft crate).
- **WASM pivot:** eframe feature flags (`native` / `web`) allow future browser deployment.

## Docs

See `docs/` for physics background and architecture notes. Start with `00-project-vision.md`.
