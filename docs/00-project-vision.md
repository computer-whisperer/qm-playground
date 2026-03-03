# Project Vision: Emergence Sim

## The Core Idea

Physics has a remarkable property: behavior at macro scales is *dramatically simpler* than
the micro-scale physics that produces it. A bouncing ball doesn't require tracking 10²³
individual molecules — a handful of numbers (position, velocity, elasticity, mass) suffice.

This simplification isn't a hack. It's a structural property of reality called **universality**:
enormous classes of micro-configurations produce identical macro behavior. Temperature doesn't
care which specific molecules are fast. Pressure doesn't care which molecules hit the wall when.

The question this project explores: **can neural networks learn these compressions automatically?**

## The Layer Stack

```
Layer 0: Quantum mechanics (ground truth)
         Solve the Schrödinger equation for small systems
         State: wavefunctions in high-dimensional space
         ↓ train neural model on QM simulation data ↓

Layer 1: Learned interatomic potentials
         Neural net predicts forces between atoms
         State: atom positions + types
         ↓ train on molecular dynamics trajectories ↓

Layer 2: Coarse-grained dynamics
         Groups of atoms become single "beads"
         State: bead positions, effective interactions
         ↓ train on coarse-grained simulations ↓

Layer 3: Continuum / mesoscale
         Fields on a grid (density, velocity, stress)
         State: discretized PDEs
         ↓ train on field simulations ↓

Layer 4: Macro-scale
         Objects with material properties
         State: rigid bodies, deformable meshes
```

Each layer compresses by orders of magnitude. Each layer's training data comes from running
the layer below on many scenarios.

## Why Start at Quantum Mechanics?

We *could* start at an even deeper level (quantum field theory, lattice QCD), and we plan
to explore that direction for educational purposes. But for the emergence pipeline, QM is
the practical sweet spot:

- It's the level from which chemistry emerges
- Chemistry is the gateway to materials, biology, everything macro
- QFT corrections are negligible for condensed matter / chemistry
- The computational challenge (exponential state space) is exactly what neural nets might compress
- There's existing validation: neural network potentials (SchNet, NequIP, MACE) already demonstrate Layer 0→1

## The Build-to-Learn Approach

This project serves dual purposes:

1. **Educational**: Build up from first-principles physics to develop real intuition
   for quantum mechanics, statistical mechanics, and the emergence of macro behavior.

2. **Research**: Explore whether stacked learned approximations can bridge the
   micro→macro gap without catastrophic error accumulation.

We expect to rip out and rework major pieces multiple times. The point is learning — both
for the human and for the neural models.

## Known Failure Modes

Where this is most likely to break down:

- **Phase transitions**: Macro behavior suddenly depends on micro details at critical
  points (freezing, fracture). Expert models need to detect when their approximation breaks.

- **Error compounding**: Even 99% accuracy per layer might produce qualitative failures
  after stacking 4-5 layers. Though universality suggests errors might wash out.

- **Training data coverage**: Layer 0 is expensive, limiting Layer 1's training data.
  Each layer's data budget constrains the next.

These are researchable questions, not known impossibilities.

## Pedagogical Progression

The simulation engine itself follows a build-to-learn progression:

1. **Single-particle QM** (1D → 3D Schrödinger equation) ← *current stage*
2. **Multi-particle QM** (entanglement, exponential scaling)
3. **Second quantization** (particles → fields)
4. **Lattice gauge theory** (discretized QFT, the actual subatomic ground truth)

Each stage builds the conceptual vocabulary needed for the next.
