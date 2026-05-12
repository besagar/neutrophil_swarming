# Setup 4 — 3D-cue swarming (placeholder)

> Originally drafted as Setup 3; renamed to Setup 4 when the 2D-radial swarm
> (same nondim groups as Setup 2, with a central trapping disk) was carved
> out into the new Setup 3. See [setup3_swarm.md](setup3_swarm.md).

**Status: under-specified.** The user has not yet provided the dynamics for
the cue field `L(x, y, z, t)` in the 3D-diffusion-above-2D-cells geometry.
Do not invent them.

## Geometry (agreed)

- Cells live on a 2D floor `z = 0`, positions `{xᵢ(t), yᵢ(t)}`.
- Cue `L` is defined in the 3D half-space `z ≥ 0` and diffuses in 3D.
- Cells sample `L|_{z=0}` and its in-plane gradient `∇₂ L|_{z=0}` to drive
  their polarization (same SDE as Setup 1/2 generalized to 2D `p = (p_x, p_y)`).

## Open questions to resolve with the user

1. What sources/sinks for `L`? Cells emit (autocrine/paracrine), boundary
   fluxes, decay rate? Is the LTB₄ relay explicit (only "emitting" cells
   produce, triggered by local `L` threshold)?
2. Boundary conditions at `z = 0` (no-flux vs. absorbing vs. mixed) and at
   `z = z_max`?
3. Discretization: full 3D PDE solve (FFT? finite difference?), or
   quasi-steady Green's-function approach with point sources at cell
   positions? The latter is much cheaper for ABM.
4. How many cells? (sets visualization tech: ≲10³ → Canvas/SVG; ≳10⁴ → WebGL.)

## Tentative scaffold to build now

- 2D agent renderer (regl or PixiJS) with N cells, each with position and
  polarization vector.
- Per-agent SDE step reusing the shared integrator from Setups 1–2.
- Pluggable `L`-provider interface: function `(x, y, t) → {L, ∇L}` so we can
  swap a stub (e.g. user-painted static field) for the real PDE solver
  later without touching agent code.
- UI for N, domain size, and a "cue source" stub.
