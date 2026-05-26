# Implementation plan — GL Motility interactive tool

This plan covers the four interactive subpages described in
[../ginzburg_landay_neutrophils.md](../ginzburg_landay_neutrophils.md).
Per-setup physics and nondimensionalization choices live in
[physics/setup1_uniform.md](physics/setup1_uniform.md),
[physics/setup2_wave.md](physics/setup2_wave.md),
[physics/setup3_swarm.md](physics/setup3_swarm.md),
[physics/setup4_swarm3d.md](physics/setup4_swarm3d.md). UI conventions in
[design/ui_conventions.md](design/ui_conventions.md).

---

## 0. Goals and non-goals

**Goals.**
- One static site, three subpages, Mathematica-`Manipulate`-style live
  exploration of the GL polarization model.
- Strict separation: physics core runs in nondimensional units; UI exposes
  dimensional and/or nondimensional knobs as the user chooses, kept in
  sync via per-setup linkage functions.
- Run locally with `python -m http.server` (or any static host). No build,
  no backend.

**Non-goals.**
- Not a publication-grade simulator. We accept simple integrators
  (Euler–Maruyama / RK4) and modest agent counts.
- Not a generic GL toolkit; it is hard-coded to this neutrophil model.
- No installation, no npm, no Webpack. CDN ES modules only.

---

## 1. Tech stack

| Concern            | Choice                          | Why                                        |
|--------------------|---------------------------------|--------------------------------------------|
| Page shell         | Vanilla HTML + ES modules       | Zero build, easy to share / version        |
| Math / RNG         | Hand-rolled in `shared/`        | Avoid heavy deps; keep numerics auditable  |
| Static plots       | Plotly.js (CDN)                 | Free-energy curves, bifurcation, summaries |
| Streaming traces   | uPlot (CDN)                     | Light, fast, no React wrapper needed       |
| 2D field render    | Canvas2D (Setup 2 cue strip)    | Trivial 1D heatmap                         |
| Swarm render       | regl or PixiJS (CDN)            | 10³–10⁵ agents at 60 fps                   |
| Heavy sim loops    | Web Worker per setup            | Keep UI thread responsive                  |
| State / reactivity | Tiny pub/sub in `shared/state.js` | No framework needed for ~30 knobs        |

If at some point we need a build (TypeScript? bundling regl?), revisit —
but only if the cost of doing without becomes obvious.

---

## 2. Repo and module layout

```
├── index.html                  # landing page, links to 3 setups
├── shared/
│   ├── rng.js                  # seeded PRNG + Gaussian (Box–Muller / Ziggurat)
│   ├── sde.js                  # Euler–Maruyama, RK4, adaptive optional
│   ├── units.js                # generic dim ↔ nondim linkage primitives
│   ├── state.js                # tiny observable store
│   ├── widgets.js              # slider / toggle / number / play-pause helpers
│   ├── plot_plotly.js          # thin Plotly wrappers
│   ├── plot_uplot.js           # thin uPlot wrappers
│   └── style.css
├── setup1/
│   ├── index.html              # page shell
│   ├── model.js                # nondim SDE (uniform L)
│   ├── nondim.js               # dim ↔ nondim mapping for setup 1
│   ├── ui.js                   # widget config + wiring
│   ├── plots.js                # F̃(p̃), trace, histogram, bifurcation
│   └── worker.js               # runs the SDE in a worker
├── setup2/
│   ├── index.html
│   ├── model.js                # nondim SDE + Gaussian wave + v(p) maps
│   ├── nondim.js
│   ├── ui.js
│   ├── plots.js                # cue+cell animation, traces, Δx̃ summary
│   └── worker.js
├── setup3/
│   ├── index.html
│   ├── agents.js               # per-cell state, batched SDE step
│   ├── L_provider.js           # pluggable cue field (stub for now)
│   ├── nondim.js
│   ├── ui.js
│   ├── render.js               # regl/PixiJS draw loop
│   └── worker.js
└── setup4/
    ├── index.html
    ├── solvers/
    │   ├── index.js            # createSolver(geometry, model) factory
    │   ├── field.js            # grid data + PIC ops: sample, accumulateSource, getRadialProfile
    │   ├── solver_m1.js        # M1 stepper (geometry-aware)
    │   ├── solver_m2.js        # M2 + per-cell R_i
    │   ├── solver_m3.js        # M3 + extracellular R PDE
    │   └── solver_m4.js        # M4 + extracellular A PDE
    ├── agents.js               # per-cell state, SDE step; calls field.sample() only
    ├── nondim.js               # dim ↔ nondim linkage (two geometry branches)
    ├── render.js               # L heatmap + radial profile + emission coloring
    ├── ui.js                   # KNOBS, Calculate button, time scrub slider
    └── worker.js               # batch loop: accumulate → step field → step agents
```

Every `model.js` exports a pure stepping function
`step(state, params, dt) → state'` that takes nondimensional inputs and
returns nondimensional outputs. UI never calls into the model with
dimensional values; it converts first via the setup's `nondim.js`.

**JSDoc typedefs.** Each setup declares `@typedef` blocks for its core
objects (e.g. `AgentState`, `FieldState`, `SimParams`) at the top of
`agents.js` / `worker.js`. No TypeScript, but JSDoc types catch dim/nondim
confusion and wrong-field-passed-to-worker bugs at zero build cost.

---

## 3. The unit system (the part most likely to bite us)

### Per-setup nondimensionalization summary

| Setup | Length | Time            | Polarization | Cue        |
|-------|--------|-----------------|--------------|------------|
| 1     | —      | `w/u²`          | `√(u/w)`     | `L_c`      |
| 2     | `σ`    | `σ/c`           | `√(u/w)`     | `L_max`    |
| 3     | `σ`    | `σ/c`           | `√(u/w)`     | `L_max`    |
| 4 (2D–2D) | `D_L/c*`, `c* = √(aσD_L/(hL_0))` | `D_L/c*²` | `√(u/w)` | `L_0` |
| 4 (2D–3D) | `D_L/c*`, `c* = aσ/L_0`           | `D_L/c*²` | `√(u/w)` | `L_0` |

(Details and derivations in `docs/physics/setup{n}.md`.)

### Knob declaration

In each `setup{n}/ui.js`, every parameter is declared once:

```js
const KNOBS = [
  {
    id: 'r0', symbol: 'r₀',
    exposure: 'dim',                  // expose dimensionally only
    dim:    { default: 1.0, min: 0, max: 10, units: 's⁻¹ nM⁻¹' },
  },
  {
    id: 'chi_t', symbol: 'χ̃',
    exposure: 'nondim',               // dimensionless number, no dim form
    nondim: { default: 1.0, min: 0, max: 5 },
  },
  {
    id: 'L', symbol: 'L',
    exposure: 'both',                 // both knobs visible, linked live
    dim:    { default: 1.0, min: 0, max: 10, units: 'nM' },
    nondim: { symbol: 'L̃', default: 1.0, min: 0, max: 10 },
    toNondim: (dim, ctx) => dim / ctx.Lc,
    toDim:    (nd,  ctx) => nd  * ctx.Lc,
  },
];
```

The `units.js` helper turns this list into:
- live-linked sliders (changing the dim slider updates the nondim readout
  and vice versa, per the linkage functions),
- a single nondimensional `params` object passed to the simulation worker.

The user is the arbiter of which knobs are dimensional vs nondimensional;
flipping a knob's `exposure` is a one-line change.

### The "context" `ctx`

The linkage functions take a `ctx` of the *currently chosen* base scales
(e.g. `{Lc, u, w, sigma, c}`). When the user changes a base scale, all
dependent linked knobs recompute. This avoids a tangle of pairwise
conversions.

---

## 4. Setup-by-setup deliverables

### 4.1 Setup 1 — uniform-L Brownian polarization

**Sim.** SDE from `physics/setup1_uniform.md`. State = scalar `p̃` (1D
polarization is enough; if user later wants 2D, generalize to `p̃ = (p̃_x, p̃_y)`
with isotropic noise — model.js should be written so this is a one-line
change).

**Plots (live, updating as worker streams chunks):**
1. `F̃(p̃; ã)` curve with current `p̃` marked as a moving dot.
2. Time trace `p̃(t̃)` (uPlot, scrolling window).
3. Running histogram / KDE of `p̃` over a configurable window.
4. Bifurcation diagram of stable/unstable extrema of `F̃` vs. `ã`, with
   current `ã` marked.

**Knobs:** `ã` (or dim `r₀, L, L_c, ν`), `D̃` (or dim `θ`), `dt̃`, RNG seed,
window length, play speed.

**Done when:** user can sweep `ã` through the first-order transition and
visually see hysteresis-like behavior in steady-state histogram; can verify
detailed balance / Boltzmann distribution `P(p̃) ∝ exp(-F̃/D̃)` at small noise.

### 4.2 Setup 2 — Gaussian running wave

**Sim.** Coupled `(p̃, x̃)` SDE/ODE from `physics/setup2_wave.md` with a
linear velocity response `ṽ = μ̃ p̃`.

**Plots:**
1. 1D animation: cue strip (Canvas2D heatmap) with cell as a marker; trail
   of last `T_trail` positions.
2. Traces of `x̃(t̃)`, `p̃(t̃)`, and `∇L̃` sampled at the cell.
3. "Wave-passage" summary: `Δx̃` after one full wave passage as a function
   of one swept knob (typically `χ̃`).
4. Toggleable comparison overlay: pure-gradient model (no `p` dynamics)
   vs. full model — to make the polarization-induced rectification
   visible at a glance.

**Knobs:** `χ̃, α, β, b̃, D̃, μ̃`, plus dim equivalents per the
linkage.

**Done when:** the GL cell shows a finite net displacement `Δx̃` after
one wave passage while the gradient-only reference cell stays near
`Δx̃ ≈ 0`.

### 4.3 Setup 3 — 2D radial swarm with central trap

**Sim.** 2D extension of Setup 2 — same nondim groups (`M, C, λ, ϑ, χ̃, μ̃`),
polarization generalized to vector `P = (P_x, P_y)`, cue is a radial Gaussian
wave launched from the center and propagating outward (constant amplitude
`M`, no `1/r` decay). Hard reflective outer boundary at `r̃ = R̃_dish`; hard
trap at `r̃ < R̃_trap` (cells frozen on entry).

**Physical regime.** Defaults are chosen in the *anti-wave* regime of the
GL response (cell net displacement *against* wave-propagation direction),
so the outward wave drives cells inward toward the trap. See
[physics/setup3_swarm.md](physics/setup3_swarm.md).

**Plots:**
1. 2D animation: disk + trap disk + wavefront ring at `r̃ = C t̃`; cells
   colored by `|P|` (free = orange, trapped = green).
2. Mean free-cell radius `⟨r̃⟩_free(t̃)`.
3. Trapped fraction `n_trap(t̃) / N`.
4. Radial density `ρ(r̃)` at the current scrub time.

**Knobs.** Wave: `M, C` (+ dim `L_max, σ, c`). Coupling: `χ̃, μ̃`. Intrinsic:
`λ, ϑ`. Geometry: `N, R̃_dish, R̃_trap`. Numerics: `dt̃, seed, play speed`.

**Implementation.**
- Agent state in `Float32Array`s (`x, y, P_x, P_y, trapped_flag`), N=1000 default.
- Full trajectory (positions + trapped flag per step) precomputed
  deterministically per seed and stored, like Setup 2 — supports full
  time scrub. Memory ~32 MB at default `N = 1000, N_steps ≈ 1500`.
- Render: Canvas2D (sufficient at N ≤ 10⁴; reach for regl only if needed).
- Single `setup3/main.js` consolidating model + UI + plots.

**Done when:** with default parameters, after one wave passage the bulk of
free cells has migrated inward and a visible majority is trapped;
`⟨r̃⟩_free(t̃)` is monotonically decreasing during the wave passage.

### 4.4 Setup 4 — Emergent-wave swarming with dynamic cue field

The cue dynamics spec is now complete. See
[physics/setup4_swarm3d.md](physics/setup4_swarm3d.md) for the full
physics and architecture, and
[physics/setup4_cue_models.md](physics/setup4_cue_models.md) for the
model catalog. Summary:

**What's new relative to Setup 3.** The cue wave is emergent (L satisfies
its own PDE; wave amplitude and speed are outputs). Two geometry variants:
2D–2D and 2D–3D (latter solved on a direct 3D thin-slab grid — no
analytical reduction). Four cue models (M1–M4) plus optional GRK2
receptor desensitization (M5). All threshold functions are Hill-smoothed;
every Hill invocation carries its own independent exponent (n_L, n_R,
n_{Lr}, n_A, n_{LA}), default 10, with only those active in the selected
model shown in the UI. UI paradigm differs from Setups 1–3:
batch compute rather than real-time (see deliverables).

**Deliverables:**

- `setup4/solvers/field.js` — shared grid data structure and PIC
  operations. Owns the `Float32Array` for L (and R or A). Exports four
  functions only; nothing else is public:
  ```js
  sample(x, y)               // bilinear interp of L at agent position
  accumulateSource(x, y, w)  // smear w units of emission onto grid nodes
                             //   using the same bilinear weights as sample()
  step(dt, params)           // advance PDE one substep (model-specific)
  getRadialProfile()         // Float32Array: azimuthal average of L vs r̃
  ```
  `agents.js` and `worker.js` call only these functions — never the raw
  grid array.

- `setup4/solvers/index.js` — `createSolver(geometry, model)` factory.
  Returns the geometry+model-appropriate solver object. Each of
  `solver_m1.js` … `solver_m4.js` exports `makeStepFn(geometry, params)`
  which returns the `step()` closure bound to the right equations and
  grid layout. Model branching lives entirely in the factory; no
  if/else chains in `worker.js`.

- `setup4/agents.js` — per-cell state in `Float32Array`s
  (`x, y, Px, Py`, plus `Ri` or `Gi` for M2/M5). Updates `emitting_i`
  each step (see `setup4_swarm3d.md` §10). Calls `field.sample()` to
  read L; never touches the grid array. Smearing is done in `worker.js`
  via `field.accumulateSource()` before each PDE step (PIC consistency).

- `setup4/render.js` — Canvas2D: L heatmap behind cells; cells colored
  by `emitting_i` (warm/bright = emitting, cool/muted = not) with |P|
  modulating brightness; optional inhibitor panel; radial profile
  `𝓛(r̃, t̃)` linked to the shared time slider.

- `setup4/worker.js` — batch loop. Per time step:
  1. For each agent: `field.accumulateSource(x, y, emission * dt)`
  2. `field.step(dt_sub, params)` (with sub-stepping for 2D-2D)
  3. `agents.step(field, dt, params)`
  Posts frames on a fixed interval; posts progress updates every ~1%:
  ```js
  // every frame:
  { type: 'frame', step, t, radialProfile: Float32Array,
    agentX: Float32Array, agentY: Float32Array,
    emitting: Uint8Array }
  // every K frames (K ≈ 10–50, user-configurable):
  { type: 'frame', ..., Lfield: Float32Array }  // full 128×128 heatmap
  { type: 'progress', pct }
  { type: 'done' }
  ```

- `setup4/ui.js` — model toggle (M1/M2/M3/M4), geometry toggle, cell
  model toggle (GL/Hill), M5 add-on toggle; **Calculate** button; time
  scrub slider. Knobs per `setup4_swarm3d.md` §9. Hill exponents
  (n_L, n_R, n_{Lr}, n_A, n_{LA}) are per-function and independently
  adjustable; see §9 "Threshold regularisation" table for which are
  active per model (e.g. M3 shows n_L and n_{Lr} only). All default 10.

**Interaction paradigm: calculate-then-explore.** Setup 4 does not update
in real time. Workflow: set parameters → press "Calculate" → simulation
runs in a Web Worker saving frames → time slider and plots become active
for scrubbing. Any parameter change invalidates the cached run.

**Validation targets:**
- M1 in 2D–2D: measured wave speed approaches `c* = √(aσD_L/(hL_0))`
  (nondim value 1); exact only in the continuum / sharp-threshold limit.
- M1 in 2D–3D: measured wave speed approaches `(2/π) c*` (nondim ≈ 0.637);
  discrete-cell corrections are expected and not a solver error.
- M2/M3/M4: wave self-extinguishes; arrest radius visible on L heatmap
  and radial profile.
- M5 toggle: measured wave speed is approximately density-independent.
- Cell drift direction matches Setup 3 anti-wave behavior when χ̃ in the
  same parameter range.

---

## 5. Performance notes

- Setups 1–2: trivially real-time on any laptop. Worker is mostly there
  for clean separation, not raw speed.
- Setup 3: budget for `N ≈ 10⁴` agents at 30 fps. Keep all per-agent state
  in `Float32Array`s; SDE step is one `for` loop over the array. Renderer
  uses instanced quads.
- Setup 4: PDE loop (128×128 grid × many substeps) is the bottleneck.
  Runs in a Web Worker as a batch job; frames saved at fixed intervals and
  replayed via time slider (no real-time rendering). 2D–3D uses
  Crank-Nicolson with pre-factorized LU to avoid the prohibitively small
  CFL from the fine z-grid. No WebGL required for the batch path.
  **Frame storage.** Every frame carries: radial profile (128×4 B = 512 B),
  agent state (4 arrays × N×4 B ≈ 8 KB at N=500), emitting flags (N B), and a
  **Uint8-compressed L heatmap** (128×128 × 1 B = 16 KB, per-frame
  max-normalized to 0–255 in the worker; absolute peak sent separately as a
  scalar). The renderer's viridis LUT only needs the 0–255 ratio, so Uint8 is
  visually indistinguishable from Float32. At 1000 frames: ~1000×(9 + 16) KB
  ≈ 25 MB — same budget as the original K-frame scheme, but the heatmap
  animates smoothly every frame instead of every K-th.

---

## 6. Validation checklist

Before declaring any subpage "done":

- [ ] Nondim equations on the page match `docs/physics/setup{n}.md` exactly.
- [ ] Every visible knob is in the `KNOBS` config; no hardcoded magic
      numbers in `model.js`.
- [ ] Reducing `dt̃` by 2× changes results by less than ~1% on a chosen
      benchmark trace.
- [ ] At zero noise, deterministic limit reproduces a known fixed
      point / trajectory analytically derivable from the equations.
- [ ] At zero noise *and* zero `χ`, Setup 2 cell stays still.
- [ ] At large noise / small barrier, Setup 1 histogram matches
      `exp(-F̃/D̃)` to eye.

---

## 7. Staging / order of work

1. **Scaffolding & shared core.** `shared/{rng,sde,units,state,widgets}.js`,
   page shell, nav. No physics yet.
2. **Setup 1.** Smallest physics surface, tests the unit system end-to-end.
3. **Setup 2.** Reuses Setup 1's SDE; adds wave + position + response maps.
4. **Setup 3.** 2D-radial swarm — reuses Setup 2's nondim groups, adds
   agent buffer + Canvas2D renderer + hard trap + radial wave.
5. **Setup 4 scaffold.** Renderer + agent loop + stub `L`. Stop here and
   wait for the cue-dynamics spec from the user.
6. **Setup 4 full.** Add real `L` dynamics once specified; tune perf.

Do not start step `n+1` until step `n` is checked against §6.

---

## 8. Decisions (resolved 2026-05-05)

- **Setup 1:** support both 1D and 2D polarization, switchable in the UI.
- **Setup 2:** single passing Gaussian only; no wave-train.
- **Setup 3 (2026-05-07):** carved out the 2D-radial-swarm subproblem from
  the old Setup 3 placeholder and made it the new Setup 3. Same nondim
  groups as Setup 2; outward Gaussian wave from center; cells in the
  anti-wave regime drift inward and are caught by a hard central trap.
  The original 3D-half-space-cue placeholder is now Setup 4.
- **Setup 4 (2026-05-12, updated 2026-05-13):** cue dynamics fully
  specified. See `docs/physics/setup4_cue_models.md` (model catalog) and
  `docs/physics/setup4_swarm3d.md` (physics + architecture). Key
  decisions: length scale `ℓ_0 = D_L/c*` where `c*` is the parametric M1
  wave speed (geometry-dependent: `c* = √(aσD_L/(hL_0))` for 2D–2D,
  `c* = aσ/L_0` for 2D–3D); four cue model variants (M1–M4) plus M5 GRK2
  add-on; 2D–3D solved on a direct 3D thin-slab grid with Crank-Nicolson —
  no Green's function approximation. Absorbing BC (𝓛=0) at dish edge;
  seed-disk IC. Batch-compute paradigm (no real-time parameter updates).
  Cell SDE and λ, ϑ, χ̃, μ̃ definitions unchanged from Setups 2–3.
- **Setup 4 architecture (2026-05-13):** solver factory pattern —
  `createSolver(geometry, model)` in `setup4/solvers/index.js` dispatches
  to model-specific stepper modules. `field.js` owns the grid and exposes
  only `sample / accumulateSource / step / getRadialProfile`; agents never
  import the raw array. Worker loop order: accumulate emissions → step PDE
  → step agents (enforces PIC consistency).
- **Setup 4 Hill exponents (2026-05-13):** each Hill function invocation
  has its own independent exponent n (n_L, n_R, n_{Lr}, n_A, n_{LA}) — no
  global n knob. All default to 10. Only the exponents active in the
  selected model are shown in the UI.
- **Setup 4 frame storage (2026-05-13):** Option A — radial profile + agent
  state saved every frame; full 128×128 L field saved every K≈10–50 frames.
  Total budget ~22 MB at 1000 frames, K=20. Heatmap scrubbing is at K-frame
  resolution; radial profile and cell animation are always full-resolution.
- **Hill alternative model (2026-05-07):** Setups 2 and 3 each expose a
  per-setup toggle between the default Ginzburg–Landau polarization SDE
  and a Hill / linear-relaxation alternative
  (`dP/dt̃ = χ̃ ∂_x̃ 𝓛 − P + √(2ϑ𝓛) η`, `dx̃/dt̃ = ṽ₀ |P|ⁿ/(1+|P|ⁿ) P̂`).
  Noise term is intentionally identical to GL (multiplicative `√(2ϑ𝓛)`).
  `χ̃, ϑ` sliders are shared between models; `μ̃, λ` are GL-only; `ṽ₀, n`
  Hill-only. F̃ panel hides in Hill mode (Setup 2). See
  [physics/hill_model.md](physics/hill_model.md).
- **Plots are always nondimensional.** No physical units on any axis,
  tick, label, tooltip, or animation coordinate. Dimensional sliders are
  *inputs* only — moving them recomputes the linked nondim parameter,
  and the sim+plots only see nondim values. (Recorded in `CLAUDE.md`.)
- **Dim ↔ nondim linkage is one-way semantically:** the nondim value is
  what the simulation uses; the dim slider is a convenient way to scale
  it. Doubling a dim knob doubles every nondim quantity proportional to
  it (via the linkage function).

## 9. Pragmatic file consolidation

The "one file per concern" layout in §2 is aspirational. For v1, each
setup may consolidate `model.js + nondim.js + ui.js + plots.js` into a
single `main.js` to ship faster. Refactor when a setup outgrows ~600
lines or when a second developer joins (neither applies yet).
