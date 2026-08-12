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
analytical reduction). Cue models M1–M4 and M6.1, plus optional GRK2
receptor desensitization (M5). M6.1 (basal adenosine density sensor) is the
model that reproduces the experimentally observed density-independent wave
speed and suppresses channelisation; see `setup4_cue_models.md` §7. All
threshold functions are Hill-smoothed;
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
  `solver_m1.js` … `solver_m4.js`, `solver_m6_1.js` exports `makeStepFn(geometry, params)`
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

- `setup4/ui.js` — model toggle (M1/M2/M3/M4/M6.1/M6.2), geometry toggle, cell
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
- **No analytical `c(σ̃)` target for M6.1/M6.2.** The mean-field pushed-front
  flux balance is not a reliable description of this system; it is no longer
  quoted anywhere (pages, plots, KPIs, catalog §7/§7b), and the sweep plots
  measured points only. Validate against *internal* consistency instead:
- M6.1: `λ_A = 0` reproduces M1 exactly; measured `𝓐(r̃)` plateau matches the
  basal tone `σ̃/Γ_A` (2D–2D) / `σ̃/√(DΓ_A)` at z=0 (2D–3D) — an exact steady
  state of the linear 𝓐 equation; raising `λ_A` bends the measured `c(σ̃)`
  away from the `λ_A = 0` sweep (compare sweeps, not formulas). Channelisation
  present at `λ_A = 0` should weaken for `λ_A > 0` with `ℓ_A = √(D/Γ_A)` of
  order the relay length.
- M6.2: `β = 0` reproduces M1 exactly; measured `𝓠(r̃)` plateau matches
  `βσ̃/γ`; the `m = 1` sweep behaves like an M6.1 sweep at `λ_A = β, Γ_A = γ`;
  the `m = 2` sweep is compared against the `m = 1` and `β = 0` sweeps for
  its shape (whether it turns over is a measurement, not a target). Full
  checklist in
  [physics/setup4_m6_2_implementation_plan.md](physics/setup4_m6_2_implementation_plan.md) §9.
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
  `c* = aσ/L_0` for 2D–3D); cue model variants M1–M4 and M6.1 plus M5 GRK2
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
- **Setup 4 M2 implementation (2026-05-26):** added per-cell intracellular
  inhibitor R̃_i alongside the existing M1. Nondim follows
  [physics/setup4_swarm3d.md](physics/setup4_swarm3d.md) §5.1/§5.3 (β̃ =
  β t_0, γ̃ = γ t_0, L̃_r = L_r/L_0, R̃ = R/R_c so R_c ≡ 1 in nondim).
  PDE structure unchanged from M1 — only the per-cell source weight gains
  an `H^-(R̃_i;1;n_R)` factor, and a new per-cell Euler-integrated ODE
  `dR̃_i/dt̃ = β̃ H^+(𝓛;L̃_r;n_{Lr}) − γ̃ R̃_i` runs alongside the SDE.
  Solver factory delegates M2 to the same M1 step functions (the gating
  is applied in the worker's emission accumulation loop). M2 uses a
  **different initiation scheme** than M1: instead of a time-limited
  bulk firing source, cells whose current position lies inside r̃ < r̃_fire
  are *force-emitted* indefinitely (their `H^+(𝓛;1;n_L)` gate is set to
  1; `H^-(R̃;1;n_R)` still applies). The wave terminates organically
  once central cells accumulate R̃_i > 1 — effective firing duration is
  set by γ̃. Cells deactivated by R̃ share the same gray color as
  never-activated cells; no third color.
- **Setup 4 M6.1 basal-adenosine model (2026-07-07):** added to the catalog
  ([physics/setup4_cue_models.md](physics/setup4_cue_models.md) §7) to
  reproduce the experimentally observed **density-independent wave speed**
  and suppress channelisation. Distinct from M4: adenosine A is produced
  *basally* (constant per-cell rate, no L-gate) with first-order decay γ_A,
  and couples by *shifting the relay threshold* `L_0 → L_0 + gA` rather than
  gating production. Nondim (intrinsic scales ℓ_0 = a/(L_0 D_L),
  A_0 = b/(hD_L)): `∂_t̃𝓛 = ∇̃²𝓛 + α Σδ̃ H^+(𝓛;1+λ𝓐;n_L)`,
  `∂_t̃𝓐 = D∇̃²𝓐 + Σδ̃ − Γ_A𝓐`, with knobs α, λ, D, Γ_A. `λ=0` ≡ M1.
  Theory (pushed/ignition-front flux balance `c²𝓛_c = ασ̃`): basal tone
  `𝓐_ss = σ̃/Γ_A` makes the threshold `𝓛_c = 1+λσ̃/Γ_A` scale *with* density,
  so `c(σ̃) = √(ασ̃/(1+λσ̃/Γ_A)) → √(αΓ_A/λ)` (density-independent) above
  `σ̃★ = Γ_A/λ`. A triggered inhibitor (M4) instead pins 𝓛_c independent of
  σ̃ and keeps `c ∝ √σ̃`. Channel suppression additionally needs the sensing
  length `ℓ_A = √(D/Γ_A)` of order the relay length.
- **Setup 4 M6.1 implementation (2026-07-08):** implemented `solver_m6_1.js`
  (2D–2D only, per catalog §10). M6.1 reuses the M1 2D–2D L stepper verbatim and
  adds a second `createField` instance for the adenosine field 𝓐 with its own
  explicit-Euler diffusion stepper (coefficient `D = D_A/D_L`, decay `Γ_A`,
  CFL `dt_sub ≤ Δx̃²/(4D)`). The solver exposes `solver.fieldA`. Worker wiring:
  (1) per-cell emission gate uses the shifted threshold
  `H⁺(𝓛_i; 1 + λ_A·𝓐_i; n_L)` (𝓐_i sampled from `fieldA`); (2) every cell
  accumulates a basal 𝓐 source `field.accumulateSource(x,y,dt)` (coefficient 1,
  no 1/h̃ — the nondim A-source coefficient is exactly 1); (3) M6.1 shares M1's
  time-limited firing-source nucleation. `emitting` flag = `𝓛 > 1 + λ_A·𝓐`
  (computed in `stepAgents`, which now takes an optional `fieldA` arg). UI: new
  M6.1 cue-model toggle option + nondim sliders λ_A (linear, 0 ≡ M1), D, Γ_A;
  the R-panel is repurposed to show 𝓐(r̃) with a dashed guide at the mean-field
  tone 𝓐_ss = σ̃/Γ_A. Symbol note: the catalog's λ is displayed as λ_A in the
  GUI to avoid a clash with the cell-side GL group λ. New page `setup4-m6.1/`
  (4c), a thin shell reusing `../setup4/ui.js` with `{ model: 'M6.1' }`.
  2D–3D M6.1 throws a descriptive error in the factory.
- **Setup 4 c(σ̃) density sweep (2026-07-08):** added a multi-run diagnostic
  for M6.1 (its own **Sweep c(σ̃)** button + dedicated Web Worker, separate from
  the single-run *Calculate*). Extracted the per-step physics into a pure,
  node-testable module `setup4/sim_core.js` (`createSim` / `advance` /
  `measureWaveSpeed`) used by BOTH the single run and the sweep, so they can't
  diverge. `measureWaveSpeed` runs a frame-less sim and fits the relay front
  speed from the outermost 𝓛 = 1 crossing over a mid-time window (past the
  firing plateau, before the wall). The sweep loops σ̃ ∈ [0.1, 8] (10 log
  points), recomputing N = round(σ̃·π·R̃²) per point, on a capped geometry
  (N_grid=80, R̃=12, t̃_max=25, robust firing) for speed; it streams points as
  each finishes. New `drawCsweep` (render.js) plots measured points on a log-σ̃
  axis against the mean-field pushed-front theory `c=√(ασ̃/(1+λ_A σ̃/Γ_A))`
  (α=1/h̃), the λ_A=0 reference `√(ασ̃)`, the saturation asymptote
  `√(αΓ_A/λ_A)`, and the crossover σ̃★=Γ_A/λ_A — all live with the current
  λ_A/Γ_A. Verified in node: measured c tracks √σ̃ for λ_A=0 and flattens
  toward the asymptote for λ_A>0.
- **Setup 4 M6.1 in 2D–3D (2026-07-09):** extended M6.1 to the 2D–3D geometry
  (previously 2D–2D only). Generalized `solver_m1.js`'s 2D–3D stepper into a
  reusable `makeSlabStepper(field2d, {N_z,h_0,alpha, getD, getDecay})`
  reaction–diffusion slab solver (`makeStepFn_2d3d` is now a thin wrapper with
  D=1, decay=Γ_L — M1's path is unchanged). `solver_m6_1.js` gains
  `makeStepFn_m6_1_2d3d`, which drives 𝓛 (D=1, Γ_L) and a full 3D adenosine field
  𝓐 (D=D_A/D_L, Γ_A) with two slab steppers; cells read 𝓐|_{z=0}. Factory
  routes M6.1+2D–3D accordingly (no longer throws). The z=0 basal tone is the
  screened surface source 𝓐_ss = σ̃/√(D Γ_A) (verified in node: 1.49 vs 1.41 at
  σ̃=1), so the 2D–3D speed law is c = (2/π)σ̃/(1+λ_A σ̃/√(D Γ_A)) → (2/π)√(D Γ_A)/λ_A
  (vs 2D–2D √(ασ̃/(1+λ_A σ̃/Γ_A))). The c(σ̃) sweep now (1) uses the SELECTED
  geometry — it was hard-coding 2D–2D — and (2) picks geometry-appropriate
  settings (2D–3D: smaller grid, stronger/briefer firing, longer window, fewer
  points, since 3D dilutes the surface stimulus and relay speeds are lower).
  `drawCsweep` switches its theory/asymptote/crossover formulas by geometry.
  Verified: 2D–3D λ_A=0 gives c∝σ̃ (linear, distinct from 2D–2D √σ̃), λ_A>0
  flattens toward (2/π)√(DΓ_A)/λ_A; low-σ̃/high-λ_A points legitimately fail to
  ignite in 3D. M6.1 page + catalog §10 updated.
- **Setup 4 sweep false-positive fix (2026-07-09):** the sweep's
  `measureWaveSpeed` originally tracked the 𝓛 = 1 contour, which — with the
  strong ignition kick and Γ_L = 0 — a purely *diffusive* spread of the fired
  blob inflates as √t, misreported by the linear fit as a finite wave speed
  (finite c at densities where the single run shows only diffusion). Rewrote it
  to track the **emission front** (90th-percentile radius of relaying cells,
  which only advances where cells cross their threshold) with two gates —
  *spread* (front clears 1.5·r_fire) and *sustained* (front hasn't receded by
  the end) — so a diffusing halo of the ignition source reads as c = 0. Speed
  is a robust two-point average over the propagation phase (handles both slow
  fronts and near-instant fills, no t_fire gate). Also fixed the single-run
  `c_eff` KPI (ui.js), which tracked the profile *peak* — pinned at r≈0 for the
  Γ_L = 0 filling front, so it read ≈0 regardless of α/σ̃ — to a front-crossing
  least-squares fit. Verified: no wave for low σ̃ at λ_A=5/Γ_A=1.72 (matches
  the ABM), √σ̃ growth preserved for λ_A=0, flattening for λ_A>0.
- **Setup 4 M6.2 defaults (2026-07-28):** β = 15, D = 7, γ = 5, m = 2,
  t̃_max = 500 (user-specified). Tone 𝓠_ss = βσ̃/γ = 3 at the default σ̃ ≈ 1, so
  with m = 2 production runs at α_eff/α = 0.1 — the throttle is firmly on —
  and ℓ_Q = √(D/γ) ≈ 1.2 is of order the relay length, the condition for 𝓠 to
  resolve local density bumps (catalog §7b). β's slider range was raised to
  0–50 to reach 15. β/D/γ/m are declared once in `params`; the page-level
  t̃_max override goes through a new `MODEL_DEFAULTS` table applied in
  `buildUI()` before the sliders are built, and writes the dim canonical
  `dim.t_max_dim` too (t̃_max is dim-canonical — setting only `params.t_max`
  would be silently reverted by `recomputeFromDim()`, breaking display ≡
  simulated). One Calculate at these defaults costs ≈ 60 s and ≈ 50 MB of
  cached frames; the throttled front needs the long window (at t̃ = 20 it has
  barely cleared the ignition disk).
- **Setup 4 channelisation order parameter Ψ (2026-07-28):** added a
  quantitative measure of the Höfer–Maini streaming instability, since
  "M6 suppresses channelisation" was previously only assessable by eye.
  Streaming = angular structure in the cell positions, so Ψ is built from the
  azimuthal Fourier modes of the cell angles (cells outside a 0.1·R̃_dish core):
  `c_m = (1/N)Σe^{imθ}`, `Ψ_m = √(max(0,(N|c_m|²−1)/(N−1)))` ∈ [0,1]. The −1
  removes the shot-noise floor `E[|c_m|²]=1/N` that a perfectly uniform swarm
  produces — without it the metric would drift with density and be useless for
  exactly the density sweeps this setup is about. Two panels: `Ψ(t̃) =
  max_{m≥2} Ψ_m` (one number for "how channelised"), and the spectrum Ψ_m at
  the scrubbed time with the 5% shot-noise level √(2/(N−1)) dashed and the
  dominant mode m* highlighted — m* IS the channel count. m=1 is excluded
  (bulk off-centre drift, not spokes) and drawn grey. Computed in the frame
  handler (O(N·m_max), no worker change); new `drawAngularSpectrum` in
  render.js. Panels added to the M1 (4a), M6.1 (4c) and M6.2 (4d) pages —
  the diagnostic is model-independent and the comparison that matters is
  across pages at equal σ̃, with M1 as the channelising baseline. Definition
  documented in physics/setup4_swarm3d.md §10. Verified in node at κ=0.4,
  μ=0.05, σ̃=1: Ψ starts at the noise level (0.047 vs 0.040) and grows to
  0.24 (M1) vs 0.17 (M6.2, β=0.2, m=2) by t̃=60.
- **Setup 4 M6.2 second petri dish (2026-07-28):** the M6.2 page now shows two
  dishes side by side — the same cells over 𝓛 (magma) and over 𝓠 (greens, a
  deliberately different hue family so the panels are never confused). The
  worker ships a Uint8-compressed `Qfield` (M6.2 only — an extra N_grid²
  bytes/frame is waste for models nothing renders it for), normalized by its
  own GRID max rather than the radial-profile max used for 𝓛: 𝓠 is nearly
  uniform, so radial normalization would flatten it to a featureless wash,
  whereas the grid max keeps per-cell hot spots and density bumps visible.
  `drawDish` gained `fieldKey` / `palette` / `haloScale` params (one routine,
  two panels, so cell markers are pixel-identical and directly comparable);
  the emission halo is damped to 0.45 on the 𝓠 panel so it does not swamp the
  field. The ⟨|ṽ_r|⟩ trace moved to its own full-width row below the dishes.
- **Setup 4: analytical `c(σ̃)` withdrawn (2026-07-28):** the mean-field
  pushed-front flux balance (`c² 𝓛_c = α_eff σ̃`, and everything derived from
  it — `√(ασ̃/(1+λ_Aσ̃/Γ_A))`, the saturation asymptote, `σ̃★`, `σ̃_peak`, and
  the 2D–3D counterparts) does **not** describe this system reliably and has
  been removed from the M6.1 and M6.2 pages, from the `c(σ̃)` sweep plot
  (`drawCsweep` now draws measured points only — no theory curve, no
  asymptote, no crossover/peak markers), from the KPI panel (`σ̃*` and
  `σ̃_peak` KPIs deleted; `tone_ss` and `ℓ_screen` stay — they are exact
  properties of the *linear* auxiliary field and were checked numerically),
  and from catalog §7/§7b, which now state the mechanism (`𝓛_c = 1 + λ_A𝓐_ss`
  for M6.1, `α_eff = α/(1+𝓠_ss^m)` for M6.2 — both exact) and stop there.
  Wave speed is treated as an OUTPUT: compare sweeps (β=0 vs β>0, m=1 vs
  m=2) rather than a sweep against a formula. Speed formulas surviving in
  older entries below are superseded by this one. Unaffected: the Dieterle
  M1 profile overlay on the radial-profile panel (a literature result for a
  different quantity) and §1's M1 speed expressions.
- **Setup 4 M6.2 implemented (2026-07-28):** the model of the previous
  entry is live as page `setup4-m6.2/` (4d). Solver: `solver_m6_1.js` was
  generalized into `solver_auxfield.js` — "M1 relay + one basally-produced
  auxiliary reaction–diffusion field" — which now serves BOTH M6.1 and M6.2
  in both geometries; the auxiliary PDE is identical (coefficients read from
  `params.D_aux`/`params.gamma_aux`), and only the per-cell L source weight
  differs, computed in `sim_core.advance`: M6.1 `H⁺(𝓛;1+λ_A𝓐;n_L)`, M6.2
  `H⁺(𝓛;1;n_L)·H⁻(𝓠;1;m)` (`hillNeg(Q,1,m_Q)` — the throttle IS a negative
  Hill). `solver.fieldA` → `solver.fieldAux` throughout. New knobs β, D, γ, m
  plus a `𝓠` initial-condition toggle (steady tone / zero — seeding avoids a
  silent t̃~1/γ transient in which the run measures unthrottled M1 physics);
  new KPIs tone_ss, σ̃★, ℓ_screen, σ̃_peak (also filled in for M6.1). `emitting`
  for M6.2 is the L gate ALONE, deliberately breaking the `w>½` truth-table
  convention — with the convention, `𝓠_ss>1` above σ̃★ would flag every cell
  silent and zero out the measured front; emission *strength* is shown as dot
  and halo opacity instead. Verified in node: β=0 reproduces M1 bit-for-bit;
  seeded 𝓠_ss=βσ̃/γ is an exact steady state; the m=1 sweep tracks an M6.1
  sweep at λ_A=β (mean-field degeneracy) until M6.1 stops igniting while M6.2
  keeps propagating; m=2 at β=0.2, γ=0.5 gives a non-monotonic c(σ̃) peaking
  near σ̃≈2 (theory 2.5). See
  [physics/setup4_m6_2_implementation_plan.md](physics/setup4_m6_2_implementation_plan.md).
- **Setup 4 sweep ignition-halo fix (2026-07-28):** `measureWaveSpeed` was
  timing the *stimulus*, not the relay. The firing source injects a 𝓛 mass
  M = s_fire·σ̃·πr̃_fire²·t̃_fire which, at Γ_L=0, is conserved — so a purely
  diffusive halo holds ⟨𝓛⟩>1 out to r̃_halo = r̃_fire√(s_fire σ̃ t̃_fire),
  pushing cells over threshold with no relay at all. Control: M6.2 with
  β=10⁴ (relay impossible) reported c̃≈0.6. Fixes: (1) timing starts outside
  the halo (r_start = max(1.5r̃_fire, 1.15r̃_halo)), with a distinct
  `unmeasurable` status when the halo would flood the dish; (2) the sweep is
  mass-targeted, s_fire = K/(σ̃·t̃_fire) (`fire_K`), so r̃_halo = r̃_fire√K is
  density-independent instead of growing ∝√σ̃ (2D–3D uses the half-space bound
  r̃_halo = (π s_fire σ̃ r̃_fire² t̃_fire)^{1/3}); (3) the timing arc is gated on
  radial span, not sample count, so a near-zero baseline can't read as a
  spuriously fast wave while a genuinely fast front still measures; (4) sweep
  geometry resized (2D–2D R̃=16/N_grid=104/t̃_max=40/K=20; 2D–3D R̃=12/K=20 with
  σ̃ capped at 3, since 3D speeds are linear in σ̃ and fill the dish); (5) points carry a status
  (ok / noignite / unresolved / unmeasurable) and are drawn distinctly, so a
  slow-but-real front is no longer plotted as a zero. This also corrects the
  M6.1 sweep, whose old settings gave r̃_halo≈12.7 in an R̃=12 dish: its
  high-σ̃ points were inflated and its high-λ_A points now honestly read
  "no ignition".
- **Setup 4 M6 → M6.1 rename + M6.2 specified (2026-07-28):** the basal-adenosine
  model is now **M6.1** everywhere (model id `'M6.1'`, `solver_m6_1.js`,
  page `setup4-m6.1/`, catalog §7) to make room for a second quorum-sensing
  variant. **M6.2 (quorum-throttled production)** is specified in catalog
  §7b: the same basal, diffusing, decaying density signal `𝓠`, but coupled
  to the *production rate* rather than the threshold —
  `∂_t̃𝓛 = ∇̃²𝓛 + α Σδ̃ H⁺(𝓛;1;n_L)H⁻(𝓠;1;m)`,
  `∂_t̃𝓠 = D∇̃²𝓠 + β Σδ̃ − γ𝓠`, with `β = b/(hD_LQ_0)`, `γ = γ_Q t_0`,
  `D = D_Q/D_L`, and a new independent Hill exponent `m`. Mean field:
  `𝓠_ss = βσ̃/γ` is uniform, so M6.2 ≡ M1 with `α → α/(1+𝓠_ss^m)` and
  `c(σ̃) = √(ασ̃/(1+(βσ̃/γ)^m))` — flat at `m=1` (degenerate with M6.1 under
  `λ_A↔β`), *decreasing* for `m>1` (peak at `σ̃_peak=(γ/β)(m−1)^{−1/m}`,
  i.e. `dc/dσ̃<0`: anti-streaming with the opposite sign, not just zero
  slope). `β=0 ≡ M1` (`m=0` gives a factor ½, not 1). Implementation plan,
  including the two traps — the `emitting` flag must **not** use the
  `w>½` truth-table convention here (it would zero the swept `c̃` above
  `σ̃★`), and `β` needs an explicit `h̃` factor in 2D–3D because `Q_0` is
  physical and cannot absorb `b/h` the way M6.1's `A_0` does — lives in
  [physics/setup4_m6_2_implementation_plan.md](physics/setup4_m6_2_implementation_plan.md).
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
- **Setup 4 target sticking (2026-07-28):** every full-swarm page (M1, M2,
  M6.1, M6.2 — they share `setup4/ui.js` + `worker.js`) gained a
  *Stick to the target* toggle (ON by default) and a nondim
  `R̃_target` slider (default 2): a circle at the dish centre standing for the
  pathogen the swarm converges on. It is an **absorbing** inner boundary — the
  first crossing projects the cell onto the circle and freezes its position
  permanently (`agents.stuck`), while its polarization SDE and its cue emission
  keep running, so no cue-model equation changes. The target disk is excluded
  from the initial cell placement (N unchanged, so σ̃ = N/(πR̃²_dish) still names
  the same knob). The sweep passes the same knobs, so `c̃(σ̃)` describes the same
  swarm as the single run. UI: teal dashed circle + teal ring on engaged cells
  in both dish panels, `on target` KPI = recruitment curve. Spec:
  [physics/setup4_swarm3d.md](physics/setup4_swarm3d.md) §2c.
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
