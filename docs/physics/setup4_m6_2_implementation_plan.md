# Setup 4 — M6.2 (quorum-throttled production): implementation plan

Physics: [setup4_cue_models.md §7b](setup4_cue_models.md). Staging:
[../PLAN.md](../PLAN.md). This document is the *how*; the catalog is the
*what*. **Status: implemented** (2026-07-28) — steps 0–9 done, with two
deviations from the original plan recorded in §12.

## 0. Rename M6 → M6.1 (DONE)

| before | after |
|--------|-------|
| model id string `'M6'` | `'M6.1'` |
| `setup4/solvers/solver_m6.js` | `setup4/solvers/solver_m6_1.js` |
| `makeStepFn_m6`, `makeStepFn_m6_2d3d` | `makeStepFn_m6_1`, `makeStepFn_m6_1_2d3d` |
| page dir `setup4-m6/` | `setup4-m6.1/` |
| nav label `4c · M6 basal adenosine` | `4c · M6.1 basal adenosine` |
| catalog §7 heading, PLAN entries | `M6.1` throughout |

Model ids are compared as strings only (`params.model === 'M6.1'`), so the
dot is harmless. `lam_A`, `D_A_nd`, `gamma_A` keep their names (they are
M6.1-specific adenosine symbols); the solver's `fieldA` became the
model-agnostic `fieldAux` in step 2.

## 1. Nondim contract (what the code must implement)

```
∂_t̃ 𝓛 = ∇̃²𝓛 + α Σ_i δ̃²(r̃−r̃_i) H⁺(𝓛_i;1;n_L) H⁻(𝓠_i;1;m) − Γ_L 𝓛
∂_t̃ 𝓠 = D ∇̃²𝓠 + β Σ_i δ̃²(r̃−r̃_i) − γ 𝓠
```

with `α = 1/h̃`, `β = b/(h D_L Q_0)`, `D = D_Q/D_L`, `γ = γ_Q t_0`, and
`H⁻(𝓠;1;m) = 1/(1+𝓠^m)` **exactly** — reuse `hillNeg(Q, 1, m_Q)` from
`agents.js`, no new math primitive.

New nondim knobs (all nondim-primary; add to `KNOBS`/`params` in `ui.js`):

| code key | symbol | default | range | note |
|----------|--------|---------|-------|------|
| `beta_Q`  | `β` | 15 | 0–50 lin | `0 ≡ M1` (the off dial — **not** `m=0`) |
| `D_Q_nd`  | `D` | 7 | 0.05–50 log | sets `ℓ_Q = √(D/γ)` |
| `gamma_Q` | `γ` | 5 | 1e-3–10 log | tone reaches 1 at σ̃ = γ/β = 1/3 |
| `m_Q`     | `m` | 2 | 0.5–8 lin, step 0.5 | throttle sharpness; `m>1` makes the drive fall with density |
| `q_ic_ss` | — | on | toggle | pre-seed `𝓠` at its steady tone (see §5) |

Defaults (set 2026-07-28) put the throttle firmly ON at the default single-run
`σ̃ ≈ 1`: tone `𝓠_ss = βσ̃/γ = 3`, so with `m = 2` production runs at
`α_eff/α = 1/(1+3²) = 0.1`. The screening length `ℓ_Q = √(D/γ) = √(7/5) ≈ 1.2`
is of order the relay length, which is the condition for 𝓠 to resolve *local*
density bumps rather than act as a global rescaling of α (catalog §7b,
anti-streaming). The throttle turn-on σ̃ = γ/β = 1/3 sits below the whole sweep
range, i.e. every swept density is in the throttled regime.

β, D, γ, m are declared once in `ui.js`'s `params`. The M6.2 *page* additionally
overrides `t_max = 500` via a `MODEL_DEFAULTS` table applied in `buildUI()`
before the sliders are built (t̃_max is a shared numerics knob, and the other
Setup-4 pages keep 50). A tenfold-throttled front is slow — at t̃ = 20 it has
only just left the ignition disk — so the long window is needed to see it
propagate. Cost at these defaults: ≈ 60 s of compute and ≈ 50 MB of cached
frames for one Calculate.

Note `t_max` is **dim-canonical**: the override must also write
`dim.t_max_dim = t_max·t₀`, or `recomputeFromDim()` overwrites it on the first
sync and the slider would show a value the simulation never used.

## 2. Solver (small change — the PDE is M6.1's)

M6.1's solver is already "M1 relay stepper + one auxiliary
reaction–diffusion field with `(D, decay)` and a PIC source". M6.2 needs
*exactly* that; the model difference lives entirely in the per-cell source
weight computed upstream in `sim_core.advance()`.

Plan: rename `solver_m6_1.js` → `solver_auxfield.js`, exporting
`makeStepFn_auxfield(field2d)` and `makeStepFn_auxfield_2d3d(field2d, p3d)`
(the current bodies, verbatim, with `𝓐` → "aux field" in comments), and let
`solvers/index.js` route **both** `'M6.1'` and `'M6.2'` to it. Rename the
returned `fieldA` → `fieldAux` (≈10 call sites across `index.js`,
`worker.js`, `sim_core.js`, `agents.js`) so the second model does not read
its quorum field through a variable named "adenosine".

Reject nothing: both geometries already work (`makeSlabStepper` handles the
2D–3D aux field with `getD`/`getDecay`).

## 3. Per-cell coupling (`sim_core.advance`)

```js
} else if (model === 'M6.2') {
  const Q_i = solver.fieldAux ? Math.max(0, solver.fieldAux.sample(xi, yi)) : 0;
  w = hillPos(L_i, 1, n_L) * hillNeg(Q_i, 1, m_Q);   // threshold UNSHIFTED
  // Basal quorum source. β is the 2D–2D group (contains 1/h); the 2D–3D
  // surface-source coefficient is β·h̃ — same h̃ bookkeeping as the L source,
  // which uses α = 1/h̃ in 2D–2D and 1 in 2D–3D.
  solver.fieldAux.accumulateSource(xi, yi, beta_Q * (is2d2d ? 1 : h_tilde) * dt);
}
```

and the firing-source branch gains `|| model === 'M6.2'` (nucleation is the
same time-limited disk as M1/M6.1).

**Do not throttle the firing source by `𝓠`.** It represents an external
stimulus (micropipette/uncaging), not cell production; throttling it would
make high-density runs fail to ignite for an artefactual reason and put
spurious zeros in the `c(σ̃)` sweep.

## 4. `emitting` flag — the one real trap

`stepAgents` sets `emitting[i]` and `measureWaveSpeed` derives the front
radius (and hence every swept `c̃`) from it. The repo's truth-table
convention is "source weight `w > 1/2`", which for M6.2 reads
`𝓛 > 1 AND 𝓠 < 1`. That convention **breaks the diagnostic here**: above
`σ̃ = γ/β` the steady tone `𝓠_ss > 1` everywhere, so every cell would be flagged
non-emitting and the sweep would report `c = 0` at exactly the densities
M6.2 exists to describe — a plausible-looking cliff that is pure
instrumentation.

Decision: for M6.2, `emitting[i] = (L_i > 1)` (the *relay* gate only), and
encode the throttle factor `H⁻(𝓠_i;1;m)` as the **opacity/brightness** of
the emitting marker in `render.js`. The physics ("still firing, but weakly")
is then visible without corrupting the front metric. Document the departure
from the `w > 1/2` convention in `agents.js` next to the branch.

## 5. Initial condition for `𝓠`

`field.reset()` zeros everything, so `𝓠` fills from 0 over `t̃ ~ 1/γ`. At
the default `γ = 0.5` that is ~2 time units (fine), but at `γ = 10⁻²` the
tone never reaches steady state within `t̃_max` and the run silently
measures the *unthrottled* M1 speed.

Add an `equilibrate 𝓠 at t=0` toggle (`q_ic_ss`, default on): seed the aux
field with `βσ̃/γ` inside the dish (2D–2D) or the screened profile
`(β h̃ σ̃/√(Dγ))·e^{−z̃√(γ/D)}` (2D–3D). Force it **on** for sweep runs.
Worth back-porting to M6.1 afterwards for parity (same silent transient).

## 6. UI (`ui.js`)

- model toggle: add `{ id: 'M6.2', label: 'M6.2 (quorum-throttled prod.)' }`.
- new `section('Cue model (M6.2)', [β, D, γ, m, IC toggle])`, shown/hidden
  by the same `m61SectionEl` mechanism (add `m62SectionEl`).
- KPIs: `𝓠_ss = βσ̃/γ` and `ℓ_Q = √(D/γ)` (geometry-switch `𝓠_ss` to
  `β h̃ σ̃/√(Dγ)`). **No speed-derived KPIs** — see §14.
- `usesFiring` gains `'M6.2'`.
- sweep `baseP`: `model: 'M6.2'` + the four knobs; sweep-done label reports
  `β` and `m` instead of `λ_A`.
- `decoratePlot('cv-rmean', …)`: axis label `\mathcal{Q}(\tilde r)`.

## 7. Plots (`render.js`)

- `drawRadialR`: third branch — guide line at `βσ̃/γ`, its own colour
  (magma-family per the dish palette; keep M6.1 blue, give M6.2 e.g. teal).
- `drawCsweep`: **superseded by §14** — no theory curve is drawn at all
  (the original plan branched analytic `c(σ̃)` curves on `params.model`).
- Emission opacity ∝ `H⁻(𝓠;1;m)` on the dish (see §4) — needs a new
  per-cell array in the frame message (`agentW`, the source weight) or
  reuse of the existing `agentR` channel, which already carries `𝓐` for
  M6.1; simplest is to make `agentR` carry `𝓠` and compute the opacity in
  the renderer from `(β,γ,m)`.

## 8. New page

`setup4-m6.2/index.html` (`4d`), a thin shell over `../setup4/ui.js` with
`{ model: 'M6.2' }`, mirroring `setup4-m6.1/index.html`: KaTeX equations of
§7b, the nondim-symbol table, the `c(σ̃)` sweep description (with the `m`
regimes), and the M6.1-vs-M6.2 discriminator paragraph. Add the nav link to
all seven existing pages plus the landing list.

## 9. Validation (PLAN §6 checklist, M6.2-specific)

1. `β = 0` reproduces M1 at the same seed (front radius vs t̃ identical to
   round-off).
2. Measured `𝓠(r̃)` plateau = `βσ̃/γ` within a few % (2D–2D) and
   `β h̃ σ̃/√(Dγ)` at `z=0` (2D–3D) — same check that caught the M6.1 3D
   tone (1.49 vs 1.41).
3. `m = 1` sweep behaves like an M6.1 sweep at `λ_A = β`, `Γ_A = γ`
   (density enters both the same way — a cross-model consistency check).
4. `m = 2` sweep is compared for *shape* against the `m = 1` and `β = 0`
   sweeps. Whether it turns over is a measurement, not a pass/fail target.
5. Grid refinement `N_grid` 80 → 160 changes `c̃` by < 5% (see §10).
6. Streaming test at fixed `σ̃`: channels present in M1, suppressed in M6.1,
   suppressed *more* in M6.2 at `m > 1` (`dc/dσ̃ < 0`).

## 10. Known issues / things that will bite

- **Grid-dependent self-tone.** In 2D the screened Green's function is
  log-divergent at the source, so a cell's *own* contribution to the `𝓠` it
  samples grows like `(β/2πD)·ln(ℓ_Q/Δx̃)` as the grid refines. M6.1 has the
  same exposure through `λ_A𝓐`, but M6.2 raises it to the power `m`, so it
  is amplified for `m > 1`. Mitigation, in order: keep `Δx̃ ≳ σ̃^{−1/2}`
  (grid coarser than the mean cell spacing — the PIC smearing then *is* the
  coarse-graining); run check §9.5 before trusting any sweep; if it fails,
  subtract the analytic self-term at sampling time.
- **`h̃` factor in 2D–3D** (§3). Easiest silent error in the whole change:
  `Q_0` is physical, so unlike M6.1's `A_0` it cannot absorb `b/h`. Getting
  it wrong rescales the tone — and so the density at which the throttle turns
  on — by `h̃`, a factor of ~10 at defaults, and the sweep will still look
  "reasonable".
- **Sweep zeros are ambiguous.** For M6.1 a zero means "did not ignite";
  for M6.2 it can also mean "ignited, front too slow to clear `1.5 r_fire`
  in `t̃_max`". `measureWaveSpeed` should distinguish them (e.g. return
  `ignited` separately from `c`, which it already computes) and the sweep
  plot should mark starved points differently from non-ignited ones.
- **Startup transient** (§5) — silently returns M1 physics at small `γ`.
- **Numerics are *easier* than M6.1**: production is bounded by `α` and
  strictly positive, so there is no runaway threshold; no new CFL beyond
  the existing `dt_sub ≤ Δx̃²/(4D)` for the aux field. Explicit Euler +
  Euler–Maruyama, unchanged.
- **`m` is a new independent Hill exponent**, so it belongs in the per-model
  exponent list in `setup4_swarm3d.md` §9 alongside `n_L, n_R, n_{Lr},
  n_A, n_{LA}` (update when implementing).

## 11. Open questions for the user

1. `m` continuous (slider) or fixed integer set {1, 2, 4}? A continuous
   slider makes the `m = 1` degeneracy with M6.1 easy to demonstrate.
2. Should `Q` also be produced by the artificial firing disk? Plan says no
   (§3) — confirm.
3. Keep M6.1 and M6.2 as separate pages (4c, 4d), or one page with a
   coupling-site toggle? Plan assumes separate pages, matching M1/M2.

## 12. What actually changed (vs this plan)

Two things surfaced during implementation:

1. **The sweep diagnostic was measuring the ignition halo, not the relay.**
   Predicted in §10 only as "sweep zeros are ambiguous"; the reality was
   worse and in the opposite direction. The firing source injects a 𝓛 mass
   `M = s_fire·σ̃·π r̃_fire²·t̃_fire`, and with `Γ_L = 0` that mass is
   conserved, so a purely diffusive halo carries `⟨𝓛⟩ > 1` out to
   `r̃_halo = r̃_fire√(s_fire σ̃ t̃_fire)` — pushing cells over threshold with
   no relay at all. Control run: M6.2 with `β = 10⁴` (production throttled to
   nothing, relay impossible) reported `c̃ ≈ 0.6`. Two fixes, both in
   `sim_core.js` / `worker.js`:
   - `measureWaveSpeed` computes `r̃_halo` and starts timing only outside it
     (`r_start = max(1.5 r̃_fire, 1.15 r̃_halo)`); if the halo would flood the
     dish it returns `status: 'unmeasurable'` instead of a number.
   - The sweep is now **mass-targeted**: `s_fire = K/(σ̃ t̃_fire)` (`fire_K`),
     so the injected mass — and `r̃_halo = r̃_fire√K` — is the same at every
     density instead of growing ∝ σ̃.
   - In 2D–3D the same mass dilutes into the half-space, so the bound is
     `r̃_halo = (π s_fire σ̃ r̃_fire² t̃_fire)^{1/3}`; the 2D formula there would
     reject every point.
   - The timing arc is gated on **radial span** (> max(0.5 r̃_fire, 0.1 r̃_max)),
     not on sample count: a near-zero baseline gives a noisy two-point average
     that reads as a spuriously *fast* wave (a slow M6.2 run once reported
     2.6× the M1 speed), while a genuinely fast front legitimately crosses the
     window in few samples.
   - Sweep geometry resized so a full window survives the halo: 2D–2D
     `R̃ = 16, N_grid = 104, t̃_max = 40, r̃_fire = 1.5, K = 20`, σ̃ ∈ [0.1, 8]
     (~30 s for 10 points). 2D–3D `R̃ = 12, N_grid = 52, K = 20` with σ̃ capped
     at **3**: 3D speeds are linear in σ̃, so the top of the 2D–2D range fills
     that dish inside one sampling interval and cannot be timed at all.
   This also corrects the **M6.1** sweep, whose high-σ̃ points were inflated
   by the same halo (its old settings gave `r̃_halo ≈ 12.7` in an `R̃ = 12`
   dish). M6.1 curves are lower and its high-λ_A points now honestly read
   "no ignition".
2. **Sweep points carry a `status`** (`ok` / `noignite` / `unresolved` /
   `unmeasurable`) and are drawn differently (filled orange / filled red /
   hollow amber / hollow gray). `unresolved` — a real but very slow front
   that never clears the halo within `t̃_max` — is common on the falling
   `m > 1` tail and must not be read as "no wave".

Not done (deliberately): back-porting the equilibrated-IC toggle to M6.1
(§5). M6.1 still starts 𝓐 at zero, so its results are unchanged by this work
apart from the sweep-metric fix above.

## 13. Measured behaviour (node, 2D–2D unless stated)

- `β = 0` reproduces M1 **bit-for-bit** (front radius identical to 6 d.p.),
  as does M6.1 at `λ_A = 0`.
- 2D–2D tone: seeded `𝓠_ss = βσ̃/γ` is an exact steady state — it holds to
  3 d.p. over `t̃ = 10`. Free (unseeded) fill-up converges to ~1.11× the
  mean-field value at `Δx̃ ≈ 0.31` (the log self-term of §10).
- 2D–3D tone at `z = 0`: 0.618 vs mean-field 0.707 (−12.5%), and M6.1 in the
  same configuration is off by the same relative factor (1.243 vs 1.414) —
  i.e. the pre-existing surface-source discretisation offset, not the β
  bookkeeping. `𝓠(M6.2, β=0.5) / 𝓐(M6.1) = 0.497` confirms β enters linearly
  and the `h̃` factor is right. Independent of `z_max` (tested to 262 ℓ₀).
- `m = 1` cross-check: at `β = λ_A = 0.4`, `γ = Γ_A = 0.5` the two models
  track each other across σ̃ = 0.25–2 (0.20/0.22, 0.38/0.36, 0.53/0.60,
  0.73/0.78). Above σ̃ ≈ 2 they part as expected from the *mechanisms*: M6.1
  stops igniting (threshold), M6.2 keeps propagating (starvation only).
- `m = 2`, `β = 0.2`, `γ = 0.5`: measured 0.27, 0.50, 0.71, **0.98, 0.93,
  0.69** over σ̃ = 0.25…8 — non-monotonic, turning over near σ̃ ≈ 2, against a
  monotonically rising M1 sweep (0.33, 0.55, 0.87, 1.28, 1.79, 2.57) on the
  same settings. Reported as measurements; no formula is fitted to them.
- Fully throttled control (`β = 10⁴`) reports `c̃ = 0` at every density.

## 14. Analytical c(σ̃) withdrawn (2026-07-28)

After implementation the mean-field pushed-front flux balance
`c² 𝓛_c = α_eff σ̃` was judged not to describe this system reliably, so
**every analytical `c(σ̃)` estimate for M6.1 and M6.2 was removed** rather
than kept with caveats:

- `drawCsweep` (render.js) plots measured points only — no theory curve, no
  `β = 0` reference curve, no saturation asymptote, no `σ̃★` / `σ̃_peak`
  markers. The only vertical line left is the current single-run σ̃.
- The `σ̃*` and `σ̃_peak` KPIs are gone. `tone_ss` and `ℓ_screen` remain —
  they are exact properties of the *linear* auxiliary field, and the tone was
  verified numerically to be an exact steady state.
- The M6.1 and M6.2 pages state the mechanism — `𝓛_c = 1 + λ_A𝓐_ss` and
  `α_eff = α/(1 + 𝓠_ss^m)`, both exact given the uniform tone — and then say
  the speed is measured, not predicted.
- Catalog §7/§7b carry an explicit "no closed-form c(σ̃)" note in place of
  the old derivations.

Consequence for §§9 and 13: comparisons are now *sweep vs sweep* (β = 0 vs
β > 0, m = 1 vs m = 2, M6.1 vs M6.2) rather than sweep vs formula. The
measured numbers in §13 stand; the theory values they were compared against
do not.

## 15. Second petri dish (2026-07-28)

The M6.2 page shows two dishes side by side: cells over 𝓛 (magma) and the
same cells over 𝓠 (greens). Implementation notes:

- `drawDish` is parameterized (`fieldKey`, `palette`, `haloScale`) rather than
  duplicated, so the cell markers in both panels come from the same code path
  and can be compared pixel for pixel.
- The worker ships `Qfield` only for M6.2, Uint8-compressed like `Lfield` but
  normalized by its own **grid max**, not the radial-profile max: 𝓠 is nearly
  uniform, and radial normalization would render it as a flat wash. The 𝓠
  panel therefore shows *structure* (where the signal piles up), not level —
  level is read off the 𝓠(r̃) profile panel.
- Emission halos are damped (`haloScale: 0.45`) on the 𝓠 panel so the warm
  bloom does not swamp the green field.
