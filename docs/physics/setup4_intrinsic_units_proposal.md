# Setup 4 — Intrinsic-units nondim proposal (deferred, 2026-05-26)

> Status: **deferred**. Implemented end-to-end on 2026-05-26 then reverted
> because the cascading numerical/UX consequences were too disruptive at
> this stage. The math is consistent and the scheme is workable; this doc
> records the full derivation, the implementation surface, and the bugs
> that surfaced, so the work can be redone deliberately if/when we
> revisit. The pre-rewrite snapshot is preserved at git tag
> `setup4-pre-intrinsic-units` and the reverted commits are
> `41ec45e..1481f7f` (5 commits) in the repo's reflog.

---

## 1. Motivation

The current Setup 4 nondim bakes the cell density σ into the unit system
through the choice of velocity scale:

  2D-3D:  c* = a σ / L_0       →  ℓ_0 = D_L / c* = D_L L_0 / (a σ)
  2D-2D:  c* = √(a σ D_L / (h L_0))  →  ℓ_0 = √(D_L h L_0 / (a σ))

so ℓ_0 depends on σ. Every time you nudge N or R_dim (which changes
σ = N/(π R_dim²) under the β=1 convention), the *unit of length itself*
shifts. Downstream consequences:

- R̃ = R_dim/ℓ_0, dx̃ = 2R̃/N_grid, and Λ, χ̃, μ̃, λ, ϑ, Γ̃_L all
  rescale together when N or R_dim moves.
- At high N (small ℓ_0), the firing footprint r̃_fire = 2 in nondim falls
  inside a single grid cell — only ~1–2 cells get ignited and the wave
  is invisible.
- σ̃ = σ·ℓ_0² is a derived quantity that moves inversely with σ in 2D-3D
  (because ℓ_0 ∝ 1/σ), contradicting the intuition that σ̃ should be a
  "dimensionless cell density."
- Dieterle wave speed is c̃ = 2/π by construction, not a measurable
  result. The actual physics-encoded ratio (cells-per-influence²) is
  hidden in the unit system.

## 2. Proposed alternative

Use the **single-cell influence radius** as the length scale — a property
of one cell + diffusion + threshold, with no reference to other cells:

  ℓ_0 ≡ ℓ_a = a / (D_L L_0)                  [µm]
  t_0 ≡ t_a = ℓ_a² / D_L = a² / (D_L³ L_0²)  [s]

Physical interpretation: solving the steady-state 3D-diffusion profile
around an isolated emitting cell, L = a/(4π D_L r), and setting L = L_0
gives r = ℓ_a / (4π) — so ℓ_a is (up to a 4π factor) the distance at
which a single cell's signal reaches the relay threshold.

Nondim diffusivity is then $\tilde D = D_L\, t_a/\ell_a^2 = 1$ identically.

Only the **three params (a, D_L, L_0)** define the unit system; σ, R, N, h,
and every cell-side parameter live INSIDE the units, so dragging N or
R_dim does not rescale ℓ_a, t_a, or any cell-side coefficient.

## 3. LTB4 PDE in the new units

### 3.1 2D-3D (bulk z>0 diffusion, surface source at z̃=0)

  ∂_{t̃} L̃ = ∇̃² L̃                                 (bulk)
  −∂_{z̃} L̃|_{z̃=0}
     = Σ_i H⁺(L̃_i;1;n_L) δ̃²(r̃ − r̃_i)               (discrete)
     = σ̃ H⁺(L̃;1;n_L)                                  (continuum)

The per-cell surface-source prefactor is exactly 1 — no 1/σ̃ factor
anywhere. σ̃ emerges only in the continuum limit Σ_i δ̃²(r̃ − r̃_i) → σ̃.

### 3.2 2D-2D (in-plane diffusion + averaged bulk source over layer height h)

  ∂_{t̃} L̃
     = ∇̃² L̃  +  (1/h̃) Σ_i H⁺(L̃_i;1;n_L) δ̃²(r̃ − r̃_i) − Γ̃_L L̃    (discrete)
     = ∇̃² L̃  +  (σ̃/h̃) H⁺(L̃;1;n_L) − Γ̃_L L̃                       (continuum)

where $\tilde h ≡ h/\ell_a$ is a **new** dimensionless layer height — it
did not appear in the σ-baked scheme because h was folded into c*.

## 4. Wave speeds — results, not unit definitions

  2D-3D continuum: c̃ = σ̃,                  c̃_Dieterle = (2/π) σ̃
  2D-2D continuum: c̃ = √(σ̃/h̃),            c̃_Dieterle = (2/π)√(σ̃/h̃)

In dim units these collapse to the same expressions as before — the
*physical* wave is unchanged; only the unit system has moved.

Front-layer width in ℓ_a units: ~ 1/σ̃ — wide at low cell density, narrow
at high density. So σ̃ acts as a regime indicator:

  σ̃ ≪ 1 : discrete-cell regime (Dieterle prediction does not apply)
  σ̃ ≫ 1 : continuum / Dieterle regime

## 5. Full nondim coefficient table

With $p_0 = \sqrt{u/w}$ as before:

| symbol | old formula (σ-baked) | new formula (intrinsic) | σ-dependent? |
|---|---|---|---|
| Λ | r_0 L_c t_old | **r_0 L_c t_a** | no |
| L_c_nd | L_c / L_0 | same | no |
| χ̃ | χ L_0 / (p_0 c*) | **χ L_0 t_a / (p_0 ℓ_a)** | no |
| μ̃ | μ p_0 / c* | **μ p_0 t_a / ℓ_a** | no |
| λ | u²/(w r_0 L_c) | same | no |
| ϑ | θ w/(u r_0) | same | no |
| Γ̃_L | Γ_L t_old | **Γ_L t_a** | no |
| σ̃ | σ ℓ_old² (∝ 1/σ in 2D-3D) | **σ ℓ_a² (true density)** | yes — emerges in continuum |
| R̃ | R/ℓ_old (∝ σ) | **R/ℓ_a** | no |
| h̃ | — | **h/ℓ_a** | no |

**Every cell-side coefficient becomes σ-independent.** Only $\tilde\sigma$
and (in 2D-2D) $\tilde h$ carry geometry/density information.

## 6. Dieterle overlay update

The OLD overlay (in OLD units) was

  L̃(ξ̃_old) = √(−2 ξ̃_old)                              (behind)
            = (π/2) e^{−2ξ̃_old/π} / √(2 ξ̃_old)         (ahead)

with $\tilde c = 2/\pi$ universally.

In new units, substituting $\xi_{old} = \tilde\sigma\,\xi_{new}$ (since
$\ell_a/\ell_{old} = \tilde\sigma$):

  L̃(ξ̃_new) = √(−2 σ̃ ξ̃_new)                            (behind)
            = (π/2) e^{−2σ̃ξ̃_new/π} / √(2 σ̃ ξ̃_new)     (ahead)

  c̃_Dieterle = (2/π) σ̃

i.e. the analytic profile **scales by σ̃** in the new units; the wave
speed is now a function of density. The buffer ±0.5/σ̃ excluded around
the front (where both branches diverge) widens at low σ̃.

## 7. Implementation surface (what would change in the code)

| file | change |
|---|---|
| **setup4/nondim.js** | rewrite `dimToNondim`: ℓ_0 = a/(D_L L_0), t_0 = ℓ_0²/D_L; output `h_tilde`; c* kept as diagnostic only |
| **setup4/worker.js** | drop 1/σ̃ from cell-source pass-through; add 1/h̃ factor in 2D-2D; firing source needs σ̃·cellPrefactor rescaling (so s_fire=1 ≡ saturated-cell-equivalent flux, preserving the OLD physical interpretation) |
| **setup4/solvers/*.js** | **no code change** — solver consumes src agnostically; comment headers updated for the new formulas |
| **setup4/agents.js** | no change |
| **setup4/render.js** | Dieterle profile takes σ̃ instead of β; analytic wave speed annotation reads c̃ = (2/π)σ̃ |
| **setup4/ui.js** | new slider defaults: Λ~0.04, χ̃~0.015, μ̃~0.001; dt up to ~10 (CFL allows large dt because Λ small); t_max up to ~10000 (wave traversal is large in t_a units); h_tilde, z_max KPIs; **R̃ becomes the nondim Geometry knob** (R_dim follows as KPI) |
| **docs/physics/setup4_swarm3d.md** | rewrite §5 (nondim + Dieterle overlay) |
| memory `project_setup4_dieterle.md`, `project_setup4_cue_models.md` | document new units |

The pipeline changes are surgical (only worker.js and nondim.js touch
the math); solver_m1.js stays bit-identical. The bulk of the diff is
slider defaults / ranges / KPI labels.

## 8. Bugs and surprises found during the deferred implementation

These are the issues that surfaced and the reasoning around each — they
should be considered carefully before trying again.

### 8.1 Per-cell source prefactor (math is correct; UX requires care)

In the OLD scheme the worker passes `h_emission · (1/σ̃_old) · dt` per
cell; in the new scheme it passes `h_emission · dt` (2D-3D) or
`h_emission · dt / h̃` (2D-2D). Both produce the same dim L profile when
all else is right. The continuum check:
$\Sigma_i w_{ik}/dx^2 \to \tilde\sigma$, so `surfSrc_continuum = 2 σ̃ H⁺/h_0`
in 2D-3D and `directSrc_continuum = σ̃ H⁺/h̃` in 2D-2D — exactly the
continuum source in the new PDE.

### 8.2 Firing source: σ̃-rescaling required

The OLD firing call `addFiringSource(r_fire, s_fire·dx²·dt)` gives
continuum surface flux = s_fire (in OLD units, where cell flux = H⁺ ≤ 1).
So OLD s_fire=1 = "saturated cell's worth of flux".

In new units, the continuum cell flux is σ̃·H⁺ ≤ σ̃ — orders of magnitude
smaller. Keeping the OLD `s_fire·dx²·dt` call makes s_fire=1 ≈ **1/σ̃
saturated cells worth** (e.g. ≈10× at σ̃=0.1), so the firing region
dominates everything and the radial profile drowns the Dieterle overlay.

**Fix:** multiply firing amount by σ̃·cellPrefactor so s_fire=1 means
"one saturated cell's worth per area" in the *current* units, matching
OLD convention's physical scale. Dim equivalent is identical between
schemes; only the nondim s_fire knob retains the same meaning.

### 8.3 z-grid truncation — the real "L is 10× Dieterle" bug

The CN z-grid has total extent
$z_{max} = h_0 \cdot (\alpha^{N_z} − 1)/(\alpha − 1)$. Defaults
$h_0=0.1, \alpha=1.1, N_z=16$ give $z_{max} \approx 3.6$ nondim.

In OLD units $t̃_{traverse} = R̃$ (≈ 10–30 for typical defaults), so the
z-diffusion spread $\sqrt{t̃_{traverse}} \approx 3–6$ was *just barely*
contained inside z_max=3.6. Coincidence — but it worked.

In new units $t̃_{traverse} = R̃/\tilde c = R̃ \cdot \pi/(2\tilde\sigma) \approx 1000$
at defaults → z-spread $\approx \sqrt{1000} \approx 32$ in ℓ_a units. But
z_max stayed at 3.6 — **10× too small**. L hits the no-flux BC at z_max,
reflects, and the column accumulates instead of leaking into half-space.
Net effect: L̃[0] at center is ~10× too large.

**Fix:** bump α_z default to 1.5 (giving z_max ≈ 131, 4× margin over
wave-spread). Slider max raised to 2.0. Add a `z_max` KPI so any future
recurrence is visible by inspection.

### 8.4 Persistent firing source perturbs the Dieterle profile

Default `t_fire = 1500` (the user's "60s" spec in new dim units) is
longer than the wave traversal time t̃_traverse ≈ 1640. So the firing
source pumps the firing region for the *entire* simulation. A sustained
surface source over a half-space gives an erfc-like buildup

  $\tilde L|_{z=0}(t) \approx 2 \tilde\sigma s_{fire} \sqrt{t/\pi}$

which at defaults integrates to ≈ 4 at center over t̃=1500 — comparable
to the Dieterle prediction itself (≈ 4.4). So even after the σ̃ rescaling
and z_max fix, sim shows ≈ Dieterle + firing-buildup ≈ 8 at center.

**Fix:** drop default t_fire to ~200 (just past the ~86 nondim ignition
time at σ̃=0.0955, s_fire=1) — robust kick, then off, wave self-sustains.
Also exclude r<r_fire from the Dieterle overlay (the firing region is
fundamentally outside the pure-Dieterle assumption).

This was the issue the user resisted most: their "60s pulse" was a
physical experimental intent, but in the new units it became
incompatible with clean Dieterle visualization. Any future revisit needs
to either (a) accept t_fire ≪ t_max as a precondition for Dieterle
match, or (b) actively decay the firing-region L̃ contribution somehow.

### 8.5 At default σ̃ ≪ 1, Dieterle is a continuum-limit prediction

With R_dim=200µm, N=3000 → σ̃ ≈ 0.0955 in new units — deep in the
discrete-cell regime. The sim is jagged (single-cell hot spots) and the
smooth analytic profile is only an order-of-magnitude reference. For a
clean visual match, push N to ~30000 (σ̃ ≈ 1) to enter continuum. This
is not a bug, but worth noting as expected behavior the user may
mistake for one.

## 9. Why this was deferred

After all four fixes (8.1–8.4) were applied, the visual Dieterle match
was still imperfect because:

- Default σ̃ ≈ 0.1 sits in the discrete regime (8.5).
- The user's preferred physical default t_fire = 60s conflicts with
  Dieterle's "trigger then self-sustain" assumption (8.3, 8.4).
- Slider re-tuning and default changes cascaded through the whole UI —
  every default value and slider range had to be re-derived. Several
  rounds of "this looks weird now" iteration followed.

In effect: the math was right and the implementation was correct, but
the cumulative UX disruption (changed defaults, changed slider semantics,
new KPIs, doc rewrite, several rounds of "why doesn't the overlay
match?") exceeded the benefit at this stage. The original σ-baked scheme
is empirically tuned to give a believable Dieterle overlay at small N,
and the user prefers to keep that.

## 10. If we revisit — checklist

1. **Re-read §8 in full** before touching any code.
2. **z_max must scale with t̃_traverse from the outset** — don't keep
   the OLD α_z=1.1 default. Pick α_z and h_0 such that $z_{max} > 5
   \sqrt{R̃/\tilde\sigma}$ for the expected defaults.
3. **Firing-source rescaling is non-negotiable** — without σ̃·cellPrefactor
   on the firing call, s_fire=1 dominates the cell sources by factor
   1/σ̃ and Dieterle never matches.
4. **t_fire ≪ t_max** is a precondition for visual Dieterle match.
   Either default to a short pulse and let advanced users override, or
   keep the long pulse but tell users explicitly the overlay won't
   match in the firing region (which is also masked from the overlay
   draw — see commit 1481f7f for the masking code).
5. **σ̃ as a regime KPI** — surface a discrete/continuum indicator so
   the user knows when the Dieterle overlay is supposed to match.
6. **Default N≥30000 (σ̃ ≈ 1)** if you want the overlay to match at
   defaults, otherwise document the discrete-regime caveat next to the
   plot.
7. **Test against the parent Dieterle script** end-to-end — pick a
   common (a, σ, R, h) set, run both the parent and the new Setup 4,
   confirm the L̃ profiles match (modulo discrete noise) at the same
   dim time.
8. **Don't rewrite the docs until the code stabilizes** — last attempt
   the cue_models doc and the swarm3d doc went through several
   inconsistent intermediate states.

## 11. Git breadcrumbs

- `setup4-pre-intrinsic-units` (tag): clean revert target. Equivalent
  to commit `6bd17cd`.
- Commits `41ec45e..1481f7f` (5 commits): the reverted implementation
  in order — nondim switch, Geometry-R̃ swap, firing rescaling,
  z_max bump, t_fire/overlay masking. Recover from reflog if needed.
- This doc itself: written on the post-revert HEAD; safe to keep
  forever as the design memo.
