# Setup 4 — Intrinsic-units nondim implementation plan (2026-05-31)

> Status: **planned, not yet implemented.** Supersedes the deferred
> proposal at [setup4_intrinsic_units_proposal.md](setup4_intrinsic_units_proposal.md)
> (kept as the postmortem of the 2026-05-26 attempt). User reviewed the
> current code surface and re-derived the nondim equations independently;
> this plan reflects the updated, locked-in choices.

---

## 1. Locked-in scheme

**Length and time:**
```
ℓ_0 = a / (L_0 D_L)
t_0 = a² / (L_0² D_L³)   ( = ℓ_0² / D_L )
```
σ drops out of the unit system entirely.

**Field normalizations:**
```
𝓛   = L / L_0
𝓛_c = L_c / L_0
𝓛_r = L_r / L_0
ℛ_i = R_i / R_0      (R_0 = inhibition threshold; renamed from R_c)
P_i = p_i / √(u/w)
```

**Nondim cell groups (five independent):**
```
λ = r·t_0·L_0                      ← uses L_0, not L_c
ν = u²·t_0 / w
κ = a·χ / (D_L²·√(u/w))
μ = μ_dim · t_0 / ℓ_0              (μ_dim is the dim slider; m treated as 1/μ_dim notationally)
ϑ = (w·L_0·t_0 / u) · θ
```
No Λ wrapper.

**Other groups:**
```
σ̃   = σ · ℓ_0²                    (true cell density now, σ-independent units)
h̃   = h / ℓ_0                     (new — only relevant in 2D-2D)
Γ_L = γ_L · t_0
Γ_R = γ_R · t_0
β   = b · t_0 / R_0
```

**Cell SDE:**
```
Ṗ_i = κ ∇𝓛(r̃_i) + λ(𝓛 − 𝓛_c)P_i + ν(|P|² − |P|⁴)P_i + √(2ϑ𝓛) ξ
ṙ̃_i = μ P_i
```
Noise has **no Λ factor** (different from current code).

**L PDE (discrete-ABM form):**
```
2D–3D:  ∂_t̃ 𝓛 = ∇̃²𝓛 + δ(z̃) Σ_i H⁺(𝓛_i;1;n_L) H⁻(ℛ_i;1;n_R) δ²(r̃ − r̃_i) − Γ_L 𝓛
2D–2D:  ∂_t̃ 𝓛 = ∇̃²𝓛 + (1/h̃) Σ_i H⁺(𝓛_i;1;n_L) H⁻(ℛ_i;1;n_R) δ²(r̃ − r̃_i) − Γ_L 𝓛
```
Per-cell prefactor: 1 in 2D-3D, 1/h̃ in 2D-2D. No 1/σ̃ anywhere.

**Inhibitor ODE (M2, per cell):**
```
∂_t̃ ℛ_i = β H⁺(𝓛(r̃_i); 𝓛_r; n_{Lr}) − Γ_R ℛ_i
```

**UI primary knobs:**
- Geometry section: R̃ primary (nondim), σ primary (dim slider), N = round(σ·π·R_dim²) derived (KPI).
- R_dim moves to KPI; no longer a primary knob.

**Firing source convention:** s_fire is the dim-equivalent factor multiplying `a` inside the firing disk. In nondim, that means the firing source is rescaled by σ̃ so that s_fire=1 corresponds to "every point in the firing disk emits as strongly as one saturated cell."

**Dieterle overlay:** kept, redrawn with σ̃ rescaling. c̃_Dieterle = (2/π)σ̃; ξ → σ̃·ξ in both branches.

---

## 2. Phase A — Math layer (no UI churn yet)

### A1. Rewrite `setup4/nondim.js`

- Replace `cstar` branch with `ell0 = a/(L_0*D_L); t0 = ell0*ell0/D_L`.
- Replace output fields:
  - Drop: `Lambda, chi_nd, mu_nd`.
  - Add: `lam, nu, kap, mu, tht` with formulas from §1.
- Keep `sigma_tilde = sigma·ell0²` (now σ-independent units; σ̃ is a true density).
- Add `h_tilde = h/ell0`.
- Keep `gamma_L, beta_R, gamma_R, L_r_nd, L_c_nd` (formulas unchanged in spirit; all rates now × t_0_new).
- Drop `cstar` from active outputs (retain as comment-only diagnostic).
- Re-tune `DIM_DEFAULTS` to put λ, ν, κ, μ, ϑ in slider-meaningful magnitudes at startup. Bench-tuning is D1's job.

### A2. Rewrite `setup4/worker.js` source-accumulation loop

- Per-cell prefactor:
  - 2D-3D: drop `inv_sigma`; per-cell weight is `w · dt`.
  - 2D-2D: pass `1/h_tilde` instead; per-cell weight is `w · dt / h_tilde`.
- Firing source: multiply amount by `sigma_tilde` (memo §8.2):
  ```js
  field.addFiringSource(r_fire, s_fire * sigma_tilde * field.dx * field.dx * dt);
  ```
- Update the comment blocks at worker.js:222-229 and 254-263 to reflect the new PDE form and firing convention.
- Param destructure: replace `Lambda, chi_nd, mu_nd` with `lam, nu, kap, mu, tht`; add `h_tilde`.

### A3. Update `setup4/agents.js` SDE assembly

- New SDE drift:
  ```
  Ṗ = κ·∇𝓛 + λ·(𝓛 − 𝓛_c)·P + ν·(|P|² − |P|⁴)·P + √(2ϑ𝓛)·ξ
  ```
- Adaptive sub-stepping criterion changes. Replace
  `drift_lin_mag = Lambda * Math.abs(L_i - L_c)`
  with
  `drift_mag = Math.max(lam * Math.abs(L_i - L_c), nu * (m2 + m4))`
  (use max of the two drift coefficients to pick `n_sub`).
- Noise: `noise_sub = √(2·tht·L_i) · √(dt_sub)` — drop the Λ factor.
- Position update: `xi += mu * px * dt_sub` (μ now means the new `mu` field).

### A4. Solver files

- `setup4/solvers/solver_m1.js` and `setup4/solvers/field.js`: no code change (the solvers consume `src[]` agnostically).
- Update header comment blocks:
  - solver_m1.js:1-24 — drop 1/σ̃ in 2D-3D nondim PDE form; add 1/h̃ in 2D-2D.
  - field.js:116-127 — rewrite the `addFiringSource` docstring to reflect σ̃-rescaled meaning.

---

## 3. Phase B — UI layer

### B1. `setup4/ui.js` params block + sliders

- **Cell-nondim section:** replace `Lambda, chi_nd, mu_nd` sliders with `lam, nu, kap, mu, tht`. Five independent nondim sliders, with re-derived ranges (D1).
- **Geometry section primary-knob flip:**
  - R̃ slider becomes primary; R_dim moves to KPI.
  - σ becomes a dim slider in Geometry.
  - N = round(σ_dim · π · R_dim²) — derived, displayed as KPI.
  - Tracked-cell-slider max follows derived N.
- **Defaults:**
  - `params.R_dish`: pick sensible ℓ_0-units value (memo suggested R̃ ≈ 10–30).
  - `t_max, dt, r_fire, t_fire`: re-derive (memo §8.4 — t_max much larger in new units; t_fire ≪ t_max for clean Dieterle).
- **`recomputeFromDim()` flow inverts:**
  - σ slider → params.sigma_tilde directly.
  - R̃ slider → R_dim KPI via R_dim = R̃ · ℓ_0.
  - N updates via N = round(σ_dim · π · R_dim²).
- Push new nondim values to sliders for λ, ν, κ, μ, ϑ.
- **KPIs (ui.js:415-424):**
  - Drop redundant entries.
  - Add: `c̃_Dieterle = (2/π)σ̃`, `z_max` (= h_0·(α^N_z−1)/(α−1)), `R_dim` (derived from R̃·ℓ_0).
  - Keep: σ̃ (color-tagged as regime indicator: ≪1 = discrete, ≳1 = continuum), N.

### B2. 2D-3D z-grid defaults (silent 10× L bug — memo §8.3)

- Bump default `alpha_z` from 1.1 → 1.5.
- Raise slider max to 2.0.
- Verify at defaults: `z_max ≈ h_0·(α^N_z − 1)/(α − 1) > 5·√(R̃/σ̃)`.
- Add `z_max` KPI for inspection (per B1).

### B3. h dim slider exposure

- h is now relevant only as the numerator of h̃ = h/ℓ_0 (1/h̃ prefactor in 2D-2D).
- Slider remains visible only in 2D-2D (current behavior).
- Default h re-tuned to give h̃ ≈ O(1).

---

## 4. Phase C — Visualization

### C1. Dieterle overlay (`setup4/render.js:317-389`)

- Replace `dieterleProfile(xi, beta)` with `dieterleProfile(xi, sigma_tilde)`:
  - `z = sigma_tilde * xi` instead of `beta * xi`.
- Wave-speed annotation: `c̃ = (2/π)σ̃` (no longer constant 2/π).
- Buffer scales as `0.5/σ̃` (front-layer width).
- Drop the `beta = (N/(πR²))/σ̃` calculation — identically 1 now.
- Discrete-regime warning: when σ̃ ≲ 0.3, render an "approximate (discrete regime)" note next to the overlay annotation.

### C2. `setup4/index.html` description block (lines 63-104)

- Replace c*-derivation paragraph with the intrinsic-units derivation: ℓ_0 = a/(L_0 D_L), t_0 = ℓ_0²/D_L. Cite the "single-cell influence radius" interpretation.
- Replace "c*=1, c=1 vs c=2/π" text with: "Measured wave speed in nondim units is c̃_M1 = σ̃ (2D-3D) or √(σ̃/h̃) (2D-2D), times the (2/π) Dieterle prefactor in 2D-3D."
- Replace the nondim-table to list λ, ν, κ, μ, ϑ (not Λ, χ̃, μ̃).
- Update the bead-plot free-energy formula:
  ```
  F(P) = −[λ(𝓛 − 𝓛_c)|P|²/2 + ν(|P|⁴/4 − |P|⁶/6)] − κ·∇𝓛·P
  ```

### C3. Bead-plot rendering (`setup4/render.js` `drawBead*`)

- Free-energy formula uses (λ, ν, κ) now.
- `beadParams` passed in ui.js:331-335 renamed: `Lambda→lam`, `lam→nu`, `chi→kap`, `L_c→L_c_nd`.
- Audit the three drawBead* functions in render.js to use the new names consistently.

---

## 5. Phase D — Calibration + verification

### D1. Pick defaults that give a believable startup state

- Goal: at default (σ_dim, a, D_L, L_0, R_dim), the values of λ, ν, κ, μ, ϑ, σ̃, R̃ land in slider-meaningful ranges.
- Acceptance:
  - Heatmap shows a propagating wave within t̃ ≤ default t_max.
  - Dieterle overlay visually within 20% of the simulated front (in continuum regime σ̃ ≳ 1).

### D2. Cross-check against a known dim case

- Pick one (a, σ, D_L, L_0, R_dim) set.
- Run the parent `la_2d3d_solver.py` and the rewritten Setup 4.
- Confirm L̃ profile vs r at matched dim time matches within discrete-noise tolerance (memo §10 #7).

### D3. Numerical-stability self-check

- Console warning when, at any cell:
  - `dt · max(λ·|𝓛 − 𝓛_c|, ν, Γ_R) > 0.3`, or
  - `z_max < 2·√(t_max)`.

---

## 6. Phase E — Docs and memory (last, per memo §10 #8)

### E1. `docs/physics/setup4_swarm3d.md`

- Rewrite §5 (nondim), §9 derived-groups, §11.3 (δ(z) cross-check now references c̃ = σ̃·2/π).

### E2. `docs/physics/setup4_intrinsic_units_proposal.md`

- Status: "deferred" → "implemented YYYY-MM-DD" with reference to implementing commits.
- Carry forward §8 lessons as a postmortem.

### E3. Memory files

- Update `project_setup4_dieterle.md` and `project_setup4_cue_models.md` with new formulas and the active-scheme statement.

### E4. `setup4/index.html` inline description

- Already in C2 — flagged here as part of user-facing docs.

---

## 7. Deliberately not in this plan

- The mass parameter m as a separate dim knob (chose to keep μ_dim, treat m as 1/μ_dim notationally).
- Reframing cell-side groups as KPI-only (chose nondim sliders).
- Dropping the firing source (chose to keep with σ̃-rescaling).
- Rewriting solver_m1.js math (confirmed unchanged — only comments).

---

## 8. Risks flagged (carried from deferred memo §8)

1. **z_max silent 10× L bug** — handled proactively in B2.
2. **Firing-source dominance without σ̃-rescaling** — handled in A2.
3. **σ̃ ≪ 1 default = discrete regime** — handled by σ̃ regime indicator + warning in C1.
4. **Default cascading UX disruption** — D1 is its own phase. Do not merge to main until D1 looks right.

---

## 9. Pre-implementation hygiene

- Create backup branch and/or git tag `setup4-pre-intrinsic-units-v2` before Phase A. The previous attempt's tag was `setup4-pre-intrinsic-units` (commit 6bd17cd).
