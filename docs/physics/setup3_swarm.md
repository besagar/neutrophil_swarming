# Setup 3 — 2D radial swarm with central trap

A "petri dish" of `N` neutrophils on a 2D disk. A single Gaussian cue wave is
launched from the center and propagates radially outward. A central trapping
disk of radius `R̃_trap` holds any cell that enters it.

This is a 2D extension of [setup2_wave.md](setup2_wave.md): the polarization
SDE is unchanged in form, only generalized from scalar `P` to a 2D vector
`P = (P_x, P_y)`, and the cue field is radially symmetric instead of planar.
**All nondimensional groups are identical to Setup 2** — `M, C, λ, ϑ, χ̃, μ̃`
keep their meanings and ranges.

## Geometry

- 2D disk `r̃ < R̃_dish` (default `R̃_dish = 10`, in units of the wave width `σ`).
- `N` cells (default `N = 1000`) initialized uniformly on the annulus
  `R̃_trap < r̃ < R̃_dish` (uniform in area: sample `r̃ = √(R̃_trap² + U·(R̃_dish² − R̃_trap²))`,
  `θ = 2πU'`).
- Hard reflective outer boundary at `r̃ = R̃_dish`.
- Hard one-way inner boundary at `r̃ = R̃_trap`: cells entering are clamped to
  `r̃ ≤ R̃_trap` thereafter and never re-emerge ("trapped" flag).

## Cue field

Sum of `n_waves` outward-propagating radial Gaussian pulses, launched at
evenly-spaced times `t̃_j = j · Δt̃_wave` for `j = 0, …, n_waves − 1`:

```
𝓛_j(r̃, t̃) = M · exp( -(r̃ - C (t̃ - t̃_j))² / 2 )    for t̃ ≥ t̃_j, else 0
𝓛(r̃, t̃)   = Σ_j 𝓛_j(r̃, t̃)
```

Each `𝓛_j` is a radially-symmetric pulse whose peak sits at radius
`r̃ = C·(t̃ − t̃_j)`. Because `𝓛` depends on position only through
`r̃ = |X̃|`, its 2D gradient is purely radial:

```
∇𝓛 = (∂_r̃ 𝓛) · r̂,    where    ∂_r̃ 𝓛 = -Σ_j 𝓛_j · (r̃ - C(t̃ - t̃_j))
```

The factor `(r̃ − C(t̃−t̃_j))` is the signed offset from the j-th wavefront
peak — positive ahead of the peak (gradient points inward toward higher
cue) and negative behind it (gradient points outward).

Amplitude `M` is constant in `r̃` (no `1/r̃` decay) so all Setup 2 nondim
groups carry over unchanged. The single-wave case (`n_waves = 1`) is the
literal Setup 2 cue extended to 2D radial geometry.

## Dynamics

For each cell `i`, with `X_i = (x̃_i, ỹ_i)`, `r̃_i = |X_i|`,
`r̂_i = X_i / r̃_i`, polarization vector `P_i = (P_xi, P_yi)`:

```
dP_i/dt̃   = χ̃ ∇𝓛(r̃_i, t̃) + (𝓛(r̃_i, t̃) - 1) P_i
            + λ (|P_i|² - |P_i|⁴) P_i
            + √(2 ϑ 𝓛(r̃_i, t̃)) η_i              (η_i: 2D iid Gaussian per dt̃)
dX_i/dt̃   = μ̃ P_i        (free cells)
dX_i/dt̃   = 0           (trapped cells)
```

Cells whose `r̃_i` first crosses below `R̃_trap` are flagged trapped. Their
polarization keeps evolving (so the GL well retains memory) but their
position is frozen at the entry point. Outer boundary: if a step would take
`r̃_i > R̃_dish`, reflect the position component of the step radially.

## Alternative model (UI toggle)

Setup 3's UI exposes a model toggle: the default GL polarization SDE above,
or a Hill / linear-relaxation alternative. See
[hill_model.md](hill_model.md) for the Hill 2D-vector equations. The trap,
dish boundary, multi-wave cue, and all aggregate plots are model-agnostic.
`χ̃` and `ϑ` sliders are shared between models; `μ̃, λ` are GL-only;
`ṽ₀, n` are Hill-only.

## Physical regime of interest

In Setup 2 the GL cell ends up displaced **in the same direction the wave
travels** at large `χ̃` and moderate `λ, ϑ`. There exists a parametric regime
(typically smaller `χ̃` and/or specific `λ, ϑ`, see Setup 2 sweep panel) in
which the cell ends up displaced **against** the wave direction. **This
"anti-wave" regime is the one of biological interest for swarming**: with an
outward-propagating wave from a central candida, anti-wave displacement
drives cells inward, where the trap captures them.

We do not flip the wave direction or the sign of `χ̃` to engineer attraction;
the inward drift is a genuine prediction of the GL model in that regime.

## Nondim groups (identical to Setup 2)

```
M  = L_max / L_c                           peak cue / threshold
C  = c / (σ r_0 L_c)                       wave speed (σ per t_0)
χ̃  = χ / (σ r_0 √(u/w))                    chemotactic strength
μ̃  = μ √(u/w) / (σ r_0 L_c)                cell motility
λ  = u² / (w r_0 L_c)                      well depth (intrinsic)
ϑ  = θ w / (u r_0)                         noise (intrinsic)
```

Length scale `σ` (wave width), time scale `t_0 = 1/(r_0 L_c)`, polarization
scale `p_0 = √(u/w)` — same as Setup 2.

## What to plot

- Live 2D animation: dish (disk), trap (faint inner disk), wavefront ring at
  `r̃ = C t̃`, cells colored by `|P|` (free = orange, trapped = green).
  Optionally short polarization arrows.
- Mean radius of free cells `⟨r̃⟩_free(t̃)`.
- Trapped fraction `n_trap(t̃) / N`.
- Radial density profile `ρ(r̃)` at the current scrub time (binned histogram).
- Time trace of the mean *radial* polarization of free cells,
  `⟨P·r̂⟩_free(t̃)` with `r̂ = X̃/r̃` the outward unit vector. This is signed:
  positive = the swarm is polarized outward (with the wave), negative = inward
  (the anti-wave regime). Preferred over `⟨|P|⟩` because the magnitude alone
  cannot distinguish the two, and the sign is the observable of interest here.
  A cell exactly at the origin has no defined `r̂` and contributes 0.
  The plot keeps `y = 0` in view with a dashed zero line.

## Sanity bounds and defaults

- Trajectory window `t̃ ∈ [0, t̃_max]`, where `t̃_max` is the integration end
  time — *not* just a plot range — entered as a typed number rather than a
  slider, since the useful range spans orders of magnitude and the cost is
  linear in it (`N_steps = ⌈t̃_max/dt̃⌉ + 1`). Default `t̃_max = 52`.
  For reference, the last wave fully clears the dish at
  `T̃_clear = (n_waves − 1) · Δt̃_wave + (R̃_dish + 8) / C`, shown as a hint
  under the field. Truncating before `T̃_clear` is fine — the tail of the last
  wave leaving the dish carries no dynamics of interest.
- `dt̃ ≤ 0.01` typical; reduce when `χ̃` or `ϑ` is large.
- `N = 1000` default; renderer (Canvas2D) supports up to ~10⁴.
- Time scrub: full trajectory (`x̃_i, ỹ_i, P_xi, P_yi, trapped_i` per step)
  is precomputed and stored; user scrubs via a time slider as in Setup 2.
  Memory is `~5 × 4 B × N × N_steps`, e.g. ~32 MB at `N = 1000`,
  `N_steps = 1500`.
