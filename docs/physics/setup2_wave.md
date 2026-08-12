# Setup 2 — Cell vs. 1D Gaussian running wave

Single cell at position `x(t)` in 1D. Cue field is a prescribed traveling
Gaussian:

```
L(x, t) = L_max · exp( -(x - c t)² / (2 σ²) )
∇L     = -L · (x - c t) / σ²
```

Polarization obeys the same SDE as Setup 1 plus the chemotactic drive
`χ ∇L`. Cell velocity is the linear response `v = μ p`.

## Dimensional equations

```
dp/dt = χ ∇L + r₀(L - L_c) p + u p³ - w p⁵ + √(2 θ L) ξ
dx/dt = μ p
```

## Nondimensionalization

Same notation as Setup 1, with one addition: distance is rescaled by the
wave width `σ`. **Time, polarization, and cue scales are identical to
Setup 1**, so `λ` and `ϑ` carry over unchanged.

```
𝓛   = L / L_c              (= M · exp(-(x̃ - C t̃)²/2)  for the wave)
P    = p / p_0,    p_0 = √(u/w)
t̃    = t / t_0,    t_0 = 1 / (r_0 L_c)        ← same as Setup 1
x̃    = x / σ                                  ← new, wave-specific
λ    = u² / (w r_0 L_c)                       ← same as Setup 1
ϑ    = θ w / (u r_0)                          ← same as Setup 1
```

Two new wave-related dimensionless groups:

```
M  = L_max / L_c                               peak cue / threshold
C  = c / (σ r_0 L_c)                           nondim wave speed (σ per t_0)
```

Two new coupling groups:

```
χ̃  = χ / (σ r_0 √(u/w))                       chemotactic strength
μ̃  = μ √(u/w) / (σ r_0 L_c)                   nondim cell motility
```

The cue field becomes

```
𝓛(x̃, t̃) = M · exp(-(x̃ - C t̃)² / 2)
∂_x̃ 𝓛   = -𝓛 · (x̃ - C t̃)
```

and the nondim equations of motion are

```
dP/dt̃   = χ̃ ∂_x̃ 𝓛 + (𝓛 - 1) P + λ (|P|² P - |P|⁴ P) + √(2 ϑ 𝓛) η
dx̃/dt̃   = μ̃ P
```

i.e. the polarization SDE has the same form as Setup 1 with the cue `𝓛`
now spatially and temporally varying, plus the chemotactic drive
`χ̃ ∂_x̃ 𝓛` and the deterministic cell motion `dx̃/dt̃ = μ̃ P`.

### Chemotactic-drive diagnostic `s`

For visualization we name the chemotactic drive term

```
s ≡ χ̃ ∂_x̃ 𝓛
```

so the polarization SDE reads `dP/dt̃ = s + (𝓛-1)P + λ(...) + noise`. As the
wave passes, `s` rises, reverses sign (at the wave crest, where `∂_x̃ 𝓛 = 0`),
then falls back — driving `P` around a loop. Plotting the parametric curve
`(P(t̃), s(t̃))` exposes this **phase trajectory in `P`–`s` space**: a pure
gradient sensor collapses onto a single-valued curve, whereas the GL cell's
polarization inertia opens the loop into a hysteresis cycle whose enclosed
area is the geometric signature of rectification (`Δx̃ ≠ 0`).

### Reduced-drive variable `q` and response ratio `U`

For a rigidly traveling cue, `∂_t̃ 𝓛 = -C ∂_x̃ 𝓛`, so at a *fixed* point in space
the chemotactic drive can be absorbed into a shifted polarization

```
q ≡ P + χ̃ 𝓛 / C
```

whose evolution carries no explicit forcing: `dq/dt̃` is just the *relaxational*
part of the SDE. The companion diagnostic is the ratio of the cell's velocity to
that relaxational rate,

```
U ≡ μ̃ P / [ (𝓛 - 1) P + λ (P² - P⁴) P ]
  = μ̃ / [ (𝓛 - 1) + λ (P² - P⁴) ]        (P cancels; used in code, no 0/0 at P = 0)
```

i.e. how much ground the cell covers per unit of relaxation of the reduced
variable. Following the cell, `d𝓛/dt̃ = ∂_x̃𝓛 · (dx̃/dt̃ - C)`, so `U = dx̃/dq`
holds exactly only when the cell's own drift is slow against the wave
(`μ̃P ≪ C`); the deviation from the measured `dx̃/dq` is itself a readout of how
much the cell is surfing the wave. `U` diverges where the bracket crosses zero (the cell sits at a
marginal point of the effective potential, `F̃'' = 0`), so the `U`–`q` plot is
drawn on a robust (percentile-clipped) `U` range. In Hill mode the relaxational
part is `-P`, so `U = ṽ₀ |P|ⁿ/(1+|P|ⁿ) · sign(P) / (-P)` and `λ, μ̃` do not enter.

### The singular curve `Σ`

Write the relaxational (autonomous internal-response) part of the drift as

```
f(P, 𝓛) ≡ (𝓛 - 1) P + λ (P³ - P⁵)        [GL]
        ≡ -P                              [Hill]
```

so the SDE reads `dP/dt̃ = s + f(P, 𝓛) + noise` and, since `q` absorbs `s`,
`q̇ = f(P, 𝓛)`. The **singular curve**

```
Σ = {(𝓛, P) : f(P, 𝓛) = 0}
```

is the instantaneous nullcline of the autonomous internal response. It is where
`U` diverges, and crossing it flips the sign of `q̇`, cutting each cycle into
alternating `q̇ > 0` and `q̇ < 0` intervals. For GL it factors exactly,

```
f = P · [𝓛 - 𝓛_Σ(P)],      𝓛_Σ(P) = 1 - λ P² + λ P⁴
```

so `Σ` has two branches — `P = 0`, and the explicit curve `𝓛 = 𝓛_Σ(P)` (no root
finding: `𝓛_Σ(0) = 1`, dipping to `1 - λ/4` at `P = ±1/√2`, back to `1` at
`P = ±1`). The factorization also gives the sign directly,
`sign(q̇) = sign(P) · sign(𝓛 - 𝓛_Σ(P))`. For Hill, `Σ` is just `P = 0` and
`sign(q̇) = -sign(P)`.

`Σ` is drawn on the `𝓛`–`P` phase portrait rather than being clipped away, with
the two sign regions faintly tinted. (Setup 1's `cv-bif` panel — extrema of `F̃`
vs `𝓛` — is the same curve: `dF̃/dP = 0` ⇔ `f = 0`.)

## What to plot

- Live 1D snapshot of `𝓛(x̃, t̃)` with the cell as a marker.
- Free-energy `F̃(P; 𝓛(x̃_cell, t̃), λ)` snapshot with current `P` marked.
- Time traces of `x̃(t̃)`, `P(t̃)`, and `∂_x̃ 𝓛` at the cell.
- On the `P(t̃)` trace, a gray background showing the **adiabatic
  (quasi-static) locus** — *every* root of `dF̃_eff/dP = 0`,

  ```
  χ̃ ∂_x̃𝓛 + (𝓛 - 1) P + λ (P³ - P⁵) = 0
  ```

  evaluated at the cell's instantaneous cue `(𝓛, ∂_x̃𝓛)`. This is a
  quintic in `P`, so the locus is **multivalued**: over the folded
  interval it has up to three stable branches (minima, `F̃''_eff > 0`,
  drawn solid) separated by unstable barriers (maxima, drawn dashed).
  The whole S-curve is drawn — branches are born and annihilate in pairs
  at saddle-node folds; the curve is *not* collapsed to a single
  continuation branch, so no artificial jump appears where a branch
  loses stability.

  The gap between `P(t̃)` and the locus is exactly the polarization
  inertia that rectifies the wave; in the adiabatic limit the trace
  rides one of the stable branches and `Δx̃ → 0`. In Hill mode there is
  no `F̃`; the analogous quasi-static balance `χ̃ ∂_x̃𝓛 − P = 0` is
  single-valued and is drawn as one solid gray curve.
- Phase trajectory in `P`–`s` space, `s = χ̃ ∂_x̃ 𝓛`.
- Response hysteresis in `𝓛`–`P` space: the parametric loop
  `(𝓛(x̃_cell, t̃), P(t̃))` — the cell's polarization plotted against the
  instantaneous cue it sees. As the wave passes, `𝓛` rises from `0` to `M`
  and back; a memoryless sensor would trace a single-valued curve, whereas
  the GL cell's polarization inertia opens the trace into a hysteresis loop.
  The loop's enclosed area is a direct visual signature of rectification:
  the response on the rising flank differs from the response on the falling
  flank, so the net drive over one passage does not cancel (`Δx̃ ≠ 0`).
  The singular curve `Σ` is overlaid on this panel, with the `q̇ > 0` / `q̇ < 0`
  regions faintly tinted.
- Phase trajectory in `U`–`q` space, `q = P + χ̃𝓛/C`, `U = μ̃P/[(𝓛-1)P + λ(P²-P⁴)P]`.
- Net displacement `Δx̃` after one wave passage as a function of `χ̃`.
- Comparison overlay: gradient-only model (no polarization inertia) vs.
  full GL model — the former should give `Δx̃ ≈ 0`.

## Alternative model (UI toggle)

Setup 2's UI exposes a model toggle: the default GL polarization SDE above,
or a Hill / linear-relaxation alternative. See
[hill_model.md](hill_model.md) for the Hill nondim equations. `χ̃` and `ϑ`
sliders are shared between models; `μ̃, λ` are GL-only; `ṽ₀, n` are
Hill-only. The free-energy panel `F̃(P)` is hidden in Hill mode (Hill
dynamics is not a gradient flow).

## Sanity bounds

- Trajectory window `t̃ ∈ [-X_init/C, +X_init/C]` so the wave starts at
  `x̃ = -X_init` and arrives at the cell (`x̃ = 0`) at `t̃ = 0`. Default
  `X_init = 8`; with `C = 1` the window is `[-8, +8]` (twice the
  wave–cell transit time).
- `dt̃ ≤ 0.01` typically; reduce when `χ̃` or `ϑ` is large.
- Box length: simulate on `x̃ ∈ [-X, X]` with `X ≳ 10` so the wave is
  effectively absent at the boundaries.
