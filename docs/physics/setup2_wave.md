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

## What to plot

- Live 1D snapshot of `𝓛(x̃, t̃)` with the cell as a marker.
- Free-energy `F̃(P; 𝓛(x̃_cell, t̃), λ)` snapshot with current `P` marked.
- Time traces of `x̃(t̃)`, `P(t̃)`, and `∂_x̃ 𝓛` at the cell.
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
