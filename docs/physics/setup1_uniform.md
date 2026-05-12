# Setup 1 — Polarization in a uniform cue

Single cell, uniform background concentration `L = L₀` (constant in space
and time, so `∇L = 0`). The cell's velocity is irrelevant here — we only
watch the polarization vector `p(t)` perform Brownian dynamics in the
6th-order GL free energy.

## Dimensional equation

From `ginzburg_landay_neutrophils.md` §2, with `∇L = 0`. Polarization is
either a scalar (1D) or 2D vector `p = (p_x, p_y)`:

```
dp_α/dt = r₀ (L - L_c) p_α + u |p|² p_α - w |p|⁴ p_α + √(2 θ L) ξ_α(t)
```

with independent white noises `⟨ξ_α(t) ξ_β(t')⟩ = δ_{αβ} δ(t - t')`.

## Nondimensionalization

Cue, polarization, and time scales:

```
𝓛 = L / L_c,   P = p / p_0,   p_0 = √(u/w),   t̃ = t / t_0,   t_0 = 1 / (r_0 L_c).
```

Two intrinsic dimensionless groups (cue-independent):

```
λ = u² / (w r_0 L_c)    — strength of the cubic+quintic relative to the linear term
ϑ = θ w / (u r_0)       — noise level
```

The nondim SDE is

```
dP_α/dt̃ = (𝓛 - 1) P_α + λ (|P|² P_α - |P|⁴ P_α) + √(2 ϑ 𝓛) η_α(t̃),
```

i.e. overdamped gradient flow `dP_α/dt̃ = -∂F̃/∂P_α + √(2 ϑ 𝓛) η_α` for
the free energy

```
F̃(P; 𝓛, λ) = -½ (𝓛 - 1) |P|² - ¼ λ |P|⁴ + ⅙ λ |P|⁶.
```

Three independent nondim knobs `{𝓛, λ, ϑ}`. Only `𝓛` carries the
physical cue concentration; `λ` and `ϑ` are intrinsic.

## What to plot

All axes nondimensional.

- Free-energy landscape `F̃(|P|; 𝓛, λ)` with current `|P|` marked.
- Time trace of `P` (1D mode) or `|P|` (2D mode).
- 2D mode only: scatter / trail of polarization vector tip in the
  `(P_x, P_y)` plane, with the radial well rings overlaid.
- Running histogram of `|P|` (steady state) with Boltzmann overlay
  `P_{ss} ∝ |P|^{d-1} exp(-F̃ / (ϑ 𝓛))`.
- Bifurcation diagram: extrema of `F̃` vs. `𝓛` (curves depend on `λ`),
  with the current `𝓛` marked.
- `⟨|P|⟩(𝓛)` from the Boltzmann distribution, both 1D and 2D.

## Sanity bounds

Bifurcation depends on the single combination `(𝓛 - 1) / λ`:

- `(𝓛 - 1) / λ > 0` (i.e. `𝓛 > 1`): polarized; `P = 0` is unstable.
- `(𝓛 - 1) / λ < -1/4` (i.e. `𝓛 < 1 - λ/4`): only `P = 0` survives.
- First-order coexistence window: `1 - λ/4 < 𝓛 < 1`. Polarized minimum at
  `|P|² = (1 + √(1 + 4(𝓛 - 1)/λ))/2`, unstable maximum at
  `|P|² = (1 - √(1 + 4(𝓛 - 1)/λ))/2`.
- Default `dt̃ = 0.01`.
