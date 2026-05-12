# Hill / linear-relaxation polarization model

An alternative to the Ginzburg–Landau (GL) polarization model used by default
in Setups 2 and 3. Both models are exposed in those setups via a UI toggle;
this doc covers the Hill model and its nondimensionalization. The doc is
shared between Setups 2 and 3 — the polarization equation is identical
(scalar in Setup 2, 2D vector in Setup 3); only the geometry of the cue
field differs.

## Dimensional equations

```
dp/dt = χ ∇L  −  ν p  +  √(2 θ L) ξ                       (polarization SDE)
dx/dt = v₀ · σ(p),    σ(p) = |p|ⁿ / (p₀ⁿ + |p|ⁿ) · p̂      (Hill velocity)
⟨ξ(t) ξ(t')⟩ = δ(t − t')                                  (white noise)
```

Differences from the GL model (`docs/physics/setup2_wave.md`):

- **No bistable well.** The polarization equation has only the linear
  relaxation term `−νp`. There is no GL 6th-order well, so no first-order
  phase transition and no free-energy `F(p)`.
- **Hill-saturating velocity.** Cell speed saturates at `v₀` for `|p| ≫ p₀`,
  with sharpness controlled by the Hill exponent `n`. The GL model uses the
  linear response `v = μ p`.
- **Same noise.** The noise term is identical to GL: multiplicative,
  `√(2 θ L) ξ`. We deliberately keep this form (rather than the additive
  `√(2 D_p) ζ` originally written for Hill in the literature) so that the
  noise slider `ϑ` is shared between the two models in the UI.

## Nondimensionalization (mirroring Setup 2)

Choose the same length scale `σ` (Gaussian wave width) and cue scale `L_c`
(threshold from GL — kept here as a user-chosen scale even though it has no
direct role in the Hill dynamics, so that the wave knob `M = L_max/L_c`
remains a free parameter and the dim ↔ nondim slider linkage matches GL
exactly). Time and polarization scales come from the model itself:

```
x̃   = x / σ
t̃   = ν t                    ← time scale 1/ν (replaces 1/(r₀ L_c))
P    = p / p₀                 ← Hill saturation polarization
𝓛   = L / L_c                (same as GL)
```

Substituting into the dimensional equations and grouping:

```
dP/dt̃   = χ̃ ∂_x̃ 𝓛  −  P  +  √(2 ϑ 𝓛) η(t̃)
dx̃/dt̃   = ṽ₀ · |P|ⁿ / (1 + |P|ⁿ) · P̂
```

with nondim groups

```
M   = L_max / L_c                       peak cue / threshold       (shared with GL)
C   = c / (σ ν)                         wave speed (σ per t₀)      (shared role; t₀ differs from GL)
χ̃   = χ L_c / (σ ν p₀)                  chemotactic coupling       (shared slider with GL)
ṽ₀  = v₀ / (σ ν)                        peak cell speed            (Hill-only)
ϑ   = θ w / (u r₀)  →  θ / ν            noise strength             (shared slider with GL; see note)
n   ∈ {1, 2, …, 10}                     Hill exponent              (Hill-only)
```

**Note on `ϑ`.** The dimensional groups for `ϑ` differ between models because
the underlying time and polarization scales differ. The slider exposes the
*nondim* value `ϑ` directly; the user does not see the dimensional
construction. The two models therefore share a common slider, with the
understanding that "the same `ϑ` value" means "the same nondim noise
strength" in both, not "the same physical `θ`".

## Comparison of nondim groups

| Quantity            | GL                                  | Hill                          | Slider behavior |
|---------------------|-------------------------------------|-------------------------------|-----------------|
| `M`, `C`            | as in Setup 2                       | as defined above              | shared          |
| `χ̃` (chemotactic)   | `χ/(σ r₀ √(u/w))`                   | `χ L_c/(σ ν p₀)`              | shared slider, same letter |
| `ϑ` (noise)         | `θ w/(u r₀)`                        | `θ/ν`                         | shared slider, same letter |
| GL motility `μ̃`    | `μ √(u/w)/(σ r₀ L_c)`               | —                             | GL-only |
| GL well `λ`         | `u² / (w r₀ L_c)`                   | —                             | GL-only |
| Hill peak speed `ṽ₀`| —                                   | `v₀/(σ ν)`                    | Hill-only |
| Hill exponent `n`   | —                                   | (already dimensionless)       | Hill-only |

## Adiabatic / no-inertia reference cell (Setup 2)

The "ghost" reference trajectory shown alongside the GL cell in Setup 2 has
a Hill analog: take the instantaneous balance `χ̃ ∂_x̃ 𝓛 − P_eq = 0`, giving
`P_eq = χ̃ ∂_x̃ 𝓛`, and feed this through the Hill velocity:

```
dx̃_ref/dt̃ = ṽ₀ · |P_eq|ⁿ / (1 + |P_eq|ⁿ) · sign(P_eq)
```

This is the gradient-only response with no relaxation lag — it should give
`Δx̃ ≈ 0` after one wave passage, isolating the rectification due to
polarization inertia.

## Setup 3 (2D radial swarm)

In Setup 3 the polarization is a 2D vector `P = (P_x, P_y)`; the SDE has the
same form per component, with `∂_x̃ 𝓛` replaced by `∇_x̃ 𝓛` and the noise
applied independently to each component. The Hill velocity is

```
dX̃/dt̃ = ṽ₀ · |P|ⁿ / (1 + |P|ⁿ) · P̂,    P̂ = P / |P|
```

with the convention that velocity is zero when `|P| = 0` (the prefactor
`|P|ⁿ/(1+|P|ⁿ)` vanishes faster than `1/|P|` for `n ≥ 1`).

The trap and dish boundary handling, multi-wave cue, and aggregate plots
(`⟨r̃⟩_free, n_trap/N, ρ(r̃)/r̃, ⟨|P|⟩_free`) are all model-agnostic and
shared between GL and Hill.
