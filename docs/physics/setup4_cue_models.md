# Setup 4 — Cue field model catalog

This document catalogs every PDE/ODE model for the LTB4 cue field L and
auxiliary inhibitor fields extracted from the literature and notes in
`Cue model summary/`. It is the source-of-truth for which models go into
Setup 4 and what their equations are.

**All derivations and model variants are reproduced here in the same
notation as the main physics spec** (`ginzburg_landay_neutrophils.md`):
L is LTB4 concentration, σ (or ρ) is surface cell density [cells/area],
a is emission rate per cell, L_c (= L_0 in many sources) is the relay
activation threshold.

---

## 1. Geometry variants

The same model equations arise in two geometries that differ in how many
spatial dimensions L occupies.

### 1a. 2D–2D geometry

Cells and L are both confined to a thin layer of height h. Sources appear
as volume densities ÷ h.

```
∂_t L = D_L ∇²_{2D} L  +  (a σ / h) f(L, ...)
```

The 2D Laplacian `∇²_{2D}` acts in the (x, y) plane.

### 1b. 2D–3D geometry

Cells sit on the floor z = 0; L diffuses in the 3D half-space z ≥ 0.
Cells are idealised as surface sources concentrated at z = 0:

```
∂_t L = D_L ∇²_{3D} L  +  a σ f(L, ...) δ(z)
```

The 3D Laplacian `∇²_{3D}` acts in (x, y, z). Cells measure L|_{z=0} and
∇_{2D} L|_{z=0}. Compared to 2D–2D the only change is the geometric factor
(1/h vs δ(z)) and the effective wave speed: for the basic relay without
inhibitor,

```
c_{2D-2D}  = √(a σ D_L / (h L_0))
c_{2D-3D}  = 2 a σ / (π L_0)
```

### Numerical treatment of 2D–3D geometry

Solve the full 3D PDE directly on a thin slab grid `N_x × N_y × N_z`.
No analytical reduction, no approximation.

**z-grid.** Use a non-uniform grid that is fine near z = 0 (where fields
are strong) and coarse at large z. An exponentially expanding spacing
works well:

```
Δz_j = h_0 · α^j,    j = 0, 1, ..., N_z - 1
```

where h_0 and α are chosen so the total slab height z_max is large enough
that the z = z_max boundary does not affect dynamics at z = 0. The
required slab height depends on how far L penetrates in z:
- Γ_L > 0: L decays exponentially with scale √(D_L/Γ_L); z_max of a few
  decay lengths suffices.
- Γ_L = 0: L decays algebraically; z_max must be verified by checking
  insensitivity of z = 0 results to z_max.

**Boundary conditions.**
- z = 0: Neumann ∂L/∂z = 0 (symmetry — we simulate z ≥ 0 by the
  image-method equivalence of the full space).
- z = z_max: Neumann ∂L/∂z = 0 (no flux out the top).

**Discretizing δ(z) — the critical step.**
The δ(z) source is implemented as a volumetric source term at z = 0 grid
nodes only, with coefficient `2/h_0` (not `1/h_0`):

```
source at z=0 node  =  f(L|_{z=0}) × 2/h_0
```

The factor 2 arises from the half-space control volume: the boundary node
at z = 0 has a control volume of thickness h_0/2 (half a cell, since the
boundary is at the edge). The δ(z) integrated over that control volume
is 1, so the discrete amplitude is 1/(h_0/2) = 2/h_0. An equivalent
derivation: integrate the PDE from z = -ε to z = +ε, use the full-space
symmetry ∂_z L|_{-ε} = -∂_z L|_{+ε}, to get the jump condition
`D_L ∂_z L|_{0^+} = f(L|_{z=0})/2`; implementing this as a Neumann flux
into the half-space is consistent with source amplitude 2/h_0 at the
first node.

For a non-uniform z-grid, h_0 is the first step size Δz_0.

**Time integration.** Crank-Nicolson (implicit) for the diffusion
operator is preferred over explicit Euler because the first z-step h_0
is much smaller than the lateral spacing, making the explicit CFL
`dt ≤ h_0²/6` prohibitively restrictive. Crank-Nicolson is
unconditionally stable for diffusion. The nonlinear source terms remain
explicit. The reference Python implementation (`la_2d3d_solver.py`) uses
this scheme with pre-factorized LU decomposition.

### 1c. Threshold smoothing (Hill regularisation)

All Heaviside step functions in the model equations are replaced by
smoothed Hill functions in the simulation:

```
H^+(x; x₀; n)  =  x^n / (x₀^n + x^n)    [activation: → 1 as x → ∞]
H^-(x; x₀; n)  =  x₀^n / (x₀^n + x^n)   [inhibition: → 0 as x → ∞]
```

Both equal 1/2 at x = x₀, and in the sharp limit n → ∞:

```
H^+(x; x₀; n)  →  Θ(x − x₀)
H^-(x; x₀; n)  →  Θ(x₀ − x)
```

Each Hill function invocation carries its own independent exponent n,
defaulting to 10. All are individually adjustable as UI knobs. In the
equations below every `Θ(x − x₀)` is written `H^+(x; x₀; n)` and
every `Θ(x₀ − x)` is written `H^-(x; x₀; n)`, with distinct n for
each invocation as listed in each model's parameter table.

---

## 2. Model M1: Basic relay, no inhibitor

**Name.** Diffusive relay (Dieterle et al. 2020).

**Geometry.** Both 2D–2D and 2D–3D. Equations given in 2D–2D form below.

```
∂_t L = D_L ∇² L  +  (a σ / h) H^+(L; L_0; n_L)  [- Γ_L L]
```

The degradation term `- Γ_L L` is optional (absent in the simplest
treatment).

**Source term.** Every cell that currently sees L ≥ L_0 emits LTB4 at
rate a (per cell). In the ABM, cells are discrete sources:

```
source(r, t) = a Σ_i H^+(L(r_i, t); L_0; n_L) δ(r - r_i)
```

**Parameters.**
```
D_L   diffusivity of LTB4   [µm²/s]
a     emission rate per cell [nM·µm²/s per cell, i.e. nM·µm³/s / (cell·µm)]
σ     cell surface density   [cells/µm²]
h     layer height (2D-2D)   [µm]
L_0   relay emission threshold [nM]   (distinct from the cell model's L_c)
Γ_L   LTB4 degradation rate  [1/s]   (optional)
n_L   Hill exponent for L activation gate  [unitless, default 10]
```

**Wave speed.** In the simplest (Γ_L = 0) traveling-wave limit the wave
speed is set by the balance of diffusion and relay production; see §1.

---

## 3. Model M2: Self-extinguishing — per-cell intracellular inhibitor R

**Name.** Self-extinguishing relay with cell-autonomous inhibitor
(Strickland, Pan, Ji, Amir, Weiner 2024 — "two-threshold model").

**Type of inhibitor.** R is a **per-cell** (intracellular) ODE variable,
not a diffusing field. Each cell i carries its own R_i ∈ [0, ∞).

**Equations (2D–2D, dimensional).**

```
∂_t L = D_L ∇² L  +  (a σ / h) H^-(R_i*; R_c; n_R) H^+(L; L_0; n_L)  -  Γ_L L

dR_i/dt = β H^+(L(r_i, t); L_r; n_{Lr})  -  γ R_i
```

where `R_i*` is the inhibitor value of the cell sitting at position r_i,
and `R_c` is the inhibition threshold (the value of R_i at which emission
is shut off). Cell i emits only when R_i < R_c (not yet inhibited) AND
the local cue exceeds the relay threshold (L > L_c). The inhibitor
accumulates whenever L exceeds a second, independent threshold L_r and
decays exponentially.

**Key feature.** The wave is self-extinguishing: the inhibitor builds up as
the wave passes, shuts off emission in the wake, and the wave collapses.
The termination radius is set by the ratio β/γ and the threshold L_r.

**Limitation noted in email.** This model seems very parameter-sensitive, and
the wave shape is a step-front (disk), not a pulse (ring), unless the
parameters are carefully tuned. It probably does not naturally reproduce
density-independent wave speed.

**Parameters.**
```
R_c      inhibitor threshold for emission shutoff  [same units as R_i]
β        inhibitor production rate                 [1/s]
γ        inhibitor degradation rate                [1/s]
L_r      second activation threshold               [nM]
n_L      Hill exponent for H^+(L; L_0) in L equation      [default 10]
n_R      Hill exponent for H^-(R_i; R_c) in L equation    [default 10]
n_{Lr}   Hill exponent for H^+(L; L_r) in R_i ODE         [default 10]
```

---

## 4. Model M3: ROS diffusing extracellular inhibitor

**Name.** ROS degradation model (proposed in email to Orion & Nitya).

**Type of inhibitor.** R is an **extracellular diffusing field**. ROS
(H₂O₂) is produced by NADPH oxidase, diffuses fast (D_ROS ≈ 2000 µm²/s),
and is converted by MPO to HOCl, which degrades LTB4 bimolecularly.

**Physical picture.** Because D_R >> D_L (ratio ≈ 10–16), R spreads
rapidly and acts nearly as a global scalar over the wave lengthscale.
This fast-diffusion limit gives robust self-extinguishing behavior.

**Equations (2D–2D, dimensional).**

```
∂_t L = D_L ∇² L  +  (a σ / h) H^+(L; L_0; n_L)  -  k L R

∂_t R = D_R ∇² R  +  (b σ / h) H^+(L; L_r; n_{Lr})  -  Γ_R R
```

The L–R interaction is bimolecular: the ROS field R degrades LTB4 at
rate k L R.

**Parameters.**
```
D_R      ROS diffusivity            [µm²/s,  D_R/D_L ≈ 10–16]
b        ROS emission rate per cell [µm³/s per cell]
Γ_R      ROS degradation rate       [1/s]
L_r      ROS activation threshold   [nM]
k        bimolecular rate constant  [1/(nM·s)]
n_L      Hill exponent for H^+(L; L_0) in L equation       [default 10]
n_{Lr}   Hill exponent for H^+(L; L_r) in R equation       [default 10]
```

**Quorum sensing.** Rescaling r → r σ^{1/2} shows the arrest radius
scales as r_arrest ~ 1/√σ (2D–2D) or r_arrest ~ 1/σ (2D–3D).

---

## 5. Model M4: Adenosine diffusing extracellular inhibitor

**Name.** Adenosine inhibitor model (proposed in email).

**Type of inhibitor.** A is an **extracellular diffusing field**. Adenosine
is co-secreted by activated cells and, once accumulated above a threshold
A_c, inhibits LTB4 biosynthesis. The inhibition acts as a production gate
(`H^-(A; A_c)` factor on the source term), not as degradation.

**Equations (dimensional; 2D–3D form shown).**

```
∂_t L = D_L ∇²_{3D} L  +  a σ H^-(A; A_c; n_A) H^+(L; L_0; n_L) δ(z)

∂_t A = D_A ∇²_{3D} A  +  b σ H^+(L; L_0; n_{LA}) δ(z)
```

In 2D–2D replace `∇²_{3D} ... δ(z)` with `∇²_{2D} ... / h`.

**Key difference from M3.** Inhibition enters as a production switch
`H^-(A; A_c; n_A)` rather than a degradation rate, and there is no degradation
term for A shown (it accumulates globally until the wave collapses).

**Parameters.**
```
D_A      adenosine diffusivity          [µm²/s,  D_A/D_L ≈ 2–5]
b        adenosine emission rate/cell   [µm³/s per cell]
A_c      adenosine inhibition threshold [nM]
n_L      Hill exponent for H^+(L; L_0) in L equation       [default 10]
n_A      Hill exponent for H^-(A; A_c) in L equation       [default 10]
n_{LA}   Hill exponent for H^+(L; L_0) in A equation       [default 10]
```

Note: in the email and source simulations, A_c and L_c happen to share
the same normalized value after nondimensionalization, but they are
distinct dimensional parameters.

**Note.** The email reports a 2D–3D simulation showing self-extinguishing
solutions with reasonable termination radius in this model.

---

## 6. Model M5: Receptor desensitization via GRK2 (per-cell ODE G_i)

**Name.** GRK2 kinetic model / receptor desensitization (Kienle et al.
2021).

**Purpose.** This is an add-on to any of M1–M4 that makes the relay
activation threshold L_c depend on the history of L exposure at each cell.
It resolves the observation that wave speed is **density-independent**,
which contradicts c ~ σ/L_c when L_c is constant.

**Type of variable.** G is a **per-cell** (intracellular) ODE variable
tracking the degree of BLT1 receptor desensitization.

**Full kinetic scheme** (receptors in states u: unbound, b: bound, d:
desensitized):

```
d_t b_i = k_b L(r_i) u_i  -  (k_u + k_d) b_i
d_t d_i = k_d b_i  -  k_r d_i
u_i + b_i + d_i = N_R      (conservation)
```

**Adiabatic limit** (k_b L, k_u, k_d >> k_r; b and u equilibrate fast):

```
d_t G_i = k_GL L(r_i)  -  γ_G G_i
```

where G_i is proportional to the desensitized receptor fraction d_i.
The activation threshold of cell i becomes:

```
L_c(i, t) = g_L G_i
```

**Quorum-sensing consequence.** A basal LTB4 level `L_basal ∝ σ` (from
random cell activity) → at steady state `G_ss ∝ σ` → `L_c ∝ σ` →
wave speed c ∝ σ/L_c = const, independent of density. This mechanism
resolves the density-independent wave speed without changing the basic
relay structure.

**Parameters.**
```
k_GL         effective GRK2 rate constant
γ_G          receptor recycling / resensitization rate
g_L          proportionality of threshold to desensitized fraction
```

---

## 7. Model M6.1: Basal adenosine — density-sensing adaptive threshold

**Name.** Basal adenosine density sensor (proposed here to reproduce the
experimentally observed *density-independent wave speed*).

**Motivation.** M1–M4 all give a wave speed that grows with cell density
(`c ~ σ` in 2D–3D, `c ~ √σ` in 2D–2D; see §1, §2). Experimentally
(Strickland et al. 2024 system) the relay speed is *almost
density-independent*. A density-dependent trigger-wave speed is also the
textbook driver of the **streaming / channelisation instability** (Höfer &
Maini, *Phys. Rev. E* **56**, 2074, 1997): a local density bump propagates
the wave faster, cells chemotax in, density rises further. Our ABM
channelises for exactly this reason; real neutrophils do not. M6.1 removes
both symptoms with one mechanism.

**Type of inhibitor.** A is an **extracellular diffusing field** produced
*basally* — every cell emits adenosine at a constant rate independent of L
— and cleared by first-order decay. This is the key contrast with M4,
where adenosine production is **LTB4-triggered** (`H^+(L; L_r)` gate). A
also couples differently: it **shifts the relay threshold** (`L_0 → L_0 +
gA`) rather than multiplicatively gating production.

**Equations (2D–2D, dimensional).**

```
∂_t L = D_L ∇² L  +  (a σ / h) H^+(L; L_0 + g A; n_L)

∂_t A = D_A ∇² A  +  (b σ / h)  -  γ_A A
```

At the ABM level cells are discrete basal sources of A (no L-gate):

```
∂_t L = D_L ∇² L + (a/h) Σ_i δ_2D(r - r_i) H^+(L; L_0 + g A; n_L)
∂_t A = D_A ∇² A + (b/h) Σ_i δ_2D(r - r_i)  -  γ_A A
```

**Nondimensionalisation.** Using the intrinsic (density-free) scales

```
ℓ_0 = a/(L_0 D_L),   t_0 = ℓ_0²/D_L = a²/(L_0² D_L³),   A_0 = b/(h D_L),
```

so that L is measured in units of L_0 (threshold → 1) and A in units of
A_0, the dimensionless groups are

```
α = a/(h D_L L_0),   λ = g b/(h D_L L_0),   D = D_A/D_L,
Γ_A = γ_A a²/(L_0² D_L³) = γ_A t_0,
```

and the dimensionless equations are

```
∂_t̃ 𝓛 = ∇̃² 𝓛 + α Σ_i δ̃_2D(r̃ - r̃_i) H^+(𝓛; 1 + λ 𝓐; n_L),
∂_t̃ 𝓐 = D ∇̃² 𝓐 + Σ_i δ̃_2D(r̃ - r̃_i)  -  Γ_A 𝓐.
```

`λ = 0` recovers M1 exactly, so **λ is a single dial from "channelises"
(M1) to "density-independent"**. (In the Setup-4 GUI this threshold-shift
coupling is displayed as `λ_A` — code symbol `lam_A` — to avoid a clash with
the cell-side GL activation group `λ = r_0 t_0 L_0`; they are unrelated.) Because the length scale ℓ_0 is intrinsic
(no σ), density enters only through the nondim areal density `σ̃ = σ ℓ_0²`
(cells per ℓ_0²) — the same symbol used throughout Setup 4; this is the
control parameter for density sweeps. (Do not write `n` for density — `n`
is reserved for Hill exponents n_L, n_R, ….)

**What the mechanism does (mean field, uniform σ̃).**
The basal source gives a *steady, spatially uniform* adenosine tone that is
present everywhere — including ahead of the front:

```
𝓐_ss = σ̃ / Γ_A                (2D–2D)
𝓐_ss = σ̃ / √(D Γ_A)           (2D–3D, at z = 0; screened surface source)
```

so the threshold every cell must cross is

```
𝓛_c = 1 + λ 𝓐_ss  ∝  σ̃    (at large σ̃),
```

i.e. the *threshold* scales with density in step with the *production*
(α σ̃), instead of the production growing alone. That is the mechanism, and
that is as far as this catalog goes.

> **No closed-form c(σ̃).** Earlier versions of this document quoted a
> mean-field pushed-front flux balance (`c² 𝓛_c = α σ̃`) and the speed law
> that follows from it. That estimate does **not** describe this system
> reliably and has been removed rather than left in place with caveats. The
> density dependence of the relay speed is an **output**: measure it with the
> Setup-4 `c(σ̃)` sweep, which plots measured points only. The same applies
> to M4-style comparisons: the qualitative contrast below (an inhibitor that
> reads *density* vs one that reads the *wave*) is a statement about what the
> threshold does, not a derivation of a speed.

Contrast M4: with LTB4-triggered production the adenosine level at the
ignition point is set by the trigger threshold L_r and does *not* scale with
σ̃, so the effective threshold stays put as density rises. Only an inhibitor
that reads **density** (basal), not the **wave** (triggered), makes 𝓛_c
track σ̃.

**Anti-streaming condition.** The Höfer–Maini instability is driven by the
*density dependence* of the trigger-wave speed, so weakening that dependence
weakens the instability. But channels are driven by *local* density
fluctuations, so A must resolve them: the adenosine screening length

```
ℓ_A = √(D / Γ_A)   (nondim)
```

should be of order the relay length (≈ 1 = ℓ_0). If ℓ_A ≫ fluctuation
width, A acts as a near-global threshold shift — it adapts the *ensemble*
threshold but leaves a local bump faster than its surroundings, so channels
can survive. Tune the tone via λ (and Γ_A) and ℓ_A via D.

**Limitation — no self-termination.** Basal-only adenosine adapts the
threshold but gives no arrest: cells behind the front keep firing (the
threshold stays at 1 + λ𝓐_ss), so the wave is a propagating/filling front,
not a finite ring. For a finite radius as well, combine M6.1 (basal arm) with
a wake-accumulating term (M4-style triggered); the two act ahead of vs.
behind the front and do not interfere.

**Relation to M5.** M6.1 (paracrine, diffusing field) and M5 (cell-autonomous
GRK2 ODE) are two implementations of the *same* principle — an effective
threshold ∝ local density. M6.1 senses density
over ℓ_A via a shared field; M5 senses it per-cell via each cell's basal-L
integral.

**Parameters.**
```
D_L   LTB4 diffusivity            [µm²/s]
a     LTB4 emission rate per cell  [as in M1]
b     adenosine emission rate/cell [basal, L-independent]
g     threshold-shift coupling     [dim of L per unit A]  → nondim λ
γ_A   adenosine decay rate         [1/s]                  → nondim Γ_A
D_A   adenosine diffusivity        [µm²/s]                → nondim D
L_0   bare relay threshold         [nM]  (nondim → 1)
n_L   Hill exponent for H^+(L; 1+λA) gate  [default 10]
```

---

## 7b. Model M6.2: Quorum-throttled LTB4 production

**Name.** Quorum sensor on the *production rate* (companion to M6.1, which
puts the same density signal on the *threshold*).

**Motivation.** M6.1 achieves a density-independent speed by letting the
relay threshold rise with density. M6.2 is the minimal alternative: the
threshold stays at `L_0`, and the density signal `Q` throttles how much
LTB4 each cell makes. Both are "quorum sensing"; they differ in *where* the
density enters the front-speed balance, and — as shown below — they are
**mean-field degenerate at m = 1** but diverge for `m ≠ 1`. Running both is
therefore a discriminating experiment, not a duplication.

**Type of inhibitor.** `Q` is an **extracellular diffusing field** produced
*basally* (every cell, constant rate `b`, no L-gate — exactly as M6.1's `A`)
and cleared by first-order decay `γ_Q`. The coupling is **multiplicative on
production**: `a_0 → a_0 / (1 + (Q/Q_0)^m)`.

**Equations (2D–2D, dimensional, mean field).**

```
∂_t L = D_L ΔL + [a_0 / (1 + (Q/Q_0)^m)] (σ/h) Θ(L - L_0)
∂_t Q = D_Q ΔQ + (b σ / h) - γ_Q Q
```

At the ABM level cells are discrete sources (the `Q`-throttle and the
`Θ`-gate are both evaluated at the emitting cell's position `r_i`):

```
∂_t L = D_L ΔL + (a_0/h) Σ_i δ_2D(r - r_i) Θ(L(r_i) - L_0) / (1 + (Q(r_i)/Q_0)^m)
∂_t Q = D_Q ΔQ + (b/h)   Σ_i δ_2D(r - r_i)                 - γ_Q Q
```

**Nondimensionalisation.** Same intrinsic (density-free) scales as M1/M6.1,
with `Q` measured in units of the *physically fixed* half-saturation
constant `Q_0` (this is the key bookkeeping difference from M6.1, where the
unit `A_0 = b/(hD_L)` was chosen to absorb the emission rate `b`; here `b`
survives as an explicit group `β`):

```
ℓ_0 = a_0/(L_0 D_L),   t_0 = ℓ_0²/D_L = a_0²/(L_0² D_L³),
𝓛 = L/L_0,   𝓠 = Q/Q_0,
α = a_0/(h D_L L_0) = 1/h̃,   β = b/(h D_L Q_0),
D = D_Q/D_L,   γ = γ_Q t_0 = γ_Q a_0²/(L_0² D_L³),
```

giving

```
∂_t̃ 𝓛 = ∇̃²𝓛 + α Σ_i δ̃_2D(r̃ - r̃_i) H^+(𝓛;1;n_L) H^-(𝓠;1;m),
∂_t̃ 𝓠 = D ∇̃²𝓠 + β Σ_i δ̃_2D(r̃ - r̃_i) - γ 𝓠,
```

where `Θ(L−L_0) → H^+(𝓛;1;n_L)` is the usual Hill regularisation of the
step (repo rule: every threshold is Hill-smoothed) and, exactly,

```
1/(1 + 𝓠^m) ≡ H^-(𝓠; 1; m).
```

So the throttle is a standard negative Hill factor with threshold 1 and its
own independent exponent `m` (code key `m_Q`; the catalog's `n_X` convention
would call it `n_Q` — the note's `m` is kept as the display symbol).

**M1 limit.** `β = 0` recovers M1 exactly. Note `m = 0` does **not** —
`1/(1+𝓠^0) = 1/2`, a uniform halving of production. The single "off" dial
is `β`.

**What the mechanism does (mean field, uniform σ̃).** The basal source is
unconditional, so `𝓠` relaxes to a spatially uniform tone everywhere —
including both ahead of and behind the front:

```
𝓠_ss = β σ̃ / γ                (2D–2D)
𝓠_ss = β h̃ σ̃ / √(D γ)        (2D–3D, at z = 0; see the geometry note below)
```

Since `𝓠` is uniform, M6.2 is *exactly M1 with a rescaled source* and an
unshifted threshold (`𝓛_c = 1`):

```
α → α_eff(σ̃) = α / (1 + 𝓠_ss^m) = α / (1 + (β σ̃ / γ)^m).
```

That statement is exact. Its consequence for the *drive* per unit area,
`α_eff σ̃`, is then elementary:

| exponent | drive `α_eff σ̃` at high density | reading |
|----------|----------------------------------|---------|
| `m < 1`  | grows like `σ̃^{1-m}`             | throttle too soft to matter |
| `m = 1`  | saturates at `α γ / β`           | drive becomes density-independent |
| `m > 1`  | *falls* like `σ̃^{1-m}`           | drive is strongest at intermediate density |

`σ̃★ = γ/β` marks where the tone reaches 1 (the throttle turns on), and the
screening length is `ℓ_Q = √(D/γ)` as in M6.1.

> **No closed-form c(σ̃).** How the front speed follows from `α_eff σ̃` is
> *not* given here. The mean-field pushed-front flux balance (`c² 𝓛_c =
> α_eff σ̃`) that earlier versions of this section used does not describe
> this system reliably, and has been removed rather than kept with caveats.
> Treat the wave speed as an output: measure it with the Setup-4 `c(σ̃)`
> sweep (which plots measured points only), and compare *sweeps* — β = 0 vs
> β > 0, m = 1 vs m = 2 — rather than a sweep against a formula.

**Relation to M6.1 (important).** Both models feed the same basal,
density-proportional tone into the relay — one at the threshold, one at the
production rate — and at `m = 1` the density enters both the same way, under
the dictionary `λ_A ↔ β`, `Γ_A ↔ γ`. Sweeping both at `m = 1` is therefore a
cross-model consistency check (the two ABMs should behave alike). They are
distinguished by:

1. **`m ≠ 1`**: only M6.2 can make the drive *fall* with density (`m > 1`).
   A measured `c(σ̃)` that turns over at high density would point to M6.2
   with `m > 1`; M6.1's threshold shift cannot produce that.
2. **Failure mode at very high density.** M6.1 fails by *threshold*: cells
   never cross `1 + λ_A𝓐` and the relay simply does not ignite. M6.2 fails
   by *starvation*: cells still cross threshold and still emit, but too
   weakly to drag the next shell over — the front slows continuously to
   zero rather than switching off. Visually: M6.1 shows a dark dish with
   few emitters; M6.2 shows many weak emitters and a stalled front.
3. **Front structure.** M6.1 changes `𝓛_c` (the level the front must reach);
   M6.2 changes the amplitude behind the front. The saturated `𝓛` level in
   the wake is `∝ α_eff` in M6.2 but unchanged (`∝ α`) in M6.1 at fixed
   `σ̃` — measurable on the `𝓛(r̃)` radial profile.

**2D–3D geometry.** Both fields diffuse in the half-space with the cells as
a `δ(z)` surface source at `z = 0`. The nondim source coefficient of a
surface source is `h̃ ×` its 2D–2D value (same relation as the L source:
`α = 1/h̃` in 2D–2D vs `1` in 2D–3D), so

```
β_3D = β h̃ = b/(D_L Q_0 ℓ_0),      𝓠_ss(z=0) = β_3D σ̃ / √(D γ).
```

The effective source is throttled by the same factor `1/(1 + 𝓠_ss^m)`
evaluated with this tone; the resulting front speed, as in 2D–2D, is left to
the sweep. The `β` slider stays the 2D–2D group of the note; the `h̃` factor is applied
internally when the geometry is 2D–3D. (M6.1 hides this factor by
redefining its arbitrary unit `A_0` per geometry — legitimate there, not
available here because `Q_0` is a physical constant.)

**Anti-streaming.** Same requirement as M6.1: `ℓ_Q = √(D/γ)` must be of
order the relay length (≈ 1 = ℓ_0) so that `Q` resolves *local* density
bumps. The hoped-for advantage over M6.1 is that with `m > 1` the drive
itself falls with density, so — if that carries through to the front speed —
a local density excess would propagate *slower* than its surroundings and
the Höfer–Maini feedback would change sign rather than merely weaken. That
is a prediction to test with the sweep, not a result.

**Limitation — no self-termination.** As in M6.1, basal-only production of
`Q` gives no arrest: the wave is a propagating/filling front, not a finite
ring (throttled but never switched off, since `H^-` is positive for all
finite `𝓠`). Combine with a triggered/wake inhibitor (M4-style) for a
finite radius.

**Parameters.**
```
D_L   LTB4 diffusivity              [µm²/s]
a_0   LTB4 emission rate per cell   [as in M1]                  → α = 1/h̃
b     quorum-signal emission rate   [basal, L-independent]      → β
Q_0   quorum half-saturation const  [conc.]  (nondim → 1)
γ_Q   quorum-signal decay rate      [1/s]                       → γ
D_Q   quorum-signal diffusivity     [µm²/s]                     → D = D_Q/D_L
L_0   relay threshold               [nM]  (nondim → 1, unshifted)
m     Hill exponent of the production throttle H^-(𝓠;1;m)   [default 2]
n_L   Hill exponent of the relay gate H^+(𝓛;1;n_L)          [default 10]
```

---

## 8. Summary table

| Model | L equation | Inhibitor field | Inhibitor type | Wave self-extinguishing | Quorum sensing |
|-------|-----------|-----------------|----------------|------------------------|----------------|
| M1 basic relay | diffusion + relay prod [- decay] | none | — | no (or trivial) | no |
| M2 per-cell R | diffusion + gated prod | R_i (per-cell ODE) | intracellular | yes | no (param-sensitive) |
| M3 ROS diffusing | diffusion + H^+(L;L_0;n_L) prod - k L R | R (PDE, D_R >> D_L) | extracellular, fast | yes | yes (r ~ 1/√σ) |
| M4 adenosine | diffusion + H^-(A;A_c;n_A)H^+(L;L_0;n_L) prod | A (PDE, D_A ~ D_L) | extracellular, slow | yes | yes (r ~ 1/σ) |
| M5 GRK2 | adds to any of above | G_i (per-cell ODE) | modifies threshold | (inherited) | yes (c ~ const) |
| M6.1 basal adenosine | diffusion + threshold-shifted relay H^+(𝓛;1+λ𝓐;n_L) | A (PDE, **basal** + decay) | extracellular, density-sensing | no (filling front) | **yes (threshold 𝓛_c ∝ σ̃, tracking production)** |
| M6.2 quorum-throttled production | diffusion + H^+(𝓛;1;n_L)·H^-(𝓠;1;m) prod, threshold unshifted | Q (PDE, **basal** + decay) | extracellular, density-sensing | no (filling front) | **yes (source α→α/(1+(βσ̃/γ)^m); drive saturates at m=1, falls for m>1)** |

Each `H^±` carries its own independent exponent n; see each model's parameter
table for the full list of n symbols.

---

## 9. Combinations

The models are not mutually exclusive:
- **M2 + M5**: self-extinguishing + density-independent speed
- **M3 + M5**: fast-diffusing ROS + density-independent speed
- **M4 + M5**: adenosine + density-independent speed
- **M6.1 + M4** (or any triggered/wake inhibitor): density-independent
  **speed** (M6.1 basal arm, ahead of front) + finite arrest **radius**
  (triggered arm, in the wake). The two arms decouple by location.
- **M6.1 vs M5**: redundant *mechanisms* for the same effect (paracrine field
  vs. cell-autonomous ODE); useful to run both and compare which better
  matches GRK-inhibitor / adenosine-deaminase (ADA) perturbation data.
- **M6.1 vs M6.2**: same basal density signal, different point of action
  (threshold vs production rate). Density enters both the same way at `m = 1`
  (`λ_A ↔ β`); distinguished by `m ≠ 1` (only M6.2 can make the drive fall
  with density), by the wake amplitude of 𝓛, and by the high-density failure
  mode (ignition failure vs front starvation). See §7b.
- **M6.2 + M6.1**: not exclusive — a single basal field could both shift the
  threshold and throttle production; not planned for the GUI (two dials for
  one measurement).

In the GUI, M5 should be implemented as an optional add-on toggle
independent of the choice of M1–M4. M6.1 is a standalone model (its basal-A
field replaces M4's triggered-A field), but its threshold-shift coupling
can also be layered onto M1–M4.

---

## 10. Supported combinations for Setup 4

For the initial implementation, we propose:

| Tier | Models | Geometry | Justification |
|------|--------|----------|---------------|
| Core | M1 | 2D–2D, 2D–3D | simplest, validates PDE solver |
| Tier 1 | M2 | 2D–2D | "old model", per-cell R, easy to add |
| Tier 1 | M3 | 2D–2D | new/preferred model, diffusing R |
| Tier 1 | M6.1 | 2D–2D, 2D–3D | density-adapted threshold; anti-streaming candidate; only M1 + one field + threshold shift. Tone σ̃/Γ_A (2D–2D), σ̃/√(D Γ_A) at z=0 (2D–3D). Speed measured, not predicted. |
| Tier 1 | M6.2 | 2D–2D, 2D–3D | quorum throttle on production; same solver as M6.1 (M1 + one basal field), differs only in how the field enters the per-cell source weight. Effective source α/(1+𝓠_ss^m); tone βσ̃/γ (2D–2D), β h̃ σ̃/√(Dγ) at z=0 (2D–3D). Speed measured, not predicted. |
| Tier 2 | M4 | 2D–3D | adenosine, biologically motivated |
| Optional | M5 | any | add-on toggle on top of M1–M4 |

M6.1 is a high-priority addition: it is the minimal change to M1 (one extra
diffusing field + a threshold shift) aimed at the density-independent
speed / channelisation mismatch between simulation and experiment.

M4 in 2D–2D is also valid but the 2D–3D version appears in the source
documents and is the experimentally relevant geometry for that model.

### Target (sticking boundary) — orthogonal to the cue model

Every full-swarm page (M1, M2, M6.1, M6.2) additionally carries a
*"stick to the target"* checkbox (ON by default) and a nondim radius
`R̃_target` (default 2): a circle at the dish centre representing the pathogen
the swarm converges on. A cell that reaches it is pinned onto the circle
permanently, but keeps polarizing and keeps emitting the cue — so the cue-model
equations of this catalog are untouched by it. Full specification:
[setup4_swarm3d.md](setup4_swarm3d.md) §2c.
