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

## 7. Summary table

| Model | L equation | Inhibitor field | Inhibitor type | Wave self-extinguishing | Quorum sensing |
|-------|-----------|-----------------|----------------|------------------------|----------------|
| M1 basic relay | diffusion + relay prod [- decay] | none | — | no (or trivial) | no |
| M2 per-cell R | diffusion + gated prod | R_i (per-cell ODE) | intracellular | yes | no (param-sensitive) |
| M3 ROS diffusing | diffusion + H^+(L;L_0;n_L) prod - k L R | R (PDE, D_R >> D_L) | extracellular, fast | yes | yes (r ~ 1/√σ) |
| M4 adenosine | diffusion + H^-(A;A_c;n_A)H^+(L;L_0;n_L) prod | A (PDE, D_A ~ D_L) | extracellular, slow | yes | yes (r ~ 1/σ) |
| M5 GRK2 | adds to any of above | G_i (per-cell ODE) | modifies threshold | (inherited) | yes (c ~ const) |

Each `H^±` carries its own independent exponent n; see each model's parameter
table for the full list of n symbols.

---

## 8. Combinations

The models are not mutually exclusive:
- **M2 + M5**: self-extinguishing + density-independent speed
- **M3 + M5**: fast-diffusing ROS + density-independent speed
- **M4 + M5**: adenosine + density-independent speed

In the GUI, M5 should be implemented as an optional add-on toggle
independent of the choice of M1–M4.

---

## 9. Supported combinations for Setup 4

For the initial implementation, we propose:

| Tier | Models | Geometry | Justification |
|------|--------|----------|---------------|
| Core | M1 | 2D–2D, 2D–3D | simplest, validates PDE solver |
| Tier 1 | M2 | 2D–2D | "old model", per-cell R, easy to add |
| Tier 1 | M3 | 2D–2D | new/preferred model, diffusing R |
| Tier 2 | M4 | 2D–3D | adenosine, biologically motivated |
| Optional | M5 | any | add-on toggle on top of M1–M4 |

M4 in 2D–2D is also valid but the 2D–3D version appears in the source
documents and is the experimentally relevant geometry for that model.
