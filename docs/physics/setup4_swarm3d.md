# Setup 4 — Emergent-wave swarming with dynamic cue field

> Replaces the earlier placeholder (which only noted the geometry and
> listed open questions). Model equations now specified; see
> [setup4_cue_models.md](setup4_cue_models.md) for the full cue-field
> model catalog.

---

## 1. What Setup 4 adds relative to Setup 3

In Setup 3 the wave is prescribed: a Gaussian pulse with fixed peak M and
speed C, launched by hand. Setup 4 makes the wave **emergent**: L satisfies
its own PDE, and the wave nucleates and propagates because activated cells
release LTB4 which triggers their neighbours. The cell polarization model
is unchanged (GL or Hill, as in Setups 2–3); the only new physics is the
`∂_t L` equation and the optional inhibitor fields.

Consequence: the wave amplitude, width, and speed are **outputs** of the
simulation, not inputs. Setup 4 therefore does not carry the wave
parameters M and C as sliders. Instead it carries the relay parameters
(emission rate Ã, diffusivity D̃_L, etc.).

---

## 2. Geometry

Two geometry variants are offered via a UI toggle:

### 2a. 2D–2D

Cells and L both live in a 2D disk of radius R_dish.
Source term ∝ 1/h (layer height).

### 2b. 2D–3D

Cells sit on the floor z = 0 of a 3D half-space; L diffuses in 3D.
The simulation solves the full 3D PDE on a thin slab grid `N_x × N_y × N_z`
with Crank-Nicolson for diffusion (explicit nonlinear source) — no
analytical reduction or approximation. Cell sources enter as a
volumetric term at the z = 0 nodes with coefficient `2/h_0` (equivalent
to a Neumann flux into the half-space; see `setup4_cue_models.md` §1b). The required slab
thickness depends on how far L penetrates in z (see `setup4_cue_models.md`
§1b for details).

---

## 3. Cell dynamics

Identical to Setup 3: a UI toggle selects between the GL and Hill
polarization models. The cue field 𝓛 and its gradient are computed by
the PDE solver rather than read from a formula; everything else is
unchanged.

**GL model:**
```
dP_i/dt̃  =  χ̃ ∇̃𝓛(r̃_i, t̃)  +  (𝓛(r̃_i, t̃) - 1) P_i
             +  λ (|P_i|² - |P_i|⁴) P_i
             +  √(2 ϑ 𝓛(r̃_i, t̃)) η_i

dx̃_i/dt̃  =  μ̃ P_i
```

**Hill model:** see [hill_model.md](hill_model.md). Shared knobs: χ̃, ϑ.
GL-only: λ, μ̃. Hill-only: ṽ₀, n.

---

## 4. Cue field models

Four cue field models are available (see `setup4_cue_models.md` for
full equations):

| Toggle label | Document label | Inhibitor | Key parameter(s) added |
|---|---|---|---|
| "Relay only" | M1 | none | [Γ̃_L]  *(M1 nondim is parameter-free by construction; c* absorbs a, σ, D_L, L_0, h)* |
| "Per-cell inhibitor" | M2 | R_i (intracellular ODE) | β̃, γ̃, L̃_r |
| "ROS (fast diffuser)" | M3 | R (extracellular PDE) | D̃_R, b̃, Γ̃_R, L̃_r, k̃ |
| "Adenosine (slow diffuser)" | M4 | A (extracellular PDE) | D̃_A, b̃, Ã_c |

An optional toggle "GRK2 desensitization" adds the per-cell threshold
modulation (M5) on top of any of the above.

The inhibitor type determines whether the simulation carries an extra
PDE grid (extracellular) or extra per-cell scalars (intracellular).

---

## 5. Nondimensionalization

### 5.1 Guiding principle: scales come from the L-wave dynamics

In Setups 2–3 the wave was prescribed, so its width σ was a free input
used as the length scale. Here the wave is emergent — σ is not an input.
The natural scales are instead set by the M1 wave speed.

**Cue unit: L_0** (relay threshold). This makes the relay threshold
`H^+(L; L_0; n_L)` become `H^+(𝓛; 1; n_L)` in nondim, so the M1 L
equation carries only the Hill exponent n_L (see below). Note `L_0 ≠ L_c` in general; the
ratio `𝓛_c = L_c / L_0` appears as a new dimensionless parameter in the
cell equation.

**Length and time scales from c\* = D_L / ℓ_0:**
```
ℓ_0 = D_L / c*,    t_0 = D_L / c*²  (= ℓ_0 / c*)
```
where c\* is the parametric M1 wave speed (numerical prefactors such as
2/π are part of the *actual* wave speed but are not included in c\*,
keeping the nondim clean):

```
2D–2D:  c*  =  √( a σ D_L / (h L_0) )       →  ℓ_0 = √( D_L h L_0 / (a σ) )
2D–3D:  c*  =  a σ / L_0                     →  ℓ_0 = D_L L_0 / (a σ)
```

These are different because 2D-2D depends on h (layer height) while
2D-3D does not. In the actual simulation the measured wave speed will
be c\* multiplied by a numerical prefactor (1 for 2D-2D and 2/π for
2D-3D, from the traveling-wave solutions).

**M1 L equation in nondim units — continuum form (one parameter: n_L):**
```
2D–2D:  ∂_t̃ 𝓛 = ∇̃²_{2D} 𝓛  +  H^+(𝓛; 1; n_L)
2D–3D:  ∂_t̃ 𝓛 = ∇̃²_{3D} 𝓛  +  H^+(𝓛; 1; n_L) δ̃(z̃)
```

**Discrete-ABM (per-cell δ-source) form** — what the simulation actually
integrates:
```
2D–2D:  ∂_t̃ 𝓛 = ∇̃²_{2D} 𝓛  +  (1/σ̃) Σ_i H^+(𝓛_i; 1; n_L) δ̃(x̃ − x̃_i)
2D–3D:  ∂_t̃ 𝓛 = ∇̃²_{3D} 𝓛  +  (1/σ̃) Σ_i H^+(𝓛_i; 1; n_L) δ̃(x̃ − x̃_i) δ̃(z̃)
```
with σ̃ = σ · ℓ_0² (dim cell density in nondim area units). Note that **N
appears only inside the sum** — the per-cell emission prefactor 1/σ̃ is
fixed by dim params and is independent of how many discrete cells N are
simulated. In the continuum limit `N → σ̃ · A_dish`, the discrete form
reduces to the continuum H^+ source. With N exceeding the dim-implied
count, total per-area emission scales with N — as for a physically denser
cell layer.

Inhibitor fields and degradation introduce additional dimensionless
parameters (see §5.3).

**M2 L equation in nondim units (per-cell inhibitor R_i).** With R̃_i = R_i / R_c
and the discrete-ABM δ-source form:
```
2D–2D: ∂_t̃ 𝓛 = ∇̃²_{2D} 𝓛
              + (1/σ̃) Σ_i H^-(R̃_i; 1; n_R) H^+(𝓛_i; 1; n_L) δ̃(x̃ − x̃_i)
              − Γ̃_L 𝓛
2D–3D: ∂_t̃ 𝓛 = ∇̃²_{3D} 𝓛
              + (1/σ̃) Σ_i H^-(R̃_i; 1; n_R) H^+(𝓛_i; 1; n_L) δ̃(x̃ − x̃_i) δ̃(z̃)
              − Γ̃_L 𝓛
```
Per-cell inhibitor ODE (same in both geometries):
```
dR̃_i / dt̃ = β̃ H^+(𝓛(r̃_i); L̃_r; n_{Lr})  −  γ̃ R̃_i
```
The PDE solver structure is identical to M1 — only the per-cell source
weight changes from `H^+(𝓛_i;1;n_L)` to the gated product
`H^-(R̃_i;1;n_R) · H^+(𝓛_i;1;n_L)`. R̃_i lives on the agents (no extra
grid), updated by explicit Euler at the agent dt̃.

### 5.2 Cell dynamics in L-wave nondim

Polarization scale: `p_0 = √(u/w)` (unchanged from Setups 1–3).

With `t̃ = t/t_0`, `x̃ = x/ℓ_0`, `𝓛 = L/L_0`, the GL SDE becomes:

```
dP_i/dt̃  =  χ̃ ∇̃𝓛  +  Λ [(𝓛 − 𝓛_c) P_i + λ (|P_i|² − |P_i|⁴) P_i
             +  √(2 ϑ 𝓛) η_i ]

dx̃_i/dt̃  =  μ̃ P_i
```

New geometry-dependent dimensionless groups:

```
Λ   =  t_0 r_0 L_c  =  r_0 L_c D_L / c*²   (wave time / cell activation time)
𝓛_c =  L_c / L_0                             (cell threshold / relay threshold)
χ̃   =  χ L_0 / (p_0 ℓ_0 / t_0)  =  χ L_0 / (p_0 c*)
μ̃   =  μ p_0 / c*
```

`Λ` is new: it encodes how fast the cell responds relative to how fast
the wave propagates. Large Λ means cells polarize quickly on the wave
timescale (adiabatic); small Λ means cells lag behind the wave.

The intrinsic cell parameters λ and ϑ are unchanged:
```
λ  =  u² / (w r_0 L_c)
ϑ  =  θ w / (u r_0)
```
These are independent of geometry and length scale.

### 5.3 Inhibitor nondim groups

**Convention:** all concentration fields (L, R, A) are normalized by L_0.
Rates are normalized by t_0. Diffusivity ratios by D_L.

| Symbol | Meaning | Formula | Model |
|---|---|---|---|
| 𝓛_c | cell threshold / relay threshold | L_c / L_0 | all |
| σ̃ | dim cell density in nondim area units | σ · ℓ_0² | all (sets per-cell discrete source prefactor 1/σ̃) |
| Γ̃_L | LTB4 degradation rate | Γ_L t_0 | M1–M4 |
| β̃ | per-cell R_i production rate | β t_0 | M2 |
| γ̃ | per-cell R_i degradation rate | γ t_0 | M2 |
| L̃_r | second activation threshold | L_r / L_0 | M2, M3 |
| D̃_R = D_R/D_L | ROS diffusivity ratio | ~10–16 | M3 |
| b̃ | inhibitor emission rate / LTB4 emission rate | b/a | M3, M4 |
| Γ̃_R | extracellular R degradation rate | Γ_R t_0 | M3 |
| k̃ | bimolecular degradation rate | k L_0 t_0 | M3 |
| D̃_A = D_A/D_L | adenosine diffusivity ratio | ~2–5 | M4 |
| Ã_c | adenosine inhibition threshold | A_c / L_0 | M4 |

The nondim A equation in M4 has the same structure as M1 (parameter-free
source), carrying only b̃, D̃_A, and Ã_c as model parameters. The ratio
b̃ = b/a is geometry-independent (geometric factors cancel in both 2D–2D
and 2D–3D when all fields are normalized by L_0).

### 5.4 Geometry switching in the UI

The two geometries share most dimensional parameters:

| Parameter | 2D–2D | 2D–3D |
|---|---|---|
| a, σ, D_L, L_0, L_c, r_0, u, w, χ, μ, θ | ✓ | ✓ |
| h (layer height) | ✓ | — |

When the user switches 2D-2D ↔ 2D-3D:
- All shared dimensional sliders keep their values.
- The h slider appears or disappears.
- ℓ_0, t_0, and all dimensionless groups (Λ, χ̃, μ̃, β̃, Γ̃_R, b̃, …) recompute
  automatically from the dimensional values via the geometry-appropriate
  formulas above.

This is exactly the same dim→nondim linkage already implemented for
Setups 1–3; no new UI mechanism is needed, only new linkage functions
in `nondim.js`.

### 5.5 Comparison with Setups 2–3

| Quantity | Setups 2–3 | Setup 4 |
|---|---|---|
| Length scale | σ (prescribed wave width) | D_L/c\* (emergent wave lengthscale) |
| Time scale | 1/(r_0 L_c) | D_L/c\*² (differs per geometry) |
| Cue unit | L_c | L_0 |
| Wave amplitude M | free input | emergent output |
| Wave speed C | free input | emergent output (≈ c\* × numerical prefactor) |
| λ, ϑ | from r_0, L_c | same formulas, same values |
| χ̃, μ̃ | σ, c in denominator | c\* in denominator |
| New parameters | — | Λ, 𝓛_c |

---

## 6. Model branching tree

```
Setup 4
│
├── Geometry
│   ├── 2D–2D (L and cells in same 2D layer)
│   └── 2D–3D (cells at z=0, L diffuses in 3D half-space → thin slab grid)
│
├── Cue model (exclusive choice)
│   ├── M1: basic relay (no inhibitor)
│   ├── M2: + per-cell inhibitor R_i  [intracellular]
│   ├── M3: + extracellular ROS field R  [fast diffuser, D_R/D_L ~ 10]
│   └── M4: + extracellular adenosine A  [slow diffuser, D_A/D_L ~ 2–5]
│
├── Cell polarization model (exclusive choice, same as Setups 2–3)
│   ├── GL (default)
│   └── Hill
│
└── Optional add-ons (independent toggles)
    └── M5: GRK2 receptor desensitization (per-cell G_i, density-independent speed)
```

---

## 7. Visual layout (view from above, same as Setup 3)

- 2D disk representing the petri dish, radius R̃_dish.
- Cells as dots colored by |P| (same color scheme as Setup 3).
- L field rendered as a heatmap behind the cells (Canvas 2D or WebGL).
- Optional: separate panel showing inhibitor field R or A when active.
- Wave self-extinguishes or propagates depending on model: no prescribed
  central trap needed (unlike Setup 3), but a central "candida" seed
  region can be included as the initial source that triggers the relay.

---

## 8. PDE solver architecture

### Grid

**2D–2D:** a 2D grid of `N_x × N_y` nodes covers the petri dish.

**2D–3D:** a 3D grid of `N_x × N_y × N_z` nodes; N_z typically 8–32
depending on how far L penetrates in z (see `setup4_cue_models.md` §1b).
The z = 0 face is where cells live and sources enter (Neumann flux BC).
The z = z_max face uses zero-flux BC.

L (and inhibitor fields R or A if present) are stored on their respective
grids. Agent positions are off-grid; they sample L by bilinear (2D) or
trilinear (3D at z=0 layer) interpolation.

### Time integration for L

**2D–2D.** Explicit Euler with diffusion substeps is sufficient: the
lateral grid spacing Δx̃ sets the CFL `dt̃_sub ≤ Δx̃²/4`, and the number
of substeps per agent SDE step is `n_sub = ceil(dt̃ / (Δx̃²/4))`.

**2D–3D.** Crank-Nicolson (implicit) for the diffusion operator, explicit
for the nonlinear source terms. The z-grid is non-uniform and fine near
z = 0; the smallest z-step h_0 << Δx̃ makes explicit Euler's CFL
`dt̃ ≤ h_0²/6` prohibitively small. Crank-Nicolson is unconditionally
stable for diffusion. The sparse linear system is pre-factorized (LU)
once before the time loop and applied each step.

The δ(z) source is discretized as a volumetric term at z = 0 nodes with
coefficient 2/h_0 — see `setup4_cue_models.md` §1b for the derivation.
The reference implementation is `waves_pde/LA model/2d_3d/la_2d3d_solver.py`.

### Inhibitor fields

- **Extracellular R or A** (M3, M4): second grid, same dimensionality as
  L-grid; same Euler + substep scheme.
- **Per-cell R_i or G_i** (M2, M5): one scalar per agent; Euler step at
  the same dt̃ as the agent SDE.

---

## 9. Knobs

### Geometry
```
R̃_dish       nondim dish radius (in units of ℓ_0)
N             number of cells (default 500, ≤ 2000 for Canvas2D)
N_grid        PDE grid points per side (default 128)
```

### Dimensional parameters (shared across geometries)
```
a      LTB4 emission rate per cell           [nM·µm³/s]
b      inhibitor emission rate per cell      [µm³/s per cell]   (M3, M4)
σ      cell surface density                  [cells/µm²]
D_L    LTB4 diffusivity                      [µm²/s]
L_0    relay threshold                       [nM]
L_c    cell polarization threshold           [nM]
r_0    polarization rate constant            [1/(s·nM)]
u, w   GL well coefficients                  [see setup1]
χ      chemotactic coefficient               [µm²/(s·nM) per polarization unit]
μ      motility coefficient                  [µm/s per polarization unit]
θ      noise amplitude                       [nM/s]
h      layer height (2D–2D only)             [µm]
```

### Derived dimensionless groups (recomputed on geometry switch)
```
Λ      wave time / cell time     r_0 L_c D_L / c*²
𝓛_c    cell/relay threshold      L_c / L_0
χ̃      chemotactic strength      χ L_0 / (p_0 c*)
μ̃      cell motility             μ p_0 / c*
λ      GL well depth             u² / (w r_0 L_c)
ϑ      noise                     θ w / (u r_0)
```

### Inhibitor parameters (visible when model M2/M3/M4 selected)

Dimensional input (shared by M3 and M4):
```
b      inhibitor emission rate per cell   [µm³/s per cell]
```

Nondimensional (shown/hidden based on active model):
```
β̃      per-cell R_i production rate   β t_0          (M2 only)
γ̃      per-cell R_i degradation rate  γ t_0          (M2 only)
L̃_r    second activation threshold    L_r / L_0      (M2, M3)
D̃_R    ROS diffusivity ratio          D_R/D_L        (M3 only, ~10–16)
b̃      inhibitor/LTB4 emission ratio  b/a            (M3, M4)
Γ̃_R    R-field degradation rate       Γ_R t_0        (M3 only)
k̃      bimolecular degradation        k L_0 t_0      (M3 only)
D̃_A    adenosine diffusivity ratio    D_A/D_L        (M4 only, ~2–5)
Ã_c    adenosine inhibition threshold A_c / L_0      (M4 only)
Γ̃_L    LTB4 degradation rate          Γ_L t_0        (all; 0 = no decay)
σ̃      dim cell density (nondim area) σ · ℓ_0²       (all; per-cell source rate is 1/σ̃)
```

### GRK2 (M5 toggle only)
```
k̃_G          GRK2 rate constant
γ̃_G          receptor recycling rate
g̃_L          threshold proportionality constant
```

### Threshold regularisation (Hill exponents)

One independent exponent per Hill function; all default to 10 (n → ∞ recovers
Heaviside step). Each is a separately adjustable UI knob.

| Exponent  | Hill function                                        | Active in  |
|-----------|------------------------------------------------------|------------|
| n_L       | H^+(𝓛; 1) — L activation gate in L equation         | M1,M2,M3,M4|
| n_R       | H^-(R̃; 1) — R inhibition gate in L equation         | M2         |
| n_{Lr}    | H^+(𝓛; L̃_r) — L activation gate in R/R_i equation  | M2, M3     |
| n_A       | H^-(Ã; Ã_c) — A inhibition gate in L equation       | M4         |
| n_{LA}    | H^+(𝓛; 1) — L activation gate in A equation         | M4         |

### Numerics
```
dt̃           time step (agent SDE; PDE uses sub-steps automatically)
seed          RNG seed
play speed   real-time simulation multiplier
```

---

## 10. What to plot

- **Live 2D view.** L-field heatmap (nondim), cells colored by emission
  state and |P| (see below), optional polarization arrows.

- **Cell emission state.** Each agent carries a boolean `emitting_i`
  updated each step. The rule is: a cell is emitting when every Hill
  factor in its source term individually exceeds 0.5 — equivalently,
  when the argument of each `H^+` exceeds its threshold and the argument
  of each `H^-` falls below its threshold:

  | Model | `emitting_i = true` when |
  |---|---|
  | M1 | 𝓛(r̃_i) > 1 |
  | M2 | (𝓛(r̃_i) > 1  OR  r̃_i < r̃_fire)  AND  R̃_i < 1 |
  | M3 | 𝓛(r̃_i) > 1 |
  | M4 | 𝓛(r̃_i) > 1  AND  Ã(r̃_i) < Ã_c |
  | M5 | 𝓛(r̃_i) > 𝓛̃_c(i) = g̃_L G̃_i  (per-cell threshold) |

  (M5 can be combined with any of M1–M4; the `𝓛 > 1` condition is
  replaced by `𝓛 > 𝓛̃_c(i)` throughout.) Emitting cells are rendered in
  a distinct warm/bright color; non-emitting cells in a muted/cool color.
  |P| modulates brightness within each group. **A cell deactivated by R̃
  (M2) is rendered with the same muted color as a never-activated cell —
  there is no third color to distinguish "shut off" from "never started".**

- **Optional second panel.** Inhibitor field R̃ or Ã heatmap (when active).

- **Radial L-profile 𝓛(r̃, t̃).** A 1D line plot of the azimuthally
  averaged cue field versus r̃, at the time t̃ shown in the petri dish
  view. The time slider is shared: scrubbing it updates both the dish
  animation and this profile simultaneously. The wave front, its
  amplitude, and the post-wave trough are all directly visible here.

- **Mean cell radius ⟨r̃⟩_free(t̃).** Inward drift diagnostic.

- **Wave speed C_eff(t̃).** Measured from L-field peak propagation (output,
  not input; useful for comparing 2D–2D vs 2D–3D and for M5 quorum sensing).

---

## 11. Implementation decisions and notes

1. **Boundary conditions for L.** Decided: **absorbing** — Dirichlet
   `𝓛 = 0` at the dish edge r̃ = R̃_dish. LTB4 is consumed at the agar
   boundary; zero-flux would let the concentration pile up at the edge
   and re-enter the domain.

2. **Initial condition — model-dependent.** Cells always start unpolarised
   (P_i = 0) and uniformly distributed in the dish; 𝓛 starts at zero
   everywhere. The wave is launched differently depending on the model:

   **M1 (relay only) — time-limited firing source.** A bulk source term
   `s_fire` is added to the L-field at every grid node inside r̃ < r̃_fire
   for t̃ < t̃_fire, then turned off. After cutoff the relay term sustains
   the wave (or it dies, depending on Γ̃_L). Three knobs: r̃_fire, t̃_fire,
   s_fire. This is the parent_solver convention.

   **M2 (per-cell inhibitor) — forced-emission disk, no time cutoff.**
   Cells whose *current* position satisfies r̃ < r̃_fire are forced to
   emit at full strength regardless of the local 𝓛 — the `H^+(𝓛_i;1;n_L)`
   gate in the source weight is overridden to 1 for these cells. The
   `H^-(R̃_i;1;n_R)` inhibitor gate still applies. R̃_i continues to be
   integrated for *all* cells (in or out of the disk). The wave terminates
   organically when R̃_i builds up in the central cells and shuts off
   their emission: the effective firing duration is set by the R-ODE
   timescale (∼1/γ̃ once R̃_i saturates near β̃/γ̃, or by the time to
   cross R̃ ≈ 1 from rest, ∼ ln(β̃/(β̃−γ̃))/γ̃). Knobs in M2: only r̃_fire;
   `t̃_fire` and `s_fire` are not used. Membership in the forced-emit
   region is checked per step against the cell's *current* position, so
   a cell that migrates in or out switches its forced state accordingly.

3. **δ(z) normalisation cross-check.** The 2/h_0 source coefficient must
   be verified against the analytic M1 wave speed. In the nondim units
   defined in §5, the continuum-limit wave speed is exactly 2/π (where
   c\* = 1 by construction). However, a discrete cell distribution
   introduces finite-density corrections, so the measured simulation
   speed will approach but not equal 2/π. This is expected and is not a
   solver error; the correction vanishes as N/R̃²_dish → ∞.

4. **Cell–grid coupling (all models).** Cells are off-grid point sources.
   The same interpolation weights must be used in both directions:
   sampling L (and R or A) at the agent position, and distributing the
   agent's emission back onto grid nodes. Bilinear (2D-2D) or trilinear
   (2D-3D, z = 0 layer) interpolation is the natural choice; using
   nearest-grid-point for smearing but bilinear for sampling (or vice
   versa) creates a systematic flux imbalance. Per-cell ODE variables
   (R_i in M2, G_i in M5) are updated using the already-interpolated
   L(r_i) and require no additional grid coupling.

5. **Cell SDE adaptive sub-stepping.** Explicit Euler–Maruyama on the GL
   polarization equation requires `dt · Λ · |𝓛 − 𝓛_c| ≲ 0.3` for stability.
   In regimes where the wave saturates strongly (e.g. low Γ̃_L, large σ̃·N),
   𝓛 can reach 10³–10⁴ and a single agent-dt step would blow `P` up to NaN
   within a handful of iterations. `agents.js stepAgents` therefore
   adaptively substeps each cell's SDE with `n_sub = ⌈Λ·|𝓛−𝓛_c|·dt/0.3⌉`
   and a per-substep noise amplitude `√(2 Λ ϑ 𝓛) · √(dt_sub)` (total
   variance over the agent step is preserved). L and ∇L are frozen across
   substeps — matches the worker's PIC accumulation cadence (one L-sample
   per agent step). Cells far from saturation use `n_sub = 1` (no overhead).

6. **Heatmap normalization uses the radial-profile maximum, not the grid
   maximum.** Discrete-ABM cells produce single-grid-cell PIC hot spots
   that can stochastically be 5–20× higher than the smooth wave amplitude.
   Normalizing the per-frame heatmap by the grid max lets one stochastic
   peak dominate the colormap, dimming every other cell frame-to-frame
   (visible to the user as cells "flickering on and off one-by-one"). The
   radial profile is azimuthally averaged within each annulus, so its max
   is stable and reflects the smooth wave amplitude; cell hot spots above
   that scale clip to saturated yellow and stay uniformly bright. True
   grid max is sent as a scalar (`Lmax`) for the KPI.

7. **Interaction paradigm: batch compute, not real-time.** The PDE
   time-loop (128² grid × many substeps) is too expensive for real-time
   parameter updates. Setup 4 therefore uses a **calculate-then-explore**
   workflow:
   - User sets all parameters via sliders.
   - User presses **"Calculate"**: the simulation runs to t̃_max in a Web
     Worker (non-blocking), saving frames at fixed intervals.
   - When done, the time slider and petri-dish view become active;
     the user scrubs t̃ through the saved trajectory.
   - Changing any parameter invalidates the cached run and re-enables
     the "Calculate" button.
   This is a deliberate departure from the real-time Manipulate style of
   Setups 1–3, which is feasible there because those solvers are cheap.
