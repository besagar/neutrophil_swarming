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

### 2c. Target (sticking boundary) — all full-swarm models

Optional inner boundary, on by default, available for **every** full-swarm cue
model (M1, M2, M6.1, M6.2) because it is a property of the cells, not of the
cue field. A circle of radius `R̃_target` (nondim, in units of ℓ_0; default 2)
at the dish centre stands for the object the swarm converges on — a candida
cluster, a bead, a sterile injury. Its dynamics are deliberately minimal:

```
stick_target = on:
    cell reaches r̃_i ≤ R̃_target  →  r̃_i is projected onto the circle
                                     (r̃_i := R̃_target r̂_i) and FROZEN
                                     for the rest of the run.
```

- The crossing is **absorbing, not reflecting**: the first crossing pins the
  cell permanently (adhesion/engagement), a one-way flag `stuck_i`.
- A stuck cell keeps its full **polarization SDE** and keeps **emitting** the
  cue exactly as any other cell (P_i still evolves, only ṙ̃_i = μP_i is
  suppressed). Killing the emission would silently change the cue model.
- Initial condition: the target disk is **excluded** from the uniform cell
  placement, so no cell starts already engaged. N is unchanged, so σ̃ still
  means N/(πR̃²_dish); the accessible area shrinks by (R̃_target/R̃_dish)²
  (0.4% at the defaults R̃_target = 2, R̃_dish = 30).
- `stick_target = off` restores the plain Setup-4 dish (reflective outer
  boundary only, uniform placement over the whole disk).
- Outer boundary (reflective at r̃ = R̃_dish) is unaffected either way.

Diagnostics: the dish view draws the circle in teal and rings each engaged cell;
a KPI reports the engaged count, i.e. the recruitment curve of the target.

---

## 3. Cell dynamics

Identical to Setup 3: a UI toggle selects between the GL and Hill
polarization models. The cue field 𝓛 and its gradient are computed by
the PDE solver rather than read from a formula; everything else is
unchanged.

**GL model:**
```
Ṗ_i = κ ∇̃𝓛(r̃_i, t̃) + λ(𝓛(r̃_i, t̃) − 𝓛_c) P_i
      + ν(|P_i|² − |P_i|⁴) P_i + √(2ϑ𝓛(r̃_i, t̃)) ξ

ṙ̃_i = μ P_i
```

Five independent nondim groups λ, ν, κ, μ, ϑ — see §5.2 for formulas.
No Λ wrapper; each coefficient is independently tunable.

**Hill model:** see [hill_model.md](hill_model.md). Shared knobs: κ, ϑ.
GL-only: λ, ν, μ. Hill-only: ṽ₀, n.

---

## 4. Cue field models

Four cue field models are available (see `setup4_cue_models.md` for
full equations):

| Toggle label | Document label | Inhibitor | Key parameter(s) added |
|---|---|---|---|
| "Relay only" | M1 | none | [Γ_L]  *(M1 source is parameter-free per cell; wave speed emerges as (2/π)σ̃)* |
| "Per-cell inhibitor" | M2 | R_i (intracellular ODE) | β, Γ_R, 𝓛_r |
| "ROS (fast diffuser)" | M3 | R (extracellular PDE) | D̃_R, b̃, k̃, 𝓛_r |
| "Adenosine (slow diffuser)" | M4 | A (extracellular PDE) | D̃_A, b̃, Ã_c |

An optional toggle "GRK2 desensitization" adds the per-cell threshold
modulation (M5) on top of any of the above.

The inhibitor type determines whether the simulation carries an extra
PDE grid (extracellular) or extra per-cell scalars (intracellular).

---

## 5. Nondimensionalization

### 5.1 Length and time scales — intrinsic-units scheme

**Cue unit: L_0** (relay threshold). This makes the relay threshold
`H^+(L; L_0; n_L)` become `H^+(𝓛; 1; n_L)` in nondim, so the M1 L
equation carries only the Hill exponent n_L (see below). Note `L_0 ≠ L_c` in general; the
ratio `𝓛_c = L_c / L_0` appears as a new dimensionless parameter in the
cell equation.

**Length and time scales (single-cell intrinsic — σ does NOT appear):**
```
ℓ_0 = a / (L_0 D_L)                  [µm]
t_0 = a² / (L_0² D_L³)  (= ℓ_0²/D_L) [s]
```

Physical interpretation: the steady-state 3D diffusion profile around a
single emitting cell is L = a/(4π D_L r); setting L = L_0 gives
r = ℓ_0/(4π). So ℓ_0 is (up to a 4π factor) the **single-cell influence
radius** — the distance at which one cell's signal reaches the relay
threshold. σ, R, N, h, and all cell-side parameters live *inside* the
units; dragging N or R_dim does not rescale ℓ_0 or t_0.

The nondim diffusivity is D̃ = D_L t_0/ℓ_0² = 1 identically.

**M1 L equation in nondim units — discrete-ABM (per-cell δ-source) form:**
```
2D–3D:  ∂_t̃ 𝓛 = ∇̃²𝓛 + δ̃(z̃) Σ_i H⁺(𝓛_i;1;n_L) δ̃²(r̃ − r̃_i) − Γ_L 𝓛
2D–2D:  ∂_t̃ 𝓛 = ∇̃²𝓛 + (1/h̃) Σ_i H⁺(𝓛_i;1;n_L) δ̃²(r̃ − r̃_i) − Γ_L 𝓛
```
Per-cell prefactor: **1** in 2D-3D, **1/h̃** in 2D-2D. No 1/σ̃ anywhere.
σ̃ emerges only in the continuum limit Σ_i δ̃²(r̃ − r̃_i) → σ̃ — it is a
true cell density (cells per ℓ_0²), not a derived unit-system quantity.

**Continuum forms (for reference):**
```
2D–3D (continuum):  ∂_t̃ 𝓛 = ∇̃²𝓛 + δ̃(z̃) σ̃ H⁺(𝓛;1;n_L) − Γ_L 𝓛
2D–2D (continuum):  ∂_t̃ 𝓛 = ∇̃²𝓛 + (σ̃/h̃) H⁺(𝓛;1;n_L) − Γ_L 𝓛
```

**M2 L equation in nondim units (per-cell inhibitor R_i).** With ℛ_i = R_i / R_0
(R_0 ≡ 1 implicit) and the discrete-ABM δ-source form:
```
2D–3D: ∂_t̃ 𝓛 = ∇̃²𝓛 + δ̃(z̃) Σ_i H⁺(𝓛_i;1;n_L) H⁻(ℛ_i;1;n_R) δ̃²(r̃ − r̃_i) − Γ_L 𝓛
2D–2D: ∂_t̃ 𝓛 = ∇̃²𝓛 + (1/h̃) Σ_i H⁺(𝓛_i;1;n_L) H⁻(ℛ_i;1;n_R) δ̃²(r̃ − r̃_i) − Γ_L 𝓛
```
Per-cell inhibitor ODE (same in both geometries):
```
∂_t̃ ℛ_i = β H⁺(𝓛(r̃_i); 𝓛_r; n_{Lr}) − Γ_R ℛ_i
```
The PDE solver structure is identical to M1 — only the per-cell source
weight changes from `H⁺(𝓛_i;1;n_L)` to the gated product
`H⁺(𝓛_i;1;n_L) · H⁻(ℛ_i;1;n_R)`. ℛ_i lives on the agents (no extra
grid), updated by explicit Euler at the agent dt̃.

### 5.2 Cell dynamics in L-wave nondim

Polarization scale: `p_0 = √(u/w)` (unchanged from Setups 1–3).

With `t̃ = t/t_0`, `x̃ = x/ℓ_0`, `𝓛 = L/L_0`, `P_i = p_i/p_0`, the GL SDE becomes:

```
Ṗ_i = κ ∇̃𝓛(r̃_i) + λ(𝓛 − 𝓛_c)P_i + ν(|P_i|² − |P_i|⁴)P_i + √(2ϑ𝓛) ξ
ṙ̃_i = μ P_i
```

Five independent geometry-independent dimensionless groups:

```
λ = r_0 · t_0 · L_0                    (linear activation; uses L_0, not L_c)
ν = u² · t_0 / w                        (GL nonlinearity / well depth)
κ = a χ / (D_L² · p_0)                  (chemotactic coupling)
μ = μ_dim · t_0 / ℓ_0                   (cell motility)
ϑ = (w L_0 t_0 / u) · θ                 (noise amplitude)
𝓛_c = L_c / L_0                         (cell threshold / relay threshold)
```

There is **no Λ wrapper** and **no geometry factor** in any of these.
λ, ν, κ, μ, ϑ are entirely determined by single-cell and molecular
parameters; σ and the geometry are absent.

### 5.3 Inhibitor nondim groups

**Convention:** all concentration fields (L, R, A) are normalized by L_0.
Rates are normalized by t_0. Diffusivity ratios by D_L.

| Symbol | Meaning | Formula | Model |
|---|---|---|---|
| 𝓛_c | cell threshold / relay threshold | L_c / L_0 | all |
| σ̃ | true cell density (cells per ℓ_0²) | σ · ℓ_0² | all; emerges in continuum limit of Σ_i δ̃² |
| h̃ | nondim layer height | h / ℓ_0 | 2D-2D only; sets 1/h̃ source prefactor |
| Γ_L | LTB4 degradation rate | γ_L · t_0 | M1–M4 |
| β | per-cell R_i production rate | b · t_0 / R_0 (R_0 ≡ 1) | M2 |
| Γ_R | per-cell R_i degradation rate | γ_R · t_0 | M2 |
| 𝓛_r | second activation threshold | L_r / L_0 | M2, M3 |
| D̃_R = D_R/D_L | ROS diffusivity ratio | ~10–16 | M3 |
| b̃ | inhibitor emission rate / LTB4 emission rate | b/a | M3, M4 |
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
- The h slider appears (2D-2D) or disappears (2D-3D).
- **ℓ_0 and t_0 are unchanged** — they depend only on a, L_0, D_L.
- **λ, ν, κ, μ, ϑ are unchanged** — all five cell-side groups are
  geometry-independent.
- Only σ̃ (source density in continuum limit) and h̃ (2D-2D source
  prefactor 1/h̃) carry the geometry/density information.

This contrasts with the old σ-baked scheme where ℓ_0 depended on σ and
geometry, so *every* group rescaled on a geometry switch.

### 5.5 Comparison with Setups 2–3

| Quantity | Setups 2–3 | Setup 4 |
|---|---|---|
| Length scale | σ (prescribed wave width) | single-cell influence ℓ_0 = a/(L_0 D_L) |
| Time scale | 1/(r_0 L_c) | ℓ_0²/D_L (geometry-independent) |
| Cue unit | L_c | L_0 |
| Wave amplitude M | free input | emergent output |
| Wave speed C | free input | emergent output: c̃_M1 = (2/π)σ̃ (2D-3D) |
| λ, ϑ | from r_0, L_c | replaced by five groups λ, ν, κ, μ, ϑ |
| χ̃, μ̃ | σ, c in denominator | removed; replaced by κ, μ (no Λ wrapper) |
| New parameters | — | λ, ν, κ, μ replace Λ; σ̃ and h̃ carry geometry/density |

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

### Target (all full-swarm models: M1, M2, M6.1, M6.2)
```
stick to the target   on/off checkbox (default ON)
R̃_target             nondim target radius (in units of ℓ_0, default 2);
                      cells that reach it adhere and stop moving (§2c)
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

### Derived dimensionless groups (geometry-independent; recomputed from dim sliders)
```
λ      linear activation          r_0 · t_0 · L_0          (uses L_0, not L_c)
ν      GL nonlinearity            u² · t_0 / w
κ      chemotactic coupling       a χ / (D_L² · p_0)
μ      cell motility              μ_dim · t_0 / ℓ_0
ϑ      noise amplitude            (w L_0 t_0 / u) · θ
𝓛_c    cell/relay threshold       L_c / L_0
σ̃      true cell density          σ · ℓ_0²   (regime indicator: ≪1 discrete, ≳1 continuum)
h̃      nondim layer height        h / ℓ_0    (2D-2D only; sets 1/h̃ source prefactor)
```

Diagnostics (KPI, not sliders):
```
ℓ_0    single-cell influence radius   a / (L_0 D_L)       [µm]
t_0    intrinsic time unit            ℓ_0² / D_L           [s]
c̃_Dieterle  analytic wave speed      (2/π) σ̃   (2D-3D M1 continuum limit)
```

### Inhibitor parameters (visible when model M2/M3/M4 selected)

Dimensional input (shared by M3 and M4):
```
b      inhibitor emission rate per cell   [µm³/s per cell]
```

Nondimensional (shown/hidden based on active model):
```
β      per-cell R_i production rate   b · t_0 / R_0 (R_0 ≡ 1)   (M2 only)
Γ_R    per-cell R_i degradation rate  γ_R · t_0                  (M2 only)
𝓛_r    second activation threshold    L_r / L_0                  (M2, M3)
D̃_R    ROS diffusivity ratio          D_R/D_L                    (M3 only, ~10–16)
b̃      inhibitor/LTB4 emission ratio  b/a                        (M3, M4)
k̃      bimolecular degradation        k L_0 t_0                  (M3 only)
D̃_A    adenosine diffusivity ratio    D_A/D_L                    (M4 only, ~2–5)
Ã_c    adenosine inhibition threshold A_c / L_0                  (M4 only)
Γ_L    LTB4 degradation rate          γ_L · t_0                  (all; 0 = no decay)
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

- **Channelisation (streaming) order parameter Ψ(t̃).** Quantifies the
  Höfer–Maini streaming instability — the spontaneous breaking of rotational
  symmetry into radial "channels"/spokes of cells. Streaming shows up as
  structure in the *angular* distribution of cell positions, so measure it
  with the azimuthal Fourier modes of the cell angles θ_i:

  ```
  c_m = (1/N) Σ_i e^{i m θ_i},        m = 1, 2, …, m_max        (m_max = 16)
  ```

  taken over the cells outside a small excluded core (r̃_i > 0.1 R̃_dish, where
  θ is meaningless and cells crowd), N = number of cells in that region.

  For N *independent, uniformly distributed* angles E[|c_m|²] = 1/N, i.e. a
  perfectly unstructured swarm still gives non-zero power purely from shot
  noise, and that floor grows as density falls. Subtracting it gives a
  density-independent order parameter on [0, 1]:

  ```
  Ψ_m = √( max(0, (N|c_m|² − 1)/(N − 1)) )
  ```

  - `Ψ_m = 0` — angularly uniform (within shot noise).
  - `Ψ_m = 1` — all cells at one angle, m-fold periodic.
  - The 5% shot-noise level is `Ψ_noise = √(2/(N−1))` (since 2N|c_m|² is
    χ²₂-distributed under the null); values below it are not significant.

  Plot two things:
  1. **Ψ(t̃) = max over m ≥ 2 of Ψ_m** — a single "how channelised is it"
     trace. m = 1 is excluded on purpose: it measures a bulk off-centre
     drift of the whole swarm, not spokes.
  2. **The spectrum Ψ_m vs m at the scrubbed time**, with the dominant mode
     m* labelled — m* *is* the channel count, so a swarm that has broken into
     six spokes peaks at m = 6.

  This is a *cell-position* diagnostic and is model-independent, so the same
  panel appears on the M1, M6.1 and M6.2 pages: M1 is the channelising
  baseline the M6 variants are meant to suppress, and comparing Ψ(t̃) across
  those pages at equal σ̃ is the actual test of that claim.

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
   be verified against the analytic M1 wave speed. In the intrinsic-units
   nondim defined in §5, the continuum-limit 2D-3D wave speed is
   c̃ = (2/π)σ̃ (where σ̃ = σ·ℓ_0² is the true cell density). The 2/h_0
   control-volume factor and the convergence of the discrete sum to the
   continuum source as N/R̃²_dish → ∞ both remain unchanged. However, a
   discrete cell distribution introduces finite-density corrections, so
   the measured simulation speed will approach but not equal (2/π)σ̃. This
   is expected and is not a solver error; the correction vanishes as
   N/R̃²_dish → ∞. The Dieterle overlay degrades at both extremes of σ̃:
   σ̃ ≪ 1 (discrete-cell regime, where the analytic front is a poor
   approximation) and σ̃·Δx̃ ≳ 0.5 (sub-grid front, resolution limit).

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
   polarization equation requires `dt · max(λ·|𝓛 − 𝓛_c|, ν) ≲ 0.3` for
   stability. In regimes where the wave saturates strongly (e.g. low Γ_L,
   large σ̃·N), 𝓛 can reach 10³–10⁴ and a single agent-dt step would blow
   `P` up to NaN within a handful of iterations. `agents.js stepAgents`
   therefore adaptively substeps each cell's SDE with
   `n_sub = ⌈max(λ·|𝓛−𝓛_c|, ν)·dt/0.3⌉`
   and a per-substep noise amplitude `√(2ϑ𝓛) · √(dt_sub)` (total
   variance over the agent step is preserved; noise has no Λ factor). L
   and ∇L are frozen across substeps — matches the worker's PIC
   accumulation cadence (one L-sample per agent step). Cells far from
   saturation use `n_sub = 1` (no overhead).

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
