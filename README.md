# Neutrophil Swarming — GL Motility Interactive Tool

An interactive browser-based simulation of a **Ginzburg–Landau (GL) model** for neutrophil polarization and chemotaxis. Manipulate-style sliders drive live plots and animations entirely in nondimensional units — no build step, no backend.

## Physics background

Classical chemotaxis models (Keller–Segel) predict zero net displacement when a traveling chemical wave passes a cell: the upward leg exactly cancels the downward leg. This tool explores a minimal fix: treating the cell's **polarization vector p** as an internal state with its own relaxation dynamics, governed by a sixth-order GL free energy:

$$\mathcal{F} = -\frac{r_0(L-L_c)}{2}p^2 - \frac{u}{4}p^4 + \frac{w}{6}p^6 - \chi\,\mathbf{p}\cdot\nabla L$$

The polarization lag gives the cell an effective "memory" that breaks the symmetry of the wave, producing net chemotactic drift and, at the collective level, a Keller–Segel-type clumping instability.

## Three setups

| Setup | Description |
|-------|-------------|
| **1 — Uniform cue** | Single cell, spatially uniform $L = L_0$. Watch the polarization perform Brownian dynamics in the GL free-energy landscape; explore the bifurcation from unpolarized to polarized states. |
| **2 — Gaussian wave** | Single cell in 1D tracking a traveling Gaussian cue pulse. Live phase-space and trajectory plots show how polarization lag produces net displacement. |
| **3 — Radial swarm** | $N \sim 10^3$ cells on a 2D disk. Repeated radial Gaussian waves launch from the center; cells move inward (against the wave) and are trapped at the origin, reproducing the anti-wave swarming phenotype. |

## Running locally

```bash
python -m http.server 8000
# open http://localhost:8000
```

No npm, no bundler. CDN ES modules only (Plotly.js, uPlot, KaTeX).

## Repo layout

```
├── index.html          # landing page and navigation
├── setup1/             # uniform-cue single-cell simulation
├── setup2/             # 1D Gaussian wave single-cell simulation
├── setup3/             # 2D radial swarm (N-agent)
├── shared/             # RNG, Euler–Maruyama integrator, DOM helpers
docs/
├── PLAN.md             # implementation plan and validation checklists
├── physics/            # per-setup nondimensionalization and equations
└── design/             # UI conventions
ginzburg_landay_neutrophils.md   # full physics specification
```

## Design principles

- **All simulations run in nondimensional units.** Dimensional sliders convert to nondim via per-setup linkage functions; plots never show unit strings.
- **Numerics are explicit and named.** Stochastic terms use Euler–Maruyama; deterministic blocks use RK4.
- **Pure simulation core.** No DOM access inside physics modules; UI layer wires sliders to state.