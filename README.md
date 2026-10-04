# Neutrophil Swarming — GL Motility Interactive Tool

An interactive browser-based simulation of a **Ginzburg–Landau (GL) model** for neutrophil polarization and chemotaxis. Manipulate-style sliders drive live plots and animations entirely in nondimensional units — no build step, no backend.

Live site: <https://besagar.github.io/neutrophil_swarming/> (GitHub Pages, served from `main` root).

## Physics background

Classical chemotaxis models (Keller–Segel) predict zero net displacement when a traveling chemical wave passes a cell: the upward leg exactly cancels the downward leg. This tool explores a minimal fix: treating the cell's **polarization vector p** as an internal state with its own relaxation dynamics, governed by a sixth-order GL free energy:

$$\mathcal{F} = -\frac{r_0(L-L_c)}{2}p^2 - \frac{u}{4}p^4 + \frac{w}{6}p^6 - \chi\,\mathbf{p}\cdot\nabla L$$

The polarization lag gives the cell an effective "memory" that breaks the symmetry of the wave, producing net chemotactic drift and, at the collective level, a Keller–Segel-type clumping instability.

## Setups

| Page | Description |
|------|-------------|
| **1 — Uniform cue** (`setup1/`) | Single cell, spatially uniform $L = L_0$. Polarization performs Brownian dynamics in the GL free-energy landscape; bifurcation from unpolarized to polarized. |
| **2 — Gaussian wave** (`setup2/`) | Single cell in 1D vs. a prescribed traveling Gaussian cue pulse. Polarization lag produces net displacement. GL or Hill polarization model. |
| **3 — Radial swarm** (`setup3/`) | $N \sim 10^3$ cells on a 2D disk; prescribed radial Gaussian waves from the center; cells move inward (against the wave) into a central trap. |
| **4a — M1 relay** (`setup4/`) | Emergent wave: cells relay-emit LTB4 above threshold; the cue field is solved as a PDE. Front speed vs. Dieterle's $(2/\pi)\tilde\sigma$. |
| **4b — M2** (`setup4/m2/`) | Self-extinguishing swarm: per-cell inhibitor $R_i$ shuts emission off in the wake. |
| **4c — M6.1** (`setup4/m6.1/`) | Basal adenosine field raises the relay threshold → density-independent wave speed. |
| **4d — M6.2** (`setup4/m6.2/`) | Basal quorum signal throttles the LTB4 production rate. |

## Running locally

```bash
python3 serve.py          # http://localhost:8000, caching disabled (use this)
# or: python3 -m http.server 8000   (but browsers will cache Setup 4's worker modules)
```

No npm, no bundler. CDN ES modules only (KaTeX, reveal.js for the deck).

## What is where

```
.
├── README.md                     # this file — the map
├── CLAUDE.md / AGENTS.md         # rules for AI agents working in this repo
├── ginzburg_landay_neutrophils.md  # source physics spec (top of the doc hierarchy)
├── index.html                    # landing page / navigation
├── serve.py                      # no-cache dev server
│
├── shared/                       # code shared by all pages
│   ├── rng.js                    #   seeded PRNG + Gaussian
│   ├── dom.js                    #   sliders/knobs, KaTeX axis-label overlays
│   ├── canvas.js                 #   Canvas2D plotting primitives (nondim axes)
│   ├── svgctx.js, svgexport.js   #   vector (SVG) export of canvas plots
│   └── style.css
├── setup1/                       # uniform cue — sim_core.js (physics+drawing, reused by slides/) + main.js (page)
├── setup2/                       # 1D running wave — main.js
├── setup3/                       # 2D radial swarm, prescribed wave — main.js
├── setup4/                       # emergent-wave swarm (dynamic cue PDE)
│   ├── index.html                #   4a page (M1)
│   ├── m2/ m6.1/ m6.2/           #   4b–4d pages: thin HTML shells over ui.js
│   ├── ui.js nondim.js render.js #   UI, dim↔nondim linkage, plotting
│   ├── worker.js sim_core.js agents.js  # Web Worker → pure physics core → agents
│   └── solvers/                  #   cue-field PDE solvers (factory, grid, M1, aux-field)
│
├── docs/
│   ├── PLAN.md                   # staging, deliverables, validation checklist (§6), decisions
│   ├── physics/                  # per-setup equations + nondimensionalization (source of truth for knobs)
│   │   ├── setup1_uniform.md  setup2_wave.md  setup3_swarm.md
│   │   ├── setup4_swarm3d.md     #   Setup 4 geometry, agents, nondim
│   │   ├── setup4_cue_models.md  #   catalog of cue models M1–M6.2
│   │   ├── hill_model.md         #   alternative polarization model (Setups 2–3)
│   │   └── figs/                 #   Setup 1 theory-check figures + the .py scripts that make them
│   ├── plans/                    # implementation plans / proposals (historical record, all implemented)
│   └── design/                   # ui_conventions.md, slides.md
│
├── slides/                       # reveal.js talk deck reusing the live applets
│   ├── deck.html                 #   the deck
│   ├── setup1-slide.{html,js}, slide.css
│   ├── mathcha/                  #   rendered Mathcha math pages (page-NN.png)
│   ├── build-math.sh             #   regenerates mathcha/ from mathcha_src/pdf_slides.pdf
│   └── mathcha_src/              #   (gitignored) local Mathcha export
│
├── paper/                        # (gitignored) paper manuscript, LaTeX
└── materials/                    # (gitignored) non-code research material, local only
    ├── literature/               #   cue-model source PDFs + notes behind setup4_cue_models.md
    ├── talks/                    #   group-meeting slide PDFs
    └── abstracts/                #   conference/poster abstracts (LaTeX + PDF)
```

### Documentation hierarchy

Edits cascade top-down; code never leads the docs:
`ginzburg_landay_neutrophils.md` → `docs/physics/*.md` → `docs/PLAN.md` → source.

## Design principles

- **All simulations run in nondimensional units.** Dimensional sliders convert to nondim via per-setup linkage functions; plots never show unit strings.
- **Numerics are explicit and named.** Stochastic terms use Euler–Maruyama; deterministic blocks use RK4.
- **Pure simulation core.** No DOM access inside physics modules; UI layer wires sliders to state.
