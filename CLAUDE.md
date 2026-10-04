# CLAUDE.md — GL Motility interactive tool

This repo builds an interactive web demonstration of a Ginzburg–Landau (GL) model
of neutrophil polarization and chemotaxis. The physical model lives in
[ginzburg_landay_neutrophils.md](ginzburg_landay_neutrophils.md). The detailed
implementation plan lives in [docs/PLAN.md](docs/PLAN.md). Read both before
making non-trivial changes.

## What this project is (and is not)

- A self-contained interactive HTML/JS tool, runnable by opening `index.html`
  in a browser. Three subpages, one per physical setup. Manipulate-style widgets
  (sliders + live plots/animations).
- **Not** a production web app. No backend, no auth, no build pipeline unless
  strictly necessary. Prefer zero-build (ES modules + CDN libs) so the user can
  serve with `python -m http.server` and iterate fast.

## Hard rules

1. **All simulation math AND all plots run in nondimensional units.** Plot
   axes, ticks, labels, animation coordinates, tooltip values are nondim,
   always — never write a unit string like `s`, `µm`, `nM` on a plot.
   Dimensional values exist only as *slider inputs*: moving a dim slider
   recomputes the linked nondim value via the unit-system linkage, and the
   simulation/plots see only the nondim value. The nondim ↔ dim mapping is
   per-setup and documented in `docs/physics/setup{1,2,3}.md`.
2. **Per-setup nondimensionalization may differ.** Do not try to unify them
   across the three subpages — pick the natural scales for each problem.
3. **Setup 3 (swarming) is under-specified.** The cue field `L(x,y,z,t)`
   equation has not been provided yet. Build the agent/visualization scaffold
   but stub `L` until the user supplies dynamics.
4. **No silent parameter changes.** Every knob the user sees corresponds to a
   single named symbol in the physics docs; if you add a knob, add it to the
   physics doc first.
5. **Numerical integrators are explicit and named.** Stochastic terms use
   Euler–Maruyama; deterministic-only blocks may use RK4. Document the scheme
   and the dt → nondim-time mapping next to the integrator.

> **TODO (open, user to decide):** two rules here look stale — Hard rule 3
> (Setup 3 cue "under-specified"; Setup 3 now has a prescribed wave, the
> dynamic cue lives in Setup 4) and Workflow rule 3 (`setup{n}/ui.js` KNOBS
> exists only in Setup 4). Revisit before relying on them.

## Repo layout

Full annotated map: [README.md](README.md) § "What is where". In short:

```
GL motility/
├── ginzburg_landay_neutrophils.md   # source physics spec (do not edit casually)
├── CLAUDE.md / AGENTS.md / README.md
├── index.html  serve.py             # landing page; no-cache dev server
├── shared/                          # RNG, widgets (dom.js), canvas plotting, SVG export
├── setup1/ setup2/ setup3/          # one folder per setup page
├── setup4/                          # emergent-wave swarm: 4a page + engine;
│   └── m2/ m6.1/ m6.2/              #   4b–4d variant pages (thin shells over ../ui.js)
├── docs/
│   ├── PLAN.md                      # implementation plan (staging source of truth)
│   ├── physics/                     # per-setup nondim + equations (+ figs/)
│   ├── plans/                       # setup-4 implementation plans/proposals (history)
│   └── design/                      # ui_conventions.md, slides.md
├── slides/                          # reveal.js deck reusing live applets
├── paper/                           # (gitignored) paper LaTeX
├── materials/                       # (gitignored) literature, talk PDFs, abstracts
└── .superpowers/                    # local skill notes
```

## Workflow

- New physics or new knob → edit `docs/physics/setup{n}.md` first, then code.
- New visual idea → sketch in `docs/design/ui_conventions.md`, then code.
- The plan in `docs/PLAN.md` is the source of truth for staging. Update it when
  scope changes; do not let code drift ahead of the plan silently.

## Workflow rules (operational)

Follow these in order whenever scope or behavior changes:

1. **Source-of-truth ordering.** Edits cascade
   `ginzburg_landay_neutrophils.md` → `docs/physics/setup{n}.md` →
   `docs/PLAN.md` → source files. Never let source files lead the docs.
2. **Per-setup work uses the validation checklist** in `docs/PLAN.md` §6
   before declaring a subpage done.
3. **Knobs are declared once.** Each user-facing parameter lives in a
   `KNOBS` config in `setup{n}/ui.js` with an `exposure: 'dim' | 'nondim' |
   'both'` field and a linkage function if linked. To re-expose a knob as
   dim instead of nondim, change one field — never duplicate.
4. **Local serving.** `python3 serve.py` from the repo root (no-cache
   server; plain `http.server` caches worker modules). No npm, no
   bundler, no build step. CDN ES modules only.
5. **`.superpowers/`.** Drop short notes here for skills (`/loop`,
   `/review`, etc.), playbooks ("how to add a setup"), and ADR-style
   decisions. Not shipped, not loaded at runtime.
6. **Memory.** User-specific facts and feedback live in the global Claude
   memory store, not in this repo. Don't write user notes into repo files.

## Math typography

- Every formula in a description, every plot axis label, every plot title
  is rendered with **KaTeX**. Never write a plot label as plain unicode
  (`p̃`, `∂_x̃ L̃`, etc.) — use `\tilde p`, `\partial_{\tilde x} \tilde L`.
- KaTeX is loaded via CDN in each setup page; descriptions use `$...$` or
  `\(...\)` and are processed by `renderMathInElement` on load.
- Plot axis labels are HTML overlays positioned over the canvas (see
  `shared/dom.js` `decoratePlot`). The canvas itself draws only ticks and
  data, never axis labels or titles.

## Style

- Plain ES modules, no bundler. Plotly.js (or uPlot for perf-critical) and
  regl/PixiJS for the swarm view, loaded via CDN.
- Keep simulation core pure (no DOM access). UI layer wires sliders to state.
- Comments only where the *why* is non-obvious (e.g. nondim choice rationale,
  numerical stability bounds).
