# Presentation slides (`slides/`)

HTML slides built on top of the live setup applets, for talks. Lives in the same
repo so slides reuse the setup simulation code directly (zero-build, relative
imports) instead of duplicating physics.

## Pieces

- **`slides/deck.html`** — the reveal.js deck (loaded from CDN, 16:9, keyboard
  nav, PDF export via `?print-pdf` / the `e` key). Each slide is a `<section>`.
  Live-applet slides embed a slide page in an `<iframe>`.
- **`slides/setup1-slide.html` + `setup1-slide.js`** — the Setup-1 slide: the
  five 1D applets in a 3×2 grid with a compact control panel (𝓛, λ, ϑ +
  play/pause/reset + KPIs) in the freed bottom-right cell.
- **`slides/slide.css`** — stage + grid + compact-panel styling.

## Single source of truth

The Setup-1 physics and Canvas2D drawing live in **`setup1/sim_core.js`**
(`createSetup1()`), imported by **both** `setup1/main.js` (the website page) and
`slides/setup1-slide.js` (the slide). No equation is written twice; the panel
sliders remain the single source of truth for the simulated params (`bind:`), so
displayed value ≡ simulated value on the slide exactly as on the site. All plot
axes stay nondimensional (hard rule in `CLAUDE.md`).

## Two sizing gotchas (already handled — don't reintroduce)

1. **No `transform: scale()` on the stage.** `autoFit()` writes
   `getBoundingClientRect().width` back as the canvas layout width; under a
   scaled ancestor that rect is the *visually-scaled* width, which feeds back and
   shrinks every plot to zero. The stage is a **real responsive 16:9 box**
   (`min(100vw, 100vh·16/9)`), and canvases size in `cqh` (fraction of the grid
   height via a container query). Reveal.js scales its slides with a transform,
   but the applet is in an `<iframe>` — a separate viewport — so that transform
   never touches the canvas measurement.
2. **Re-fit on resize.** If `autoFit()`'s first measure happens while the stage
   is 0-sized (a reveal iframe sized after `load`, pre-fullscreen), the plot
   locks collapsed. `setup1-slide.js` runs a `ResizeObserver` on the grid that
   clears the inline `canvas.style.{width,height}` on any resize, so the CSS
   reapplies and the next frame re-measures against the real box.

## Mathcha slides (`slides html/`, gitignored)

`slides html/` is a Mathcha document export (one long scrolling page of
handwritten math + SVG diagrams, ~8.5 MB with bundled fonts) — **not** a
paginated deck, so it doesn't drop into reveal cleanly. Keep Mathcha as the
authoring tool for static math; to bring a derivation into the deck, export that
diagram from Mathcha as PNG/SVG and place it on a slide as an `<img>` (see the
placeholder slide in `deck.html`). Don't rebuild the math in HTML.

## Adding a setup slide

Extract that setup's physics+drawing into a `setupN/sim_core.js` factory (as
Setup 1 does), then add `slides/setupN-slide.{html,js}` reusing it, and a
`<section>` iframe in `deck.html`.
