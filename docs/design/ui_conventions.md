# UI conventions

## Layout

Each setup page has the same shell:

```
┌─────────────────────────────────────────────────────────┐
│ Top nav: [Setup 1] [Setup 2] [Setup 3]   title         │
├──────────────────────┬──────────────────────────────────┤
│                      │                                  │
│   Plots / animation  │   Controls (sliders, toggles)    │
│   (left, ~2/3 width) │   (right, ~1/3 width)            │
│                      │                                  │
│                      │   ─ Nondimensional knobs ─       │
│                      │   ─ Dimensional knobs ─          │
│                      │   ─ Numerical knobs ─            │
│                      │                                  │
│                      │   [▶ Play] [⏸] [⟳ Reset]        │
└──────────────────────┴──────────────────────────────────┘
```

## Knob conventions

Every knob is declared once in a JS config object with:

```js
{
  symbol: 'χ',                // the physics symbol
  exposure: 'both',           // 'dim' | 'nondim' | 'both'
  dim:    { min, max, step, default, units: 'µm² s⁻¹ /(nM/µm)' },
  nondim: { min, max, step, default, symbol: 'χ̃' },
  link:   (dim, units) => nondim,   // computed both ways from current units
}
```

Sliders are logarithmic by default for any knob spanning >1 decade.

## Plots

- Plotly.js for static-ish plots (free-energy landscape, bifurcation, summary).
- uPlot for streaming time traces (cheap, smooth).
- regl / PixiJS only for the swarm view (Setup 3).

## Color / theme

- Light theme by default, dark-theme toggle later.
- Use a sequential colormap (viridis) for any field rendering, never rainbow.
