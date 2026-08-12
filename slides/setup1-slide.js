// Setup 1 presentation slide — 16:9 layout + compact control panel.
// Reuses the exact same physics/drawing as the website page via sim_core.js:
// no duplicated equations, and the on-panel sliders remain the single source
// of truth for the simulated params (bind:) so displayed value ≡ simulated value.

import { makeSlider, makeButtonRow, section, decoratePlot } from '../shared/dom.js';
import { createSetup1 } from '../setup1/sim_core.js';

const sim = createSetup1({ seed: 42 });
const { params } = sim;

const STEPS_PER_FRAME_BASE = 50;
let stepAccum = 0;
let running = true;

// ─── compact control panel (𝓛, λ, ϑ, sim speed + transport) ───────────────
const controlsEl = document.getElementById('controls');

const sL     = makeSlider({ id: 'L',     symbol: '\\mathcal{L}',       bind: [params, 'L'],     min: 0,    max: 5,   step: 0.01, fmt: v => v.toFixed(2) });
const sLam   = makeSlider({ id: 'lam',   symbol: '\\lambda',           bind: [params, 'lam'],   min: 0.01, max: 10,  log: true,  fmt: v => v.toPrecision(3) });
const sTht   = makeSlider({ id: 'tht',   symbol: '\\vartheta',         bind: [params, 'tht'],   min: 1e-4, max: 1,   log: true,  fmt: v => v.toExponential(2) });
const sSpeed = makeSlider({ id: 'speed', symbol: '\\text{sim speed}',  bind: [params, 'speed'], min: 0.01, max: 100, log: true,  fmt: v => `${v.toPrecision(2)}×` });

const buttons = makeButtonRow([
  { label: '⏸', onClick() {
      running = !running;
      buttons.refs['⏸'].textContent = running ? '⏸' : '▶';
    } },
  { label: '⟳', ghost: true, onClick: () => sim.reset() },
]);

controlsEl.appendChild(section('controls', [sL.el, sLam.el, sTht.el, sSpeed.el, buttons.el]));

// ─── KaTeX axis labels (same conventions as the website page) ─────────────
function decorateAll() {
  decoratePlot('cv-F',     { titleTex: '\\text{free energy } F(P)',         xLabelTex: 'P', yLabelTex: 'F' });
  decoratePlot('cv-bif',   { titleTex: '\\text{extrema of } F',            xLabelTex: '\\mathcal{L}', yLabelTex: 'P_{\\text{ext}}' });
  decoratePlot('cv-meanp', { titleTex: '\\langle|P|\\rangle(\\mathcal{L})', xLabelTex: '\\mathcal{L}', yLabelTex: '\\langle|P|\\rangle' });
  decoratePlot('cv-trace', { titleTex: '\\text{trace } P(t)',              xLabelTex: 't', yLabelTex: 'P' });
  decoratePlot('cv-hist',  { titleTex: '\\text{histogram of } |P|',         xLabelTex: '|P|', yLabelTex: '\\mathcal{P}' });
}
if (window.katex) decorateAll();
else window.addEventListener('load', decorateAll);

// ─── animation loop ───────────────────────────────────────────────────────
const cvF = document.getElementById('cv-F');
const cvBif = document.getElementById('cv-bif');
const cvTrace = document.getElementById('cv-trace');
const cvHist = document.getElementById('cv-hist');
const cvMeanp = document.getElementById('cv-meanp');
const plots = [cvF, cvBif, cvTrace, cvHist, cvMeanp];

// autoFit() writes the measured size to canvas.style.{width,height} and caches
// it; if that first measure happened while the stage was still 0-sized (e.g. a
// reveal.js iframe sized after load, or before fullscreen), every plot stays
// locked to a collapsed size. Whenever the grid actually resizes, clear those
// inline sizes so the CSS width/height (cqh) reapplies and the next frame's
// autoFit re-measures against the real box.
const grid = document.querySelector('.stage-grid');
if (window.ResizeObserver && grid) {
  const ro = new ResizeObserver(() => {
    for (const cv of plots) { cv.style.width = ''; cv.style.height = ''; }
  });
  ro.observe(grid);
}

function frame() {
  if (running) {
    stepAccum += STEPS_PER_FRAME_BASE * params.speed;
    const N = Math.floor(stepAccum);
    stepAccum -= N;
    sim.step(N);
  }
  sim.drawF(cvF);
  sim.drawBif(cvBif);
  sim.drawMeanP(cvMeanp, { show2D: false });   // 1D slide — no 2D branch
  sim.drawTrace(cvTrace);
  sim.drawHist(cvHist);
  requestAnimationFrame(frame);
}

// Start only once layout is final. autoFit() measures the canvas on its first
// draw and caches that size; starting the loop before the grid has laid out
// would lock every plot to a collapsed width. Waiting for `load` avoids that.
function start() {
  sim.reset();
  requestAnimationFrame(frame);
}
if (document.readyState === 'complete') start();
else window.addEventListener('load', start);
