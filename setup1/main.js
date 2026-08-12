// Setup 1 — Polarization in a uniform cue (website page).
// All physics + drawing lives in ./sim_core.js (shared with slides/). This
// file owns only the website's control panel, layout wiring, and RAF loop.

import { el, makeSlider, makeToggle, makeButtonRow, makeKpis, section, detailsSection, decoratePlot } from '../shared/dom.js';
import { attachSvgExports } from '../shared/svgexport.js';
import { createSetup1 } from './sim_core.js';

const sim = createSetup1({ seed: 42 });
const { params, dim, linkage } = sim;
const { LfromDim, lambdaFromDim, thetaFromDim } = linkage;

const STEPS_PER_FRAME_BASE = 50;
let stepAccum = 0;
let running = true;

// ─── controls ────────────────────────────────────────────────────────────
const controlsEl = document.getElementById('controls');
const kpis = makeKpis([
  { id: 'L',    label: '𝓛' },
  { id: 'lam',  label: 'λ' },
  { id: 'tht',  label: 'ϑ' },
  { id: 'p',    label: '|P|' },
]);

const sL    = makeSlider({ id: 'L',   symbol: '\\mathcal{L}', bind: [params, 'L'],   min: 0,    max: 5, step: 0.01, fmt: v => v.toFixed(2) });
const sLam  = makeSlider({ id: 'lam', symbol: '\\lambda',     bind: [params, 'lam'], min: 0.01, max: 10, log: true, fmt: v => v.toPrecision(3) });
const sTht  = makeSlider({ id: 'tht', symbol: '\\vartheta',   bind: [params, 'tht'], min: 1e-4, max: 1,  log: true, fmt: v => v.toExponential(2) });
let applyingDim = false;
function refreshKpis() {
  kpis.set('L',   params.L.toFixed(3));
  kpis.set('lam', params.lam.toPrecision(3));
  kpis.set('tht', params.tht.toExponential(2));
}
sL.onChange(()   => refreshKpis());
sLam.onChange(() => refreshKpis());
sTht.onChange(() => {});

function pushAllNondim() {
  applyingDim = true;
  sL.set(LfromDim()); sLam.set(lambdaFromDim()); sTht.set(thetaFromDim());
  applyingDim = false;
}

// Dim sliders (collapsed under a spoiler in the panel).
const linkedL  = () => `→ 𝓛 = ${LfromDim().toFixed(3)}`;
const linkedR0 = () => `→ λ = ${lambdaFromDim().toPrecision(3)},  ϑ = ${thetaFromDim().toExponential(2)}`;
const linkedTh = () => `→ ϑ = ${thetaFromDim().toExponential(2)}`;

const dimL  = makeSlider({ id: 'L_dim',  symbol: 'L',       bind: [dim, 'L'],     min: 0,    max: 5,  step: 0.01, units: '[L_c]', fmt: v => v.toFixed(2), linkedLabel: linkedL });
const dimR0 = makeSlider({ id: 'r0',     symbol: 'r_{0}',   bind: [dim, 'r0'],    min: 0.1,  max: 5,  step: 0.01, fmt: v => v.toFixed(2), linkedLabel: linkedR0 });
const dimTh = makeSlider({ id: 'theta',  symbol: '\\theta', bind: [dim, 'theta'], min: 1e-4, max: 1,  log: true,  fmt: v => v.toExponential(2), linkedLabel: linkedTh });

function refreshDimReadouts() {
  dimL.setLinkedText(linkedL()); dimR0.setLinkedText(linkedR0()); dimTh.setLinkedText(linkedTh());
}
dimL.onChange(()  => { pushAllNondim(); refreshDimReadouts(); refreshKpis(); });
dimR0.onChange(() => { pushAllNondim(); refreshDimReadouts(); refreshKpis(); });
dimTh.onChange(() => { pushAllNondim(); refreshDimReadouts(); refreshKpis(); });

const sDt    = makeSlider({ id: 'dt',    symbol: 'd\\tilde{t}',         bind: [params, 'dt'],    min: 1e-4, max: 0.05, log: true, fmt: v => v.toExponential(2) });
const sSpeed = makeSlider({ id: 'speed', symbol: '\\text{sim speed}',   bind: [params, 'speed'], min: 0.01, max: 100,  log: true, fmt: v => `${v.toPrecision(2)}×` });
const sSeed  = makeSlider({ id: 'seed',  symbol: '\\text{seed}',        value: 42,               min: 1,    max: 9999, step: 1,   transform: Math.round, fmt: v => v.toFixed(0) });
sDt.onChange(() => {});
sSpeed.onChange(() => {});
sSeed.onChange(v => sim.setSeed(v));

const modeToggle = makeToggle({
  label: 'polarization dimension',
  options: [{ id: '1d', label: '1D' }, { id: '2d', label: '2D' }],
  value: '1d',
  onChange: v => {
    sim.setMode(v);
    document.getElementById('cell-2d').style.display = v === '2d' ? '' : 'none';
  },
});

const buttons = makeButtonRow([
  { label: '⏸  pause', onClick() {
      running = !running;
      buttons.refs['⏸  pause'].textContent = running ? '⏸  pause' : '▶  play';
    } },
  { label: '⟳  reset', ghost: true, onClick: () => sim.reset() },
]);

controlsEl.appendChild(section('mode', [modeToggle.el]));
controlsEl.appendChild(section('nondim parameters (drive simulation)', [sL.el, sLam.el, sTht.el]));
controlsEl.appendChild(detailsSection('dim sliders (push linked nondim)', [dimL.el, dimR0.el, dimTh.el]));
controlsEl.appendChild(section('numerics', [sDt.el, sSpeed.el, sSeed.el, buttons.el]));
controlsEl.appendChild(kpis.el);
controlsEl.appendChild(el('div', { class: 'note' }, [
  'Bifurcation: ',
  el('span', { style: { color: '#2b6cb0', fontWeight: 600 } }, 'solid blue'),
  ' = stable extrema; ',
  el('span', { style: { color: '#999', fontWeight: 600 } }, 'dashed gray'),
  ' = unstable. Vertical line marks current 𝓛. ',
  'Histogram green dashed: Boltzmann P_ss ∝ exp(-F̃/(ϑ𝓛)) (× |P| in 2D), peak-normalized. ',
  '⟨|P|⟩ vs 𝓛: blue solid = 1D, orange dashed = 2D (with radial Jacobian).',
]));

sim.setSeed(sSeed.value);
refreshKpis();
refreshDimReadouts();

// ─── decorate plots with KaTeX axis labels ─────────────────────────────
function decorateAll() {
  decoratePlot('cv-F',     { titleTex: '\\text{free energy}',                xLabelTex: 'P', yLabelTex: '\\tilde F' });
  decoratePlot('cv-bif',   { titleTex: '\\text{extrema of } \\tilde F',      xLabelTex: '\\mathcal{L}', yLabelTex: 'P_{\\text{ext}}' });
  decoratePlot('cv-trace', { titleTex: '\\text{trace (1D: signed }P\\text{; 2D: }|P|\\text{)}',
                             xLabelTex: '\\tilde t', yLabelTex: 'P' });
  decoratePlot('cv-hist',  { titleTex: '\\text{histogram of } |P|',          xLabelTex: '|P|', yLabelTex: '\\tilde P' });
  decoratePlot('cv-meanp', { titleTex: '\\langle|P|\\rangle(\\mathcal{L}) \\text{ — Boltzmann (blue=1D, orange dashed=2D)}',
                             xLabelTex: '\\mathcal{L}', yLabelTex: '\\langle|P|\\rangle' });
  decoratePlot('cv-2d',    { titleTex: '\\text{polarization tip (2D)}',      xLabelTex: 'P_x', yLabelTex: 'P_y' });
  attachSvgExports({
    'cv-F':     () => sim.drawF(cvF),
    'cv-bif':   () => sim.drawBif(cvBif),
    'cv-trace': () => sim.drawTrace(cvTrace),
    'cv-hist':  () => sim.drawHist(cvHist),
    'cv-meanp': () => sim.drawMeanP(cvMeanp),
    'cv-2d':    () => sim.draw2D(cv2d),
  }, { prefix: 'setup1' });
}
if (window.katex) decorateAll();
else window.addEventListener('load', decorateAll);

// ─── animation loop ─────────────────────────────────────────────────────
const cvF = document.getElementById('cv-F');
const cvBif = document.getElementById('cv-bif');
const cvTrace = document.getElementById('cv-trace');
const cvHist = document.getElementById('cv-hist');
const cvMeanp = document.getElementById('cv-meanp');
const cv2d = document.getElementById('cv-2d');

function frame() {
  if (running) {
    stepAccum += STEPS_PER_FRAME_BASE * params.speed;
    const N = Math.floor(stepAccum);
    stepAccum -= N;
    sim.step(N);
  }
  sim.drawF(cvF);
  sim.drawTrace(cvTrace);
  sim.drawHist(cvHist);
  sim.drawBif(cvBif);
  sim.drawMeanP(cvMeanp);
  if (params.mode === '2d') sim.draw2D(cv2d);
  refreshKpis();
  kpis.set('p', sim.magnitude().toFixed(3));
  requestAnimationFrame(frame);
}
sim.reset();
refreshKpis();
requestAnimationFrame(frame);
