// Setup 1 — Polarization in a uniform cue.
// Per docs/physics/setup1_uniform.md:
//   dP_α/dt̃ = (𝓛 - 1) P_α + λ (|P|² P_α - |P|⁴ P_α) + √(2 ϑ 𝓛) η_α
// Three independent nondim knobs: 𝓛, λ, ϑ.

import { makeRng } from '../shared/rng.js';
import { el, makeSlider, makeToggle, makeButtonRow, makeKpis, section, detailsSection, decoratePlot } from '../shared/dom.js';
import { autoFit, makeAxis, drawFrame, strokePath, dot, clipPlot } from '../shared/canvas.js';

// ─── nondim ↔ dim linkage ────────────────────────────────────────────────
// Dim params: u, w, Lc, r0, theta, L. Derived nondim:
//   𝓛 = L / Lc
//   λ = u² / (w r0 Lc)
//   ϑ = θ w  / (u r0)
const dim = { u: 1.0, w: 1.0, Lc: 1.0, r0: 1.0, theta: 0.05, L: 1.5 };
function LfromDim()      { return dim.L / dim.Lc; }
function lambdaFromDim() { return (dim.u * dim.u) / (dim.w * dim.r0 * dim.Lc); }
function thetaFromDim()  { return (dim.theta * dim.w) / (dim.u * dim.r0); }

// ─── simulation parameters ──────────────────────────────────────────────
const params = {
  L:    1.5,    // 𝓛
  lam:  1.0,    // λ
  tht:  0.05,   // ϑ
  dt:   0.01,
  mode: '1d',
  speed: 1.0,
};
const STEPS_PER_FRAME_BASE = 50;
let stepAccum = 0;

const rng = makeRng(42);
let p = [0.05, 0.0];
let t = 0;
let running = true;

const TRACE_N = 800;
const trace = { t: new Float64Array(TRACE_N), x: new Float64Array(TRACE_N), y: new Float64Array(TRACE_N), n: 0 };
const HIST_BINS = 40;
let HIST_PMAX = 1.6;
let hist = new Float64Array(HIST_BINS);
let histTotal = 0;
const TIP_TRAIL_N = 600;
const tip = { x: new Float64Array(TIP_TRAIL_N), y: new Float64Array(TIP_TRAIL_N), n: 0 };

function reset() {
  p[0] = 0.05 * (rng.uniform() * 2 - 1);
  p[1] = 0.05 * (rng.uniform() * 2 - 1);
  t = 0; trace.n = 0; tip.n = 0;
  hist = new Float64Array(HIST_BINS); histTotal = 0;
}

function autoExpandHistRange() {
  const m = Math.sqrt(p[0] * p[0] + p[1] * p[1]);
  if (m > 0.9 * HIST_PMAX) {
    const newMax = Math.max(HIST_PMAX, 1.2 * m);
    if (newMax !== HIST_PMAX) {
      HIST_PMAX = newMax;
      hist = new Float64Array(HIST_BINS); histTotal = 0;
    }
  }
}

function stepOnce() {
  const { L, lam, tht, dt, mode } = params;
  const sigmaSqrt = Math.sqrt(2 * tht * Math.max(L, 0) * dt);
  const px = p[0], py = mode === '2d' ? p[1] : 0;
  const m2 = px * px + py * py;
  const m4 = m2 * m2;
  // dP_α = [(𝓛-1) + λ(|P|² - |P|⁴)] P_α dt + noise
  const driftFactor = (L - 1) + lam * (m2 - m4);
  p[0] = px + driftFactor * px * dt + sigmaSqrt * rng.gauss();
  p[1] = mode === '2d' ? py + driftFactor * py * dt + sigmaSqrt * rng.gauss() : 0;
  t += dt;

  const i = trace.n % TRACE_N;
  trace.t[i] = t; trace.x[i] = p[0]; trace.y[i] = p[1];
  trace.n++;

  if (mode === '2d') {
    const j = tip.n % TIP_TRAIL_N;
    tip.x[j] = p[0]; tip.y[j] = p[1];
    tip.n++;
  }

  autoExpandHistRange();
  const mag = Math.sqrt(p[0] * p[0] + p[1] * p[1]);
  if (mag < HIST_PMAX) {
    const b = Math.min(HIST_BINS - 1, Math.floor(mag / HIST_PMAX * HIST_BINS));
    hist[b] += 1;
    histTotal++;
  }
}

// ─── F̃ helpers ──────────────────────────────────────────────────────────
function Fof(P, L, lam) {
  const P2 = P * P;
  return -0.5 * (L - 1) * P2 - 0.25 * lam * P2 * P2 + (1/6) * lam * P2 * P2 * P2;
}

function outerExtremum(L, lam) {
  if (lam <= 0) return 0;
  const u = (L - 1) / lam;
  const disc = 1 + 4 * u;
  if (disc < 0) return 0;
  const m = (1 + Math.sqrt(disc)) / 2;
  return m >= 0 ? Math.sqrt(m) : 0;
}
function innerExtremum(L, lam) {
  if (lam <= 0) return 0;
  const u = (L - 1) / lam;
  const disc = 1 + 4 * u;
  if (disc < 0) return -1;
  const m = (1 - Math.sqrt(disc)) / 2;
  return m >= 0 ? Math.sqrt(m) : -1;
}

// ─── plotting ───────────────────────────────────────────────────────────
function drawF() {
  const cv = document.getElementById('cv-F');
  const ctx = autoFit(cv);
  const w = cv.clientWidth, h = cv.clientHeight;
  const { L, lam } = params;
  const pOuter = outerExtremum(L, lam);
  const pmax = pOuter > 0 ? Math.max(1.0, 1.3 * pOuter) : 1.0;

  const N = 240;
  const xs = new Float64Array(2 * N + 1);
  const ys = new Float64Array(2 * N + 1);
  let yMin = Infinity, yMax = -Infinity;
  for (let i = -N; i <= N; i++) {
    const pp = (i / N) * pmax;
    const F = Fof(pp, L, lam);
    xs[i + N] = pp; ys[i + N] = F;
    if (F < yMin) yMin = F; if (F > yMax) yMax = F;
  }
  const pad = 0.15 * (yMax - yMin || 1e-3);
  const ax = makeAxis({ xMin: -pmax, xMax: pmax, yMin: yMin - pad, yMax: yMax + pad, w, h });
  drawFrame(ctx, ax);
  clipPlot(ctx, ax);
  strokePath(ctx, ax, xs, ys, { color: '#2b6cb0', width: 1.8 });
  const pos = params.mode === '2d' ? Math.sqrt(p[0] * p[0] + p[1] * p[1]) : p[0];
  dot(ctx, ax, pos, Fof(pos, L, lam), 5, '#b34700');
  ctx.restore();
}

function drawTrace() {
  const cv = document.getElementById('cv-trace');
  const ctx = autoFit(cv);
  const w = cv.clientWidth, h = cv.clientHeight;
  const n = Math.min(trace.n, TRACE_N);
  const frameOpts = { showTickLabelsX: false, showGridX: false, showGridY: false };
  if (n < 2) {
    const ax = makeAxis({ xMin: 0, xMax: 1, yMin: -1.5, yMax: 1.5, w, h });
    drawFrame(ctx, ax, frameOpts);
    return;
  }
  const start = trace.n >= TRACE_N ? trace.n % TRACE_N : 0;
  const ts = new Float64Array(n);
  const ys = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const k = (start + i) % TRACE_N;
    ts[i] = trace.t[k];
    if (params.mode === '2d') {
      const px = trace.x[k], py = trace.y[k];
      ys[i] = Math.sqrt(px * px + py * py);
    } else {
      ys[i] = trace.x[k];
    }
  }
  const tMin = ts[0], tMax = ts[n - 1];
  let pMin = Infinity, pMax = -Infinity;
  for (let i = 0; i < n; i++) {
    if (ys[i] < pMin) pMin = ys[i]; if (ys[i] > pMax) pMax = ys[i];
  }
  if (params.mode === '1d') {
    const m = Math.max(Math.abs(pMin), Math.abs(pMax), 0.5);
    pMin = -m; pMax = m;
  } else {
    pMin = 0; pMax = Math.max(pMax, 0.5);
  }
  const pad = 0.1 * (pMax - pMin || 1);
  const ax = makeAxis({ xMin: tMin, xMax: tMax, yMin: pMin - pad, yMax: pMax + pad, w, h });
  drawFrame(ctx, ax, frameOpts);
  clipPlot(ctx, ax);
  strokePath(ctx, ax, ts, ys, { color: '#2b6cb0' });
  ctx.restore();
}

function drawHist() {
  const cv = document.getElementById('cv-hist');
  const ctx = autoFit(cv);
  const w = cv.clientWidth, h = cv.clientHeight;
  let maxC = 1;
  for (const c of hist) if (c > maxC) maxC = c;
  const ax = makeAxis({ xMin: 0, xMax: HIST_PMAX, yMin: 0, yMax: 1.05, w, h });
  drawFrame(ctx, ax);
  clipPlot(ctx, ax);
  ctx.fillStyle = 'rgba(43,108,176,0.6)';
  for (let b = 0; b < HIST_BINS; b++) {
    const lo = b / HIST_BINS * HIST_PMAX;
    const hi = (b + 1) / HIST_BINS * HIST_PMAX;
    const px0 = ax.xToPx(lo), px1 = ax.xToPx(hi);
    const fy = ax.yToPx(hist[b] / maxC);
    const fy0 = ax.yToPx(0);
    ctx.fillRect(px0, fy, px1 - px0 - 1, fy0 - fy);
  }

  const Deff = params.tht * Math.max(params.L, 1e-9);
  if (Deff > 1e-6) {
    const N = 200;
    const xs = new Float64Array(N + 1);
    const ys = new Float64Array(N + 1);
    let m = 0;
    for (let i = 0; i <= N; i++) {
      const pp = i / N * HIST_PMAX;
      const F = Fof(pp, params.L, params.lam);
      let weight = Math.exp(-F / Deff);
      if (params.mode === '2d') weight *= pp;
      xs[i] = pp; ys[i] = weight;
      if (weight > m) m = weight;
    }
    if (m > 0) for (let i = 0; i <= N; i++) ys[i] /= m;
    strokePath(ctx, ax, xs, ys, { color: '#2f7a3a', width: 1.5, dash: [4, 3] });
  }
  ctx.restore();
}

function drawBif() {
  const cv = document.getElementById('cv-bif');
  const ctx = autoFit(cv);
  const w = cv.clientWidth, h = cv.clientHeight;
  // x-axis: 𝓛. Curves depend on λ.
  const Lcur = params.L, lam = params.lam;
  const Lmin = 0, Lmax = Math.max(3, 1.2 * Lcur);
  // y-range from outer extremum at the largest 𝓛 in view
  const yMax = Math.max(1.5, 1.3 * outerExtremum(Lmax, lam));
  const ax = makeAxis({ xMin: Lmin, xMax: Lmax, yMin: -yMax, yMax: yMax, w, h });
  drawFrame(ctx, ax);
  clipPlot(ctx, ax);

  // P=0 branch: stable for 𝓛 < 1, unstable for 𝓛 > 1.
  if (Lmin < 1) strokePath(ctx, ax,
    new Float64Array([Lmin, Math.min(1, Lmax)]), new Float64Array([0, 0]),
    { color: '#2b6cb0', width: 2 });
  if (Lmax > 1) strokePath(ctx, ax,
    new Float64Array([Math.max(1, Lmin), Lmax]), new Float64Array([0, 0]),
    { color: '#999', width: 1.6, dash: [5, 4] });

  // Outer extrema (stable minima): solid blue. Inner extrema (unstable
  // maxima): dashed gray.
  const N = 220;
  const Louter = [], yOpos = [], yOneg = [];
  const Linner = [], yIpos = [], yIneg = [];
  for (let i = 0; i <= N; i++) {
    const Lt = Lmin + (Lmax - Lmin) * i / N;
    const o = outerExtremum(Lt, lam);
    const inn = innerExtremum(Lt, lam);
    if (o > 0) { Louter.push(Lt); yOpos.push(o); yOneg.push(-o); }
    if (inn > 0) { Linner.push(Lt); yIpos.push(inn); yIneg.push(-inn); }
  }
  if (Louter.length > 1) {
    const xs = Float64Array.from(Louter);
    strokePath(ctx, ax, xs, Float64Array.from(yOpos), { color: '#2b6cb0', width: 2 });
    strokePath(ctx, ax, xs, Float64Array.from(yOneg), { color: '#2b6cb0', width: 2 });
  }
  if (Linner.length > 1) {
    const xs = Float64Array.from(Linner);
    strokePath(ctx, ax, xs, Float64Array.from(yIpos), { color: '#999', width: 1.6, dash: [5, 4] });
    strokePath(ctx, ax, xs, Float64Array.from(yIneg), { color: '#999', width: 1.6, dash: [5, 4] });
  }
  ctx.restore();

  // 𝓛 cursor
  ctx.save();
  ctx.strokeStyle = 'rgba(0,0,0,0.5)';
  ctx.setLineDash([2, 2]);
  ctx.beginPath();
  ctx.moveTo(ax.xToPx(Lcur), ax.padT);
  ctx.lineTo(ax.xToPx(Lcur), ax.padT + ax.plotH);
  ctx.stroke();
  ctx.restore();

  const pos = params.mode === '2d' ? Math.sqrt(p[0] * p[0] + p[1] * p[1]) : p[0];
  if (Math.abs(pos) <= yMax && Lcur >= Lmin && Lcur <= Lmax) {
    dot(ctx, ax, Lcur, pos, 4, '#b34700');
  }
}

// ─── ⟨|P|⟩ vs 𝓛 (numerical Boltzmann integral) ─────────────────────────
const MEANP_NL = 121;
const MEANP_NP = 400;
const MEANP_LMAX = 5;
const MEANP_PMAX = 3.5;
let meanPCache = null;        // { lam, tht, Ls, m1, m2 }

function recomputeMeanP() {
  const Ls = new Float64Array(MEANP_NL);
  const m1 = new Float64Array(MEANP_NL);
  const m2 = new Float64Array(MEANP_NL);
  const dp = MEANP_PMAX / MEANP_NP;
  const ps = new Float64Array(MEANP_NP + 1);
  for (let j = 0; j <= MEANP_NP; j++) ps[j] = j * dp;
  for (let i = 0; i < MEANP_NL; i++) {
    const Lt = (i / (MEANP_NL - 1)) * MEANP_LMAX;
    Ls[i] = Lt;
    const Deff = params.tht * Math.max(Lt, 1e-12);
    const Fs = new Float64Array(MEANP_NP + 1);
    let Fmin = Infinity;
    for (let j = 0; j <= MEANP_NP; j++) {
      Fs[j] = Fof(ps[j], Lt, params.lam);
      if (Fs[j] < Fmin) Fmin = Fs[j];
    }
    let Z1 = 0, N1 = 0, Z2 = 0, N2 = 0;
    for (let j = 0; j <= MEANP_NP; j++) {
      const w = Math.exp(-(Fs[j] - Fmin) / Deff);
      const trapz = (j === 0 || j === MEANP_NP) ? 0.5 : 1;
      const wT = w * trapz * dp;
      Z1 += wT;
      N1 += ps[j] * wT;
      Z2 += ps[j] * wT;
      N2 += ps[j] * ps[j] * wT;
    }
    m1[i] = Z1 > 0 ? N1 / Z1 : 0;
    m2[i] = Z2 > 0 ? N2 / Z2 : 0;
  }
  meanPCache = { lam: params.lam, tht: params.tht, Ls, m1, m2 };
}

function drawMeanP() {
  if (!meanPCache || meanPCache.lam !== params.lam || meanPCache.tht !== params.tht) {
    recomputeMeanP();
  }
  const cv = document.getElementById('cv-meanp');
  const ctx = autoFit(cv);
  const w = cv.clientWidth, h = cv.clientHeight;
  const { Ls, m1, m2 } = meanPCache;
  let yMax = 0;
  for (let i = 0; i < Ls.length; i++) {
    if (m1[i] > yMax) yMax = m1[i];
    if (m2[i] > yMax) yMax = m2[i];
  }
  yMax = Math.max(yMax * 1.15, 0.2);
  const ax = makeAxis({ xMin: 0, xMax: MEANP_LMAX, yMin: 0, yMax, w, h });
  drawFrame(ctx, ax);
  clipPlot(ctx, ax);
  strokePath(ctx, ax, Ls, m1, { color: '#2b6cb0', width: 2 });
  strokePath(ctx, ax, Ls, m2, { color: '#b34700', width: 2, dash: [5, 3] });
  ctx.restore();
  ctx.save();
  ctx.strokeStyle = 'rgba(0,0,0,0.5)';
  ctx.setLineDash([2, 2]);
  ctx.beginPath();
  ctx.moveTo(ax.xToPx(params.L), ax.padT);
  ctx.lineTo(ax.xToPx(params.L), ax.padT + ax.plotH);
  ctx.stroke();
  ctx.restore();
}

function draw2D() {
  const cv = document.getElementById('cv-2d');
  const ctx = autoFit(cv);
  const w = cv.clientWidth, h = cv.clientHeight;
  const pOuter = outerExtremum(params.L, params.lam);
  const R = Math.max(1.5, 1.3 * pOuter);
  const ax = makeAxis({ xMin: -R, xMax: R, yMin: -R, yMax: R, w, h, aspect: 1 });
  drawFrame(ctx, ax);
  clipPlot(ctx, ax);

  const cx = ax.xToPx(0), cy = ax.yToPx(0);
  const pxPerUnit = (ax.xToPx(1) - cx);
  if (pOuter > 0) {
    ctx.save();
    ctx.strokeStyle = '#2b6cb0';
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    ctx.arc(cx, cy, Math.abs(pxPerUnit * pOuter), 0, 2 * Math.PI);
    ctx.stroke();
    ctx.restore();
  }
  const pInner = innerExtremum(params.L, params.lam);
  if (pInner > 0) {
    ctx.save();
    ctx.strokeStyle = '#999';
    ctx.setLineDash([5, 4]);
    ctx.lineWidth = 1.4;
    ctx.beginPath();
    ctx.arc(cx, cy, Math.abs(pxPerUnit * pInner), 0, 2 * Math.PI);
    ctx.stroke();
    ctx.restore();
  }

  const n = Math.min(tip.n, TIP_TRAIL_N);
  if (n > 1) {
    const start = tip.n >= TIP_TRAIL_N ? tip.n % TIP_TRAIL_N : 0;
    ctx.save();
    ctx.strokeStyle = 'rgba(43,108,176,0.4)';
    ctx.lineWidth = 1;
    ctx.beginPath();
    for (let i = 0; i < n; i++) {
      const k = (start + i) % TIP_TRAIL_N;
      const x = ax.xToPx(tip.x[k]), y = ax.yToPx(tip.y[k]);
      if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    ctx.stroke();
    ctx.restore();
  }
  ctx.save();
  ctx.strokeStyle = '#b34700';
  ctx.lineWidth = 1.5;
  ctx.beginPath();
  ctx.moveTo(ax.xToPx(0), ax.yToPx(0));
  ctx.lineTo(ax.xToPx(p[0]), ax.yToPx(p[1]));
  ctx.stroke();
  ctx.restore();
  dot(ctx, ax, p[0], p[1], 5, '#b34700');
  ctx.restore();
}

// ─── controls ────────────────────────────────────────────────────────────
const controlsEl = document.getElementById('controls');
const kpis = makeKpis([
  { id: 'L',    label: '𝓛' },
  { id: 'lam',  label: 'λ' },
  { id: 'tht',  label: 'ϑ' },
  { id: 'p',    label: '|P|' },
]);

const sL    = makeSlider({ id: 'L',   symbol: '𝓛', value: params.L,   min: 0,    max: 5, step: 0.01, fmt: v => v.toFixed(2) });
const sLam  = makeSlider({ id: 'lam', symbol: 'λ', value: params.lam, min: 0.01, max: 10, log: true, fmt: v => v.toPrecision(3) });
const sTht  = makeSlider({ id: 'tht', symbol: 'ϑ', value: params.tht, min: 1e-4, max: 1,  log: true, fmt: v => v.toExponential(2) });
let applyingDim = false;
function refreshKpis() {
  kpis.set('L',   params.L.toFixed(3));
  kpis.set('lam', params.lam.toPrecision(3));
  kpis.set('tht', params.tht.toExponential(2));
}
sL.onChange(v   => { params.L   = v; refreshKpis(); });
sLam.onChange(v => { params.lam = v; refreshKpis(); });
sTht.onChange(v => { params.tht = v; });

function pushAllNondim() {
  applyingDim = true;
  sL.set(LfromDim()); sLam.set(lambdaFromDim()); sTht.set(thetaFromDim());
  applyingDim = false;
}

// Dim sliders (collapsed under a spoiler in the panel).
const linkedL  = () => `→ 𝓛 = ${LfromDim().toFixed(3)}`;
const linkedR0 = () => `→ λ = ${lambdaFromDim().toPrecision(3)},  ϑ = ${thetaFromDim().toExponential(2)}`;
const linkedTh = () => `→ ϑ = ${thetaFromDim().toExponential(2)}`;

const dimL  = makeSlider({ id: 'L_dim',  symbol: 'L',  value: dim.L,     min: 0,    max: 5,  step: 0.01, units: '[L_c]', fmt: v => v.toFixed(2), linkedLabel: linkedL });
const dimR0 = makeSlider({ id: 'r0',     symbol: 'r₀', value: dim.r0,    min: 0.1,  max: 5,  step: 0.01, fmt: v => v.toFixed(2), linkedLabel: linkedR0 });
const dimTh = makeSlider({ id: 'theta',  symbol: 'θ',  value: dim.theta, min: 1e-4, max: 1,  log: true,  fmt: v => v.toExponential(2), linkedLabel: linkedTh });

function refreshDimReadouts() {
  dimL.setLinkedText(linkedL()); dimR0.setLinkedText(linkedR0()); dimTh.setLinkedText(linkedTh());
}
dimL.onChange(v  => { dim.L     = v; pushAllNondim(); refreshDimReadouts(); refreshKpis(); });
dimR0.onChange(v => { dim.r0    = v; pushAllNondim(); refreshDimReadouts(); refreshKpis(); });
dimTh.onChange(v => { dim.theta = v; pushAllNondim(); refreshDimReadouts(); refreshKpis(); });

const sDt    = makeSlider({ id: 'dt',    symbol: 'dt̃',        value: params.dt,    min: 1e-4, max: 0.05, log: true, fmt: v => v.toExponential(2) });
const sSpeed = makeSlider({ id: 'speed', symbol: 'sim speed', value: params.speed, min: 0.01, max: 100,  log: true, fmt: v => `${v.toPrecision(2)}×` });
const sSeed  = makeSlider({ id: 'seed',  symbol: 'seed',      value: 42,           min: 1,    max: 9999, step: 1,   fmt: v => v.toFixed(0) });
sDt.onChange(v => { params.dt = v; });
sSpeed.onChange(v => { params.speed = v; });
sSeed.onChange(v => { rng.seed(Math.round(v)); });

const modeToggle = makeToggle({
  label: 'polarization dimension',
  options: [{ id: '1d', label: '1D' }, { id: '2d', label: '2D' }],
  value: '1d',
  onChange: v => {
    params.mode = v;
    document.getElementById('cell-2d').style.display = v === '2d' ? '' : 'none';
    if (v === '1d') p[1] = 0;
    tip.n = 0;
  },
});

const buttons = makeButtonRow([
  { label: '⏸  pause', onClick() {
      running = !running;
      buttons.refs['⏸  pause'].textContent = running ? '⏸  pause' : '▶  play';
    } },
  { label: '⟳  reset', ghost: true, onClick: reset },
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

// ─── sync from possibly-restored slider values ─────────────────────────
// makeSlider may load saved values from localStorage that differ from the
// hardcoded `params`/`dim` defaults; without a sync, the simulation would
// run with stale defaults until the user nudges every slider.
params.L     = sL.value;
params.lam   = sLam.value;
params.tht   = sTht.value;
params.dt    = sDt.value;
params.speed = sSpeed.value;
dim.L     = dimL.value;
dim.r0    = dimR0.value;
dim.theta = dimTh.value;
rng.seed(Math.round(sSeed.value));
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
}
if (window.katex) decorateAll();
else window.addEventListener('load', decorateAll);

// ─── animation loop ─────────────────────────────────────────────────────
function frame() {
  if (running) {
    stepAccum += STEPS_PER_FRAME_BASE * params.speed;
    const N = Math.floor(stepAccum);
    stepAccum -= N;
    for (let i = 0; i < N; i++) stepOnce();
  }
  drawF();
  drawTrace();
  drawHist();
  drawBif();
  drawMeanP();
  if (params.mode === '2d') draw2D();
  refreshKpis();
  kpis.set('p', Math.sqrt(p[0] * p[0] + p[1] * p[1]).toFixed(3));
  requestAnimationFrame(frame);
}
reset();
refreshKpis();
requestAnimationFrame(frame);
