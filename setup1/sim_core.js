// Setup 1 — Polarization in a uniform cue: reusable simulation + drawing core.
// Per docs/physics/setup1_uniform.md:
//   dP_α/dt̃ = (𝓛 - 1) P_α + λ (|P|² P_α - |P|⁴ P_α) + √(2 ϑ 𝓛) η_α
// Three independent nondim knobs: 𝓛, λ, ϑ.
//
// This module holds ALL of the setup-1 physics and Canvas2D drawing so that
// both the website page (setup1/main.js) and the presentation slide
// (slides/setup1-slide.js) share a single source of truth — no duplicated
// physics, no slider/simulator drift. Callers own only their layout + controls
// and drive the returned object's step()/draw*(canvas) methods from their own
// animation loop.

import { makeRng } from '../shared/rng.js';
import { autoFit, makeAxis, drawFrame, strokePath, dot, clipPlot } from '../shared/canvas.js';

export function createSetup1({ seed = 42 } = {}) {
  // ─── nondim ↔ dim linkage ──────────────────────────────────────────────
  // Dim params: u, w, Lc, r0, theta, L. Derived nondim:
  //   𝓛 = L / Lc,  λ = u² / (w r0 Lc),  ϑ = θ w / (u r0)
  const dim = { u: 1.0, w: 1.0, Lc: 1.0, r0: 1.0, theta: 0.05, L: 1.5 };
  const linkage = {
    LfromDim()      { return dim.L / dim.Lc; },
    lambdaFromDim() { return (dim.u * dim.u) / (dim.w * dim.r0 * dim.Lc); },
    thetaFromDim()  { return (dim.theta * dim.w) / (dim.u * dim.r0); },
  };

  // ─── simulation parameters ────────────────────────────────────────────
  const params = { L: 1.5, lam: 1.0, tht: 0.05, dt: 0.01, mode: '1d', speed: 1.0 };

  const rng = makeRng(seed);
  let p = [0.05, 0.0];
  let t = 0;

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

  // Advance the simulation by `nSteps` explicit Euler–Maruyama steps.
  function step(nSteps) {
    for (let i = 0; i < nSteps; i++) stepOnce();
  }

  // ─── F̃ helpers ────────────────────────────────────────────────────────
  function Fof(P, L, lam) {
    const P2 = P * P;
    return -0.5 * (L - 1) * P2 - 0.25 * lam * P2 * P2 + (1 / 6) * lam * P2 * P2 * P2;
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

  function magnitude() { return Math.sqrt(p[0] * p[0] + p[1] * p[1]); }

  // ─── plotting (each takes the target <canvas> element) ─────────────────
  function drawF(cv) {
    if (!cv) return;
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
    const pos = params.mode === '2d' ? magnitude() : p[0];
    dot(ctx, ax, pos, Fof(pos, L, lam), 5, '#b34700');
    ctx.restore();
  }

  function drawTrace(cv) {
    if (!cv) return;
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

  function drawHist(cv) {
    if (!cv) return;
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

  function drawBif(cv) {
    if (!cv) return;
    const ctx = autoFit(cv);
    const w = cv.clientWidth, h = cv.clientHeight;
    // x-axis: 𝓛. Curves depend on λ.
    const Lcur = params.L, lam = params.lam;
    const Lmin = 0, Lmax = Math.max(3, 1.2 * Lcur);
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

    const pos = params.mode === '2d' ? magnitude() : p[0];
    if (Math.abs(pos) <= yMax && Lcur >= Lmin && Lcur <= Lmax) {
      dot(ctx, ax, Lcur, pos, 4, '#b34700');
    }
  }

  // ─── ⟨|P|⟩ vs 𝓛 (numerical Boltzmann integral) ────────────────────────
  const MEANP_NL = 121;
  const MEANP_NP = 400;
  const MEANP_LMAX = 5;
  const MEANP_PMAX = 3.5;
  let meanPCache = null; // { lam, tht, Ls, m1, m2 }

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

  // opts.show2D === false hides the 2D (radial-Jacobian) curve — used by the
  // 1D-only presentation slide, where a 2D branch would be meaningless.
  function drawMeanP(cv, opts = {}) {
    if (!cv) return;
    const show2D = opts.show2D !== false;
    if (!meanPCache || meanPCache.lam !== params.lam || meanPCache.tht !== params.tht) {
      recomputeMeanP();
    }
    const ctx = autoFit(cv);
    const w = cv.clientWidth, h = cv.clientHeight;
    const { Ls, m1, m2 } = meanPCache;
    let yMax = 0;
    for (let i = 0; i < Ls.length; i++) {
      if (m1[i] > yMax) yMax = m1[i];
      if (show2D && m2[i] > yMax) yMax = m2[i];
    }
    yMax = Math.max(yMax * 1.15, 0.2);
    const ax = makeAxis({ xMin: 0, xMax: MEANP_LMAX, yMin: 0, yMax, w, h });
    drawFrame(ctx, ax);
    clipPlot(ctx, ax);
    strokePath(ctx, ax, Ls, m1, { color: '#2b6cb0', width: 2 });
    if (show2D) strokePath(ctx, ax, Ls, m2, { color: '#b34700', width: 2, dash: [5, 3] });
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

  function draw2D(cv) {
    if (!cv) return;
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

  function setMode(mode) {
    params.mode = mode;
    if (mode === '1d') p[1] = 0;
    tip.n = 0;
  }

  return {
    params, dim, rng, linkage,
    get p() { return p; },
    magnitude,
    step, reset, setMode,
    setSeed(s) { rng.seed(s); },
    drawF, drawBif, drawTrace, drawHist, drawMeanP, draw2D,
  };
}
