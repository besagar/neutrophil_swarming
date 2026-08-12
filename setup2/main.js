// Setup 2 — Single cell vs. 1D Gaussian running wave.
// Per docs/physics/setup2_wave.md, in nondim units (same notation as Setup 1
// plus distance rescaled by the wave width σ):
//   𝓛(x̃, t̃) = M · exp(-(x̃ - C t̃)² / 2)
//   ∂_x̃ 𝓛   = -𝓛 (x̃ - C t̃)
//   dP/dt̃   = χ̃ ∂_x̃ 𝓛 + (𝓛 - 1) P + λ (P² - P⁴) P + √(2 ϑ 𝓛) η
//   dx̃/dt̃   = μ̃ P
// Time, polarization and intrinsic groups (λ, ϑ) are identical to Setup 1.
//
// Hill / linear-relaxation alternative (docs/physics/hill_model.md):
//   dP/dt̃   = χ̃ ∂_x̃ 𝓛 − P + √(2 ϑ 𝓛) η   (replaces bistable well)
//   dx̃/dt̃   = ṽ₀ · |P|ⁿ/(1+|P|ⁿ) · sign(P)  (replaces linear μ̃ P)
// χ̃ and ϑ are SHARED between the two models (same slider, same value).
//
// The full trajectory on t̃ ∈ [-X_init/C, +X_init/C] is precomputed
// (deterministic up to the seed); the user scrubs through it via a time
// slider. Recomputed when any sim parameter changes.

import { makeRng } from '../shared/rng.js';
import { el, makeSlider, makeToggle, makeButtonRow, makeKpis, section, detailsSection, decoratePlot } from '../shared/dom.js';
import { autoFit, makeAxis, drawFrame, strokePath, dot, clipPlot } from '../shared/canvas.js';
import { attachSvgExports } from '../shared/svgexport.js';

const X_INIT = 8;          // wave starts at x̃ = -X_INIT, ends at +X_INIT
const X_HALF = 12;         // visualization x̃ range

// ─── nondim parameters ──────────────────────────────────────────────────
const params = {
  // wave
  M: 2.0, C: 1.5,
  // GL-only intrinsic well
  lam: 1.0,
  // shared: noise, chemotactic coupling
  tht: 1e-4, chi: 0.75,
  // GL-only motility
  mu: 0.05,
  // Hill-only: peak speed, Hill exponent
  v0: 0.5, nHill: 4,
  // numerics + UI state
  dt: 0.01, speed: 1.0, showRef: true, seed: 7,
  model: 'GL',   // 'GL' | 'Hill'
};

// Adaptive trajectory window: C → [-X_INIT/C, +X_INIT/C].
function trajTMin() { return -X_INIT / Math.max(params.C, 1e-6); }
function trajTMax() { return +X_INIT / Math.max(params.C, 1e-6); }

// ─── dim ↔ nondim linkage ───────────────────────────────────────────────
// Dim quantities: L_max, σ, c (wave-specific). Linkage:
//   M  = L_max / L_c
//   C  = c / (σ r_0 L_c)
//   χ̃  = χ  / (σ r_0 √(u/w))
//   μ̃  = μ √(u/w) / (σ r_0 L_c)
// Doubling L_max doubles M.
// Doubling c doubles C.
// Doubling σ halves C, χ̃, μ̃.
const dim = {
  Lmax: 1.0, sigma: 1.0, c: 1.0,
  cal: { Lmax: 1.0, sigma: 1.0, c: 1.0,
         M: 2.0, C: 1.5, chi: 0.75, mu: 0.05 },
};
function recomputeFromDim() {
  const sLmax  = dim.Lmax  / dim.cal.Lmax;
  const sC     = dim.c     / dim.cal.c;
  const sSigma = dim.sigma / dim.cal.sigma;
  params.M     = dim.cal.M     * sLmax;
  params.C = dim.cal.C * sC / sSigma;
  params.chi   = dim.cal.chi   / sSigma;
  params.mu    = dim.cal.mu    / sSigma;
}
function recalibrate() {
  if (applyingDim) return;
  dim.cal = {
    Lmax: dim.Lmax, sigma: dim.sigma, c: dim.c,
    M: params.M, C: params.C, chi: params.chi, mu: params.mu,
  };
}
let applyingDim = false;

// ─── trajectory storage ─────────────────────────────────────────────────
let traj = null;             // { ts, xs, ps, gs, ss, Ls, xR, vcArr, ranges, n }
let trajDirty = true;
let sweepDirty = true;
function markDirty() { trajDirty = true; sweepDirty = true; }

function L_at(xx, tt) { const u = xx - params.C * tt; return params.M * Math.exp(-0.5 * u * u); }
function dxL_at(xx, tt) { const u = xx - params.C * tt; return -params.M * Math.exp(-0.5 * u * u) * u; }

// ─── model-specific drift and velocity functions ─────────────────────────
// GL polarization drift: χ̃ ∂_x̃𝓛 + (𝓛-1)P + λ(P²-P⁴)P
function driftP_GL(p, x, t, C, M, chi, lam) {
  const u = x - C * t;
  const L = M * Math.exp(-0.5 * u * u);
  const dxL = -L * u;
  const m2 = p * p, m4 = m2 * m2;
  return chi * dxL + (L - 1) * p + lam * (m2 - m4) * p;
}

// Hill polarization drift: χ̃ ∂_x̃𝓛 − P  (linear relaxation, no bistable well)
function driftP_Hill(p, x, t, C, M, chi) {
  const u = x - C * t;
  const L = M * Math.exp(-0.5 * u * u);
  const dxL = -L * u;
  return chi * dxL - p;
}

// GL cell velocity: dx̃/dt̃ = μ̃ P
function vOfP_GL(p, mu) { return mu * p; }

// Hill cell velocity: dx̃/dt̃ = ṽ₀ · |P|ⁿ/(1+|P|ⁿ) · sign(P)
function vOfP_Hill(p, v0, n) {
  const a = Math.abs(p);
  if (a < 1e-12) return 0;
  const an = Math.pow(a, n);
  return v0 * (an / (1 + an)) * Math.sign(p);
}

// ALL roots of dF̃_eff/dP = 0, i.e. of the GL drift
//   g(P) = χ̃∂_x̃𝓛 + (𝓛-1)P + λ(P³-P⁵),
// sorted ascending, each tagged stable (g' < 0 → minimum of F̃_eff) or unstable
// (g' > 0 → barrier). Up to five roots, so the adiabatic locus is multivalued in
// the folded region — we keep every branch rather than picking one, so the plot
// shows the full S-curve instead of a jump at the fold.
// Grid scan for sign changes + bisection refine: cheap at ~1e3 samples and immune
// to Newton runaway near g' ≈ 0 (exactly where the folds are).
function adiabaticRootsGL(L, dxL, chi, lam, out, stabOut) {
  const g  = P => chi * dxL + (L - 1) * P + lam * (P * P * P - P * P * P * P * P);
  const gp = P => (L - 1) + lam * (3 * P * P - 5 * P * P * P * P);
  // Bracket: outer extremum of the symmetric well sets the scale; pad for the tilt.
  const R = Math.max(2, 1.5 * outerExtremum(Math.max(L, 1), lam) + 1);
  const NG = 600;
  let m = 0;
  let Pa = -R, ga = g(Pa);
  for (let k = 1; k <= NG && m < 5; k++) {
    const Pb = -R + 2 * R * k / NG, gb = g(Pb);
    if ((ga <= 0 && gb > 0) || (ga >= 0 && gb < 0)) {
      let lo = Pa, hi = Pb, glo = ga;
      for (let it = 0; it < 40; it++) {
        const mid = 0.5 * (lo + hi), gm = g(mid);
        if ((glo <= 0) === (gm <= 0)) { lo = mid; glo = gm; } else hi = mid;
      }
      const root = 0.5 * (lo + hi);
      out[m] = root;
      stabOut[m] = gp(root) < 0 ? 1 : 0;
      m++;
    }
    Pa = Pb; ga = gb;
  }
  return m;
}

// Chain the per-time root sets into continuous branches of the adiabatic locus.
// Greedy nearest-neighbour matching between consecutive time slices: a root joins
// the branch whose last value it is closest to (within TOL); otherwise it starts a
// new branch. Branches are born and die at the folds, which is precisely the shape
// we want to draw. TOL is in P units — roots drift slowly between samples, but a
// newly born pair appears a finite distance away, so a modest tolerance separates
// "same branch" from "new branch" without tuning.
const BRANCH_TOL = 0.25;
function buildAdiabaticBranches(ts, rootVals, rootStab, rootCnt, N) {
  const branches = [];
  let active = [];   // {ts:[], ps:[], stab:[]}
  for (let i = 0; i < N; i++) {
    const m = rootCnt[i], base = i * 5;
    const next = [];
    const claimed = new Uint8Array(active.length);
    for (let j = 0; j < m; j++) {
      const P = rootVals[base + j];
      let bi = -1, bd = BRANCH_TOL;
      for (let a = 0; a < active.length; a++) {
        if (claimed[a]) continue;
        const d = Math.abs(P - active[a].ps[active[a].ps.length - 1]);
        if (d < bd) { bd = d; bi = a; }
      }
      const br = bi >= 0 ? (claimed[bi] = 1, active[bi]) : { ts: [], ps: [], stab: [] };
      br.ts.push(ts[i]); br.ps.push(P); br.stab.push(rootStab[base + j]);
      next.push(br);
    }
    // Branches not extended this step have ended (fold annihilation).
    for (let a = 0; a < active.length; a++)
      if (!claimed[a] && next.indexOf(active[a]) < 0) branches.push(active[a]);
    active = next;
  }
  for (const br of active) branches.push(br);
  return branches.filter(b => b.ps.length > 1);
}

function recomputeTrajectory() {
  const rng = makeRng(params.seed);
  const dt = params.dt;
  const C = params.C, M = params.M, chi = params.chi, mu = params.mu;
  const lam = params.lam, tht = params.tht;
  const v0 = params.v0, n = params.nHill;
  const model = params.model;
  const tMin = trajTMin(), tMax = trajTMax();
  const N = Math.ceil((tMax - tMin) / dt) + 1;
  const ts = new Float64Array(N);
  const xs = new Float64Array(N);
  const ps = new Float64Array(N);
  const gs = new Float64Array(N);
  const Ls = new Float64Array(N);   // 𝓛 at the cell (for the 𝓛–s phase plot)
  const xR = new Float64Array(N);
  let p = 0, x = 0, xRef = 0, t = tMin;
  // Initial values (inline fused exp for cell gradient)
  { const u0 = x - C * t; const L0 = M * Math.exp(-0.5 * u0 * u0);
    ts[0] = t; xs[0] = x; ps[0] = p; gs[0] = -L0 * u0; Ls[0] = L0; xR[0] = xRef; }
  for (let i = 1; i < N; i++) {
    // Single exp for cell position (was two separate L_at + dxL_at calls with the same u)
    const u = x - C * t;
    const L = M * Math.exp(-0.5 * u * u);
    const dxL = -L * u;
    // Euler–Maruyama integration. Drift dispatches on model.
    const drift = model === 'GL'
      ? chi * dxL + (L - 1) * p + lam * (p * p - p * p * p * p) * p
      : chi * dxL - p;
    const noise = Math.sqrt(2 * tht * Math.max(L, 0) * dt) * rng.gauss();
    gs[i] = dxL;  // store gradient at pre-step position; 1-step lag invisible at dt ≤ 0.05
    Ls[i] = L;    // 𝓛 at the same (pre-step) position, consistent with gs
    p += drift * dt + noise;
    x += (model === 'GL' ? vOfP_GL(p, mu) : vOfP_Hill(p, v0, n)) * dt;
    if (params.showRef) {
      // Adiabatic reference: P_eq = χ̃ ∂_x̃𝓛 / (1 − 𝓛)  [GL],
      //                  or  P_eq = χ̃ ∂_x̃𝓛               [Hill, balance χ̃∂𝓛 − P = 0].
      const uR = xRef - C * t;
      const Lr = M * Math.exp(-0.5 * uR * uR);
      if (model === 'GL') {
        const denom = Math.max(Math.abs(Lr - 1), 0.05);
        xRef += mu * (chi * Lr * uR / denom) * dt;
      } else {
        // Hill: P_eq = χ̃ · (-𝓛 · uR), then Hill velocity
        const PeqHill = chi * (-Lr * uR);
        xRef += vOfP_Hill(PeqHill, v0, n) * dt;
      }
    }
    t += dt;
    ts[i] = t; xs[i] = x; ps[i] = p; xR[i] = xRef;
  }

  // Chemotactic-drive diagnostic s = χ̃ ∂_x̃𝓛 (the linear forcing of the P-SDE).
  // gs already holds ∂_x̃𝓛 at the cell, so s = χ̃ · gs. Plotted parametrically
  // against P as the "phase trajectory" in P–s space.
  const ss = new Float64Array(N);
  for (let i = 0; i < N; i++) ss[i] = chi * gs[i];

  // Reduced-drive variable q = P + χ̃𝓛/C and response ratio U = μ̃P / (relaxational drift),
  // i.e. U ≈ dx̃/dq (see docs/physics/setup2_wave.md). For GL the P in numerator and
  // denominator cancels analytically — use the cancelled form so P = 0 is not 0/0.
  // U still diverges where the bracket crosses zero (marginal point of F̃_eff); the
  // U-range below is therefore a robust percentile range, and the plot clips.
  const qs = new Float64Array(N);
  const Us = new Float64Array(N);
  const Cnz0 = Math.max(C, 1e-9);
  const EPS = 1e-6;
  for (let i = 0; i < N; i++) {
    qs[i] = ps[i] + chi * Ls[i] / Cnz0;
    let den, num;
    if (model === 'GL') {
      const p2 = ps[i] * ps[i];
      den = (Ls[i] - 1) + lam * (p2 - p2 * p2);
      num = mu;
    } else {
      den = -ps[i];
      num = vOfP_Hill(ps[i], v0, n);
    }
    if (Math.abs(den) < EPS) den = den < 0 ? -EPS : EPS;
    Us[i] = num / den;
  }

  // Adiabatic locus: every root of dF̃_eff/dP = 0 at the cell's instantaneous cue,
  // i.e. the full multivalued equilibrium curve (stable branches + the unstable
  // barrier that connects them through the folds). In Hill mode there is no F̃;
  // the analogous quasi-static balance χ̃∂_x̃𝓛 − P = 0 is single-valued.
  let adBranches, adLo = Infinity, adHi = -Infinity;
  if (model === 'GL') {
    const rootVals = new Float64Array(N * 5);
    const rootStab = new Uint8Array(N * 5);
    const rootCnt  = new Uint8Array(N);
    const tmpV = new Float64Array(5), tmpS = new Uint8Array(5);
    for (let i = 0; i < N; i++) {
      const m = adiabaticRootsGL(Ls[i], gs[i], chi, lam, tmpV, tmpS);
      rootCnt[i] = m;
      for (let j = 0; j < m; j++) {
        rootVals[i * 5 + j] = tmpV[j];
        rootStab[i * 5 + j] = tmpS[j];
        if (tmpV[j] < adLo) adLo = tmpV[j];
        if (tmpV[j] > adHi) adHi = tmpV[j];
      }
    }
    adBranches = buildAdiabaticBranches(ts, rootVals, rootStab, rootCnt, N);
  } else {
    adBranches = [{ ts: Array.from(ts), ps: Array.from(ss), stab: new Array(N).fill(1) }];
    for (let i = 0; i < N; i++) { if (ss[i] < adLo) adLo = ss[i]; if (ss[i] > adHi) adHi = ss[i]; }
  }

  // Precompute vcArr: GL → μ̃·P/C; Hill → vOfP_Hill(P)/C
  const Cnz = Math.max(C, 1e-9);
  const vcArr = new Float64Array(N);
  if (model === 'GL') {
    const vcScale = mu / Cnz;
    for (let i = 0; i < N; i++) vcArr[i] = vcScale * ps[i];
  } else {
    for (let i = 0; i < N; i++) vcArr[i] = vOfP_Hill(ps[i], v0, n) / Cnz;
  }

  // Precompute per-channel value ranges; avoids O(n) scan inside each draw call every frame.
  function chanRange(arr) {
    let lo = Infinity, hi = -Infinity;
    for (let i = 0; i < N; i++) { if (arr[i] < lo) lo = arr[i]; if (arr[i] > hi) hi = arr[i]; }
    return { lo, hi };
  }
  // U can diverge (denominator crosses zero), so its axis range comes from
  // percentiles rather than min/max, otherwise one spike flattens the whole loop.
  function robustRange(arr, frac = 0.02) {
    // Subsample before sorting: dt can be 1e-4 (N ≈ 1.6e5) and this runs on every
    // recompute, i.e. potentially every frame while a slider is dragged.
    const MAXS = 8192;
    const stride = Math.max(1, Math.ceil(arr.length / MAXS));
    const m = Math.ceil(arr.length / stride);
    const s = new Float64Array(m);
    for (let i = 0, j = 0; i < arr.length; i += stride, j++) s[j] = arr[i];
    s.sort();   // typed-array sort is numeric
    if (m === 0) return { lo: -1, hi: 1 };
    const k = Math.floor(frac * (m - 1));
    return { lo: s[k], hi: s[m - 1 - k] };
  }
  const ranges = {
    xs: chanRange(xs), ps: chanRange(ps), gs: chanRange(gs),
    xR: chanRange(xR), vcArr: chanRange(vcArr), ss: chanRange(ss), Ls: chanRange(Ls),
    qs: chanRange(qs), Us: robustRange(Us),
    adiabatic: { lo: adLo, hi: adHi },
  };

  traj = { ts, xs, ps, gs, ss, Ls, qs, Us, xR, vcArr, adBranches, ranges, n: N };
  trajDirty = false;
  // After recompute, clamp current scrub time into the new window and resync slider range.
  sTime.setMinMax(tMin, tMax);
  if (currentTime < tMin) currentTime = tMin;
  if (currentTime > tMax) currentTime = tMax;
  sTime.set(currentTime);
}

// ─── current scrub time ─────────────────────────────────────────────────
let currentTime = -X_INIT;   // initial window with C = 1
let playing = false;
let lastFrameMs = 0;

function indexAt(t) {
  if (!traj) return 0;
  const tMin = trajTMin(), tMax = trajTMax();
  const dt = (tMax - tMin) / (traj.n - 1);
  const i = Math.round((t - tMin) / dt);
  return Math.max(0, Math.min(traj.n - 1, i));
}

// ─── plotting ───────────────────────────────────────────────────────────
function drawWave() {
  const cv = document.getElementById('cv-wave');
  const ctx = autoFit(cv);
  const w = cv.clientWidth, h = cv.clientHeight;
  const ax = makeAxis({ xMin: -X_HALF, xMax: X_HALF, yMin: -0.05, yMax: Math.max(1.1, params.M * 1.15), w, h });
  drawFrame(ctx, ax);
  clipPlot(ctx, ax);
  const N = 400;
  const xs = new Float64Array(N + 1);
  const ys = new Float64Array(N + 1);
  for (let i = 0; i <= N; i++) {
    const xx = -X_HALF + (2 * X_HALF) * i / N;
    xs[i] = xx; ys[i] = L_at(xx, currentTime);
  }
  strokePath(ctx, ax, xs, ys, { color: '#2b6cb0', width: 1.6 });
  // 𝓛 = 1 reference line (the polarization threshold)
  strokePath(ctx, ax, new Float64Array([-X_HALF, X_HALF]), new Float64Array([1, 1]),
    { color: 'rgba(0,0,0,0.3)', width: 1, dash: [3, 3] });
  if (!traj) { ctx.restore(); return; }
  const i = indexAt(currentTime);
  const x = traj.xs[i], xRef = traj.xR[i];
  if (x >= -X_HALF && x <= X_HALF) dot(ctx, ax, x, L_at(x, currentTime), 6, '#b34700');
  if (params.showRef && xRef >= -X_HALF && xRef <= X_HALF)
    dot(ctx, ax, xRef, L_at(xRef, currentTime), 5, 'rgba(47,122,58,0.8)');
  ctx.restore();
}

// Symmetric GL free energy: F̃(P; 𝓛, λ) = -½(𝓛-1)P² - ¼λP⁴ + ⅙λP⁶
function Fof(P, L, lam) {
  const P2 = P * P;
  return -0.5 * (L - 1) * P2 - 0.25 * lam * P2 * P2 + (1/6) * lam * P2 * P2 * P2;
}
// Effective free energy including gradient coupling: F̃_eff = F̃ - χ̃ · ∂_x̃𝓛 · P
// The linear tilt term breaks P → -P symmetry; the full drift is -dF̃_eff/dP.
function FofEff(P, L, dxL, chi, lam) {
  return Fof(P, L, lam) - chi * dxL * P;
}
// Location of the outer (non-zero) extremum of the symmetric F̃, used to set pmax.
function outerExtremum(L, lam) {
  if (lam <= 0) return 0;
  const u = (L - 1) / lam;
  const disc = 1 + 4 * u;
  if (disc < 0) return 0;
  const m = (1 + Math.sqrt(disc)) / 2;
  return m >= 0 ? Math.sqrt(m) : 0;
}

function drawF() {
  if (params.model !== 'GL') return;   // hidden in Hill mode via applyModelVisibility
  const cv = document.getElementById('cv-F');
  const ctx = autoFit(cv);
  const w = cv.clientWidth, h = cv.clientHeight;
  if (!traj) {
    const ax = makeAxis({ xMin: -1, xMax: 1, yMin: -1, yMax: 1, w, h });
    drawFrame(ctx, ax);
    return;
  }
  const i = indexAt(currentTime);
  const xCell = traj.xs[i];
  const Lcell = L_at(xCell, currentTime);
  const dxLcell = dxL_at(xCell, currentTime);
  const lam = params.lam, chi = params.chi;

  // pmax: based on peak-cue well (constant as wave passes → axis stays stable),
  // expanded to always show the current bead position.
  const pOuterPeak = outerExtremum(params.M, lam);
  const pmax = Math.max(1.5, 1.5 * pOuterPeak, Math.abs(traj.ps[i]) * 1.2);

  const NN = 240;
  const xs   = new Float64Array(2 * NN + 1);
  const ysEff = new Float64Array(2 * NN + 1);  // tilted F̃_eff at current (𝓛, ∂_x̃𝓛)
  const ysRef = new Float64Array(2 * NN + 1);  // symmetric F̃ at 𝓛 = M (reference)
  let yMin = Infinity, yMax = -Infinity;
  for (let j = -NN; j <= NN; j++) {
    const pp = (j / NN) * pmax;
    xs[j + NN] = pp;
    ysEff[j + NN] = FofEff(pp, Lcell, dxLcell, chi, lam);
    ysRef[j + NN] = Fof(pp, params.M, lam);
    if (ysEff[j + NN] < yMin) yMin = ysEff[j + NN];
    if (ysEff[j + NN] > yMax) yMax = ysEff[j + NN];
    if (ysRef[j + NN] < yMin) yMin = ysRef[j + NN];
    if (ysRef[j + NN] > yMax) yMax = ysRef[j + NN];
  }
  const pad = 0.12 * (yMax - yMin || 1e-3);
  const ax = makeAxis({ xMin: -pmax, xMax: pmax, yMin: yMin - pad, yMax: yMax + pad, w, h });
  drawFrame(ctx, ax);
  clipPlot(ctx, ax);

  // Dashed: symmetric F̃ at 𝓛 = M — deepest untilted well, constant reference
  strokePath(ctx, ax, xs, ysRef, { color: 'rgba(43,108,176,0.25)', width: 1.2, dash: [3, 3] });
  // Solid: F̃_eff = F̃ - χ̃ · ∂_x̃𝓛 · P — actual landscape including gradient tilt
  strokePath(ctx, ax, xs, ysEff, { color: '#2b6cb0', width: 1.8 });

  const Pcur = traj.ps[i];
  dot(ctx, ax, Pcur, FofEff(Pcur, Lcell, dxLcell, chi, lam), 5, '#b34700');
  ctx.restore();
}

// Gray background: the adiabatic locus dF̃_eff/dP = 0. Each branch is split into
// runs of constant stability so minima draw solid and the connecting barrier draws
// dashed; the boundary sample is repeated in both runs so the S-curve stays closed
// through the fold.
function drawAdiabatic(ctx, ax, branches) {
  const STABLE   = { color: 'rgba(0,0,0,0.30)', width: 2.5 };
  const UNSTABLE = { color: 'rgba(0,0,0,0.18)', width: 1.4, dash: [4, 3] };
  for (const br of branches) {
    let s = 0;
    for (let e = 1; e <= br.ps.length; e++) {
      if (e < br.ps.length && br.stab[e] === br.stab[s]) continue;
      const end = Math.min(e + 1, br.ps.length);   // overlap one sample
      strokePath(ctx, ax, br.ts.slice(s, end), br.ps.slice(s, end),
                 br.stab[s] ? STABLE : UNSTABLE);
      s = e;
    }
  }
}

function drawTrajChannel(canvasId, channel, opts = {}) {
  const cv = document.getElementById(canvasId);
  const ctx = autoFit(cv);
  const w = cv.clientWidth, h = cv.clientHeight;
  const tMin = trajTMin(), tMax = trajTMax();
  if (!traj) {
    const ax = makeAxis({ xMin: tMin, xMax: tMax, yMin: -1, yMax: 1, w, h });
    drawFrame(ctx, ax);
    return;
  }
  const ts = traj.ts, ys = traj[channel];
  // y-range adapts to the primary channel only; extraChannel draws as dashed reference
  // but does not expand the range.
  let { lo, hi } = traj.ranges[channel] || { lo: -1, hi: 1 };
  if (!isFinite(lo)) { lo = -1; hi = 1; }
  // The adiabatic locus may extend the range, but only by a bounded margin: its
  // outer branches sit at |P| ~ ((𝓛-1)/λ)^{1/4}, which diverges as λ → 0 and would
  // otherwise hijack the y-axis — making the P trace appear to shrink by 5× across
  // a λ sweep that in fact changes P by ~15%. The axis stays owned by the trace;
  // locus branches outside the margin are clipped (clipPlot is already active).
  if (opts.adiabatic && traj.ranges.adiabatic) {
    const r = traj.ranges.adiabatic;
    if (isFinite(r.lo)) {
      const margin = 0.5 * (hi - lo || 1);
      lo = Math.max(lo - margin, Math.min(lo, r.lo));
      hi = Math.min(hi + margin, Math.max(hi, r.hi));
    }
  }
  const pad = 0.1 * (hi - lo || 1);
  const ax = makeAxis({ xMin: tMin, xMax: tMax, yMin: lo - pad, yMax: hi + pad, w, h });
  drawFrame(ctx, ax);
  clipPlot(ctx, ax);
  if (opts.zeroLine)
    strokePath(ctx, ax, new Float64Array([tMin, tMax]), new Float64Array([0, 0]),
      { color: 'rgba(0,0,0,0.3)', width: 1, dash: [3, 3] });
  if (opts.diagLine) {
    // wave centerline x̃ = C t̃
    strokePath(ctx, ax, new Float64Array([tMin, tMax]),
      new Float64Array([params.C * tMin, params.C * tMax]),
      { color: 'rgba(43,108,176,0.4)', width: 1, dash: [4, 3] });
  }
  if (opts.adiabatic && traj.adBranches) drawAdiabatic(ctx, ax, traj.adBranches);
  strokePath(ctx, ax, ts, ys, { color: opts.color || '#2b6cb0', width: 1.5 });
  if (opts.extraChannel && params.showRef)
    strokePath(ctx, ax, ts, traj[opts.extraChannel], { color: '#2f7a3a', width: 1.2, dash: [3, 3] });
  const i = indexAt(currentTime);
  ctx.restore();
  ctx.save();
  ctx.strokeStyle = 'rgba(0,0,0,0.5)';
  ctx.setLineDash([2, 2]);
  ctx.beginPath();
  ctx.moveTo(ax.xToPx(currentTime), ax.padT);
  ctx.lineTo(ax.xToPx(currentTime), ax.padT + ax.plotH);
  ctx.stroke();
  ctx.restore();
  dot(ctx, ax, currentTime, ys[i], 4, '#b34700');
  if (opts.extraChannel && params.showRef)
    dot(ctx, ax, currentTime, traj[opts.extraChannel][i], 3, '#2f7a3a');
}

// ─── ⟨Δx̃⟩ vs {C, χ̃} sweeps ────────────────────────────────────────────
// Each sweep point averages N_SWEEP_REPS stochastic trajectories
// (Euler-Maruyama noise on top of RK4 drift). The adiabatic reference cell
// is deterministic and needs only one integration per point.
const N_SWEEP_REPS = 10;

// Single cell passage. RK4 drift + E-M noise. Returns final x̃.
// model: 'GL' | 'Hill'   v0, nH: Hill-only (ignored in GL mode)
function integrateSinglePass(chi, M, C, lam, mu, tht, seed, model, v0, nH) {
  const tMin = -X_INIT / Math.max(C, 1e-6), tMax = -tMin;
  const dt = Math.max(0.003, Math.min(0.02, (tMax - tMin) / 2000));
  const rng = makeRng(seed);
  const driftP = model === 'GL'
    ? (p, x, t) => {
        const u = x - C * t;
        const L = M * Math.exp(-0.5 * u * u);
        const m2 = p * p;
        return chi * (-L * u) + (L - 1) * p + lam * (m2 - m2 * m2) * p;
      }
    : (p, x, t) => {
        const u = x - C * t;
        const L = M * Math.exp(-0.5 * u * u);
        return chi * (-L * u) - p;
      };
  const vP = model === 'GL'
    ? p => mu * p
    : p => vOfP_Hill(p, v0, nH);
  let p = 0, x = 0, t = tMin;
  while (t < tMax) {
    const h = Math.min(dt, tMax - t);
    const k1p = driftP(p, x, t),                        k1x = vP(p);
    const k2p = driftP(p+h/2*k1p, x+h/2*k1x, t+h/2),   k2x = vP(p+h/2*k1p);
    const k3p = driftP(p+h/2*k2p, x+h/2*k2x, t+h/2),   k3x = vP(p+h/2*k2p);
    const k4p = driftP(p+h*k3p,   x+h*k3x,   t+h),      k4x = vP(p+h*k3p);
    // E-M noise evaluated at start-of-step position
    const u0 = x - C * t;
    const L0 = M * Math.exp(-0.5 * u0 * u0);
    const noise = Math.sqrt(2 * tht * Math.max(L0, 0) * h) * rng.gauss();
    p += h/6 * (k1p + 2*k2p + 2*k3p + k4p) + noise;
    x += h/6 * (k1x + 2*k2x + 2*k3x + k4x);
    t += h;
  }
  return x;
}

// Adiabatic reference cell (pure-gradient, no inertia). Deterministic, 1 run. Returns final x̃.
// GL:   P_eq = χ̃ ∂_x̃𝓛 / (1-𝓛), velocity = μ̃ P_eq
// Hill: P_eq = χ̃ ∂_x̃𝓛          (balance χ̃∂𝓛 − P = 0), velocity = vOfP_Hill(P_eq)
function integrateRefPass(chi, M, C, mu, model, v0, nH) {
  const tMin = -X_INIT / Math.max(C, 1e-6), tMax = -tMin;
  const dt = Math.max(0.003, Math.min(0.02, (tMax - tMin) / 2000));
  const dxRef = model === 'GL'
    ? (x, t) => {
        const u = x - C * t;
        const L = M * Math.exp(-0.5 * u * u);
        const denom = Math.max(Math.abs(L - 1), 0.05);
        return mu * chi * L * u / denom;
      }
    : (x, t) => {
        const u = x - C * t;
        const L = M * Math.exp(-0.5 * u * u);
        const Peq = chi * (-L * u);   // P_eq from balance χ̃∂𝓛 − P = 0
        return vOfP_Hill(Peq, v0, nH);
      };
  let xr = 0, t = tMin;
  while (t < tMax) {
    const h = Math.min(dt, tMax - t);
    const k1 = dxRef(xr, t);
    const k2 = dxRef(xr+h/2*k1, t+h/2);
    const k3 = dxRef(xr+h/2*k2, t+h/2);
    const k4 = dxRef(xr+h*k3,   t+h);
    xr += h/6 * (k1 + 2*k2 + 2*k3 + k4);
    t += h;
  }
  return xr;
}

let sweepCacheC = null, sweepCacheChi = null;
function runSweepC() {
  const N = 81;
  const Cs = new Float64Array(N), dxFull = new Float64Array(N), dxRef = new Float64Array(N);
  // Linear sweep over wave speed; Cmin > 0 to keep the trajectory window finite.
  const Cmin = 0.05, Cmax = 5;
  const { chi, M, lam, mu, tht, seed, model, v0, nHill } = params;
  for (let k = 0; k < N; k++) {
    const C = Cmin + (Cmax - Cmin) * k / (N - 1);
    Cs[k] = C;
    dxRef[k] = integrateRefPass(chi, M, C, mu, model, v0, nHill);
    let sum = 0;
    for (let rep = 0; rep < N_SWEEP_REPS; rep++)
      sum += integrateSinglePass(chi, M, C, lam, mu, tht, seed + rep * 7919 + k * 83, model, v0, nHill);
    dxFull[k] = sum / N_SWEEP_REPS;
  }
  sweepCacheC = { Cs, dxFull, dxRef };
}
function runSweepChi() {
  const N = 61;
  const chiMax = Math.max(3, params.chi * 1.5);
  const chis = new Float64Array(N), dxFull = new Float64Array(N), dxRef = new Float64Array(N);
  const { M, C, lam, mu, tht, seed, model, v0, nHill } = params;
  for (let k = 0; k < N; k++) {
    const chi = (chiMax * k) / (N - 1);
    chis[k] = chi;
    dxRef[k] = integrateRefPass(chi, M, C, mu, model, v0, nHill);
    let sum = 0;
    for (let rep = 0; rep < N_SWEEP_REPS; rep++)
      sum += integrateSinglePass(chi, M, C, lam, mu, tht, seed + rep * 7919 + k * 83, model, v0, nHill);
    dxFull[k] = sum / N_SWEEP_REPS;
  }
  sweepCacheChi = { chis, dxFull, dxRef };
}
function runSweep() { runSweepC(); runSweepChi(); }

// Generic sweep renderer. xs is the swept parameter array, ysFull/ysRef the
// signed Δx̃ outputs. cursorX = current parameter value to mark.
function drawGenericSweep(canvasId, xs, ysFull, ysRef, cursorX, opts = {}) {
  const cv = document.getElementById(canvasId);
  const ctx = autoFit(cv);
  const w = cv.clientWidth, h = cv.clientHeight;
  if (!xs) {
    const ax = makeAxis({ xMin: opts.xMin ?? 0, xMax: opts.xMax ?? 1,
                          yMin: opts.yMin ?? -1, yMax: opts.yMax ?? 1, w, h, logX: !!opts.logX });
    drawFrame(ctx, ax);
    return;
  }
  // y-range: opts.yMin/yMax override; otherwise adapt to GL-cell data only
  // (reference curve doesn't expand the range).
  let yMin, yMax;
  if (opts.yMin != null && opts.yMax != null) {
    yMin = opts.yMin; yMax = opts.yMax;
  } else {
    let lo = Infinity, hi = -Infinity;
    for (let i = 0; i < xs.length; i++) {
      if (ysFull[i] < lo) lo = ysFull[i]; if (ysFull[i] > hi) hi = ysFull[i];
    }
    if (!isFinite(lo)) { lo = -0.1; hi = 0.1; }
    const pad = 0.1 * (hi - lo || 0.1);
    yMin = lo - pad; yMax = hi + pad;
  }
  const ax = makeAxis({ xMin: opts.xMin ?? xs[0], xMax: opts.xMax ?? xs[xs.length - 1],
                        yMin, yMax, w, h, logX: !!opts.logX });
  drawFrame(ctx, ax);
  clipPlot(ctx, ax);
  strokePath(ctx, ax, new Float64Array([ax.xMin, ax.xMax]), new Float64Array([0, 0]),
    { color: 'rgba(0,0,0,0.3)', width: 1, dash: [3, 3] });
  strokePath(ctx, ax, xs, ysFull, { color: '#b34700', width: 2 });
  strokePath(ctx, ax, xs, ysRef,  { color: '#2f7a3a', width: 1.5, dash: [4, 3] });
  ctx.restore();
  ctx.save();
  ctx.strokeStyle = 'rgba(0,0,0,0.5)';
  ctx.setLineDash([2, 2]);
  ctx.beginPath();
  ctx.moveTo(ax.xToPx(cursorX), ax.padT);
  ctx.lineTo(ax.xToPx(cursorX), ax.padT + ax.plotH);
  ctx.stroke();
  ctx.restore();
}

function drawSweepC() {
  const c = sweepCacheC;
  drawGenericSweep('cv-sweep-c',
    c && c.Cs, c && c.dxFull, c && c.dxRef, params.C,
    { xMin: 0.05, xMax: 5, yMin: -1, yMax: 1 });
}
function drawSweepChi() {
  const c = sweepCacheChi;
  drawGenericSweep('cv-sweep-chi',
    c && c.chis, c && c.dxFull, c && c.dxRef, params.chi,
    { xMin: 0, xMax: c ? c.chis[c.chis.length - 1] : Math.max(3, params.chi * 1.5) });
}

function drawVcell() {
  const cv = document.getElementById('cv-vcell');
  const ctx = autoFit(cv);
  const w = cv.clientWidth, h = cv.clientHeight;
  const tMin = trajTMin(), tMax = trajTMax();
  if (!traj) {
    const ax = makeAxis({ xMin: tMin, xMax: tMax, yMin: -1.5, yMax: 1.5, w, h });
    drawFrame(ctx, ax);
    return;
  }
  // Use precomputed vcArr = (v_cell)/C and its range (both cached on traj at recompute time).
  // y-range adapts to the data; reference lines (y=0, y=1) may fall outside and be clipped.
  const vcArr = traj.vcArr;
  let { lo, hi } = traj.ranges.vcArr || { lo: -1, hi: 1 };
  if (!isFinite(lo)) { lo = -1; hi = 1; }
  const pad = 0.1 * (hi - lo || 0.1);
  const ax = makeAxis({ xMin: tMin, xMax: tMax, yMin: lo - pad, yMax: hi + pad, w, h });
  drawFrame(ctx, ax);
  clipPlot(ctx, ax);
  // y = 0 baseline
  strokePath(ctx, ax, new Float64Array([tMin, tMax]), new Float64Array([0, 0]),
    { color: 'rgba(0,0,0,0.3)', width: 1, dash: [3, 3] });
  // y = 1 reference (cell matches wave speed)
  strokePath(ctx, ax, new Float64Array([tMin, tMax]), new Float64Array([1, 1]),
    { color: 'rgba(43,108,176,0.5)', width: 1.2, dash: [5, 3] });
  strokePath(ctx, ax, traj.ts, vcArr, { color: '#b34700', width: 1.8 });
  const i = indexAt(currentTime);
  ctx.restore();
  ctx.save();
  ctx.strokeStyle = 'rgba(0,0,0,0.5)';
  ctx.setLineDash([2, 2]);
  ctx.beginPath();
  ctx.moveTo(ax.xToPx(currentTime), ax.padT);
  ctx.lineTo(ax.xToPx(currentTime), ax.padT + ax.plotH);
  ctx.stroke();
  ctx.restore();
  dot(ctx, ax, currentTime, vcArr[i], 4, '#b34700');
}

// ─── singular curve Σ = {(𝓛, P) : f(P, 𝓛) = 0} ──────────────────────────
// f is the autonomous internal response — the relaxational part of the P-SDE,
// i.e. the whole drift minus the chemotactic drive s = χ̃∂_x̃𝓛:
//     GL:   f = (𝓛 - 1)P + λ(P³ - P⁵) = P · [𝓛 - 𝓛_Σ(P)],  𝓛_Σ(P) = 1 - λP² + λP⁴
//     Hill: f = -P
// Because q = P + χ̃𝓛/C absorbs the drive, q̇ = f (for a cell that is slow
// against the wave). So Σ is exactly the q̇ = 0 nullcline: crossing it flips the
// sign of q̇ and cuts the cycle into q̇ > 0 and q̇ < 0 intervals. It is also where
// U = μ̃P/f diverges. The GL branch is explicit in P (no root finding needed),
// and sign(f) = sign(P)·sign(𝓛 - 𝓛_Σ(P)) makes the sign shading a per-row split.
const L_SIGMA = P => 1 - params.lam * P * P * (1 - P * P);
const TINT_POS = 'rgba(43,108,176,0.07)';   // q̇ > 0
const TINT_NEG = 'rgba(179,71,0,0.07)';     // q̇ < 0

function drawSingularLP(ctx, ax) {
  const isGL = params.model === 'GL';
  // Sign shading: for a fixed P (one pixel row) the row splits at 𝓛 = 𝓛_Σ(P),
  // with sign(f) = sign(P) to the right of the split and -sign(P) to the left.
  ctx.save();
  const y0 = ax.padT, y1 = ax.padT + ax.plotH, x0 = ax.padL, x1 = ax.padL + ax.plotW;
  for (let py = y0; py < y1; py++) {
    const P = ax.yMin + (y1 - py - 0.5) / ax.plotH * (ax.yMax - ax.yMin);
    if (P === 0) continue;
    // Split abscissa in px; Hill has no 𝓛-dependence, so the whole row is one sign.
    const xs = isGL ? Math.max(x0, Math.min(x1, ax.xToPx(L_SIGMA(P)))) : x0;
    // GL: 𝓛 > 𝓛_Σ ⇒ sign(f) = sign(P). Hill: f = -P, no 𝓛-dependence, sign flipped.
    const fPos = isGL ? (P > 0) : (P < 0);
    const right = fPos ? TINT_POS : TINT_NEG;
    const left  = fPos ? TINT_NEG : TINT_POS;
    if (isGL && xs > x0) { ctx.fillStyle = left;  ctx.fillRect(x0, py, xs - x0, 1); }
    if (xs < x1)         { ctx.fillStyle = right; ctx.fillRect(xs, py, x1 - xs, 1); }
  }
  ctx.restore();
  // Σ, branch 1: P = 0 (present in both models).
  strokePath(ctx, ax, new Float64Array([ax.xMin, ax.xMax]), new Float64Array([0, 0]),
    { color: '#4a5568', width: 1.6 });
  if (!isGL) return;
  // Σ, branch 2: 𝓛 = 𝓛_Σ(P), parametrized by P over the visible range.
  const NS = 200;
  const Lx = new Float64Array(NS), Py = new Float64Array(NS);
  for (let k = 0; k < NS; k++) {
    const P = ax.yMin + (k / (NS - 1)) * (ax.yMax - ax.yMin);
    Py[k] = P; Lx[k] = L_SIGMA(P);
  }
  strokePath(ctx, ax, Lx, Py, { color: '#4a5568', width: 1.6 });
}

// Generic phase-trajectory plot: the parametric curve (xChan(t̃), yChan(t̃)).
// The full loop is drawn faint; the portion already traversed (t̃ ≤ currentTime)
// is solid, with the current state as a bead. For P–s space the enclosed area is
// the rectification signature; the 𝓛–P view shows the response hysteresis of the
// cell's polarization against the cue it currently sees.
// opts.includeOriginX/Y force the corresponding zero into the visible range.
function drawPhase(canvasId, xChan, yChan, opts = {}) {
  const cv = document.getElementById(canvasId);
  const ctx = autoFit(cv);
  const w = cv.clientWidth, h = cv.clientHeight;
  if (!traj) {
    const ax = makeAxis({ xMin: -1, xMax: 1, yMin: -1, yMax: 1, w, h });
    drawFrame(ctx, ax);
    return;
  }
  const xa = traj[xChan], ya = traj[yChan];
  let xr = traj.ranges[xChan], yr = traj.ranges[yChan];
  let xLo = xr.lo, xHi = xr.hi, yLo = yr.lo, yHi = yr.hi;
  if (!isFinite(xLo)) { xLo = -1; xHi = 1; }
  if (!isFinite(yLo)) { yLo = -1; yHi = 1; }
  if (opts.includeOriginX !== false) { xLo = Math.min(xLo, 0); xHi = Math.max(xHi, 0); }
  if (opts.includeOriginY !== false) { yLo = Math.min(yLo, 0); yHi = Math.max(yHi, 0); }
  const xPad = 0.1 * (xHi - xLo || 1), yPad = 0.1 * (yHi - yLo || 1);
  const ax = makeAxis({ xMin: xLo - xPad, xMax: xHi + xPad, yMin: yLo - yPad, yMax: yHi + yPad, w, h });
  drawFrame(ctx, ax);
  clipPlot(ctx, ax);
  // Σ and the q̇-sign shading go underneath the trajectory. Σ's P = 0 branch
  // replaces the horizontal zero reference line, so skip that one here.
  if (opts.singular) drawSingularLP(ctx, ax);
  // Zero reference lines.
  if (!opts.singular)
    strokePath(ctx, ax, new Float64Array([ax.xMin, ax.xMax]), new Float64Array([0, 0]),
      { color: 'rgba(0,0,0,0.3)', width: 1, dash: [3, 3] });
  strokePath(ctx, ax, new Float64Array([0, 0]), new Float64Array([ax.yMin, ax.yMax]),
    { color: 'rgba(0,0,0,0.3)', width: 1, dash: [3, 3] });
  // Full loop, faint.
  strokePath(ctx, ax, xa, ya, { color: 'rgba(43,108,176,0.35)', width: 1.2 });
  // Traversed portion (up to currentTime), solid.
  const i = indexAt(currentTime);
  if (i > 0)
    strokePath(ctx, ax, xa.subarray(0, i + 1), ya.subarray(0, i + 1),
      { color: '#b34700', width: 1.8 });
  ctx.restore();
  dot(ctx, ax, xa[i], ya[i], 5, '#b34700');
}

// ─── controls ───────────────────────────────────────────────────────────
const controlsEl = document.getElementById('controls');
const kpis = makeKpis([
  { id: 'x',  label: 'x̃ (cell)' },
  { id: 'xR', label: 'x̃ (ref)' },
  { id: 'p',  label: 'P' },
  { id: 't',  label: 't̃' },
]);

const sM     = makeSlider({ id: 'M',     symbol: 'M',                   bind: [params, 'M'],     min: 0,    max: 5,   step: 0.01, fmt: v => v.toFixed(2) });
const sC     = makeSlider({ id: 'C',     symbol: 'C',                   bind: [params, 'C'],     min: 0,    max: 10,  step: 0.01, fmt: v => v.toFixed(2) });
const sChi   = makeSlider({ id: 'chi',   symbol: '\\tilde{\\chi}',      bind: [params, 'chi'],   min: 0,    max: 5,   step: 0.01, fmt: v => v.toFixed(2) });
const sMu    = makeSlider({ id: 'mu',    symbol: '\\tilde{\\mu}',       bind: [params, 'mu'],    min: 0,    max: 3,   step: 0.01, fmt: v => v.toFixed(2) });
const sLam   = makeSlider({ id: 'lam',   symbol: '\\lambda',            bind: [params, 'lam'],   min: 0.01, max: 10,  log: true, fmt: v => v.toPrecision(3) });
const sTht   = makeSlider({ id: 'tht',   symbol: '\\vartheta',          bind: [params, 'tht'],   min: 1e-4, max: 1,   log: true, fmt: v => v.toExponential(2) });
// Hill-only sliders
const sV0    = makeSlider({ id: 'v0',    symbol: '\\tilde{v}_{0}',      bind: [params, 'v0'],    min: 0,    max: 3,   step: 0.01, fmt: v => v.toFixed(2) });
const sNHill = makeSlider({ id: 'nHill', symbol: 'n',                   bind: [params, 'nHill'], min: 1,    max: 10,  step: 1,    transform: Math.round, fmt: v => v.toFixed(0) });

sM.onChange(()    => { recalibrate(); markDirty(); });
sC.onChange(()    => { recalibrate(); markDirty(); });
sChi.onChange(()  => { recalibrate(); markDirty(); });
sMu.onChange(()   => { recalibrate(); markDirty(); });
sLam.onChange(()  => markDirty());
// ϑ: trajectory updates in real time; sweep recomputes only on release (averaging is expensive).
sTht.onChange(()  => { trajDirty = true; });
sTht.onRelease(() => { sweepDirty = true; });
sV0.onChange(()   => markDirty());
sNHill.onChange(() => markDirty());

const linkedReadout = () =>
  `→ M=${params.M.toFixed(2)}, C=${params.C.toPrecision(3)}, ` +
  `χ̃=${params.chi.toFixed(2)}, μ̃=${params.mu.toFixed(2)}`;
function pushAllNondimSliders() {
  applyingDim = true;
  sM.set(params.M); sC.set(params.C); sChi.set(params.chi); sMu.set(params.mu);
  applyingDim = false;
}

const sLmax  = makeSlider({ id: 'Lmax',  symbol: 'L_{\\max}', bind: [dim, 'Lmax'],  min: 0.1, max: 4, step: 0.01, fmt: v => v.toFixed(2), linkedLabel: () => linkedReadout() });
const sSigma = makeSlider({ id: 'sigma', symbol: '\\sigma',   bind: [dim, 'sigma'], min: 0.1, max: 4, step: 0.01, fmt: v => v.toFixed(2), linkedLabel: () => linkedReadout() });
const sCwave = makeSlider({ id: 'c',     symbol: 'c',         bind: [dim, 'c'],     min: 0.1, max: 4, step: 0.01, fmt: v => v.toFixed(2), linkedLabel: () => linkedReadout() });
function refreshDimReadouts() {
  sLmax.setLinkedText(linkedReadout()); sSigma.setLinkedText(linkedReadout()); sCwave.setLinkedText(linkedReadout());
}
sLmax.onChange(()  => { recomputeFromDim(); pushAllNondimSliders(); refreshDimReadouts(); markDirty(); });
sSigma.onChange(() => { recomputeFromDim(); pushAllNondimSliders(); refreshDimReadouts(); markDirty(); });
sCwave.onChange(() => { recomputeFromDim(); pushAllNondimSliders(); refreshDimReadouts(); markDirty(); });

const sDt    = makeSlider({ id: 'dt',    symbol: 'd\\tilde{t}',       bind: [params, 'dt'],    min: 1e-4, max: 0.05, log: true, fmt: v => v.toExponential(2) });
const sSpeed = makeSlider({ id: 'speed', symbol: '\\text{play speed}', bind: [params, 'speed'], min: 0.01, max: 100,  log: true, fmt: v => `${v.toPrecision(2)}×` });
const sSeed  = makeSlider({ id: 'seed',  symbol: '\\text{seed}',       bind: [params, 'seed'],  min: 1,    max: 9999, step: 1,   transform: Math.round, fmt: v => v.toFixed(0) });
sDt.onChange(()    => markDirty());
sSpeed.onChange(() => {});
sSeed.onChange(()  => markDirty());

// THE time slider — the primary way the user navigates the trajectory.
const sTime = makeSlider({
  id: 'time', symbol: '\\tilde{t}', value: -X_INIT,
  min: -X_INIT, max: +X_INIT, step: 0.01,
  fmt: v => v.toFixed(2),
});
sTime.onChange(v => { currentTime = v; });

const refToggle = makeToggle({
  label: 'reference cell (pure-gradient, no inertia)',
  options: [{ id: 'on', label: 'show' }, { id: 'off', label: 'hide' }],
  value: 'on',
  onChange: v => { params.showRef = (v === 'on'); markDirty(); },
});

// ─── model toggle ────────────────────────────────────────────────────────
const MODEL_STORAGE_KEY = 'gl-motility:setup2:model';
function loadModelChoice() {
  try { return localStorage.getItem(MODEL_STORAGE_KEY) || 'GL'; } catch { return 'GL'; }
}
function saveModelChoice(v) {
  try { localStorage.setItem(MODEL_STORAGE_KEY, v); } catch {}
}
// Apply model choice BEFORE creating sliders so localStorage-restored model is used.
params.model = loadModelChoice();

// Per-model defaults applied to all sliders (shared and model-specific) on
// model toggle. Shared sliders (M, C, χ̃, ϑ) snap to the new model's preferred
// regime; model-specific sliders are restored to their own defaults.
const GL_DEFAULTS   = { M: 2.0, C: 1.5, chi: 0.75, mu: 0.05, lam: 1.0, tht: 1e-4 };
const HILL_DEFAULTS = { M: 2.0, C: 5.0, chi: 5.0,  v0: 0.5,  nHill: 4, tht: 1e-4 };

function applyModelDefaults(modelName) {
  const d = modelName === 'GL' ? GL_DEFAULTS : HILL_DEFAULTS;
  if ('M'     in d) sM.set(d.M);
  if ('C'     in d) sC.set(d.C);
  if ('chi'   in d) sChi.set(d.chi);
  if ('tht'   in d) sTht.set(d.tht);
  if ('mu'    in d) sMu.set(d.mu);
  if ('lam'   in d) sLam.set(d.lam);
  if ('v0'    in d) sV0.set(d.v0);
  if ('nHill' in d) sNHill.set(d.nHill);
}

const modelToggle = makeToggle({
  label: 'model',
  options: [{ id: 'GL', label: 'Ginzburg–Landau' }, { id: 'Hill', label: 'Hill / linear-relaxation' }],
  value: params.model,
  onChange: v => {
    params.model = v;
    saveModelChoice(v);
    applyModelDefaults(v);
    applyModelVisibility();
    markDirty();
  },
});

const buttons = makeButtonRow([
  { label: '▶  play', onClick() {
      playing = !playing;
      if (playing && currentTime >= trajTMax() - 1e-6) {
        currentTime = trajTMin(); sTime.set(currentTime);
      }
      buttons.refs['▶  play'].textContent = playing ? '⏸  pause' : '▶  play';
    } },
  { label: '⟳  reset', ghost: true, onClick: () => {
      currentTime = trajTMin(); sTime.set(currentTime);
    } },
]);

// ─── section layout ──────────────────────────────────────────────────────
// Sections for GL-only and Hill-only parameters are always in the DOM;
// applyModelVisibility() hides the inactive sections.
const secGLVelocity   = section('GL velocity (nondim)',   [sMu.el]);
const secHillVelocity = section('Hill velocity (nondim)', [sV0.el, sNHill.el]);
const secGLWell       = section('GL well (nondim)',       [sLam.el]);
const secNoise        = section('noise (nondim)',         [sTht.el]);

controlsEl.appendChild(section('model', [modelToggle.el]));
controlsEl.appendChild(section('time scrub', [sTime.el, buttons.el]));
controlsEl.appendChild(section('wave (nondim)', [sM.el, sC.el]));
controlsEl.appendChild(section('coupling (nondim)', [sChi.el]));
controlsEl.appendChild(secGLVelocity);
controlsEl.appendChild(secHillVelocity);
controlsEl.appendChild(secGLWell);
controlsEl.appendChild(secNoise);
controlsEl.appendChild(detailsSection('dim sliders (push linked nondim)', [sLmax.el, sSigma.el, sCwave.el]));
controlsEl.appendChild(section('numerics & playback', [sDt.el, sSpeed.el, sSeed.el, refToggle.el]));
controlsEl.appendChild(kpis.el);
controlsEl.appendChild(el('div', { class: 'note' }, [
  'Trajectory window adapts: t̃ ∈ [-X/C, +X/C] with X = ' + X_INIT + '. ',
  'Wave reaches the cell at t̃ = 0 in every case. ',
  'Trajectory is recomputed deterministically (for the given seed) on any sim-parameter change. ',
  'F̃ panel (GL only): solid = current well at 𝓛(x̃_cell, t̃); faint dashed = well at peak cue 𝓛 = M.',
]));

// ─── model visibility ────────────────────────────────────────────────────
// Called after decorateAll() and on every model toggle.
function applyModelVisibility() {
  const isGL = params.model === 'GL';
  secGLVelocity.style.display   = isGL ? '' : 'none';
  secHillVelocity.style.display = isGL ? 'none' : '';
  secGLWell.style.display       = isGL ? '' : 'none';
  // F̃ panel: find the .plot-wrap wrapping cv-F (created lazily by decoratePlot).
  const cvF = document.getElementById('cv-F');
  if (cvF) {
    const wrap = cvF.closest('.plot-wrap') || cvF.parentNode;
    if (wrap) wrap.style.display = isGL ? '' : 'none';
  }
}

// All slider-driven params/dim fields are kept in sync by makeSlider's `bind:`
// (write-through on construction + every change). Initialize calibration to
// current dim/nondim state so subsequent dim slider moves scale relative to
// the loaded state, not the hardcoded one.
recalibrate();
// Adjust time slider min/max for the loaded C.
sTime.setMinMax(trajTMin(), trajTMax());
currentTime = Math.min(trajTMax(), Math.max(trajTMin(), sTime.value));
sTime.set(currentTime);

function decorateAll() {
  decoratePlot('cv-wave',  { titleTex: '\\text{wave } \\mathcal{L}(\\tilde x, \\tilde t)',
                             xLabelTex: '\\tilde x', yLabelTex: '\\mathcal{L}' });
  decoratePlot('cv-F',     { titleTex: '\\tilde F_{\\text{eff}} = \\tilde F - \\tilde\\chi\\,\\partial_{\\tilde x}\\mathcal{L}\\cdot P\\;\\text{(dashed: }\\mathcal{L}=M\\text{)}',
                             xLabelTex: 'P', yLabelTex: '\\tilde F_{\\text{eff}}' });
  decoratePlot('cv-grad',  { titleTex: '\\text{gradient at cell}',
                             xLabelTex: '\\tilde t', yLabelTex: '\\partial_{\\tilde x} \\mathcal{L}' });
  decoratePlot('cv-x',     { titleTex: '\\text{position (orange = cell, green dashed = ref)}',
                             xLabelTex: '\\tilde t', yLabelTex: '\\tilde x' });
  decoratePlot('cv-p',     { titleTex: '\\text{polarization at cell}\\;\\text{(gray: }\\mathrm{d}\\tilde F_{\\text{eff}}/\\mathrm{d}P=0\\text{, solid = stable, dashed = unstable)}',
                             xLabelTex: '\\tilde t', yLabelTex: 'P' });
  decoratePlot('cv-sweep-c',   { titleTex: '\\langle\\Delta\\tilde x\\rangle \\text{ vs wave speed}',
                                 xLabelTex: 'C', yLabelTex: '\\langle\\Delta\\tilde x\\rangle' });
  decoratePlot('cv-sweep-chi', { titleTex: '\\langle\\Delta\\tilde x\\rangle \\text{ vs chemotactic coupling}',
                                 xLabelTex: '\\tilde\\chi', yLabelTex: '\\langle\\Delta\\tilde x\\rangle' });
  decoratePlot('cv-vcell',     { titleTex: '\\text{cell velocity / wave speed (blue dashed = 1)}',
                                 xLabelTex: '\\tilde t', yLabelTex: '\\tilde v_{\\text{cell}} / C' });
  decoratePlot('cv-phase',     { titleTex: '\\text{phase trajectory in } P\\text{--}s',
                                 xLabelTex: 'P', yLabelTex: 's = \\tilde\\chi\\,\\partial_{\\tilde x}\\mathcal{L}' });
  decoratePlot('cv-phase-L',   { titleTex: '\\mathcal{L}\\text{--}P\\text{ hysteresis; gray: }\\Sigma\\;(\\dot q = 0)',
                                 xLabelTex: '\\mathcal{L}', yLabelTex: 'P' });
  decoratePlot('cv-phase-U',   { titleTex: 'q\\text{--}U\\text{ phase trajectory }(U = \\tilde\\mu P / f)',
                                 xLabelTex: 'q = P + \\tilde\\chi\\mathcal{L}/C',
                                 yLabelTex: 'U' });
  attachSvgExports({
    'cv-wave':      () => drawWave(),
    'cv-F':         () => drawF(),
    'cv-grad':      () => drawTrajChannel('cv-grad', 'gs', { color: '#2b6cb0', zeroLine: true }),
    'cv-x':         () => drawTrajChannel('cv-x',    'xs', { color: '#b34700', diagLine: true, extraChannel: 'xR' }),
    'cv-p':         () => drawTrajChannel('cv-p',    'ps', { color: '#2b6cb0', adiabatic: true }),
    'cv-sweep-c':   () => drawSweepC(),
    'cv-sweep-chi': () => drawSweepChi(),
    'cv-vcell':     () => drawVcell(),
    'cv-phase':     () => drawPhase('cv-phase',   'ps', 'ss'),
    'cv-phase-L':   () => drawPhase('cv-phase-L', 'Ls', 'ps', { singular: true }),
    'cv-phase-U':   () => drawPhase('cv-phase-U', 'qs', 'Us', { includeOriginY: false }),
  }, { prefix: 'setup2' });
  // Apply model visibility after plot-wraps exist.
  applyModelVisibility();
}
if (window.katex) decorateAll();
else window.addEventListener('load', decorateAll);

// Apply visibility immediately (before decorateAll; F wrap may not exist yet,
// but that's fine — applyModelVisibility re-runs inside decorateAll).
applyModelVisibility();

function frame(nowMs) {
  if (trajDirty) recomputeTrajectory();
  if (sweepDirty) { runSweep(); sweepDirty = false; }
  if (playing) {
    const dtReal = lastFrameMs ? (nowMs - lastFrameMs) / 1000 : 0;
    currentTime += params.speed * dtReal;
    if (currentTime >= trajTMax()) {
      currentTime = trajTMax();
      playing = false;
      buttons.refs['▶  play'].textContent = '▶  play';
    }
    sTime.set(currentTime);
  }
  lastFrameMs = nowMs;

  drawWave();
  drawF();
  drawTrajChannel('cv-grad', 'gs', { color: '#2b6cb0', zeroLine: true });
  drawTrajChannel('cv-x',    'xs', { color: '#b34700', diagLine: true, extraChannel: 'xR' });
  drawTrajChannel('cv-p',    'ps', { color: '#2b6cb0', adiabatic: true });
  drawSweepC();
  drawSweepChi();
  drawVcell();
  drawPhase('cv-phase',   'ps', 'ss');
  drawPhase('cv-phase-L', 'Ls', 'ps', { singular: true });
  drawPhase('cv-phase-U', 'qs', 'Us', { includeOriginY: false });

  if (traj) {
    const i = indexAt(currentTime);
    kpis.set('x',  traj.xs[i].toFixed(3));
    kpis.set('xR', traj.xR[i].toFixed(3));
    kpis.set('p',  traj.ps[i].toFixed(3));
    kpis.set('t',  currentTime.toFixed(2));
  }
  requestAnimationFrame(frame);
}
requestAnimationFrame(frame);
