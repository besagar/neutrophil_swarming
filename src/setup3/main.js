// Setup 3 — 2D radial swarm with central trap.
// Physics: docs/physics/setup3_swarm.md (identical nondim groups to Setup 2).
// Cue: 𝓛(r̃, t̃) = M · exp(-(r̃ - C t̃)² / 2),  ∇𝓛 = -𝓛 · (r̃ - C t̃) · r̂
// GL SDE (Euler–Maruyama):
//   dP_xi = [χ̃ dxL + (𝓛-1) P_xi + λ(|P|²-|P|⁴) P_xi] dt + √(2 ϑ 𝓛 dt) η_xi
//   dP_yi = [χ̃ dyL + (𝓛-1) P_yi + λ(|P|²-|P|⁴) P_yi] dt + √(2 ϑ 𝓛 dt) η_yi
//   dx_i = μ̃ P_xi dt,  dy_i = μ̃ P_yi dt  (free cells; trapped: frozen)
// Hill SDE (Euler–Maruyama; docs/physics/hill_model.md):
//   dP_xi = [χ̃ dxL − P_xi] dt + √(2 ϑ 𝓛 dt) η_xi   (linear relaxation, no well)
//   dP_yi = [χ̃ dyL − P_yi] dt + √(2 ϑ 𝓛 dt) η_yi
//   dx_i = ṽ₀ |P|ⁿ/(1+|P|ⁿ) P̂_x dt,  dy_i = ṽ₀ |P|ⁿ/(1+|P|ⁿ) P̂_y dt
// Boundary: r > R_dish → reflect radially; r < R_trap → trapped, freeze position.

import { makeRng } from '../shared/rng.js';
import { el, makeSlider, makeToggle, makeButtonRow, makeKpis, section, detailsSection, decoratePlot } from '../shared/dom.js';
import { autoFit, makeAxis, drawFrame, strokePath, dot, clipPlot } from '../shared/canvas.js';

// ─── nondim parameters ──────────────────────────────────────────────────────
const params = {
  model: 'GL',               // 'GL' | 'Hill'
  N: 1000,
  R_dish: 10, R_trap: 1,
  M: 2.0, C: 1.5,
  chi: 0.75, mu: 0.05,
  lam: 1.0, tht: 1e-4,
  v0: 0.5, nHill: 4,         // Hill-only: peak speed ṽ₀, Hill exponent n
  nWaves: 1, dtWave: 5.0,    // evenly spaced launches at t̃ = 0, Δt, 2Δt, …
  dt: 0.01, speed: 1.0, seed: 7,
};

// ─── model localStorage persistence ──────────────────────────────────────────
const MODEL_KEY = 'gl-motility:setup3:model';
function saveModelChoice() {
  try { localStorage.setItem(MODEL_KEY, params.model); } catch {}
}
function loadModelChoice() {
  try { return localStorage.getItem(MODEL_KEY) || 'GL'; } catch { return 'GL'; }
}
params.model = loadModelChoice();

// Trajectory window: last wave launches at (nWaves-1)·Δt; needs (R_dish+8)/C
// after that to fully leave the dish.
function trajTEnd() {
  const lastLaunch = Math.max(0, params.nWaves - 1) * params.dtWave;
  return lastLaunch + (params.R_dish + 8) / Math.max(params.C, 1e-6);
}

// ─── dim ↔ nondim linkage (verbatim from Setup 2) ───────────────────────────
// L_max → M = (Lmax/cal.Lmax) · cal.M
// c     → C = (c/cal.c) / (σ/cal.σ) · cal.C
// σ     → C, χ̃, μ̃ all scale by 1/(σ/cal.σ)
const dim = {
  Lmax: 1.0, sigma: 1.0, c: 1.0,
  cal: { Lmax: 1.0, sigma: 1.0, c: 1.0,
         M: 2.0, C: 1.5, chi: 0.75, mu: 0.05 },
};
function recomputeFromDim() {
  const sLmax  = dim.Lmax  / dim.cal.Lmax;
  const sC     = dim.c     / dim.cal.c;
  const sSigma = dim.sigma / dim.cal.sigma;
  params.M   = dim.cal.M   * sLmax;
  params.C   = dim.cal.C   * sC / sSigma;
  params.chi = dim.cal.chi / sSigma;
  params.mu  = dim.cal.mu  / sSigma;
}
function recalibrate() {
  if (applyingDim) return;
  dim.cal = {
    Lmax: dim.Lmax, sigma: dim.sigma, c: dim.c,
    M: params.M, C: params.C, chi: params.chi, mu: params.mu,
  };
}
let applyingDim = false;

// ─── trajectory storage ─────────────────────────────────────────────────────
// xs[k*N+i], ys[k*N+i]: positions; trapped[k*N+i]: 0|1.
// Aggregates per step: rmean[k], nTrap[k], pmean[k].
let traj = null;
let trajDirty = true;
function markDirty() { trajDirty = true; }

function recomputeTrajectory() {
  const { N, R_dish, R_trap, M, C, chi, mu, lam, tht, v0, nHill, dt, nWaves, dtWave } = params;
  const model_GL = params.model === 'GL';
  // Precompute wave launch times.
  const waveTs = new Float64Array(Math.max(1, nWaves));
  for (let j = 0; j < waveTs.length; j++) waveTs[j] = j * dtWave;
  const T_end = trajTEnd();
  const N_steps = Math.ceil(T_end / dt) + 1;

  if (N * N_steps > 5e7) console.warn(`Setup 3: allocating ${N * N_steps} cells×steps — consider reducing N or dt.`);

  const xs      = new Float32Array(N_steps * N);
  const ys      = new Float32Array(N_steps * N);
  const trapped = new Uint8Array(N_steps * N);
  const rmean   = new Float64Array(N_steps);
  const nTrap   = new Float64Array(N_steps);
  const pmean   = new Float64Array(N_steps);

  const rng = makeRng(params.seed);

  // Initialize cells uniformly on annulus R_trap < r < R_dish (uniform in area).
  const px = new Float32Array(N);
  const py = new Float32Array(N);
  const Px = new Float32Array(N);
  const Py = new Float32Array(N);
  const tr = new Uint8Array(N);   // current trapped flags

  const Rd2 = R_dish * R_dish, Rt2 = R_trap * R_trap;
  for (let i = 0; i < N; i++) {
    const r = Math.sqrt(Rt2 + rng.uniform() * (Rd2 - Rt2));
    const th = 2 * Math.PI * rng.uniform();
    px[i] = r * Math.cos(th);
    py[i] = r * Math.sin(th);
    // Px[i] = Py[i] = 0 (default Float32Array)
    // tr[i] = 0 (default Uint8Array)
  }

  // Store step 0.
  xs.set(px, 0); ys.set(py, 0);
  // trapped step 0 all zero (default)
  { let rsum = 0, psum = 0, nfree = 0;
    for (let i = 0; i < N; i++) {
      const r = Math.hypot(px[i], py[i]);
      rsum += r; psum += Math.hypot(Px[i], Py[i]); nfree++;
    }
    rmean[0] = nfree ? rsum / nfree : 0;
    nTrap[0] = 0;
    pmean[0] = nfree ? psum / nfree : 0;
  }

  // Integration loop (Euler–Maruyama; dt → nondim-time step = dt̃).
  for (let k = 1; k < N_steps; k++) {
    const t = (k - 1) * dt;   // pre-step time (step k advances from t to t+dt)
    const sqNoise = Math.sqrt(2 * tht * dt);  // noise prefactor (𝓛 multiplied inline)

    let rsum = 0, psum = 0, nfree = 0, ntrap = 0;
    for (let i = 0; i < N; i++) {
      const xi = px[i], yi = py[i];
      const r = Math.hypot(xi, yi);

      // Multi-wave cue: 𝓛 = Σ_j M·exp(-(r - C·(t - t_j))²/2) for t ≥ t_j.
      // ∂_r 𝓛 = Σ_j -𝓛_j·(r - C·(t - t_j)). Inline radial-to-cartesian via r̂.
      let L = 0, gradR = 0;
      for (let j = 0; j < nWaves; j++) {
        const tj = waveTs[j];
        if (t < tj) continue;
        const uu = r - C * (t - tj);
        const Lj = M * Math.exp(-0.5 * uu * uu);
        L += Lj;
        gradR += -Lj * uu;
      }
      const ginv = r > 1e-12 ? gradR / r : 0;
      const dxL = ginv * xi, dyL = ginv * yi;

      // Polarization drift coefficient: GL uses bistable well, Hill uses −1 (linear relaxation).
      const m2 = Px[i] * Px[i] + Py[i] * Py[i];
      const m4 = m2 * m2;
      const driftCoef = model_GL ? ((L - 1) + lam * (m2 - m4)) : (-1);
      const sqL = Math.sqrt(Math.max(L, 0));

      if (tr[i]) {
        // Trapped: freeze position, continue polarization.
        Px[i] += (chi * dxL + driftCoef * Px[i]) * dt + sqNoise * sqL * rng.gauss();
        Py[i] += (chi * dyL + driftCoef * Py[i]) * dt + sqNoise * sqL * rng.gauss();
        ntrap++;
      } else {
        // Free cell.
        Px[i] += (chi * dxL + driftCoef * Px[i]) * dt + sqNoise * sqL * rng.gauss();
        Py[i] += (chi * dyL + driftCoef * Py[i]) * dt + sqNoise * sqL * rng.gauss();

        // Position update: GL = linear velocity μ̃P; Hill = Hill-saturating velocity.
        let nx, ny;
        if (model_GL) {
          nx = xi + mu * Px[i] * dt;
          ny = yi + mu * Py[i] * dt;
        } else {
          const m = Math.hypot(Px[i], Py[i]);
          if (m < 1e-12) { nx = xi; ny = yi; }
          else {
            const mn = Math.pow(m, nHill);
            const factor = v0 * mn / (1 + mn);   // saturates to v0 as |P| → ∞
            nx = xi + factor * Px[i] / m * dt;
            ny = yi + factor * Py[i] / m * dt;
          }
        }

        // Outer boundary: reflect radially if r > R_dish.
        const rn = Math.hypot(nx, ny);
        if (rn > R_dish) {
          const scale = (2 * R_dish - rn) / rn;
          // If still outside after reflection, clamp to R_dish.
          const rr = 2 * R_dish - rn;
          if (rr < 0) { nx = (nx / rn) * R_dish; ny = (ny / rn) * R_dish; }
          else { nx *= scale; ny *= scale; }
        }

        // Trap boundary: if r < R_trap, mark trapped and freeze.
        const rn2 = Math.hypot(nx, ny);
        if (rn2 < R_trap) {
          tr[i] = 1;
          // Freeze at entry point (keep position where it crossed).
          px[i] = nx; py[i] = ny;
          ntrap++;
        } else {
          px[i] = nx; py[i] = ny;
          const rmag = Math.hypot(Px[i], Py[i]);
          rsum += rn2; psum += rmag; nfree++;
        }
      }
    }

    // Store.
    const base = k * N;
    for (let i = 0; i < N; i++) {
      xs[base + i] = px[i];
      ys[base + i] = py[i];
      trapped[base + i] = tr[i];
    }
    rmean[k] = nfree ? rsum / nfree : 0;
    nTrap[k] = ntrap;
    pmean[k] = nfree ? psum / nfree : 0;
  }

  // Precompute aggregate ranges for y-axis scaling.
  function arrRange(a) {
    let lo = Infinity, hi = -Infinity;
    for (let i = 0; i < a.length; i++) { if (a[i] < lo) lo = a[i]; if (a[i] > hi) hi = a[i]; }
    return { lo, hi };
  }
  const ranges = {
    rmean:   arrRange(rmean),
    nTrap:   arrRange(nTrap),
    pmean:   arrRange(pmean),
  };

  traj = { xs, ys, trapped, rmean, nTrap, pmean, ranges, N_steps, N };
  trajDirty = false;

  // Sync time slider range.
  sTime.setMinMax(0, T_end);
  if (currentTime > T_end) currentTime = T_end;
  if (currentTime < 0) currentTime = 0;
  sTime.set(currentTime);
}

// ─── scrub time ─────────────────────────────────────────────────────────────
let currentTime = 0;
let playing = false;
let lastFrameMs = 0;

function indexAt(t) {
  if (!traj) return 0;
  const T_end = trajTEnd();
  const i = Math.round(t / (T_end / (traj.N_steps - 1)));
  return Math.max(0, Math.min(traj.N_steps - 1, i));
}

// ─── drawing ─────────────────────────────────────────────────────────────────

function drawDish() {
  const cv = document.getElementById('cv-dish');
  const ctx = autoFit(cv);
  const w = cv.clientWidth, h = cv.clientHeight;
  const R = params.R_dish;
  const ax = makeAxis({
    xMin: -R * 1.05, xMax: R * 1.05,
    yMin: -R * 1.05, yMax: R * 1.05,
    w, h, aspect: 1,
  });
  drawFrame(ctx, ax, { showGridX: false, showGridY: false, showTickLabelsX: false, showTickLabelsY: false });
  clipPlot(ctx, ax);

  // Faint trap disk fill.
  const txc = ax.xToPx(0), tyc = ax.yToPx(0);
  const trapPxR = ax.xToPx(params.R_trap) - ax.xToPx(0);
  ctx.beginPath();
  ctx.arc(txc, tyc, Math.abs(trapPxR), 0, 2 * Math.PI);
  ctx.fillStyle = 'rgba(43,108,176,0.10)';
  ctx.fill();
  ctx.strokeStyle = 'rgba(43,108,176,0.6)';
  ctx.lineWidth = 1.2;
  ctx.stroke();

  // Outer dish boundary.
  const dishPxR = ax.xToPx(R) - ax.xToPx(0);
  ctx.beginPath();
  ctx.arc(txc, tyc, Math.abs(dishPxR), 0, 2 * Math.PI);
  ctx.strokeStyle = 'rgba(100,100,100,0.5)';
  ctx.lineWidth = 1.5;
  ctx.stroke();

  // Wavefront rings: one per launched wave at r = C·(t - t_j). ±1σ faint
  // companions suggest Gaussian width (σ = 1 in nondim).
  for (let j = 0; j < params.nWaves; j++) {
    const tj = j * params.dtWave;
    if (currentTime < tj) continue;
    const wfR = params.C * (currentTime - tj);
    if (wfR <= 0 || wfR > R * 1.1) continue;
    for (const dr of [-1, 1]) {
      const rr = wfR + dr;
      if (rr > 0 && rr < R * 1.1) {
        const rPx = Math.abs(ax.xToPx(rr) - ax.xToPx(0));
        ctx.beginPath();
        ctx.arc(txc, tyc, rPx, 0, 2 * Math.PI);
        ctx.strokeStyle = 'rgba(43,108,176,0.18)';
        ctx.lineWidth = 1;
        ctx.setLineDash([]);
        ctx.stroke();
      }
    }
    const wfPxR = Math.abs(ax.xToPx(wfR) - ax.xToPx(0));
    ctx.beginPath();
    ctx.arc(txc, tyc, wfPxR, 0, 2 * Math.PI);
    ctx.strokeStyle = 'rgba(43,108,176,0.75)';
    ctx.lineWidth = 1.8;
    ctx.setLineDash([]);
    ctx.stroke();
  }

  if (!traj) { ctx.restore(); return; }

  const k = indexAt(currentTime);
  const base = k * traj.N;
  const N = traj.N;

  // Draw cells: batch free and trapped separately for performance.
  const CELL_R = 2.5;

  // Free cells (orange).
  ctx.beginPath();
  for (let i = 0; i < N; i++) {
    if (traj.trapped[base + i]) continue;
    const cx = ax.xToPx(traj.xs[base + i]);
    const cy = ax.yToPx(traj.ys[base + i]);
    ctx.moveTo(cx + CELL_R, cy);
    ctx.arc(cx, cy, CELL_R, 0, 2 * Math.PI);
  }
  ctx.fillStyle = '#b34700';
  ctx.fill();

  // Trapped cells (green).
  ctx.beginPath();
  for (let i = 0; i < N; i++) {
    if (!traj.trapped[base + i]) continue;
    const cx = ax.xToPx(traj.xs[base + i]);
    const cy = ax.yToPx(traj.ys[base + i]);
    ctx.moveTo(cx + CELL_R, cy);
    ctx.arc(cx, cy, CELL_R, 0, 2 * Math.PI);
  }
  ctx.fillStyle = '#2f7a3a';
  ctx.fill();

  ctx.restore();
}

// Generic time-trace for aggregate channels.
function drawTimeTrace(canvasId, arr, ranges, color) {
  const cv = document.getElementById(canvasId);
  const ctx = autoFit(cv);
  const w = cv.clientWidth, h = cv.clientHeight;
  const T_end = trajTEnd();
  if (!traj) {
    const ax = makeAxis({ xMin: 0, xMax: T_end, yMin: 0, yMax: 1, w, h });
    drawFrame(ctx, ax);
    return;
  }
  let { lo, hi } = ranges;
  if (!isFinite(lo)) { lo = 0; hi = 1; }
  const pad = 0.1 * (hi - lo || 1);
  const ax = makeAxis({ xMin: 0, xMax: T_end, yMin: lo - pad, yMax: hi + pad, w, h });
  drawFrame(ctx, ax);
  clipPlot(ctx, ax);

  const N_steps = traj.N_steps;
  // ts is rebuilt from the stored traj time grid (same for all channels).
  const ts = traj._ts || (() => {
    const a = new Float64Array(N_steps);
    const dtt = T_end / (N_steps - 1);
    for (let k = 0; k < N_steps; k++) a[k] = k * dtt;
    traj._ts = a;
    return a;
  })();
  strokePath(ctx, ax, ts, arr, { color, width: 1.5 });

  ctx.restore();
  // Scrub cursor.
  ctx.save();
  ctx.strokeStyle = 'rgba(0,0,0,0.5)';
  ctx.setLineDash([2, 2]);
  ctx.beginPath();
  ctx.moveTo(ax.xToPx(currentTime), ax.padT);
  ctx.lineTo(ax.xToPx(currentTime), ax.padT + ax.plotH);
  ctx.stroke();
  ctx.restore();
  // Current value dot.
  const k = indexAt(currentTime);
  dot(ctx, ax, currentTime, arr[k], 4, color);
}

function drawRmean()   { if (traj) drawTimeTrace('cv-rmean',   traj.rmean, traj.ranges.rmean,   '#b34700'); else drawTimePlaceholder('cv-rmean'); }
function drawTrapped() {
  if (!traj) { drawTimePlaceholder('cv-trapped'); return; }
  // Normalize nTrap to fraction; cached on traj. Adaptive y-range from data
  // (floor at 0; pads computed inside drawTimeTrace via ranges).
  if (!traj._trapFrac) {
    const a = new Float64Array(traj.N_steps);
    let hi = 0;
    for (let k = 0; k < traj.N_steps; k++) { a[k] = traj.nTrap[k] / traj.N; if (a[k] > hi) hi = a[k]; }
    traj._trapFrac = a;
    traj._trapFracRange = { lo: 0, hi: Math.max(hi, 1e-3) };
  }
  drawTimeTrace('cv-trapped', traj._trapFrac, traj._trapFracRange, '#2f7a3a');
}
function drawPmean()   { if (traj) drawTimeTrace('cv-pmean',   traj.pmean, traj.ranges.pmean,   '#2b6cb0'); else drawTimePlaceholder('cv-pmean'); }

function drawTimePlaceholder(id) {
  const cv = document.getElementById(id);
  const ctx = autoFit(cv);
  const w = cv.clientWidth, h = cv.clientHeight;
  const ax = makeAxis({ xMin: 0, xMax: trajTEnd(), yMin: 0, yMax: 1, w, h });
  drawFrame(ctx, ax);
}

function drawRho() {
  const cv = document.getElementById('cv-rho');
  const ctx = autoFit(cv);
  const w = cv.clientWidth, h = cv.clientHeight;
  if (!traj) {
    const ax = makeAxis({ xMin: 0, xMax: params.R_dish, yMin: 0, yMax: 1, w, h });
    drawFrame(ctx, ax);
    return;
  }
  const k = indexAt(currentTime);
  const base = k * traj.N;
  const N = traj.N;
  const R = params.R_dish;
  const nBins = 30;
  const binW = R / nBins;
  const counts = new Float64Array(nBins);
  for (let i = 0; i < N; i++) {
    const r = Math.hypot(traj.xs[base + i], traj.ys[base + i]);
    const b = Math.min(nBins - 1, Math.floor(r / binW));
    counts[b]++;
  }
  // Plot ρ(r)/r ∝ surface density σ(r) = ρ(r)/(2π r). Initial uniform-in-area
  // distribution gives counts[b] ∝ r_b, so dividing by r_b yields a flat
  // baseline at t̃ = 0. Constant 1/(2π·N·binW) is dropped (visualization).
  const binXs = new Float64Array(nBins);
  const sigma = new Float64Array(nBins);
  let hi = 0;
  for (let b = 0; b < nBins; b++) {
    const rc = (b + 0.5) * binW;
    binXs[b] = rc;
    sigma[b] = counts[b] / (rc * N);
    if (sigma[b] > hi) hi = sigma[b];
  }
  const pad = 0.1 * (hi || 0.01);
  const ax = makeAxis({ xMin: 0, xMax: R, yMin: 0, yMax: hi + pad, w, h });
  drawFrame(ctx, ax);
  clipPlot(ctx, ax);
  strokePath(ctx, ax, binXs, sigma, { color: '#b34700', width: 1.8 });

  // Vertical dashed line at r = R_trap.
  ctx.restore();
  clipPlot(ctx, ax);
  ctx.strokeStyle = 'rgba(43,108,176,0.7)';
  ctx.lineWidth = 1.2;
  ctx.setLineDash([4, 3]);
  ctx.beginPath();
  ctx.moveTo(ax.xToPx(params.R_trap), ax.padT);
  ctx.lineTo(ax.xToPx(params.R_trap), ax.padT + ax.plotH);
  ctx.stroke();
  ctx.setLineDash([]);
  ctx.restore();
}

// ─── controls ────────────────────────────────────────────────────────────────
const controlsEl = document.getElementById('controls');
const kpis = makeKpis([
  { id: 't',       label: 't̃' },
  { id: 'ntrap',   label: 'n_trap/N' },
  { id: 'rmeanFr', label: '⟨r̃⟩_free' },
  { id: 'pmeanFr', label: '⟨|P|⟩_free' },
]);

const sTime  = makeSlider({ id: 'time3',  symbol: 't̃',         value: 0,           min: 0,    max: 18,  step: 0.01, fmt: v => v.toFixed(2) });
const sM     = makeSlider({ id: 'M3',     symbol: 'M',          value: params.M,    min: 0,    max: 5,   step: 0.01, fmt: v => v.toFixed(2) });
const sC     = makeSlider({ id: 'C3',     symbol: 'C',          value: params.C,    min: 0,    max: 10,  step: 0.01, fmt: v => v.toFixed(2) });
const sChi   = makeSlider({ id: 'chi3',   symbol: 'χ̃',         value: params.chi,  min: 0,    max: 5,   step: 0.01, fmt: v => v.toFixed(2) });
const sMu    = makeSlider({ id: 'mu3',    symbol: 'μ̃',         value: params.mu,   min: 0,    max: 3,   step: 0.01, fmt: v => v.toFixed(2) });
const sLam   = makeSlider({ id: 'lam3',   symbol: 'λ',          value: params.lam,  min: 0.01, max: 10,  log: true,  fmt: v => v.toPrecision(3) });
const sTht   = makeSlider({ id: 'tht3',   symbol: 'ϑ',          value: params.tht,  min: 1e-4, max: 1,   log: true,  fmt: v => v.toExponential(2) });
const sV0    = makeSlider({ id: 'v03',    symbol: 'ṽ_0',        value: params.v0,   min: 0,    max: 3,   step: 0.01, fmt: v => v.toFixed(2) });
const sNHill = makeSlider({ id: 'nHill3', symbol: 'n',          value: params.nHill, min: 1,   max: 10,  step: 1,    fmt: v => v.toFixed(0) });
const sNWaves = makeSlider({ id: 'nWaves3', symbol: 'n_waves', value: params.nWaves, min: 1, max: 10, step: 1, fmt: v => v.toFixed(0) });
const sDtWave = makeSlider({ id: 'dtWave3', symbol: 'Δt_wave', value: params.dtWave, min: 0.5, max: 30, step: 0.1, fmt: v => v.toFixed(1) });
const sN     = makeSlider({ id: 'N3',     symbol: 'N',          value: params.N,    min: 100,  max: 5000, step: 100, fmt: v => v.toFixed(0) });
const sRDish = makeSlider({ id: 'Rdish3', symbol: 'R̃_dish',    value: params.R_dish, min: 4,  max: 20,  step: 0.5,  fmt: v => v.toFixed(1) });
const sRTrap = makeSlider({ id: 'Rtrap3', symbol: 'R̃_trap',    value: params.R_trap, min: 0.1, max: 5,  step: 0.1,  fmt: v => v.toFixed(1) });
const sDt    = makeSlider({ id: 'dt3',    symbol: 'dt̃',         value: params.dt,   min: 1e-4, max: 0.05, log: true, fmt: v => v.toExponential(2) });
const sSpeed = makeSlider({ id: 'speed3', symbol: 'play speed', value: params.speed, min: 0.01, max: 100, log: true, fmt: v => `${v.toPrecision(2)}×` });
const sSeed  = makeSlider({ id: 'seed3',  symbol: 'seed',       value: params.seed,  min: 1,    max: 9999, step: 1,  fmt: v => v.toFixed(0) });

sM.onChange(v     => { params.M   = v; recalibrate(); markDirty(); });
sC.onChange(v     => { params.C   = v; recalibrate(); markDirty(); });
sChi.onChange(v   => { params.chi = v; recalibrate(); markDirty(); });
sMu.onChange(v    => { params.mu  = v; recalibrate(); markDirty(); });
sLam.onChange(v   => { params.lam = v; markDirty(); });
sTht.onChange(v   => { params.tht = v; trajDirty = true; });
sV0.onChange(v    => { params.v0    = v; markDirty(); });
sNHill.onChange(v => { params.nHill = Math.round(v); markDirty(); });
sNWaves.onChange(v => { params.nWaves = Math.round(v); markDirty(); });
sDtWave.onChange(v => { params.dtWave = v; markDirty(); });
sN.onChange(v     => { params.N   = Math.round(v); markDirty(); });
sRDish.onChange(v => { params.R_dish = v; markDirty(); });
sRTrap.onChange(v => { params.R_trap = v; markDirty(); });
sDt.onChange(v    => { params.dt    = v; markDirty(); });
sSpeed.onChange(v => { params.speed = v; });
sSeed.onChange(v  => { params.seed  = Math.round(v); markDirty(); });
sTime.onChange(v  => { currentTime = v; });

const linkedReadout = () =>
  `→ M=${params.M.toFixed(2)}, C=${params.C.toPrecision(3)}, ` +
  `χ̃=${params.chi.toFixed(2)}, μ̃=${params.mu.toFixed(2)}`;
function pushAllNondimSliders() {
  applyingDim = true;
  sM.set(params.M); sC.set(params.C); sChi.set(params.chi); sMu.set(params.mu);
  applyingDim = false;
}

const sLmax  = makeSlider({ id: 'Lmax3',  symbol: 'L_max', value: dim.Lmax,  min: 0.1, max: 4, step: 0.01, fmt: v => v.toFixed(2), linkedLabel: () => linkedReadout() });
const sSigma = makeSlider({ id: 'sigma3', symbol: 'σ',     value: dim.sigma, min: 0.1, max: 4, step: 0.01, fmt: v => v.toFixed(2), linkedLabel: () => linkedReadout() });
const sCwave = makeSlider({ id: 'c3',     symbol: 'c',     value: dim.c,     min: 0.1, max: 4, step: 0.01, fmt: v => v.toFixed(2), linkedLabel: () => linkedReadout() });
function refreshDimReadouts() {
  sLmax.setLinkedText(linkedReadout()); sSigma.setLinkedText(linkedReadout()); sCwave.setLinkedText(linkedReadout());
}
sLmax.onChange(v  => { dim.Lmax  = v; recomputeFromDim(); pushAllNondimSliders(); refreshDimReadouts(); markDirty(); });
sSigma.onChange(v => { dim.sigma = v; recomputeFromDim(); pushAllNondimSliders(); refreshDimReadouts(); markDirty(); });
sCwave.onChange(v => { dim.c     = v; recomputeFromDim(); pushAllNondimSliders(); refreshDimReadouts(); markDirty(); });

const buttons = makeButtonRow([
  { label: '▶  play', onClick() {
      playing = !playing;
      if (playing && currentTime >= trajTEnd() - 1e-6) {
        currentTime = 0; sTime.set(currentTime);
      }
      buttons.refs['▶  play'].textContent = playing ? '⏸  pause' : '▶  play';
    } },
  { label: '⟳  reset', ghost: true, onClick: () => {
      currentTime = 0; sTime.set(currentTime);
    } },
]);

// ─── model toggle (top of controls) ──────────────────────────────────────────
// Per-model defaults applied on toggle to shared (M, C, χ̃, ϑ) and
// model-specific sliders. Geometry / wave-train / numerics sliders are not
// touched.
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
    saveModelChoice();
    applyModelDefaults(v);
    applyModelVisibility();
    markDirty();
  },
});
controlsEl.appendChild(modelToggle.el);

controlsEl.appendChild(section('time scrub', [sTime.el, buttons.el]));
controlsEl.appendChild(section('wave (nondim)', [sM.el, sC.el, sNWaves.el, sDtWave.el]));
controlsEl.appendChild(section('coupling (nondim)', [sChi.el]));
const secGLvel   = section('GL velocity (nondim)',   [sMu.el]);
const secHillVel = section('Hill velocity (nondim)', [sV0.el, sNHill.el]);
const secGLwell  = section('GL well (nondim)',        [sLam.el]);
controlsEl.appendChild(secGLvel);
controlsEl.appendChild(secHillVel);
controlsEl.appendChild(secGLwell);
controlsEl.appendChild(section('noise (nondim)', [sTht.el]));
controlsEl.appendChild(section('geometry', [sN.el, sRDish.el, sRTrap.el]));
controlsEl.appendChild(detailsSection('dim sliders (push linked nondim)', [sLmax.el, sSigma.el, sCwave.el]));
controlsEl.appendChild(section('numerics & playback', [sDt.el, sSpeed.el, sSeed.el]));
controlsEl.appendChild(kpis.el);
controlsEl.appendChild(el('div', { class: 'note' }, [
  'Trajectory window: t̃ ∈ [0, (R̃_dish + 8)/C]. ',
  'Full trajectory precomputed on any sim-parameter change; time slider scrubs the precomputed result. ',
  'GL anti-wave regime: at moderate χ̃ and λ ≈ 1, cells drift inward as the outward wave passes. ',
  'Hill mode: linear-relaxation SDE with saturating velocity; no bistable well. ',
  'Orange = free cells, green = trapped. Blue ring = wavefront at r̃ = C t̃.',
]));

// ─── model visibility ─────────────────────────────────────────────────────────
function applyModelVisibility() {
  const isHill = params.model === 'Hill';
  secGLvel.style.display   = isHill ? 'none' : '';
  secHillVel.style.display = isHill ? '' : 'none';
  secGLwell.style.display  = isHill ? 'none' : '';
}
applyModelVisibility();

// ─── sync from possibly-restored slider values ───────────────────────────────
params.M      = sM.value;
params.C      = sC.value;
params.chi    = sChi.value;
params.mu     = sMu.value;
params.lam    = sLam.value;
params.tht    = sTht.value;
params.v0     = sV0.value;
params.nHill  = Math.round(sNHill.value);
params.nWaves = Math.round(sNWaves.value);
params.dtWave = sDtWave.value;
params.N      = Math.round(sN.value);
params.R_dish = sRDish.value;
params.R_trap = sRTrap.value;
params.dt     = sDt.value;
params.speed  = sSpeed.value;
params.seed   = Math.round(sSeed.value);
dim.Lmax  = sLmax.value;
dim.sigma = sSigma.value;
dim.c     = sCwave.value;
recalibrate();
sTime.setMinMax(0, trajTEnd());
currentTime = Math.min(trajTEnd(), Math.max(0, sTime.value));
sTime.set(currentTime);

// ─── KaTeX decorations ───────────────────────────────────────────────────────
function decorateAll() {
  decoratePlot('cv-dish',    { titleTex: '\\text{2D radial swarm}' });
  decoratePlot('cv-rmean',   { titleTex: '\\text{mean radius of free cells}',
                               xLabelTex: '\\tilde t',
                               yLabelTex: '\\langle\\tilde r\\rangle_{\\text{free}}' });
  decoratePlot('cv-trapped', { titleTex: '\\text{trapped fraction}',
                               xLabelTex: '\\tilde t',
                               yLabelTex: 'n_{\\text{trap}}/N' });
  decoratePlot('cv-rho',     { titleTex: '\\rho(\\tilde r)/\\tilde r\\,\\propto\\,\\sigma(\\tilde r)\\text{ (current }\\tilde t\\text{)}',
                               xLabelTex: '\\tilde r',
                               yLabelTex: '\\rho(\\tilde r)/\\tilde r' });
  decoratePlot('cv-pmean',   { titleTex: '\\text{mean polarization magnitude (free cells)}',
                               xLabelTex: '\\tilde t',
                               yLabelTex: '\\langle|P|\\rangle_{\\text{free}}' });
}
if (window.katex) decorateAll();
else window.addEventListener('load', decorateAll);

// ─── main render loop ────────────────────────────────────────────────────────
function frame(nowMs) {
  if (trajDirty) recomputeTrajectory();

  if (playing) {
    const dtReal = lastFrameMs ? (nowMs - lastFrameMs) / 1000 : 0;
    currentTime += params.speed * dtReal;
    const T_end = trajTEnd();
    if (currentTime >= T_end) {
      currentTime = T_end;
      playing = false;
      buttons.refs['▶  play'].textContent = '▶  play';
    }
    sTime.set(currentTime);
  }
  lastFrameMs = nowMs;

  drawDish();
  drawRmean();
  drawTrapped();
  drawRho();
  drawPmean();

  if (traj) {
    const k = indexAt(currentTime);
    kpis.set('t',       currentTime.toFixed(2));
    kpis.set('ntrap',   (traj.nTrap[k] / traj.N).toFixed(3));
    kpis.set('rmeanFr', traj.rmean[k].toFixed(3));
    kpis.set('pmeanFr', traj.pmean[k].toFixed(3));
  }

  requestAnimationFrame(frame);
}
requestAnimationFrame(frame);
