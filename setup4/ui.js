// setup4/ui.js
// UI assembly for Setup 4: KNOBS, Calculate button, time scrub slider, KPIs.
//
// Interaction paradigm: calculate-then-explore.
//   1. User sets parameters via sliders.
//   2. User presses "Calculate" → worker runs, saves frames.
//   3. Time slider becomes active; scrubbing updates heatmap + radial profile.
//   4. Any parameter change invalidates the run and re-enables Calculate.
//
// KNOBS follows the setup3 pattern: each parameter declared once with
// exposure: 'dim' | 'nondim' | 'both'.
// Dim sliders live in a collapsed <details> section; moving them recomputes
// nondim sliders via recomputeFromDim(). recalibrate() saves current nondim
// values as the new dim→nondim calibration reference.

import { el, makeSlider, makeToggle, makeButtonRow, makeKpis, section,
         detailsSection, decoratePlot } from '../shared/dom.js';
import { dimToNondim, DIM_DEFAULTS } from './nondim.js';
import { drawDish, drawRadialProfile, drawMeanRadius, drawBead3D, drawBead2D, drawBead1D } from './render.js';

// ─── Nondim params (simulation state) ──────────────────────────────────────
const params = {
  geometry: '2d2d',   // '2d2d' | '2d3d'
  model:    'M1',
  // N is now DERIVED from σ̃·πR̃² (β = 1 by construction); the user-facing
  // cell-density knob is σ (dim), exposed in the Geometry section.
  N:        500,
  N_grid:   128,
  R_dish:   10,   // ~ 20 µm in 2D-2D, ~ 4 µm in 2D-3D with experimental σ_dim.
                  // For experimental-scale 2D-2D dish (200 µm), push slider to R̃≈100.
                  // 2D-3D experimental scale (R̃≈500) needs much heavier numerics.
  // Cell nondim.
  Lambda:   1,
  L_c_nd:   1,
  chi_nd:   0.75,
  mu_nd:    0.05,
  lam:      1,
  tht:      1e-4,
  // Cue M1 nondim.
  n_L:      10,
  gamma_L:  0,
  sigma_tilde: 0.02,  // σ̃ = σ·ℓ_0²; set by dim linkage. Per-cell emission = 1/σ̃.
  // Time-limited firing source (parent_solver convention).
  r_fire:   2.0,
  t_fire:   5.0,
  s_fire:   1.0,
  // Numerics. Per-geometry t̃_max stored so toggling restores the active value.
  dt:       0.01,
  t_max:        30,    // active value (matches initial geometry)
  t_max_2d2d:   30,
  t_max_2d3d:   30,
  seed:     7,
  // 2D-3D only.
  N_z:      16,
  h_0:      0.1,
  alpha_z:  1.1,  // gives z_max ≈ 3.6 ℓ₀ — sufficient for the 2D-3D wave penetration depth
};

// ─── Dim state + linkage ────────────────────────────────────────────────────
// Direct dim→nondim computation (no scale-from-calibration pattern). When the
// user moves a dim slider, nondim sliders are recomputed from the formula in
// nondim.js. Editing nondim sliders directly works but the value is overwritten
// on the next dim-slider change — this is intentional for v1; revisit if it
// becomes annoying.
const dim = Object.assign({}, DIM_DEFAULTS);
// R_dim and N are the primary Geometry inputs (dim). σ is *derived*:
//   σ_dim = N / (π R_dim²)
// so β ≡ N / (σ̃ · πR̃²) = 1 by physical consistency. Both R_dim and N are
// fully independent experimental knobs; σ is whatever the dish ends up at.
dim.R_dim = 200;   // µm — dish radius (default ~ neutrophil drop scale)
dim.N     = 3000;  // cells — total seeded

let applyingDim = false;

function recomputeFromDim() {
  if (applyingDim) return;
  applyingDim = true;
  // σ_dim is derived from the two Geometry knobs (R_dim, N). β = 1 always.
  if (dim.R_dim > 0) {
    dim.sigma = dim.N / (Math.PI * dim.R_dim * dim.R_dim);
  }
  const nd = dimToNondim(dim, params.geometry);
  params.ell0 = nd.ell0;
  // R̃ = R_dim / ℓ₀ (ℓ₀ already accounts for σ via the geometry-specific c*).
  if (isFinite(nd.ell0) && nd.ell0 > 0) {
    params.R_dish = dim.R_dim / nd.ell0;
  }
  params.N = dim.N;
  sTracked && sTracked.setMinMax(0, Math.max(0, params.N - 1));
  if (kpis) {
    kpis.set('N', String(params.N));
    kpis.set('sigma_dim', dim.sigma.toExponential(2));
  }
  // Push nondim values to sliders (if sliders exist).
  sLambda && sLambda.set(nd.Lambda);
  sLcNd   && sLcNd.set(nd.L_c_nd);
  sChiNd  && sChiNd.set(nd.chi_nd);
  sMuNd   && sMuNd.set(nd.mu_nd);
  sLam    && sLam.set(nd.lam);
  sTht    && sTht.set(nd.tht);
  sGammaL && sGammaL.set(nd.gamma_L);
  // Update params.
  params.Lambda      = nd.Lambda;
  params.L_c_nd      = nd.L_c_nd;
  params.chi_nd      = nd.chi_nd;
  params.mu_nd       = nd.mu_nd;
  params.lam         = nd.lam;
  params.tht         = nd.tht;
  params.gamma_L     = nd.gamma_L;
  params.sigma_tilde = nd.sigma_tilde;
  if (kpis) kpis.set('sigma_tilde', nd.sigma_tilde.toExponential(2));
  applyingDim = false;
  markDirty();
}

// Editing a nondim slider directly is allowed; nothing to recalibrate in the
// direct-formula scheme (kept as a no-op so onChange handlers stay parallel
// to setup3's idiom in case we add cal-scaling later).
function recalibrate() { /* no-op, see comment above `dim` */ }

// ─── Frame storage ───────────────────────────────────────────────────────────
// Every frame carries a Uint8-compressed Lfield (16 KB at N_grid=128), so no
// nearest-snapshot fallback is needed — heatmap animates smoothly.
/** @type {Array<{step, t, radialProfile, agentX, agentY, emitting, Lfield, Lmax}>} */
let frames = [];

let runDirty = true;
function markDirty() { runDirty = true; enableCalculate(); }

// ─── Time scrub ─────────────────────────────────────────────────────────────
let currentFrameIdx = 0;
let currentTime = 0;
let trackedCellIdx = 0;   // which cell to follow in the bead-plot

function getFrameAt(t) {
  if (frames.length === 0) return null;
  // Binary search for nearest frame by time.
  let lo = 0, hi = frames.length - 1;
  while (lo < hi) {
    const mid = (lo + hi + 1) >> 1;
    if (frames[mid].t <= t) lo = mid; else hi = mid - 1;
  }
  return frames[lo];
}

// ─── Worker ─────────────────────────────────────────────────────────────────
let worker = null;
let workerBusy = false;

function startWorker() {
  // Terminate any in-flight worker before swapping in a new one.
  // Without terminate(), the old worker's onmessage closure still resolves the
  // *current* (newly-reset) `frames` array, leaking stale frames into the new run.
  if (worker) {
    try { worker.terminate(); } catch {}
    worker = null;
  }
  worker = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
  workerBusy = true;
  frames = [];
  currentFrameIdx = 0;

  // Build params for worker.
  const wp = {
    geometry: params.geometry,
    model:    params.model,
    N:        params.N,
    N_grid:   params.N_grid,
    R_dish:   params.R_dish,
    Lambda:   params.Lambda,
    L_c_nd:   params.L_c_nd,
    chi_nd:   params.chi_nd,
    mu_nd:    params.mu_nd,
    lam:      params.lam,
    tht:      params.tht,
    n_L:          params.n_L,
    gamma_L:      params.gamma_L,
    sigma_tilde:  params.sigma_tilde,
    r_fire:   params.r_fire,
    t_fire:   params.t_fire,
    s_fire:   params.s_fire,
    dt:       params.dt,
    t_max:    params.t_max,
    seed:     params.seed,
    N_z:      params.N_z,
    h_0:      params.h_0,
    alpha_z:  params.alpha_z,
  };

  worker.postMessage({ type: 'run', params: wp });

  worker.onmessage = function (e) {
    const msg = e.data;

    if (msg.type === 'frame') {
      frames.push({
        step: msg.step, t: msg.t,
        radialProfile: msg.radialProfile,
        agentX: msg.agentX, agentY: msg.agentY,
        Px:     msg.agentPx, Py: msg.agentPy,
        Gx:     msg.agentGx, Gy: msg.agentGy,   // ∇𝓛 sampled at each cell
        emitting: msg.emitting,
        Lfield: msg.Lfield,
        Lmax:   msg.Lmax,
      });
      // Update display to latest.
      currentFrameIdx = frames.length - 1;
      currentTime = msg.t;
      sTime && sTime.set(msg.t);
      redraw();
    }

    if (msg.type === 'progress') {
      const pct = (msg.pct * 100).toFixed(0);
      kpis && kpis.set('status', `Running… ${pct}%`);
      kpis && kpis.set('t', msg.t.toFixed(2));
    }

    if (msg.type === 'done') {
      workerBusy = false;
      runDirty = false;
      kpis && kpis.set('status', 'Done');
      sTime && sTime.setMinMax(0, params.t_max);
      sTime && sTime.set(params.t_max);
      drawMeanRadius('cv-rmean', frames, { R_dish: params.R_dish });
      updateWaveSpeed();
    }

    if (msg.type === 'error') {
      workerBusy = false;
      kpis && kpis.set('status', 'Error: ' + msg.message);
      console.error('Setup4 worker error:', msg.message);
    }
  };

  worker.onerror = function (e) {
    workerBusy = false;
    kpis && kpis.set('status', 'Worker error: ' + e.message);
    console.error('Setup4 worker onerror:', e);
  };

  kpis && kpis.set('status', 'Running…');
  disableCalculate();
}

// ─── KPI: crude wave speed estimate ─────────────────────────────────────────
function updateWaveSpeed() {
  if (frames.length < 4) return;
  // Find frame at ~10% and ~80% of t_max; track peak of radial profile.
  const f1 = getFrameAt(params.t_max * 0.1);
  const f2 = getFrameAt(params.t_max * 0.8);
  if (!f1 || !f2 || !f1.radialProfile || !f2.radialProfile) return;

  const peakR = (prof) => {
    let pkIdx = 0, pkVal = 0;
    for (let k = 0; k < prof.length; k++) {
      if (prof[k] > pkVal) { pkVal = prof[k]; pkIdx = k; }
    }
    const dx = (2 * params.R_dish) / (params.N_grid - 1);
    return pkIdx * dx;
  };

  const r1 = peakR(f1.radialProfile);
  const r2 = peakR(f2.radialProfile);
  const dt = f2.t - f1.t;
  if (dt < 1e-9) return;
  const cEff = (r2 - r1) / dt;
  kpis && kpis.set('ceff', cEff.toFixed(3));
}

// ─── Redraw ──────────────────────────────────────────────────────────────────
function redraw() {
  const f = (frames.length > 0) ? frames[currentFrameIdx] : null;
  if (!f) return;

  drawDish('cv-dish', f, {
    N_grid: params.N_grid, R_dish: params.R_dish, t: f.t,
    trackedCellIdx,
  });
  drawRadialProfile('cv-profile', f, {
    N_grid: params.N_grid, R_dish: params.R_dish,
    model: params.model, geometry: params.geometry,
    N_cells: params.N, sigma_tilde: params.sigma_tilde,
  });
  const beadParams = {
    Lambda: params.Lambda, L_c: params.L_c_nd, lam: params.lam,
    chi: params.chi_nd,
    R_dish: params.R_dish, N_grid: params.N_grid,
  };
  drawBead3D('cv-bead-3d', frames, currentFrameIdx, trackedCellIdx, beadParams);
  drawBead2D('cv-bead-2d', frames, currentFrameIdx, trackedCellIdx, beadParams);
  drawBead1D('cv-bead-1d', frames, currentFrameIdx, trackedCellIdx, beadParams);

  // KPIs.
  if (kpis) {
    kpis.set('t', f.t.toFixed(2));
    // Absolute peak L is the grid max stored on the frame; fall back to radial.
    let Lmax = f.Lmax || 0;
    if (!Lmax && f.radialProfile) {
      for (let k = 0; k < f.radialProfile.length; k++) {
        if (f.radialProfile[k] > Lmax) Lmax = f.radialProfile[k];
      }
    }
    kpis.set('Lmax', Lmax.toFixed(2));
    if (f.emitting) {
      let nE = 0;
      for (let k = 0; k < f.emitting.length; k++) nE += f.emitting[k];
      kpis.set('nemit', `${nE} / ${f.emitting.length}`);
    }
  }
}

// ─── Slider references (populated during build) ──────────────────────────────
let sLambda, sLcNd, sChiNd, sMuNd, sLam, sTht, sGammaL, sTime, sTmax, sH, sA, sTracked;
let kpis;
let btnCalc, btnReset;

function enableCalculate()  { if (btnCalc) { btnCalc.disabled = false; btnCalc.textContent = 'Calculate'; } }
function disableCalculate() { if (btnCalc) { btnCalc.disabled = true;  btnCalc.textContent = 'Running…'; } }

function resetRun() {
  // Terminate any in-flight worker and wipe accumulated state.
  if (worker) {
    try { worker.terminate(); } catch {}
    worker = null;
  }
  workerBusy = false;
  frames = [];
  currentFrameIdx = 0;
  runDirty = true;
  if (sTime) sTime.set(0);
  if (kpis) {
    kpis.set('status', 'Press Calculate');
    kpis.set('t',      '–');
    kpis.set('Lmax',   '–');
    kpis.set('nemit',  '–');
    kpis.set('ceff',   '–');
  }
  // Blank the canvases by drawing an empty frame.
  ['cv-dish', 'cv-profile', 'cv-rmean'].forEach(id => {
    const c = document.getElementById(id);
    if (c) { const ctx = c.getContext('2d'); ctx.clearRect(0, 0, c.width, c.height); }
  });
  enableCalculate();
}

// ─── Build UI ────────────────────────────────────────────────────────────────
export function buildUI(containerId) {
  const container = document.getElementById(containerId);
  if (!container) { console.error('setup4/ui: container not found:', containerId); return; }

  // ── Calculate + Reset buttons ──
  btnCalc = el('button', { type: 'button',
    style: { flex: '2', fontWeight: '700', fontSize: '1.05em' } }, 'Calculate');
  btnCalc.addEventListener('click', () => {
    if (runDirty || !workerBusy) startWorker();
  });
  btnReset = el('button', { type: 'button',
    style: { flex: '1', fontWeight: '600' } }, 'Reset');
  btnReset.addEventListener('click', resetRun);
  const btnRow = el('div', { style: { display: 'flex', gap: '6px', marginBottom: '8px' } });
  btnRow.appendChild(btnCalc);
  btnRow.appendChild(btnReset);
  container.appendChild(btnRow);

  // ── KPIs ──
  kpis = makeKpis([
    { id: 'status', label: 'Status', init: 'Press Calculate' },
    { id: 't',      label: 't̃',    init: '–' },
    { id: 'Lmax',   label: 'max 𝓛', init: '–' },
    { id: 'nemit',  label: 'emitting', init: '–' },
    { id: 'ceff',   label: 'c_eff', init: '–' },
    { id: 'sigma_tilde', label: 'σ̃', init: '–' },
    { id: 'sigma_dim',   label: 'σ', init: '–' },
    { id: 'N',     label: 'N',     init: '–' },
  ]);
  container.appendChild(kpis.el);

  // ── Time scrub slider ──
  sTime = makeSlider({
    id: 's4-time', symbol: 't̃', value: 0, min: 0, max: params.t_max, step: 'any',
    fmt: v => v.toFixed(2),
  });
  sTime.onChange(v => {
    currentTime = v;
    currentFrameIdx = 0;
    if (frames.length === 0) return;
    let best = 0;
    for (let i = 0; i < frames.length; i++) {
      if (frames[i].t <= v) best = i; else break;
    }
    currentFrameIdx = best;
    redraw();
  });
  container.appendChild(section('Time scrub', [sTime.el]));

  // ── Geometry toggle ──
  const geoToggle = makeToggle({
    label: 'Geometry',
    options: [{ id: '2d2d', label: '2D–2D' }, { id: '2d3d', label: '2D–3D' }],
    value: params.geometry,
    onChange(v) {
      // Save current t̃_max into the geometry we're leaving, restore for the new one.
      if (params.geometry === '2d2d') params.t_max_2d2d = params.t_max;
      else                            params.t_max_2d3d = params.t_max;
      params.geometry = v;
      const newTmax = (v === '2d2d') ? params.t_max_2d2d : params.t_max_2d3d;
      params.t_max = newTmax;
      sTmax && sTmax.set(newTmax);
      sTime && sTime.setMinMax(0, newTmax);
      // Show/hide geometry-specific sliders. h is the 2D-2D layer height;
      // it doesn't appear in the 2D-3D PDE.
      sH.el.style.display     = (v === '2d2d') ? '' : 'none';
      nzRow.style.display     = (v === '2d3d') ? '' : 'none';
      h0Row.style.display     = (v === '2d3d') ? '' : 'none';
      alphaRow.style.display  = (v === '2d3d') ? '' : 'none';
      recomputeFromDim();
      markDirty();
    },
  });

  // ── Geometry section ──
  // R_dim and N are the two free dim inputs. σ_dim is derived (β = 1 always).
  // R̃ = R_dim/ℓ₀ is derived downstream in recomputeFromDim.
  const sN_geom = makeSlider({
    id: 's4-N-geom', symbol: 'N', value: dim.N,
    min: 1, max: 1e5, step: null, log: true,
    fmt: v => String(Math.max(1, Math.round(v))),
  });
  sN_geom.onChange(v => {
    if (applyingDim) return;
    dim.N = Math.max(1, Math.round(v));
    recomputeFromDim();
  });

  const sRdim_geom = makeSlider({
    id: 's4-Rdim-geom', symbol: 'R_{\\text{dish}}', value: dim.R_dim,
    min: 1, max: 500, step: 1,
    units: 'µm',
  });
  sRdim_geom.onChange(v => {
    if (applyingDim) return;
    dim.R_dim = v;
    recomputeFromDim();
  });

  const sNgrid = makeSlider({ id: 's4-Ngrid', symbol: 'N_{\\text{grid}}', value: params.N_grid, min: 32, max: 256, step: 1, fmt: v => String(Math.round(v)) });
  sNgrid.onChange(v => { params.N_grid = Math.round(v); markDirty(); });

  container.appendChild(section('Geometry', [geoToggle.el, sN_geom.el, sRdim_geom.el, sNgrid.el]));

  // ── Cell model toggle ──
  const cellModelToggle = makeToggle({
    label: 'Cell model',
    options: [{ id: 'GL', label: 'GL' }, { id: 'Hill', label: 'Hill (soon)' }],
    value: 'GL',
    onChange(v) {
      if (v === 'Hill') { alert('Hill model not yet implemented for Setup 4.'); cellModelToggle.set('GL'); return; }
      markDirty();
    },
  });
  container.appendChild(section('Cell model', [cellModelToggle.el]));

  // ── Cue (M1) nondim section ──
  sGammaL = makeSlider({ id: 's4-gammaL', symbol: '\\tilde{\\Gamma}_L', value: params.gamma_L, min: 0, max: 5, step: 0.01 });
  sGammaL.onChange(v => { params.gamma_L = v; recalibrate(); markDirty(); });

  const snL = makeSlider({ id: 's4-nL', symbol: 'n_L', value: params.n_L, min: 1, max: 50, step: 1, fmt: v => String(Math.round(v)) });
  snL.onChange(v => { params.n_L = Math.round(v); markDirty(); });

  container.appendChild(section('Cue model (M1)', [sGammaL.el, snL.el]));

  // ── Cell nondim section ──
  sLambda = makeSlider({ id: 's4-Lambda', symbol: '\\Lambda', value: params.Lambda, min: 0.01, max: 20, step: 0.01 });
  sLambda.onChange(v => { params.Lambda = v; recalibrate(); markDirty(); });

  sLcNd = makeSlider({ id: 's4-Lc', symbol: '\\tilde{\\mathcal{L}}_c', value: params.L_c_nd, min: 0.1, max: 5, step: 0.05 });
  sLcNd.onChange(v => { params.L_c_nd = v; recalibrate(); markDirty(); });

  sChiNd = makeSlider({ id: 's4-chi', symbol: '\\tilde{\\chi}', value: params.chi_nd, min: 0, max: 5, step: 0.05 });
  sChiNd.onChange(v => { params.chi_nd = v; recalibrate(); markDirty(); });

  sMuNd = makeSlider({ id: 's4-mu', symbol: '\\tilde{\\mu}', value: params.mu_nd, min: 0, max: 1, step: 0.005 });
  sMuNd.onChange(v => { params.mu_nd = v; recalibrate(); markDirty(); });

  sLam = makeSlider({ id: 's4-lam', symbol: '\\lambda', value: params.lam, min: 0.01, max: 10, step: 0.01 });
  sLam.onChange(v => { params.lam = v; recalibrate(); markDirty(); });

  sTht = makeSlider({ id: 's4-tht', symbol: '\\vartheta', value: params.tht, min: 0, max: 0.01, step: 1e-5, log: false });
  sTht.onChange(v => { params.tht = v; recalibrate(); markDirty(); });

  container.appendChild(section('Cell nondim', [sLambda.el, sLcNd.el, sChiNd.el, sMuNd.el, sLam.el, sTht.el]));

  // ── IC section ──
  // Time-limited firing source (replaces the static seed IC). Same trigger
  // mechanism as parent_solver / Dieterle-limit script — adds an extra
  // source term inside r̃ < r_fire for t̃ < t_fire, with magnitude s_fire
  // (the canonical Dieterle wave at c̃ = 2/π corresponds to s_fire = 1).
  const sRFire = makeSlider({ id: 's4-rfire', symbol: '\\tilde{r}_{\\text{fire}}', value: params.r_fire, min: 0.1, max: 30, step: 0.1 });
  sRFire.onChange(v => {
    const clamped = Math.min(v, params.R_dish);
    params.r_fire = clamped;
    if (clamped !== v) sRFire.set(clamped);
    markDirty();
  });

  const sTFire = makeSlider({ id: 's4-tfire', symbol: '\\tilde{t}_{\\text{fire}}', value: params.t_fire, min: 0, max: 100, step: 0.1 });
  sTFire.onChange(v => { params.t_fire = v; markDirty(); });

  const sSFire = makeSlider({ id: 's4-sfire', symbol: 's_{\\text{fire}}', value: params.s_fire, min: 0, max: 20, step: 0.05 });
  sSFire.onChange(v => { params.s_fire = v; markDirty(); });

  container.appendChild(section('Initial stimulus (firing source)', [sRFire.el, sTFire.el, sSFire.el]));

  // ── Numerics section ──
  const sDt = makeSlider({ id: 's4-dt', symbol: 'd\\tilde{t}', value: params.dt, min: 0.001, max: 0.1, step: 0.001, log: false });
  sDt.onChange(v => { params.dt = v; markDirty(); });

  sTmax = makeSlider({ id: 's4-tmax', symbol: '\\tilde{t}_{\\max}', value: params.t_max, min: 0.1, max: 1000, step: 0.1 });
  sTmax.onChange(v => {
    params.t_max = v;
    // Keep the per-geometry slot in sync so a toggle round-trip is idempotent.
    if (params.geometry === '2d2d') params.t_max_2d2d = v;
    else                            params.t_max_2d3d = v;
    sTime.setMinMax(0, v);
    markDirty();
  });

  const sSeed = makeSlider({ id: 's4-seed', symbol: '\\text{seed}', value: params.seed, min: 1, max: 999, step: 1, fmt: v => String(Math.round(v)) });
  sSeed.onChange(v => { params.seed = Math.round(v); markDirty(); });

  container.appendChild(section('Numerics', [sDt.el, sTmax.el, sSeed.el]));

  // ── Tracked cell (for bead-in-free-energy view) ──
  // sTracked must be declared at function scope so recomputeFromDim can
  // update its max when N (derived from σ̃·πR̃²) changes.
  sTracked = makeSlider({
    id: 's4-tracked', symbol: '\\text{tracked cell}',
    value: trackedCellIdx, min: 0, max: Math.max(0, params.N - 1), step: 1,
    fmt: v => String(Math.round(v)),
  });
  sTracked.onChange(v => { trackedCellIdx = Math.round(v); redraw(); });
  container.appendChild(section('Bead view', [sTracked.el]));

  // ── 2D-3D extra sliders (shown/hidden by geometry toggle) ──
  const sNz = makeSlider({ id: 's4-Nz', symbol: 'N_z', value: params.N_z, min: 4, max: 64, step: 1, fmt: v => String(Math.round(v)) });
  sNz.onChange(v => { params.N_z = Math.round(v); markDirty(); });

  const sH0 = makeSlider({ id: 's4-h0', symbol: 'h_0', value: params.h_0, min: 0.01, max: 1, step: 0.01 });
  sH0.onChange(v => { params.h_0 = v; markDirty(); });

  const sAlpha = makeSlider({ id: 's4-alpha', symbol: '\\alpha_z', value: 1.1, min: 1.05, max: 1.5, step: 0.05 });
  sAlpha.onChange(v => { params.alpha_z = v; markDirty(); });

  const nzRow    = sNz.el;
  const h0Row    = sH0.el;
  const alphaRow = sAlpha.el;

  // h slider for 2D-2D dim section.
  sH = makeSlider({
    id: 's4-h', symbol: 'h', value: dim.h, min: 1e-3, max: 1e6, step: null, log: true,
    units: 'µm',
  });
  sH.onChange(v => { dim.h = v; recomputeFromDim(); });
  attachDimDesc(sH, 'Layer thickness (2D-2D only)');

  container.appendChild(detailsSection('2D-3D grid (z-direction)', [nzRow, h0Row, alphaRow]));

  // ── Dimensional inputs (collapsed) ──
  // Helper: attach a small description text above the slider's label row.
  function attachDimDesc(s, name) {
    if (!name) return;
    const desc = el('div', {
      class: 'dim-desc',
      style: { fontSize: '0.8em', color: '#888', margin: '6px 0 -2px 0', fontStyle: 'italic' },
    }, name);
    s.el.insertBefore(desc, s.el.firstChild);
  }
  const makeDimSlider = (id, sym, key, min, max, step, logScale, name, unit) => {
    const s = makeSlider({
      id: `s4-dim-${id}`, symbol: sym, value: dim[key],
      min, max, step, log: !!logScale,
      units: unit || '',
    });
    s.onChange(v => { dim[key] = v; recomputeFromDim(); });
    attachDimDesc(s, name);
    return s;
  };

  // R_dim and N live in the Geometry section (top of UI); σ_dim is derived
  // and shown as a KPI. No bidirectional dim mirrors here — those caused
  // confusion about which knob was "free".
  sA = makeSlider({
    id: 's4-dim-a', symbol: 'a', value: dim.a, min: 1, max: 1e9, step: null, log: true,
    units: 'nM·µm³/s',
  });
  sA.onChange(v => { dim.a = v; recomputeFromDim(); });
  attachDimDesc(sA, 'Per-cell LTB4 emission strength');

  const sDL    = makeDimSlider('DL',    'D_L',      'D_L',     1,     1000, null, true,
                              'LTB4 diffusion constant', 'µm²/s');
  const sL0    = makeDimSlider('L0',    'L_0',      'L_0',     0.01,  100,  null, true,
                              'LTB4 Hill threshold (emission)', 'nM');
  const sLcDim = makeDimSlider('Lc',    'L_c',      'L_c',     0.01,  100,  null, true,
                              'GL polarization threshold', 'nM');
  const sR0    = makeDimSlider('r0',    'r_0',      'r_0',     0.001, 10,   null, true,
                              'GL relaxation rate prefactor', '1/(s·nM)');
  const sU     = makeDimSlider('u',     'u',        'u',       0.001, 10,   null, true,
                              'GL nonlinear coefficient u', '1/s');
  const sW     = makeDimSlider('w',     'w',        'w',       0.001, 10,   null, true,
                              'GL nonlinear coefficient w', '1/s');
  const sChiD  = makeDimSlider('chi',   '\\chi',    'chi',     0.001, 100,  null, true,
                              'Chemotactic susceptibility', 'µm²/(s·nM)');
  const sMuD   = makeDimSlider('mu',    '\\mu',     'mu',      0.001, 10,   null, true,
                              'Self-propulsion speed per polarization', 'µm/s');
  const sTheta = makeDimSlider('theta', '\\theta',  'theta',   1e-7,  0.1,  null, true,
                              'Noise amplitude', 'nM/s');
  const sGLdec = makeDimSlider('GamL',  '\\Gamma_L','Gamma_L', 0,     10,   0.01, false,
                              'LTB4 decay rate', '1/s');

  // R_dim, N live in Geometry (primary knobs); σ_dim is derived (KPI only).
  const dimChildren = [sA.el, sDL.el, sL0.el, sLcDim.el, sR0.el,
                       sU.el, sW.el, sChiD.el, sMuD.el, sTheta.el, sGLdec.el, sH.el];
  container.appendChild(detailsSection('Dimensional inputs', dimChildren));

  // ── Initial visibility based on geometry ──
  nzRow.style.display    = (params.geometry === '2d3d') ? '' : 'none';
  h0Row.style.display    = (params.geometry === '2d3d') ? '' : 'none';
  alphaRow.style.display = (params.geometry === '2d3d') ? '' : 'none';
  sH.el.style.display = (params.geometry === '2d2d') ? '' : 'none';

  // Initial dim → nondim sync so params.R_dish, σ̃, etc. are populated.
  recomputeFromDim();

  // ── Decorate plots with KaTeX axis labels ──
  // Called after DOMContentLoaded; KaTeX may not be ready yet. We defer.
  function decoratePlots() {
    decoratePlot('cv-dish',    { titleTex: '\\mathcal{L}(\\tilde{x},\\tilde{y})' });
    decoratePlot('cv-profile', {
      titleTex:  '\\mathcal{L}(\\tilde{r})',
      xLabelTex: '\\tilde{r}',
      yLabelTex: '\\mathcal{L}',
    });
    decoratePlot('cv-rmean', {
      titleTex:  '\\langle\\tilde{r}\\rangle_{\\text{free}}(\\tilde{t})',
      xLabelTex: '\\tilde{t}',
      yLabelTex: '\\langle\\tilde{r}\\rangle',
    });
    decoratePlot('cv-bead-3d', { titleTex: 'F(P_x, P_y)\\ \\text{(with chemotactic tilt)}' });
    decoratePlot('cv-bead-2d', {
      titleTex:  '(P_x, P_y)\\ \\text{phase plane}',
      xLabelTex: 'P_x', yLabelTex: 'P_y',
    });
    decoratePlot('cv-bead-1d', {
      titleTex:  'F(s)\\ \\text{along}\\ \\hat{\\nabla}\\mathcal{L}',
      xLabelTex: 's = P\\!\\cdot\\!\\hat{\\nabla}\\mathcal{L}',
      yLabelTex: 'F',
    });
  }

  if (window.katex) {
    decoratePlots();
  } else {
    window.addEventListener('katex-ready', decoratePlots);
    // Fallback: try after short delay in case event already fired.
    setTimeout(decoratePlots, 300);
  }
}
