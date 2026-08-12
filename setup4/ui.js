// setup4/ui.js
// UI assembly for Setup 4: KNOBS, Calculate button, time scrub slider, KPIs.
//
// Interaction paradigm: calculate-then-explore.
//   1. User sets parameters via sliders.
//   2. User presses "Calculate" → worker runs, saves frames.
//   3. Time slider becomes active; scrubbing updates heatmap + radial profile.
//   4. Any parameter change invalidates the run and re-enables Calculate.
//
// Intrinsic-units nondim scheme (plan §1):
//   ℓ_0 = a / (L_0 D_L)   [σ-independent length scale]
//   t_0 = ℓ_0² / D_L      [σ-independent time scale]
// Primary nondim groups: λ, ν, κ, μ, ϑ (cell-side), σ̃ = σ·ℓ_0² (density),
//   h̃ = h/ℓ_0 (2D-2D only). No Λ, no χ̃, no μ̃ in this scheme.
//
// Geometry primary knobs (plan §B1):
//   - R̃_dish: nondim dish radius (user-facing slider).
//   - σ_dim: dimensional cell density [cells/µm²] (dim slider in Geometry section).
//   - N = round(σ_dim · π · R_dim²) with R_dim = R̃·ℓ_0 — derived, shown as KPI.
//   - β = 1 by construction (N = round(σ·π·R_dim²)).
//
// KNOBS pattern: each parameter declared once; exposure: 'dim'|'nondim'|'both'.
// Moving a dim slider triggers recomputeFromDim() which updates all nondim sliders.

import { el, makeSlider, makeToggle, makeButtonRow, makeKpis, section,
         detailsSection, decoratePlot } from '../shared/dom.js';
import { dimToNondim, DIM_DEFAULTS } from './nondim.js';
import { drawDish, drawRadialProfile, drawRadialR, drawTimeSeries,
         drawBead3D, drawBead2D, drawBead1D, drawCsweep,
         drawAngularSpectrum } from './render.js';

// One-shot: clear any persisted gear-popover override of N_grid's default.
// Built-in default is 128 (see params below); previous sessions may have saved
// a different default (e.g. 1024) before the bind: refactor wired display ≡
// simulated. Runs once per browser per page-scope; subsequent gear edits
// persist normally.
(function migrateNgridDefault() {
  try {
    const flag = 'gl-motility:setup4:s4-Ngrid:migrated-default-2026-06';
    if (localStorage.getItem(flag)) return;
    const key = 'gl-motility:setup4:s4-Ngrid';
    const raw = localStorage.getItem(key);
    if (raw) {
      const s = JSON.parse(raw);
      if (s && s.default != null) { delete s.default; localStorage.setItem(key, JSON.stringify(s)); }
    }
    localStorage.setItem(flag, '1');
  } catch {}
})();

// ─── Nondim params (simulation state) ──────────────────────────────────────
const params = {
  geometry: '2d2d',   // '2d2d' | '2d3d'
  model:    'M1',
  // N is DERIVED from σ_dim · π · R_dim²  (β = 1 by construction; see recomputeFromDim).
  N:        500,
  N_grid:   128,
  // R̃_dish: nondim geometry knob. Both R̃ and R_dim are slider knobs;
  // dim is canonical (params.R_dish = dim.R_dim / ℓ_0 inside recomputeFromDim).
  // Default 30 = 300 µm / 10 µm = R_dim_default / ℓ_0_default.
  R_dish:   30,
  // Five independent cell-side nondim groups (σ-free, geometry-free).
  // See docs/physics/setup4_intrinsic_units_implementation_plan.md §1.
  lam:      1,        // λ = r_0 · t_0 · L_0
  nu:       1,        // ν = u² · t_0 / w
  kap:      0.075,    // κ = a χ / (D_L² · p_0)
  mu:       0.005,    // μ = μ_dim · t_0 / ℓ_0
  tht:      1e-4,     // ϑ = (w L_0 t_0 / u) · θ
  L_c_nd:   1,        // 𝓛_c = L_c / L_0
  // L-PDE nondim groups.
  n_L:      10,
  gamma_L:  0,        // Γ_L = γ_L · t_0
  sigma_tilde: 1.0,   // σ̃ = σ·ℓ_0² (TRUE cell density; set by dim linkage via σ_dim)
  h_tilde:  1,        // h̃ = h/ℓ_0 (set by dim linkage; only used in 2D-2D source prefactor)
  // M2 per-cell inhibitor nondim groups.
  beta_R:   0.075,    // β = Beta · t_0  (R_0 ≡ 1 implicit)
  gamma_R:  1e-4,     // Γ_R = Gamma_R · t_0
  L_r_nd:   0.01,     // 𝓛_r = L_r / L_0
  n_R:      10,
  n_Lr:     10,
  // M6.1 basal-adenosine density sensor (nondim-primary groups; catalog §7).
  //   λ_A threshold-shift coupling (λ_A = 0 ≡ M1); D = D_A/D_L; Γ_A decay.
  //   Defaults put the basal tone 𝓐_ss = σ̃/Γ_A at ≈ 2 for the default σ̃ ≈ 1,
  //   i.e. a threshold shift of order the bare threshold, with screening
  //   length ℓ_A = √(D/Γ_A) ≈ 1.4 (order the relay length). n_L is reused for
  //   the shifted relay gate (M6.1 adds no new Hill exponent).
  lam_A:    1.0,      // λ_A
  D_A_nd:   1.0,      // D = D_A / D_L
  gamma_A:  0.5,      // Γ_A = γ_A · t_0
  // M6.2 quorum-throttled production (nondim-primary groups; catalog §7b).
  //   β production-throttle emission rate (β = 0 ≡ M1 — NOT m = 0, which would
  //   halve production); D = D_Q/D_L; γ decay; m the throttle Hill exponent.
  //   Defaults: 𝓠_ss = βσ̃/γ = 3σ̃ — the throttle is firmly on at the default
  //   σ̃ ≈ 1 (m = 2 ⇒ production × 1/10) — and the screening length
  //   ℓ_Q = √(D/γ) = √(7/5) ≈ 1.2 is of order the relay length, so 𝓠 resolves
  //   LOCAL density bumps instead of acting as a global rescaling of α.
  //   The M6.2 page pairs these with a long t̃_max (MODEL_DEFAULTS): a heavily
  //   throttled front is slow.
  beta_Q:   15,       // β = b/(h D_L Q_0)   [2D–2D group; ×h̃ in 2D–3D]
  D_Q_nd:   7,        // D = D_Q / D_L
  gamma_Q:  5,        // γ = γ_Q · t_0
  m_Q:      2,        // m — Hill exponent of H⁻(𝓠;1;m)
  q_ic_ss:  true,     // seed 𝓠 at its steady tone at t̃=0 (avoids a 1/γ transient)
  // Target ("stick to the target"): a circle of radius R̃_target at the dish
  // centre standing for the pathogen (candida cluster / sterile injury). Cells
  // that reach it adhere and stop moving; they keep polarizing and emitting.
  // Applies to every full-swarm cue model (M1, M2, M6.1, M6.2).
  stick_target: true,
  R_target: 2.0,
  // Time-limited firing source.
  // s_fire = 1 → "every point in the firing disk emits as strongly as one
  // saturated cell" (source amount rescaled by σ̃ in worker, plan §A2).
  r_fire:   5.0,
  t_fire:   5.0,
  s_fire:   1.0,
  // Numerics. Per-geometry t̃_max stored so toggling restores the active value.
  // t_max = 50 covers typical R̃≈20 wave traversal at default σ̃ within a
  // tighter window than the conservative 300 originally used in calibration.
  dt:       0.01,
  t_max:        50,
  seed:     7,
  // 2D-3D z-grid. alpha_z = 1.5 (was 1.1 — silent 10× L bug B2; gives
  // z_max ≈ h_0·(α^N_z−1)/(α−1) ≈ 0.1·(1.5^16−1)/0.5 ≈ 23 ℓ_0 at defaults,
  // well above the recommended z_max > 5·√(R̃/σ̃) ≈ 32 at R̃=20, σ̃=0.1).
  N_z:      16,
  h_0:      0.1,
  alpha_z:  1.5,
};

// ─── Per-page (locked-model) default overrides ───────────────────────────────
// `params` above holds the defaults shared by every Setup-4 page. A page that
// locks a cue model can override a few of them here — one table, so a default
// is still declared exactly once per (model, knob) and never duplicated into
// the page HTML. Applied in buildUI() before the sliders are constructed.
const MODEL_DEFAULTS = {
  // M6.2: a strongly-throttling, fast-and-short-ranged quorum signal
  // (ℓ_Q = √(D/γ) = √(7/5) ≈ 1.2, i.e. of order the relay length, so 𝓠
  // resolves LOCAL density bumps rather than acting as a global rescaling),
  // with the tone 𝓠_ss = βσ̃/γ = 3 σ̃ well above 1 at the default σ̃ — the
  // throttle is firmly on. The long t̃_max = 500 is needed because a heavily
  // throttled front is slow: at t̃ = 50 it has barely left the ignition disk.
  // β, D, γ, m are declared once in `params` above (15 / 7 / 5 / 2 — the M6.2
  // defaults); only t̃_max is overridden per page, since it is a shared
  // numerics knob and the other Setup-4 pages keep 50.
  'M6.2': { t_max: 500 },
};

// ─── Dim state + linkage ────────────────────────────────────────────────────
// Intrinsic-units flow:
//   dim is the canonical source of truth (incl. dim.R_dim, dim.sigma).
//   Both Geometry-section sliders (R̃, σ̃) and the collapsed dim panel mirrors
//   (R_dim, σ_dim) are views — editing either back-computes through
//   recomputeFromDim() to keep them in sync.
//   ℓ_0 = a/(L_0·D_L) is a pure function of (a, L_0, D_L) — no σ dependence.
//   R̃ = R_dim / ℓ_0,  σ̃ = σ · ℓ_0²,  N = round(σ · π · R_dim²) (β=1 by const.)
const dim = Object.assign({}, DIM_DEFAULTS);
dim.R_dim = 300;  // µm — petri dish radius (canonical dim value)
// Dim canonicals for r_fire and t_max (mirror the nondim sim-config knobs).
// Defaults satisfy r̃=dim/ℓ₀ and t̃=dim/t₀ at the default ℓ₀=10 µm, t₀=1 s:
//   r_fire_dim = 50 µm →  r̃_fire = 5
//   t_max_dim  = 50 s  →  t̃_max  = 50
// Per-geometry t_max stored in seconds so each mode remembers its last value
// across geometry toggles (replaces the old params.t_max_2d2d/2d3d).
dim.r_fire_dim    = 50;
dim.t_max_dim     = 50;
dim.t_max_dim_2d2d = 50;
dim.t_max_dim_2d3d = 50;

let applyingDim = false;

function recomputeFromDim() {
  if (applyingDim) return;
  applyingDim = true;
  const nd = dimToNondim(dim, params.geometry);
  params.ell0 = nd.ell0;
  // R_dim is DERIVED: R_dim = R̃_dish · ℓ_0 (R̃ is the primary nondim slider).
  // R_dish (nondim) is now DERIVED from dim.R_dim (canonical); both are sliders.
  params.R_dish = dim.R_dim / nd.ell0;
  // Enforce r_fire ≤ R_dish on the dim side (so the constraint survives unit
  // changes — e.g. if user shrinks R_dim below r_fire_dim, clamp it down).
  if (dim.r_fire_dim > dim.R_dim) dim.r_fire_dim = dim.R_dim;
  params.r_fire = dim.r_fire_dim / nd.ell0;
  params.t_max  = dim.t_max_dim  / nd.t0;
  // N derived: N = round(σ_dim · π · R_dim²). β = 1 by construction.
  params.N = Math.max(1, Math.round(dim.sigma * Math.PI * dim.R_dim * dim.R_dim));
  sTracked && sTracked.setMinMax(0, Math.max(0, params.N - 1));
  // Push canonical values to BOTH Geometry sliders (R_dish, sigma_tilde) and
  // dim mirrors (R_dim, sigma). applyingDim guard prevents recursive onChange.
  sRdish     && sRdish.set(params.R_dish);
  sSigmaGeom && sSigmaGeom.set(nd.sigma_tilde);
  sHtilde    && sHtilde.set(nd.h_tilde);
  sRdimDim   && sRdimDim.set(dim.R_dim);
  sSigmaDim  && sSigmaDim.set(dim.sigma);
  // r_fire and t_max: push to both nondim sliders and their dim mirrors.
  sRFire     && sRFire.set(params.r_fire);
  sTmax      && sTmax.set(params.t_max);
  sRfireDim  && sRfireDim.set(dim.r_fire_dim);
  sTmaxDim   && sTmaxDim.set(dim.t_max_dim);
  sTime      && sTime.setMinMax(0, params.t_max);
  if (kpis) {
    kpis.set('N', String(params.N));
    kpis.set('sigma_dim', dim.sigma.toExponential(2));
    kpis.set('R_dim', dim.R_dim.toFixed(1) + ' µm');
    // c̃_Dieterle = (2/π)σ̃ for M1 + 2D-3D (plan §C1).
    if (params.model === 'M1' && params.geometry === '2d3d') {
      const cD = (2 / Math.PI) * nd.sigma_tilde;
      kpis.set('cD', cD.toFixed(3));
    } else {
      kpis.set('cD', '–');
    }
    // z_max KPI for 2D-3D (plan §B2 — silent 10× L bug check).
    if (params.geometry === '2d3d') {
      const zMax = params.h_0 * (Math.pow(params.alpha_z, params.N_z) - 1) /
                   (params.alpha_z - 1);
      kpis.set('zmax', zMax.toFixed(1));
    } else {
      kpis.set('zmax', '–');
    }
    // α = 1/h̃ = a/(h·D_L·L_0): the 2D-2D relay source prefactor (M6.1 speed
    // scale). Not an independent knob — it's fixed by the layer height h̃.
    if (params.geometry === '2d2d') {
      kpis.set('alpha', (1 / nd.h_tilde).toPrecision(3));
    } else {
      kpis.set('alpha', '–');
    }
    // M6.x auxiliary-field readouts: the uniform steady tone (𝓐_ss / 𝓠_ss) and
    // the screening length √(D/decay). In 2D–3D the tone is the screened
    // surface value ∝ σ̃/√(D·decay), and M6.2's source amplitude carries an
    // extra h̃ (catalog §7b). Both are exact properties of the linear auxiliary
    // PDE and were checked against the solver; neither is turned into a wave
    // speed here.
    const is3dK = (params.geometry === '2d3d');
    if (params.model === 'M6.1') {
      const gA = Math.max(params.gamma_A, 1e-12);
      const kA = is3dK ? 1 / Math.sqrt(Math.max(params.D_A_nd * gA, 1e-12)) : 1 / gA;
      kpis.set('tone', (kA * nd.sigma_tilde).toPrecision(3));
      kpis.set('ellQ', Math.sqrt(params.D_A_nd / gA).toPrecision(3));
    } else if (params.model === 'M6.2') {
      const gQ  = Math.max(params.gamma_Q, 1e-12);
      const amp = params.beta_Q * (is3dK ? nd.h_tilde : 1);
      const kQ  = is3dK ? amp / Math.sqrt(Math.max(params.D_Q_nd * gQ, 1e-12)) : amp / gQ;
      kpis.set('tone', (kQ * nd.sigma_tilde).toPrecision(3));
      kpis.set('ellQ', Math.sqrt(params.D_Q_nd / gQ).toPrecision(3));
    } else {
      kpis.set('tone', '–'); kpis.set('ellQ', '–');
    }
  }
  // Push linked nondim values to sliders (if sliders exist).
  sLam    && sLam.set(nd.lam);
  sNu     && sNu.set(nd.nu);
  sKap    && sKap.set(nd.kap);
  sMuNd   && sMuNd.set(nd.mu);
  sTht    && sTht.set(nd.tht);
  sLcNd   && sLcNd.set(nd.L_c_nd);
  sGammaL && sGammaL.set(nd.gamma_L);
  sBetaR  && sBetaR.set(nd.beta_R);
  sGammaR && sGammaR.set(nd.gamma_R);
  sLrNd   && sLrNd.set(nd.L_r_nd);
  // Update params with new nondim values.
  params.lam         = nd.lam;
  params.nu          = nd.nu;
  params.kap         = nd.kap;
  params.mu          = nd.mu;
  params.tht         = nd.tht;
  params.L_c_nd      = nd.L_c_nd;
  params.gamma_L     = nd.gamma_L;
  params.sigma_tilde = nd.sigma_tilde;
  params.h_tilde     = nd.h_tilde;
  params.beta_R      = nd.beta_R;
  params.gamma_R     = nd.gamma_R;
  params.L_r_nd      = nd.L_r_nd;
  if (kpis) kpis.set('sigma_tilde', nd.sigma_tilde.toExponential(2));
  applyingDim = false;
  markDirty();
}

// No-op: in the direct-formula scheme there is no calibration reference to
// update. Kept so onChange handlers that call recalibrate() remain harmless.
function recalibrate() { /* no-op */ }

// ─── Frame storage ───────────────────────────────────────────────────────────
/** @type {Array<{step, t, radialProfile, agentX, agentY, emitting, Lfield, Lmax}>} */
let frames = [];

let runDirty = true;
function markDirty() {
  runDirty = true; stopPlayback(); enableCalculate();
  // Keep the c(σ̃) panel's current-σ̃ marker live (M6.x pages only).
  updateSweepPlot();
}

// ─── Time scrub ─────────────────────────────────────────────────────────────
let currentFrameIdx = 0;
let currentTime = 0;
let trackedCellIdx = 0;

function getFrameAt(t) {
  if (frames.length === 0) return null;
  let lo = 0, hi = frames.length - 1;
  while (lo < hi) {
    const mid = (lo + hi + 1) >> 1;
    if (frames[mid].t <= t) lo = mid; else hi = mid - 1;
  }
  return frames[lo];
}

// ─── Playback (replay cached frames as a real-experiment timelapse) ──────────
// The video plays `playSpeedup`× faster than the real (dimensional) experiment:
//   video duration = (dimensional t_max in seconds) / playSpeedup.
// Default 120 = "30 min of dynamics in 15 s" (1800 s / 15 s). This is NOT an FPS
// knob — it is the dim-time → video-time compression factor the user asked for.
// We capture the run's t_max (nondim) and its dimensional duration (seconds) at
// completion so playback timing is unaffected by later slider edits.
let playRAF        = null;   // requestAnimationFrame handle (null ⇒ paused/stopped)
let playStartWall  = 0;      // performance.now() at (re)start, ms
let playStartTNd   = 0;      // nondim time the current play segment started from
let playSpeedup    = 120;    // dim-seconds of dynamics per second of video
let runTmaxNd      = 0;      // params.t_max of the completed run (nondim)
let runDurationSec = 0;      // dim.t_max_dim of the completed run (real seconds)

function curDurationSec() { return runDurationSec > 0 ? runDurationSec : dim.t_max_dim; }

function updatePlayInfo() {
  if (!playInfo) return;
  const d   = curDurationSec();
  const vid = playSpeedup > 0 ? d / playSpeedup : 0;
  playInfo.textContent = `${d.toFixed(0)} s dynamics → ${vid.toFixed(1)} s video`;
}

function stopPlayback() {
  if (playRAF !== null) { cancelAnimationFrame(playRAF); playRAF = null; }
  if (btnPlay) btnPlay.textContent = '▶ Play';
}

function playTick(now) {
  if (playRAF === null) return;
  const videoTotalSec = curDurationSec() / playSpeedup;       // wall-seconds for full run
  const elapsed = (now - playStartWall) / 1000;
  const tNd = (videoTotalSec > 0)
    ? playStartTNd + (elapsed / videoTotalSec) * runTmaxNd
    : runTmaxNd;
  if (tNd >= runTmaxNd) {
    sTime.set(runTmaxNd);   // snap to final frame
    stopPlayback();
    return;
  }
  sTime.set(tNd);           // fires sTime.onChange → frame pick + redraw
  playRAF = requestAnimationFrame(playTick);
}

function togglePlay() {
  if (playRAF !== null) { stopPlayback(); return; }
  if (frames.length === 0 || runTmaxNd <= 0) return;   // nothing computed yet
  // Restart from the beginning if we're parked at (or past) the end.
  if (currentTime >= runTmaxNd - 1e-9) sTime.set(0);
  playStartTNd  = currentTime;
  playStartWall = performance.now();
  if (btnPlay) btnPlay.textContent = '⏸ Pause';
  playRAF = requestAnimationFrame(playTick);
}

// ─── c(σ̃) density sweep (M6.1 only) ────────────────────────────────────────────
// The sweep is a SEPARATE multi-run action (its own button + worker): for each
// σ̃ it runs a full, frame-less sim and measures the front speed c̃. It uses a
// capped geometry (smaller dish/grid, robust firing) so ~10 runs stay fast; the
// measured relay speed is a propagation property, independent of those details.
// The panel plots the MEASURED points only — no analytical c(σ̃) overlay.
// Per-geometry sweep settings. 2D–3D uses a smaller grid but a longer window
// (the 3D solver is costly and 2D–3D relay speeds are lower, so slow fronts
// need more time to clear the firing disk).
// fire_K: the mass-targeted ignition constant. The worker sets
// s_fire = K/(σ̃·t̃_fire) per point, so the INJECTED 𝓛 mass — and hence the
// ignition-halo radius r̃_halo = r̃_fire√K — is the same at every density. With
// the old fixed s_fire the halo grew ∝√σ̃ and eventually flooded the dish, and
// measureWaveSpeed timed that diffusive spread instead of the relay (verified:
// a fully throttled M6.2 run, relay impossible, reported c ≈ 0.6). K is capped
// by the geometry: measurement needs 1.15·r̃_fire√K < 0.75·(0.9 R̃_dish).
// Sizing rule: the timing window runs from r_start = 1.15·r̃_fire√K out to
// 0.8·(0.9 R̃_dish), so the dish must be comfortably larger than the halo.
const SWEEP_GEOM_2D2D = { N_grid: 104, R_dish: 16, t_max: 40, r_fire: 1.5, t_fire: 5, s_fire: 3,
                          fire_K: 20, sigma_min: 0.1, sigma_max: 8, nPoints: 10 };
// 2D–3D dilutes the surface firing into the half-space and its relay speeds are
// lower, so it needs a MUCH stronger, briefer ignition kick and a longer window
// — and stays slower to compute, so use fewer points.
// 2D–3D speeds are LINEAR in σ̃ ((2/π)σ̃), so the top of the 2D–2D σ̃ range
// would fill this dish within one sampling interval and be unmeasurable; the
// range is capped at 3 instead of 8. (In 3D the halo dilutes into the
// half-space, so r̃_halo = (π s_fire σ̃ r̃_fire² t̃_fire)^{1/3} — see sim_core.)
const SWEEP_GEOM_2D3D = { N_grid: 52, R_dish: 12, t_max: 40, r_fire: 2, t_fire: 4, s_fire: 15,
                          fire_K: 20, sigma_min: 0.1, sigma_max: 3,
                          N_z: 9, h_0: 0.13, alpha_z: 1.5, nPoints: 8 };
/** Sweep σ̃ range for the currently selected geometry. */
function sweepRange() {
  const G = (params.geometry === '2d3d') ? SWEEP_GEOM_2D3D : SWEEP_GEOM_2D2D;
  return { min: G.sigma_min, max: G.sigma_max };
}
let sweepWorker = null;
let sweepBusy   = false;
let sweepPoints = [];        // [{sigma, c}] accumulated from the running sweep
let btnSweep, sweepInfo;

function sweepSigmas(nPoints) {
  const { min, max } = sweepRange();
  const out = [];
  for (let i = 0; i < nPoints; i++) {
    out.push(min * Math.pow(max / min, i / (nPoints - 1)));
  }
  return out;
}

// Redraw the c(σ̃) panel (measured points only — no analytical curve; see
// drawCsweep in render.js).
function updateSweepPlot() {
  if (!document.getElementById('cv-csweep')) return;   // only on the M6.1 page
  drawCsweep('cv-csweep', {
    points: sweepPoints,
    sigma_min: sweepRange().min, sigma_max: sweepRange().max,
    sigma_current: params.sigma_tilde,
  });
}

function startSweep() {
  if (!document.getElementById('cv-csweep')) return;
  if (sweepBusy) return;
  if (sweepWorker) { try { sweepWorker.terminate(); } catch {} sweepWorker = null; }
  sweepBusy = true;
  sweepPoints = [];
  updateSweepPlot();

  // Base params: current cell/cue nondim groups + M6.1 knobs, with capped
  // (geometry-appropriate) sweep geometry and robust firing. σ̃ / N are set
  // per-point in the worker. Uses the CURRENTLY SELECTED geometry so the
  // sweep reflects 2D–2D (c ∝ √σ̃) vs 2D–3D (c ∝ σ̃) physics.
  const is3d = (params.geometry === '2d3d');
  const G = is3d ? SWEEP_GEOM_2D3D : SWEEP_GEOM_2D2D;
  const baseP = {
    geometry: params.geometry, model: params.model,
    N_grid: G.N_grid, R_dish: G.R_dish, t_max: G.t_max,
    r_fire: G.r_fire, t_fire: G.t_fire, s_fire: G.s_fire, fire_K: G.fire_K,
    lam: params.lam, nu: params.nu, kap: params.kap, mu: params.mu, tht: params.tht,
    L_c_nd: params.L_c_nd, n_L: params.n_L, gamma_L: params.gamma_L,
    h_tilde: params.h_tilde,
    lam_A: params.lam_A, D_A_nd: params.D_A_nd, gamma_A: params.gamma_A,
    beta_Q: params.beta_Q, D_Q_nd: params.D_Q_nd, gamma_Q: params.gamma_Q,
    m_Q: params.m_Q, q_ic_ss: params.q_ic_ss,
    // The target is part of the cell dynamics, so the sweep must see it too —
    // otherwise the swept c̃(σ̃) would describe a different swarm than the run.
    stick_target: params.stick_target, R_target: params.R_target,
    dt: params.dt, seed: params.seed,
    ...(is3d ? { N_z: G.N_z, h_0: G.h_0, alpha_z: G.alpha_z } : {}),
  };
  const sigmas = sweepSigmas(G.nPoints);

  sweepWorker = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
  sweepWorker.onmessage = (e) => {
    const m = e.data;
    if (m.type === 'sweep-progress') {
      sweepPoints.push({ sigma: m.sigma, c: m.c, status: m.status });
      if (sweepInfo) sweepInfo.textContent = `sweeping… ${m.i}/${m.n}  (σ̃=${m.sigma.toFixed(2)}, c̃=${m.c.toFixed(3)})`;
      updateSweepPlot();
    } else if (m.type === 'sweep-done') {
      sweepPoints = m.points;
      sweepBusy = false;
      if (sweepInfo) {
        const nUnres = m.points.filter(q => q.status === 'unresolved').length;
        const nUnmeas = m.points.filter(q => q.status === 'unmeasurable').length;
        const tag = (params.model === 'M6.1')
          ? `λ_A=${params.lam_A.toFixed(2)}`
          : `β=${params.beta_Q.toFixed(2)}, m=${params.m_Q}`;
        sweepInfo.textContent = `done — ${m.points.length} points at ${tag}` +
          (nUnres  ? ` · ${nUnres} too slow to resolve` : '') +
          (nUnmeas ? ` · ${nUnmeas} unmeasurable (halo)` : '');
      }
      if (btnSweep) { btnSweep.disabled = false; btnSweep.textContent = 'Sweep c(σ̃)'; }
      updateSweepPlot();
      try { sweepWorker.terminate(); } catch {} sweepWorker = null;
    } else if (m.type === 'error') {
      sweepBusy = false;
      if (sweepInfo) sweepInfo.textContent = 'error: ' + m.message;
      if (btnSweep) { btnSweep.disabled = false; btnSweep.textContent = 'Sweep c(σ̃)'; }
    }
  };
  sweepWorker.onerror = (e) => {
    sweepBusy = false;
    if (sweepInfo) sweepInfo.textContent = 'worker error: ' + e.message;
    if (btnSweep) { btnSweep.disabled = false; btnSweep.textContent = 'Sweep c(σ̃)'; }
  };
  sweepWorker.postMessage({ type: 'sweep', params: baseP, sigmas });
  if (btnSweep) { btnSweep.disabled = true; btnSweep.textContent = 'Sweeping…'; }
  if (sweepInfo) sweepInfo.textContent =
    `sweeping σ̃ ∈ [${sweepRange().min}, ${sweepRange().max}] (${params.geometry}${is3d ? ', slower — 3D solver' : ''}) …`;
}

// ─── Worker ─────────────────────────────────────────────────────────────────
let worker = null;
let workerBusy = false;

function startWorker() {
  if (worker) {
    try { worker.terminate(); } catch {}
    worker = null;
  }
  worker = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
  workerBusy = true;
  stopPlayback();   // a fresh run invalidates the cached frames being played
  frames = [];
  currentFrameIdx = 0;

  // Build params for worker (intrinsic-units groups only; no Lambda/chi_nd/mu_nd).
  const wp = {
    geometry:     params.geometry,
    model:        params.model,
    N:            params.N,
    N_grid:       params.N_grid,
    R_dish:       params.R_dish,
    lam:          params.lam,
    nu:           params.nu,
    kap:          params.kap,
    mu:           params.mu,
    tht:          params.tht,
    L_c_nd:       params.L_c_nd,
    n_L:          params.n_L,
    gamma_L:      params.gamma_L,
    sigma_tilde:  params.sigma_tilde,
    h_tilde:      params.h_tilde,
    beta_R:       params.beta_R,
    gamma_R:      params.gamma_R,
    L_r_nd:       params.L_r_nd,
    n_R:          params.n_R,
    n_Lr:         params.n_Lr,
    lam_A:        params.lam_A,
    D_A_nd:       params.D_A_nd,
    gamma_A:      params.gamma_A,
    beta_Q:       params.beta_Q,
    D_Q_nd:       params.D_Q_nd,
    gamma_Q:      params.gamma_Q,
    m_Q:          params.m_Q,
    q_ic_ss:      params.q_ic_ss,
    stick_target: params.stick_target,
    R_target:     params.R_target,
    r_fire:       params.r_fire,
    t_fire:       params.t_fire,
    s_fire:       params.s_fire,
    dt:           params.dt,
    t_max:        params.t_max,
    seed:         params.seed,
    N_z:          params.N_z,
    h_0:          params.h_0,
    alpha_z:      params.alpha_z,
  };

  worker.postMessage({ type: 'run', params: wp });

  worker.onmessage = function (e) {
    const msg = e.data;

    if (msg.type === 'frame') {
      // Precompute per-frame swarm-averaged diagnostics (cheap O(N)).
      //   meanAbsVr = ⟨|μ P·r̂|⟩  (radial speed, cells at r ≈ 0 excluded)
      //   meanAbsP  = ⟨|P|⟩
      const Nc = msg.agentX.length;
      const muV = params.mu;   // μ (intrinsic-units motility)
      let sumVr = 0, sumP = 0, cntVr = 0;
      const EPS_R = 1e-6;
      for (let i = 0; i < Nc; i++) {
        const xi = msg.agentX[i], yi = msg.agentY[i];
        const pxi = msg.agentPx[i], pyi = msg.agentPy[i];
        sumP += Math.hypot(pxi, pyi);
        const r = Math.hypot(xi, yi);
        if (r > EPS_R) {
          sumVr += Math.abs(muV * (pxi * xi + pyi * yi) / r);
          cntVr++;
        }
      }
      const meanAbsVr = cntVr > 0 ? sumVr / cntVr : 0;
      const meanAbsP  = Nc   > 0 ? sumP  / Nc   : 0;

      // Channelisation (streaming) order parameter — see
      // docs/physics/setup4_swarm3d.md §10. Azimuthal Fourier modes of the
      // cell angles, shot-noise corrected so the value is density-independent:
      //   c_m = (1/N)Σ e^{imθ},  Ψ_m = √(max(0,(N|c_m|²−1)/(N−1))) ∈ [0,1]
      // Ψ = max over m ≥ 2 (m = 1 is a bulk off-centre drift, not spokes);
      // the dominant mode m* is the channel count.
      const M_MAX = 16;
      const rCore = 0.1 * params.R_dish;   // exclude the core, where θ is noise
      const chanPsi = new Float32Array(M_MAX + 1);
      let chanPsiMax = 0, chanMstar = 0, chanNoise = 0;
      {
        const cRe = new Float64Array(M_MAX + 1), cIm = new Float64Array(M_MAX + 1);
        let Nang = 0;
        for (let i = 0; i < Nc; i++) {
          const xi = msg.agentX[i], yi = msg.agentY[i];
          if (Math.hypot(xi, yi) < rCore) continue;
          const th = Math.atan2(yi, xi);
          Nang++;
          for (let m = 1; m <= M_MAX; m++) {
            cRe[m] += Math.cos(m * th);
            cIm[m] += Math.sin(m * th);
          }
        }
        if (Nang > 1) {
          chanNoise = Math.sqrt(2 / (Nang - 1));   // 5% shot-noise level
          for (let m = 1; m <= M_MAX; m++) {
            const P = (cRe[m] * cRe[m] + cIm[m] * cIm[m]) / (Nang * Nang);
            const excess = (Nang * P - 1) / (Nang - 1);
            chanPsi[m] = Math.sqrt(Math.max(0, excess));
            if (m >= 2 && chanPsi[m] > chanPsiMax) { chanPsiMax = chanPsi[m]; chanMstar = m; }
          }
        }
      }

      frames.push({
        step: msg.step, t: msg.t,
        radialProfile: msg.radialProfile,
        agentX: msg.agentX, agentY: msg.agentY,
        Px:     msg.agentPx, Py: msg.agentPy,
        Gx:     msg.agentGx, Gy: msg.agentGy,
        emitting: msg.emitting,
        stuck:  msg.stuck,     // 1 = engaged with the target circle
        agentR: msg.agentR,
        Lfield: msg.Lfield,
        Lmax:   msg.Lmax,
        Qfield: msg.Qfield,   // M6.2 only (second petri dish); undefined otherwise
        Qmax:   msg.Qmax,
        meanAbsVr, meanAbsP,
        chanPsi, chanPsiMax, chanMstar, chanNoise,
      });
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
      updateWaveSpeed();
      // Capture the run's timing so playback is independent of later slider edits.
      runTmaxNd      = params.t_max;
      runDurationSec = dim.t_max_dim;
      updatePlayInfo();
      // Re-enable the button so the run can be replayed/recomputed; otherwise
      // it stays stuck on the disabled "Running…" label until a slider changes.
      enableCalculate();
    }

    if (msg.type === 'error') {
      workerBusy = false;
      kpis && kpis.set('status', 'Error: ' + msg.message);
      console.error('Setup4 worker error:', msg.message);
      enableCalculate();
    }
  };

  worker.onerror = function (e) {
    workerBusy = false;
    kpis && kpis.set('status', 'Worker error: ' + e.message);
    console.error('Setup4 worker onerror:', e);
    enableCalculate();
  };

  kpis && kpis.set('status', 'Running…');
  disableCalculate();
}

// ─── KPI: wave-speed estimate ────────────────────────────────────────────────
// Tracks the FRONT (outermost radius where the azimuthally-averaged 𝓛 crosses
// 1), NOT the peak. With Γ_L = 0 the relay is a runaway filling front whose
// peak sits at r≈0 forever — peak-tracking gives c≈0 regardless of the actual
// speed. The front-crossing radius grows at the true wave speed, so this
// reflects α (=1/h̃), σ̃ and λ_A. Least-squares slope over the frames whose
// front lies strictly between the firing disk and the wall.
function updateWaveSpeed() {
  if (frames.length < 4) return;
  const dx   = (2 * params.R_dish) / (params.N_grid - 1);
  const rMin = Math.max((params.r_fire || 0) * 1.1, dx * 2);
  const rMax = 0.9 * params.R_dish;
  const frontR = (prof) => {
    let rf = 0;
    for (let k = 0; k < prof.length; k++) if (prof[k] >= 1) rf = k * dx;
    return rf;
  };

  const ts = [], rs = [];
  for (const f of frames) {
    if (!f.radialProfile) continue;
    const rf = frontR(f.radialProfile);
    if (rf > rMin && rf < rMax) { ts.push(f.t); rs.push(rf); }
  }
  if (ts.length < 2) { kpis && kpis.set('ceff', '0.000'); return; }

  let st = 0, sr = 0, stt = 0, str = 0;
  const n = ts.length;
  for (let i = 0; i < n; i++) { st += ts[i]; sr += rs[i]; stt += ts[i] * ts[i]; str += ts[i] * rs[i]; }
  const denom = n * stt - st * st;
  const c = Math.abs(denom) < 1e-12 ? 0 : (n * str - st * sr) / denom;
  kpis && kpis.set('ceff', (c > 0 ? c : 0).toFixed(3));
}

// ─── Redraw ──────────────────────────────────────────────────────────────────
function redraw() {
  const f = (frames.length > 0) ? frames[currentFrameIdx] : null;
  if (!f) return;

  drawDish('cv-dish', f, {
    N_grid: params.N_grid, R_dish: params.R_dish, t: f.t,
    trackedCellIdx,
    stick_target: params.stick_target, R_target: params.R_target,
    // M6.2: frame.agentR carries 𝓠_i, so the renderer can dim each emitting
    // cell by its throttle factor H⁻(𝓠;1;m) — "firing, but weakly".
    model: params.model, m_Q: params.m_Q,
  });
  // Second petri dish (M6.2 page only): the same cells over the quorum field
  // 𝓠 instead of 𝓛, in a green palette so the two panels never read as the
  // same quantity. Identical cell markers, so the panels can be compared
  // side by side: bright cells on the left are firing, and their brightness
  // there is already the throttle factor set by 𝓠 shown on the right.
  drawDish('cv-dishQ', f, {
    N_grid: params.N_grid, R_dish: params.R_dish, t: f.t,
    trackedCellIdx,
    stick_target: params.stick_target, R_target: params.R_target,
    model: params.model, m_Q: params.m_Q,
    fieldKey: 'Qfield', palette: 'greens', haloScale: 0.45,
  });
  drawRadialProfile('cv-profile', f, {
    N_grid: params.N_grid, R_dish: params.R_dish,
    model: params.model, geometry: params.geometry,
    N_cells: params.N, sigma_tilde: params.sigma_tilde,
    L_r_nd: params.L_r_nd,
  });
  drawRadialR('cv-rmean', f, {
    R_dish: params.R_dish, model: params.model, geometry: params.geometry,
    // Guide line: M6.1 → basal tone 𝓐_ss = σ̃/Γ_A; M6.2 → 𝓠_ss = βσ̃/γ
    // (β·h̃ in 2D–3D, where the tone is the screened σ̃·β/√(Dγ)).
    sigma_tilde: params.sigma_tilde, gamma_A: params.gamma_A,
    beta_Q: params.beta_Q, gamma_Q: params.gamma_Q, D_Q_nd: params.D_Q_nd,
    h_tilde: params.h_tilde,
  });
  // Channelisation diagnostics (pages that include the canvases).
  drawTimeSeries('cv-chan', frames, 'chanPsiMax',
                 { color: '#2b6cb0', currentT: f.t, yMinAuto: 0.05 });
  drawAngularSpectrum('cv-chanspec', f);
  drawTimeSeries('cv-vr',   frames, 'meanAbsVr', { color: '#a06030', currentT: f.t });
  drawTimeSeries('cv-pmag', frames, 'meanAbsP',  { color: '#306060', currentT: f.t });

  // Bead-plot params use the new intrinsic-units cell-side groups.
  // lam→λ, nu→ν, kap→κ (plan §C3: Lambda→lam, lam→nu, chi→kap).
  const beadParams = {
    lam: params.lam, nu: params.nu, kap: params.kap, L_c: params.L_c_nd,
    R_dish: params.R_dish, N_grid: params.N_grid,
  };
  drawBead3D('cv-bead-3d', frames, currentFrameIdx, trackedCellIdx, beadParams);
  drawBead2D('cv-bead-2d', frames, currentFrameIdx, trackedCellIdx, beadParams);
  drawBead1D('cv-bead-1d', frames, currentFrameIdx, trackedCellIdx, beadParams);

  if (kpis) {
    kpis.set('t', f.t.toFixed(2));
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
    if (params.stick_target && f.stuck) {
      let nS = 0;
      for (let k = 0; k < f.stuck.length; k++) nS += f.stuck[k];
      kpis.set('stuck', `${nS} / ${f.stuck.length}`);
    } else {
      kpis.set('stuck', '–');
    }
  }
}

// ─── Slider references (populated during buildUI) ─────────────────────────────
// New scheme: sLam, sNu, sKap replace old sLambda, sChiNd (removed entirely).
// sMuNd and sTht are kept but now drive μ and ϑ (not the old μ̃ and ϑ).
let sLam, sNu, sKap, sMuNd, sLcNd, sTht;
let sGammaL, sBetaR, sGammaR, sLrNd;
let sLamA, sDA, sGammaA;
let sBetaQ, sDQ, sGammaQ, smQ, qIcToggle;
let sTime, sTmax, sH, sHtilde, sA, sTracked;
// Playback controls (Play button + timelapse-speed slider + info readout).
let btnPlay, sPlaySpeed, playInfo;
// Dim mirror sliders for Geometry knobs (declared at file scope so recompute
// can push values to them).
let sRdimDim, sSigmaDim;
// Nondim Geometry sliders (declared at file scope so recompute pushes work too).
let sRdish, sSigmaGeom;
let sTFire, sSFire, sRFire;
// Target (sticking boundary) controls.
let sRtarget, stickToggle;
// Dim mirrors for r_fire and t_max (canonical-side knobs).
let sRfireDim, sTmaxDim;
let m2SectionEl, m61SectionEl, m62SectionEl;
let kpis;
let btnCalc, btnReset;

function enableCalculate()  { if (btnCalc) { btnCalc.disabled = false; btnCalc.textContent = 'Calculate'; } }
function disableCalculate() { if (btnCalc) { btnCalc.disabled = true;  btnCalc.textContent = 'Running…'; } }

function resetRun() {
  if (worker) {
    try { worker.terminate(); } catch {}
    worker = null;
  }
  if (sweepWorker) {
    try { sweepWorker.terminate(); } catch {}
    sweepWorker = null;
    sweepBusy = false;
    if (btnSweep) { btnSweep.disabled = false; btnSweep.textContent = 'Sweep c(σ̃)'; }
  }
  workerBusy = false;
  stopPlayback();
  runTmaxNd = 0;
  runDurationSec = 0;
  frames = [];
  currentFrameIdx = 0;
  runDirty = true;
  if (sTime) sTime.set(0);
  updatePlayInfo();
  if (kpis) {
    kpis.set('status', 'Press Calculate');
    kpis.set('t',      '–');
    kpis.set('Lmax',   '–');
    kpis.set('nemit',  '–');
    kpis.set('stuck',  '–');
    kpis.set('ceff',   '–');
    kpis.set('cD',     '–');
    kpis.set('zmax',   '–');
  }
  ['cv-dish', 'cv-profile', 'cv-rmean', 'cv-vr', 'cv-pmag'].forEach(id => {
    const c = document.getElementById(id);
    if (c) { const ctx = c.getContext('2d'); ctx.clearRect(0, 0, c.width, c.height); }
  });
  enableCalculate();
}

// ─── Build UI ────────────────────────────────────────────────────────────────
export function buildUI(containerId, opts = {}) {
  const container = document.getElementById(containerId);
  if (!container) { console.error('setup4/ui: container not found:', containerId); return; }

  // Lock the cue model when the host page specifies one (M1 page vs M2 page).
  // When locked, the M1/M2 toggle is omitted from the controls panel.
  const lockedModel = opts.model || null;
  if (lockedModel) params.model = lockedModel;

  // Per-page default overrides. Applied BEFORE any slider is built, so the
  // sliders pick these up as their initial (and gear-popover "reset") values.
  // t̃_max is dim-canonical, so it has to be pushed through the dim mirror as
  // well or recomputeFromDim() would immediately overwrite it (the
  // display ≡ simulated invariant).
  const md = MODEL_DEFAULTS[lockedModel];
  if (md) {
    for (const [k, v] of Object.entries(md)) {
      if (k === 't_max') {
        const nd0 = dimToNondim(dim, params.geometry);
        params.t_max        = v;
        dim.t_max_dim       = v * nd0.t0;
        dim.t_max_dim_2d2d  = dim.t_max_dim;
        dim.t_max_dim_2d3d  = dim.t_max_dim;
      } else {
        params[k] = v;
      }
    }
  }

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

  // ── c(σ̃) density-sweep button (M6.x only) ──
  // Separate multi-run action: measures wave speed across a range of densities
  // to reveal the density-independent (M6.1) / non-monotonic (M6.2) regime.
  if (params.model === 'M6.1' || params.model === 'M6.2') {
    btnSweep = el('button', { type: 'button',
      style: { width: '100%', fontWeight: '600', marginBottom: '4px' } }, 'Sweep c(σ̃)');
    btnSweep.addEventListener('click', startSweep);
    sweepInfo = el('div', { class: 'linked-readout',
      style: { fontSize: '0.8em', opacity: '0.8', marginBottom: '8px' } },
      (params.model === 'M6.1')
        ? `run a c(σ̃) sweep at the current λ_A = ${params.lam_A.toFixed(2)}`
        : `run a c(σ̃) sweep at the current β = ${params.beta_Q.toFixed(2)}, m = ${params.m_Q}`);
    container.appendChild(btnSweep);
    container.appendChild(sweepInfo);
  }

  // ── KPIs ──
  // R_dim: derived from R̃·ℓ_0 (plan §B1). N: derived from σ_dim·π·R_dim².
  // cD: c̃_Dieterle = (2/π)σ̃ for M1+2D-3D. zmax: for B2 z-grid check.
  kpis = makeKpis([
    { id: 'status',      label: 'Status',        init: 'Press Calculate' },
    { id: 't',           label: 't̃',        init: '–' },
    { id: 'Lmax',        label: 'max ℓ',     init: '–' },
    { id: 'nemit',       label: 'emitting',       init: '–' },
    { id: 'stuck',       label: 'on target',      init: '–' },
    { id: 'ceff',        label: 'c_eff',          init: '–' },
    { id: 'cD',          label: 'c̃_Dieterle', init: '–' },
    { id: 'sigma_tilde', label: 'σ̃',   init: '–' },
    { id: 'sigma_dim',   label: 'σ',         init: '–' },
    { id: 'R_dim',       label: 'R_dim',          init: '–' },
    { id: 'N',           label: 'N',              init: '–' },
    { id: 'zmax',        label: 'z_max',          init: '–' },
    { id: 'alpha',       label: 'α = 1/h̃',   init: '–' },
    // M6.x density-sensor readouts (blank for other models). These are
    // properties of the LINEAR auxiliary field only — its uniform steady level
    // and its screening length — both verified numerically. No speed estimate
    // is derived from them.
    { id: 'tone',        label: 'tone_ss',        init: '–' },
    { id: 'ellQ',        label: 'ℓ_screen',   init: '–' },
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

  // ── Playback: Play/Pause button + timelapse-speed slider ──
  // Speed is the dim-time → video-time compression factor (NOT FPS): the video
  // plays `playSpeedup`× faster than real experiment time. Default 120 shows
  // "30 min of dynamics in 15 s" (1800 s / 15 s). See playback block above.
  btnPlay = el('button', { type: 'button',
    style: { fontWeight: '700', whiteSpace: 'nowrap' } }, '▶ Play');
  btnPlay.addEventListener('click', togglePlay);

  sPlaySpeed = makeSlider({
    id: 's4-playspeed', symbol: 'speed', value: 120, min: 1, max: 2000, step: 1,
    units: 'real-time', fmt: v => `${v.toFixed(0)}×`,
  });
  sPlaySpeed.el.style.flex = '1';
  playSpeedup = sPlaySpeed.value;
  sPlaySpeed.onChange(v => {
    playSpeedup = v;
    // Re-anchor an in-flight playback at the current position so changing speed
    // mid-play rescales the rest of the video instead of jumping.
    if (playRAF !== null) { playStartTNd = currentTime; playStartWall = performance.now(); }
    updatePlayInfo();
  });

  playInfo = el('div', { class: 'linked-readout',
    style: { marginTop: '2px', opacity: '0.8' } }, '');
  updatePlayInfo();

  const playRow = el('div', {
    style: { display: 'flex', gap: '8px', alignItems: 'center', marginBottom: '6px' } });
  playRow.appendChild(btnPlay);
  playRow.appendChild(sPlaySpeed.el);

  container.appendChild(section('Time scrub', [playRow, playInfo, sTime.el]));

  // ── Geometry toggle ──
  const geoToggle = makeToggle({
    label: 'Geometry',
    options: [{ id: '2d2d', label: '2D–2D' }, { id: '2d3d', label: '2D–3D' }],
    value: params.geometry,
    onChange(v) {
      // Per-geometry t_max preference is stored on the dim side (seconds).
      if (params.geometry === '2d2d') dim.t_max_dim_2d2d = dim.t_max_dim;
      else                            dim.t_max_dim_2d3d = dim.t_max_dim;
      params.geometry = v;
      dim.t_max_dim = (v === '2d2d') ? dim.t_max_dim_2d2d : dim.t_max_dim_2d3d;
      // h slider visible only in 2D-2D (h̃ = h/ℓ_0 is the 2D-2D source prefactor).
      sH.el.style.display     = (v === '2d2d') ? '' : 'none';
      sHtilde.el.style.display = (v === '2d2d') ? '' : 'none';
      nzRow.style.display     = (v === '2d3d') ? '' : 'none';
      h0Row.style.display     = (v === '2d3d') ? '' : 'none';
      alphaRow.style.display  = (v === '2d3d') ? '' : 'none';
      recomputeFromDim();
      markDirty();
    },
  });

  // ── Geometry section (nondim) ──
  // Both R̃_dish and σ̃ are slider knobs; dim mirrors (R_dim, σ_dim) live in
  // the collapsed Dimensional inputs panel. dim is canonical — editing either
  // side back-computes through recomputeFromDim().
  sRdish = makeSlider({
    id: 's4-Rdish', symbol: '\\tilde{R}_{\\text{dish}}',
    bind: [params, 'R_dish'], min: 0.5, max: 500, step: null, log: true,
  });
  sRdish.onChange(v => {
    if (applyingDim) return;
    // Back-compute canonical dim.R_dim from R_tilde and current ell_0.
    const ell0 = dim.a / (dim.L_0 * dim.D_L);
    dim.R_dim = v * ell0;
    recomputeFromDim();
  });

  // sigma_tilde (nondim) slider in Geometry; dim mirror lives in the dim panel.
  sSigmaGeom = makeSlider({
    id: 's4-sigma-geom', symbol: '\\tilde{\\sigma}',
    bind: [params, 'sigma_tilde'], min: 1e-4, max: 100, step: null, log: true,
  });
  sSigmaGeom.onChange(v => {
    if (applyingDim) return;
    // Back-compute canonical dim.sigma from sigma_tilde and current ell_0.
    const ell0 = dim.a / (dim.L_0 * dim.D_L);
    dim.sigma = v / (ell0 * ell0);
    recomputeFromDim();
  });

  const sNgrid = makeSlider({
    id: 's4-Ngrid', symbol: 'N_{\\text{grid}}', bind: [params, 'N_grid'],
    min: 32, max: 1024, step: null, log: true,
    transform: Math.round, fmt: v => String(v),
  });
  sNgrid.onChange(() => markDirty());

  // Nondim layer height h̃ = h/ℓ_0 (2D-2D only). This is the knob behind the
  // relay source prefactor α = 1/h̃ (a/(h·D_L·L_0)) — the overall speed scale
  // in M1/M6.1. Mirror of the dimensional h slider (in Dimensional inputs);
  // editing either back-computes through recomputeFromDim (like R̃/σ̃).
  sHtilde = makeSlider({
    id: 's4-htilde', symbol: '\\tilde{h}\\ (\\alpha = 1/\\tilde{h})',
    bind: [params, 'h_tilde'], min: 0.01, max: 100, step: null, log: true,
  });
  sHtilde.onChange(v => {
    if (applyingDim) return;
    const ell0 = dim.a / (dim.L_0 * dim.D_L);
    dim.h = v * ell0;   // back-compute canonical dim.h
    recomputeFromDim();
    markDirty();
  });

  container.appendChild(section('Geometry', [geoToggle.el, sRdish.el, sSigmaGeom.el, sHtilde.el, sNgrid.el]));

  // ── Target (sticking boundary) ──
  // Nondim-only knob (exposure: 'nondim'): R̃_target is a geometric radius in
  // units of ℓ_0, read straight by the agent stepper. Cells that reach the
  // circle adhere to it permanently — the pathogen the swarm converges on.
  sRtarget = makeSlider({
    id: 's4-Rtarget', symbol: '\\tilde{R}_{\\text{target}}',
    bind: [params, 'R_target'], min: 0.1, max: 30, step: 0.1,
  });
  sRtarget.onChange(() => markDirty());

  stickToggle = makeToggle({
    label: 'Stick to the target',
    options: [{ id: 'on', label: 'on' }, { id: 'off', label: 'off' }],
    value: params.stick_target ? 'on' : 'off',
    onChange(v) {
      params.stick_target = (v === 'on');
      sRtarget.el.style.display = params.stick_target ? '' : 'none';
      markDirty();
    },
  });
  sRtarget.el.style.display = params.stick_target ? '' : 'none';
  container.appendChild(section('Target', [stickToggle.el, sRtarget.el]));

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

  // ── Cue model toggle + M1/M2 nondim sections ──
  const cueModelToggle = makeToggle({
    label: 'Cue model',
    options: [
      { id: 'M1', label: 'M1 (relay only)' },
      { id: 'M2', label: 'M2 (+ per-cell R)' },
      { id: 'M6.1', label: 'M6.1 (basal adenosine)' },
      { id: 'M6.2', label: 'M6.2 (quorum-throttled prod.)' },
    ],
    value: params.model,
    onChange(v) {
      params.model = v;
      updateModelVisibility();
      // Refresh Dieterle KPI which is model-dependent.
      recomputeFromDim();
      markDirty();
    },
  });

  sGammaL = makeSlider({ id: 's4-gammaL', symbol: '\\tilde{\\Gamma}_L', bind: [params, 'gamma_L'], min: 0, max: 5, step: 0.01 });
  sGammaL.onChange(() => { recalibrate(); markDirty(); });

  const snL = makeSlider({ id: 's4-nL', symbol: 'n_L', bind: [params, 'n_L'], min: 1, max: 50, step: 1, transform: Math.round, fmt: v => String(v) });
  snL.onChange(() => markDirty());

  const cueModelChildren = lockedModel
    ? [sGammaL.el, snL.el]
    : [cueModelToggle.el, sGammaL.el, snL.el];
  container.appendChild(section('Cue model', cueModelChildren));

  // M2-specific nondim sliders.
  sBetaR = makeSlider({ id: 's4-betaR', symbol: '\\tilde{\\beta}', bind: [params, 'beta_R'], min: 0, max: 1, step: 0.001 });
  sBetaR.onChange(() => { recalibrate(); markDirty(); });

  sGammaR = makeSlider({ id: 's4-gammaR', symbol: '\\tilde{\\gamma}', bind: [params, 'gamma_R'], min: 1e-4, max: 1, step: null, log: true });
  sGammaR.onChange(() => { recalibrate(); markDirty(); });

  sLrNd = makeSlider({ id: 's4-Lr', symbol: '\\tilde{L}_r', bind: [params, 'L_r_nd'], min: 1e-4, max: 10, step: null, log: true });
  sLrNd.onChange(() => { recalibrate(); markDirty(); });

  const snR = makeSlider({ id: 's4-nR', symbol: 'n_R', bind: [params, 'n_R'], min: 1, max: 50, step: 1, transform: Math.round, fmt: v => String(v) });
  snR.onChange(() => markDirty());

  const snLr = makeSlider({ id: 's4-nLr', symbol: 'n_{Lr}', bind: [params, 'n_Lr'], min: 1, max: 50, step: 1, transform: Math.round, fmt: v => String(v) });
  snLr.onChange(() => markDirty());

  m2SectionEl = section('Cue model (M2)', [sBetaR.el, sGammaR.el, sLrNd.el, snR.el, snLr.el]);
  container.appendChild(m2SectionEl);

  // M6.1-specific nondim sliders (basal-adenosine density sensor).
  //   λ_A threshold-shift coupling (λ_A = 0 recovers M1 exactly).
  //   D = D_A/D_L adenosine diffusivity ratio.
  //   Γ_A adenosine decay rate; screening length ℓ_A = √(D/Γ_A).
  // λ_A linear (not log) so λ_A = 0 — the exact M1-recovery case — is reachable.
  sLamA = makeSlider({ id: 's4-lamA', symbol: '\\lambda_A', bind: [params, 'lam_A'], min: 0, max: 20, step: 0.05 });
  sLamA.onChange(() => markDirty());

  sDA = makeSlider({ id: 's4-DA', symbol: 'D = D_A/D_L', bind: [params, 'D_A_nd'], min: 0.05, max: 50, step: null, log: true });
  sDA.onChange(() => markDirty());

  sGammaA = makeSlider({ id: 's4-gammaA', symbol: '\\tilde{\\Gamma}_A', bind: [params, 'gamma_A'], min: 1e-3, max: 10, step: null, log: true });
  sGammaA.onChange(() => markDirty());

  m61SectionEl = section('Cue model (M6.1)', [sLamA.el, sDA.el, sGammaA.el]);
  container.appendChild(m61SectionEl);

  // M6.2-specific nondim sliders (quorum-throttled production).
  //   β emission rate of the quorum signal (β = 0 recovers M1 exactly).
  //   D = D_Q/D_L, γ decay → tone 𝓠_ss = βσ̃/γ, screening ℓ_Q = √(D/γ).
  //   m throttle exponent: it sets how sharply production falls with density
  //   (the effective source is α/(1+𝓠_ss^m)); the resulting wave speed is an
  //   OUTPUT — run the sweep. See catalog §7b.
  // β linear (not log) so β = 0 — the exact M1-recovery case — is reachable.
  sBetaQ = makeSlider({ id: 's4-betaQ', symbol: '\\beta', bind: [params, 'beta_Q'], min: 0, max: 50, step: 0.05 });
  sBetaQ.onChange(() => { updateSweepPlot(); markDirty(); });

  sDQ = makeSlider({ id: 's4-DQ', symbol: 'D = D_Q/D_L', bind: [params, 'D_Q_nd'], min: 0.05, max: 50, step: null, log: true });
  sDQ.onChange(() => { updateSweepPlot(); markDirty(); });

  sGammaQ = makeSlider({ id: 's4-gammaQ', symbol: '\\tilde{\\gamma}', bind: [params, 'gamma_Q'], min: 1e-3, max: 10, step: null, log: true });
  sGammaQ.onChange(() => { updateSweepPlot(); markDirty(); });

  smQ = makeSlider({ id: 's4-mQ', symbol: 'm', bind: [params, 'm_Q'], min: 0.5, max: 8, step: 0.5 });
  smQ.onChange(() => { updateSweepPlot(); markDirty(); });

  qIcToggle = makeToggle({
    label: '𝓠 initial condition',
    options: [
      { id: 'ss',   label: 'steady tone' },
      { id: 'zero', label: 'zero' },
    ],
    value: params.q_ic_ss ? 'ss' : 'zero',
    onChange(v) { params.q_ic_ss = (v === 'ss'); markDirty(); },
  });

  m62SectionEl = section('Cue model (M6.2)',
    [sBetaQ.el, sDQ.el, sGammaQ.el, smQ.el, qIcToggle.el]);
  container.appendChild(m62SectionEl);

  // ── Cell nondim section ──
  // Five independent groups replacing old Λ, χ̃, μ̃ scheme:
  //   λ (activation), ν (GL nonlinearity), κ (chemotactic coupling),
  //   μ (motility), ϑ (noise), plus 𝓛_c threshold.
  sLam = makeSlider({ id: 's4-lam', symbol: '\\lambda', bind: [params, 'lam'], min: 0.001, max: 50, step: null, log: true });
  sLam.onChange(() => markDirty());

  sNu = makeSlider({ id: 's4-nu', symbol: '\\nu', bind: [params, 'nu'], min: 0.001, max: 50, step: null, log: true });
  sNu.onChange(() => markDirty());

  sKap = makeSlider({ id: 's4-kap', symbol: '\\kappa', bind: [params, 'kap'], min: 0, max: 5, step: 0.005 });
  sKap.onChange(() => markDirty());

  sMuNd = makeSlider({ id: 's4-mu', symbol: '\\mu', bind: [params, 'mu'], min: 1e-5, max: 1, step: null, log: true });
  sMuNd.onChange(() => markDirty());

  sLcNd = makeSlider({ id: 's4-Lc', symbol: '\\tilde{\\mathcal{L}}_c', bind: [params, 'L_c_nd'], min: 0.1, max: 5, step: 0.05 });
  sLcNd.onChange(() => markDirty());

  sTht = makeSlider({ id: 's4-tht', symbol: '\\vartheta', bind: [params, 'tht'], min: 1e-7, max: 0.01, step: null, log: true });
  sTht.onChange(() => markDirty());

  container.appendChild(section('Cell nondim', [sLam.el, sNu.el, sKap.el, sMuNd.el, sLcNd.el, sTht.el]));

  // ── IC section ──
  // Firing source: s_fire = 1 ≡ "every point in the firing disk emits as
  // strongly as one saturated cell" (σ̃-rescaled in worker; plan §A2).
  // t_fire default 5.0 (≪ t_max=300 — wave nucleates and runs on its own).
  // r_fire is dim-canonical: nondim slider back-computes dim.r_fire_dim via the
  // current ℓ₀, then recomputeFromDim derives params.r_fire and pushes to both
  // sliders. r_fire ≤ R_dish constraint lives on the dim side (see recomputeFromDim).
  sRFire = makeSlider({ id: 's4-rfire', symbol: '\\tilde{r}_{\\text{fire}}', bind: [params, 'r_fire'], min: 0.1, max: 30, step: 0.1 });
  sRFire.onChange(v => {
    if (applyingDim) return;
    const ell0 = dim.a / (dim.L_0 * dim.D_L);
    dim.r_fire_dim = v * ell0;
    recomputeFromDim();
    markDirty();
  });

  // t_fire upper bound raised to 200 (new t_max = 300; wave can be fired later).
  sTFire = makeSlider({ id: 's4-tfire', symbol: '\\tilde{t}_{\\text{fire}}', bind: [params, 't_fire'], min: 0, max: 200, step: 0.1 });
  sTFire.onChange(() => markDirty());

  sSFire = makeSlider({ id: 's4-sfire', symbol: 's_{\\text{fire}}', bind: [params, 's_fire'], min: 0, max: 20, step: 0.05 });
  sSFire.onChange(() => markDirty());

  container.appendChild(section('Initial stimulus (firing source)', [sRFire.el, sTFire.el, sSFire.el]));

  // ── Numerics section ──
  const sDt = makeSlider({ id: 's4-dt', symbol: 'd\\tilde{t}', bind: [params, 'dt'], min: 0.001, max: 0.1, step: 0.001 });
  sDt.onChange(() => markDirty());

  // t_max is dim-canonical: nondim slider back-computes dim.t_max_dim via t₀,
  // recomputeFromDim derives params.t_max. Per-geometry preference stored in seconds.
  sTmax = makeSlider({ id: 's4-tmax', symbol: '\\tilde{t}_{\\max}', bind: [params, 't_max'], min: 0.1, max: 5000, step: 0.5 });
  sTmax.onChange(v => {
    if (applyingDim) return;
    const nd = dimToNondim(dim, params.geometry);
    dim.t_max_dim = v * nd.t0;
    if (params.geometry === '2d2d') dim.t_max_dim_2d2d = dim.t_max_dim;
    else                            dim.t_max_dim_2d3d = dim.t_max_dim;
    recomputeFromDim();
    markDirty();
  });

  const sSeed = makeSlider({ id: 's4-seed', symbol: '\\text{seed}', bind: [params, 'seed'], min: 1, max: 999, step: 1, transform: Math.round, fmt: v => String(v) });
  sSeed.onChange(() => markDirty());

  container.appendChild(section('Numerics', [sDt.el, sTmax.el, sSeed.el]));

  // ── Tracked cell (for bead-in-free-energy view) ──
  // trackedCellIdx is a module-scoped let (not a params field), so bind isn't
  // applicable; transform: Math.round canonicalizes the value at the source.
  sTracked = makeSlider({
    id: 's4-tracked', symbol: '\\text{tracked cell}',
    value: trackedCellIdx, min: 0, max: Math.max(0, params.N - 1), step: 1,
    transform: Math.round, fmt: v => String(v),
  });
  sTracked.onChange(v => { trackedCellIdx = v; redraw(); });
  container.appendChild(section('Bead view', [sTracked.el]));

  // ── 2D-3D z-grid sliders ──
  // alpha_z default 1.5 (plan §B2): at N_z=16, h_0=0.1 this gives
  //   z_max = h_0·(1.5^16−1)/0.5 ≈ 23 ℓ_0, resolving the exponential
  //   z-column comfortably. Old default 1.1 caused the silent 10× L bug (§B2).
  function updateZmaxKpi() {
    if (kpis && params.geometry === '2d3d') {
      const zMax = params.h_0 * (Math.pow(params.alpha_z, params.N_z) - 1) /
                   (params.alpha_z - 1);
      kpis.set('zmax', zMax.toFixed(1));
    }
  }

  const sNz = makeSlider({ id: 's4-Nz', symbol: 'N_z', bind: [params, 'N_z'], min: 4, max: 64, step: 1, transform: Math.round, fmt: v => String(v) });
  sNz.onChange(() => { updateZmaxKpi(); markDirty(); });

  const sH0 = makeSlider({ id: 's4-h0', symbol: 'h_0', bind: [params, 'h_0'], min: 0.01, max: 1, step: 0.01 });
  sH0.onChange(() => { updateZmaxKpi(); markDirty(); });

  // alpha_z range extended to 2.0 (plan §B2).
  const sAlpha = makeSlider({ id: 's4-alpha', symbol: '\\alpha_z', bind: [params, 'alpha_z'], min: 1.05, max: 2.0, step: 0.05 });
  sAlpha.onChange(() => { updateZmaxKpi(); markDirty(); });

  const nzRow    = sNz.el;
  const h0Row    = sH0.el;
  const alphaRow = sAlpha.el;

  // h dim slider for 2D-2D: h̃ = h/ℓ_0 is the source prefactor (plan §B3).
  // Visible only in 2D-2D. Default h = 10 µm → h̃ = 1 at ℓ_0 = 10 µm.
  sH = makeSlider({
    id: 's4-h', symbol: 'h', bind: [dim, 'h'], min: 1e-3, max: 1e6, step: null, log: true,
    units: 'µm',
  });
  sH.onChange(() => recomputeFromDim());
  attachDimDesc(sH, 'Layer thickness (2D-2D only); sets h̃ = h/ℓ₀');

  container.appendChild(detailsSection('2D-3D grid (z-direction)', [nzRow, h0Row, alphaRow]));

  // ── Dimensional inputs (collapsed) ──
  // R_dim and σ_dim are dim mirrors of the Geometry-section nondim knobs;
  // both sides edit `dim` (the canonical source); recomputeFromDim() keeps
  // them in sync. N is still a KPI (derived from σ_dim·π·R_dim²).
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
      id: `s4-dim-${id}`, symbol: sym, bind: [dim, key],
      min, max, step, log: !!logScale,
      units: unit || '',
    });
    s.onChange(() => recomputeFromDim());
    attachDimDesc(s, name);
    return s;
  };

  // Dim mirrors of the Geometry-section knobs (R_dim, σ_dim). Editing these
  // updates `dim` directly; recomputeFromDim pushes back to the nondim sliders.
  sRdimDim = makeSlider({
    id: 's4-dim-Rdim', symbol: 'R_{\\text{dish}}', bind: [dim, 'R_dim'],
    min: 1, max: 5000, step: null, log: true, units: 'µm',
  });
  sRdimDim.onChange(() => { if (!applyingDim) recomputeFromDim(); });
  attachDimDesc(sRdimDim, 'Petri dish radius (dim mirror of R̃)');

  sSigmaDim = makeSlider({
    id: 's4-dim-sigma', symbol: '\\sigma', bind: [dim, 'sigma'],
    min: 1e-6, max: 10, step: null, log: true, units: 'cells/µm²',
  });
  sSigmaDim.onChange(() => { if (!applyingDim) recomputeFromDim(); });
  attachDimDesc(sSigmaDim, 'Cell surface density (dim mirror of σ̃ = σ·ℓ₀²)');

  sA = makeSlider({
    id: 's4-dim-a', symbol: 'a', bind: [dim, 'a'], min: 1, max: 1e9, step: null, log: true,
    units: 'nM·µm³/s',
  });
  sA.onChange(() => recomputeFromDim());
  attachDimDesc(sA, 'Per-cell LTB4 emission strength; sets ℓ₀ = a/(L₀ D_L)');

  const sDL    = makeDimSlider('DL',    'D_L',      'D_L',     1,     1000, null, true,
                              'LTB4 diffusion constant', 'µm²/s');
  const sL0    = makeDimSlider('L0',    'L_0',      'L_0',     0.01,  100,  null, true,
                              'LTB4 Hill threshold (emission)', 'nM');
  const sLcDim = makeDimSlider('Lc',    'L_c',      'L_c',     0.01,  100,  null, true,
                              'GL polarization threshold', 'nM');
  const sR0    = makeDimSlider('r0',    'r_0',      'r_0',     0.001, 10,   null, true,
                              'GL relaxation rate prefactor (sets λ = r₀ t₀ L₀)', '1/(s·nM)');
  const sU     = makeDimSlider('u',     'u',        'u',       0.001, 10,   null, true,
                              'GL nonlinear coefficient u (sets ν = u² t₀/w)', '1/s');
  const sW     = makeDimSlider('w',     'w',        'w',       0.001, 10,   null, true,
                              'GL nonlinear coefficient w', '1/s');
  const sChiD  = makeDimSlider('chi',   '\\chi',    'chi',     0.001, 100,  null, true,
                              'Chemotactic susceptibility (sets κ = aχ/(D_L² p₀))', 'µm²/(s·nM)');
  const sMuD   = makeDimSlider('mu',    '\\mu',     'mu',      0.001, 10,   null, true,
                              'Cell motility (sets μ = μ_dim t₀/ℓ₀)', 'µm/s');
  const sTheta = makeDimSlider('theta', '\\theta',  'theta',   1e-7,  0.1,  null, true,
                              'Noise amplitude (sets ϑ = w L₀ t₀ θ/u)', 'nM/s');
  const sGLdec = makeDimSlider('GamL',  '\\Gamma_L','Gamma_L', 0,     10,   0.01, false,
                              'LTB4 decay rate (sets Γ_L = γ_L t₀)', '1/s');

  // M2 dim sliders.
  const sBetaDim   = makeDimSlider('Beta',   '\\beta',   'Beta',   0.001, 100, null, true,
                                  'Per-cell R production rate', '1/s');
  const sGammaRDim = makeDimSlider('GammaR', '\\Gamma_R','Gamma_R',0.001, 100, null, true,
                                  'Per-cell R degradation rate', '1/s');
  const sLrDim     = makeDimSlider('Lr',     'L_r',      'L_r',    0.001, 100, null, true,
                                  'Second activation threshold (R-ODE)', 'nM');

  // r_fire and t_max dim mirrors (canonical-side). r_fire constrained ≤ R_dim
  // inside recomputeFromDim. t_max is per-geometry (dim.t_max_dim_<geom>).
  sRfireDim = makeSlider({
    id: 's4-dim-rfire', symbol: 'r_{\\text{fire}}', bind: [dim, 'r_fire_dim'],
    min: 0.1, max: 10000, step: null, log: true, units: 'µm',
  });
  sRfireDim.onChange(() => { if (!applyingDim) recomputeFromDim(); });
  attachDimDesc(sRfireDim, 'Firing source radius (dim mirror of r̃_fire = r/ℓ₀; ≤ R_dish)');

  sTmaxDim = makeSlider({
    id: 's4-dim-tmax', symbol: 't_{\\max}', bind: [dim, 't_max_dim'],
    min: 0.01, max: 1e6, step: null, log: true, units: 's',
  });
  sTmaxDim.onChange(() => {
    if (applyingDim) return;
    if (params.geometry === '2d2d') dim.t_max_dim_2d2d = dim.t_max_dim;
    else                            dim.t_max_dim_2d3d = dim.t_max_dim;
    recomputeFromDim();
    markDirty();
  });
  attachDimDesc(sTmaxDim, 'Simulation duration (dim mirror of t̃_max = t/t₀)');

  // σ_dim is in the Geometry section (primary knob). R_dim and N are KPIs.
  // h is listed last (2D-2D only; visibility toggled by geometry).
  const dimChildren = [sRdimDim.el, sSigmaDim.el,
                       sA.el, sDL.el, sL0.el, sLcDim.el, sR0.el,
                       sU.el, sW.el, sChiD.el, sMuD.el, sTheta.el, sGLdec.el,
                       sBetaDim.el, sGammaRDim.el, sLrDim.el,
                       sRfireDim.el, sTmaxDim.el, sH.el];
  container.appendChild(detailsSection('Dimensional inputs', dimChildren));

  // ── Model visibility ──
  function updateModelVisibility() {
    const m = params.model;
    m2SectionEl.style.display = (m === 'M2') ? '' : 'none';
    m61SectionEl.style.display = (m === 'M6.1') ? '' : 'none';
    m62SectionEl.style.display = (m === 'M6.2') ? '' : 'none';
    // Firing source (initiation) is used by M1/M6.1/M6.2; M2 self-nucleates.
    const usesFiring = (m === 'M1' || m === 'M6.1' || m === 'M6.2');
    sTFire.el.style.display = usesFiring ? '' : 'none';
    sSFire.el.style.display = usesFiring ? '' : 'none';
  }

  // ── Initial visibility based on geometry ──
  nzRow.style.display    = (params.geometry === '2d3d') ? '' : 'none';
  h0Row.style.display    = (params.geometry === '2d3d') ? '' : 'none';
  alphaRow.style.display = (params.geometry === '2d3d') ? '' : 'none';
  sH.el.style.display    = (params.geometry === '2d2d') ? '' : 'none';
  sHtilde.el.style.display = (params.geometry === '2d2d') ? '' : 'none';

  updateModelVisibility();

  // Initial dim → nondim sync: populates all derived params (R_dim, N, σ̃, etc.)
  // and sets KPIs. Must run after all slider declarations.
  recomputeFromDim();

  // ── Decorate plots with KaTeX axis labels ──
  function decoratePlots() {
    decoratePlot('cv-dish',    { titleTex: '\\mathcal{L}(\\tilde{x},\\tilde{y})' });
    decoratePlot('cv-dishQ',   { titleTex: '\\mathcal{Q}(\\tilde{x},\\tilde{y})' });
    decoratePlot('cv-chan', {
      titleTex:  '\\Psi(\\tilde{t})\\ \\text{(channelisation)}',
      xLabelTex: '\\tilde{t}',
      yLabelTex: '\\Psi = \\max_{m\\ge2}\\Psi_m',
    });
    decoratePlot('cv-chanspec', {
      titleTex:  '\\Psi_m\\ \\text{(angular spectrum)}',
      xLabelTex: 'm',
      yLabelTex: '\\Psi_m',
    });
    decoratePlot('cv-profile', {
      titleTex:  '\\mathcal{L}(\\tilde{r})',
      xLabelTex: '\\tilde{r}',
      yLabelTex: '\\mathcal{L}',
    });
    const rPanelTex = (params.model === 'M6.1')
      ? { titleTex: '\\mathcal{A}(\\tilde{r})', xLabelTex: '\\tilde{r}', yLabelTex: '\\mathcal{A}' }
      : (params.model === 'M6.2')
      ? { titleTex: '\\mathcal{Q}(\\tilde{r})', xLabelTex: '\\tilde{r}', yLabelTex: '\\mathcal{Q}' }
      : { titleTex: '\\tilde{R}(\\tilde{r})',   xLabelTex: '\\tilde{r}', yLabelTex: '\\tilde{R}' };
    decoratePlot('cv-rmean', rPanelTex);
    decoratePlot('cv-vr', {
      titleTex:  '\\langle|\\tilde{v}_r|\\rangle(\\tilde{t})',
      xLabelTex: '\\tilde{t}',
      yLabelTex: '\\langle|\\tilde{v}_r|\\rangle',
    });
    decoratePlot('cv-pmag', {
      titleTex:  '\\langle|\\mathbf{P}|\\rangle(\\tilde{t})',
      xLabelTex: '\\tilde{t}',
      yLabelTex: '\\langle|\\mathbf{P}|\\rangle',
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
    if (document.getElementById('cv-csweep')) {
      decoratePlot('cv-csweep', {
        titleTex:  '\\tilde{c}(\\tilde{\\sigma})',
        xLabelTex: '\\tilde{\\sigma}',
        yLabelTex: '\\tilde{c}',
      });
    }
  }

  if (window.katex) {
    decoratePlots();
  } else {
    window.addEventListener('katex-ready', decoratePlots);
    setTimeout(decoratePlots, 300);
  }

  // Initial c(σ̃) panel: show the theory curves immediately (M6.1 page only);
  // measured points fill in after the user runs a sweep.
  updateSweepPlot();
}
