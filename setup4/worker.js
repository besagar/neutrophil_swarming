// setup4/worker.js  (Web Worker entry point, ES module)
// Batch loop: for each time step:
//   1. Accumulate cell emissions → field (PIC consistency)
//   2. Step field PDE (with internal substepping)
//   3. Step agents (SDE)
// Posts frames to main thread; posts progress every ~1%.
//
// Message protocol:
//   main → worker: { type: 'run', params: SimParams }
//                  { type: 'cancel' }
//   worker → main: { type: 'progress', pct, step, t }
//                  { type: 'frame', step, t, radialProfile: Float32Array,
//                    agentX: Float32Array, agentY: Float32Array,
//                    emitting: Uint8Array,
//                    Lfield: Uint8Array (per-frame max-normalized 0-255),
//                    Lmax: number (absolute peak for KPI / radial scale) }
//                  { type: 'done' }
//                  { type: 'error', message }

import { createField }  from './solvers/field.js';
import { createSolver } from './solvers/index.js';
import { initAgents, stepAgents, hillEmission } from './agents.js';
import { makeRng }      from '../shared/rng.js';

// L heatmap is compressed to Uint8 (per-frame max-normalized to 0-255) and
// sent EVERY frame. The renderer only uses the normalized 0-1 ratio for the
// viridis colormap, so 256 levels is visually indistinguishable from Float32.
// Memory per frame: 128² × 1 B = 16 KB; over 1000 frames ≈ 16 MB.

let cancelled = false;

self.onmessage = function (e) {
  const msg = e.data;
  if (msg.type === 'cancel') {
    cancelled = true;
    return;
  }
  if (msg.type === 'run') {
    cancelled = false;
    try {
      runSimulation(msg.params);
    } catch (err) {
      self.postMessage({ type: 'error', message: String(err) });
    }
  }
};

/**
 * @param {SimParams} p
 */
function runSimulation(p) {
  // Unpack params with defaults.
  const {
    geometry  = '2d2d',
    model     = 'M1',
    N         = 500,
    N_grid    = 128,
    R_dish    = 10,
    Lambda    = 1,
    L_c_nd    = 1,
    chi_nd    = 0.75,
    mu_nd     = 0.05,
    lam       = 1,
    tht       = 1e-4,
    n_L       = 10,
    gamma_L   = 0,
    sigma_tilde = 0.02,   // σ̃ = σ·ℓ_0². Per-cell emission rate is 1/σ̃ (fixed by dim params, NOT by N).
    // Time-limited firing source (parent_solver convention; matches
    // la_2d3d_solver.py: firing_radius=2, firing_duration=5, firing_strength=1).
    // Adds 2·s_fire·δ(z̃) inside r̃<r_fire for t̃<t_fire (2D-3D), or a flat
    // s_fire bulk source over the same disk (2D-2D). After t̃≥t_fire the
    // wave self-sustains via the relay term.
    r_fire    = 2.0,
    t_fire    = 5.0,
    s_fire    = 1.0,
    dt        = 0.01,
    t_max     = 50,
    seed      = 7,
    N_z       = 16,
    h_0       = 0.1,
    alpha_z   = 1.1,
  } = p;

  const n_steps   = Math.ceil(t_max / dt);
  // Target ~500 output frames per run (Uint8 heatmap = 16 KB, agent state
  // ≈ 8–16 KB at N≤2000 → ~12-16 MB total at 500 frames). The cap matters
  // mainly for long runs; short runs (n_steps < 500) emit every step.
  const n_frames  = Math.min(500, n_steps);
  const frameEvery = Math.max(1, Math.floor(n_steps / n_frames));
  const progressEvery = Math.max(1, Math.floor(n_steps / 100));

  // Validate.
  if (n_steps > 1e6) {
    self.postMessage({ type: 'error', message: `Too many steps: ${n_steps}. Reduce t_max or increase dt.` });
    return;
  }

  // Build field and solver.
  const field  = createField(N_grid, R_dish);
  const solver = createSolver(
    geometry, model, field,
    { N_z, h_0, alpha: alpha_z }
  );

  // Initialize field IC: zero everywhere. The wave is kicked off by the
  // firing source applied in the time loop (see below).
  solver.reset();

  // Initialize agents.
  const agents = initAgents(N, R_dish, seed);
  const rng    = makeRng(seed + 1);  // separate RNG stream from IC placement

  // Params object passed to solver and agents.
  const simParams = {
    Lambda, L_c: L_c_nd, chi: chi_nd, mu: mu_nd, lam, tht,
    n_L, gamma_L, dt, R_dish,
  };

  let frameCount = 0;

  // ─── main loop ──────────────────────────────────────────────────────────
  for (let step = 0; step <= n_steps; step++) {
    if (cancelled) return;

    const t = step * dt;

    // Emit frame before first step (step 0) and at regular intervals.
    const isFrame = (step === 0) || (step % frameEvery === 0) || (step === n_steps);

    if (isFrame) {
      const radialProfile = field.getRadialProfile();
      // Copy agent arrays (transferable for efficiency).
      const agentX    = new Float32Array(agents.x);
      const agentY    = new Float32Array(agents.y);
      const agentPx   = new Float32Array(agents.Px);
      const agentPy   = new Float32Array(agents.Py);
      const emitting  = new Uint8Array(agents.emitting);
      // Cell-sampled ∇𝓛 — used by the bead-plot's chemotactic-tilt term.
      const agentGx   = new Float32Array(agents.N);
      const agentGy   = new Float32Array(agents.N);
      for (let i = 0; i < agents.N; i++) {
        const g = field.sampleGrad(agents.x[i], agents.y[i]);
        agentGx[i] = g.gx;
        agentGy[i] = g.gy;
      }

      // Compress L field to Uint8. Normalize by the *radial-profile max*
      // (the smooth azimuthally-averaged field), NOT the grid max.
      //
      // Why: cells emit via discrete PIC δ-sources, creating tiny grid-cell
      // hot spots that are typically 5-20× the smooth wave amplitude. Using
      // the grid max as the normalizer lets one stochastically-hot cell
      // dominate the entire colormap, dimming all other cells frame-to-frame
      // — a flickering "cells stop shining one-by-one" visual.
      // The radial profile averages cell hot spots together across each
      // annulus, so its max is a stable proxy for the smooth wave amplitude.
      // Cells brighter than that just clip to 255 (saturated yellow), giving
      // a stable bright-dot appearance independent of neighbours.
      const Lf = field.getLfield();
      let lmaxGrid = 0;
      for (let k = 0; k < Lf.length; k++) if (Lf[k] > lmaxGrid) lmaxGrid = Lf[k];
      let lmaxRadial = 0;
      for (let k = 0; k < radialProfile.length; k++) {
        if (radialProfile[k] > lmaxRadial) lmaxRadial = radialProfile[k];
      }
      // Scale of bright-yellow: smooth wave amplitude with a small headroom
      // multiplier so the wave peak doesn't pin at saturated yellow either.
      const lmaxNorm = Math.max(lmaxRadial * 1.2, lmaxGrid * 1e-3, 1e-6);
      const inv = 255 / lmaxNorm;
      const Lfield = new Uint8Array(Lf.length);
      for (let k = 0; k < Lf.length; k++) {
        const v = Lf[k] * inv;
        Lfield[k] = v > 255 ? 255 : (v < 0 ? 0 : (v | 0));
      }

      const msg = {
        type: 'frame', step, t, frameCount,
        radialProfile,
        agentX, agentY, agentPx, agentPy, agentGx, agentGy, emitting,
        Lfield,             // Uint8, normalized by smooth-wave scale
        Lmax: lmaxGrid,     // true grid maximum (for KPI)
        LmaxNorm: lmaxNorm, // the actual normalizer used (Uint8 byte = L × 255 / LmaxNorm)
      };

      self.postMessage(msg, [
        radialProfile.buffer, agentX.buffer, agentY.buffer,
        agentPx.buffer, agentPy.buffer, agentGx.buffer, agentGy.buffer,
        emitting.buffer, Lfield.buffer,
      ]);
      frameCount++;
    }

    // Progress update.
    if (step > 0 && step % progressEvery === 0) {
      self.postMessage({ type: 'progress', pct: step / n_steps, step, t });
    }

    if (step === n_steps) break;  // don't advance past t_max

    // ── PIC loop order (enforces consistency): ──────────────────────────
    // 1. Accumulate cell emissions onto grid (before field steps).
    //
    // Discrete-ABM nondim source (per-cell δ-source form):
    //   ∂_t̃ 𝓛 = ∇̃²𝓛 + (1/σ̃) · Σ_i H⁺(𝓛_i;1;n_L) δ̃(x̃ − x̃_i)  [+ δ̃(z̃) in 2D-3D]
    // N enters ONLY through the sum. The per-cell prefactor 1/σ̃ is fixed by
    // dim params (σ̃ = σ_dim · ℓ_0²), independent of how many discrete cells N
    // we simulate. Doubling N doubles total emission, as in real biology.
    const inv_sigma = 1 / sigma_tilde;

    for (let i = 0; i < N; i++) {
      const L_i = Math.max(0, field.sample(agents.x[i], agents.y[i]));
      const h   = hillEmission(L_i, n_L);  // H⁺(𝓛_i; 1; n_L)
      // Per-cell contribution to grid via PIC weights; solver divides by dx²
      // (2D-2D) or by dx²·h_0/2 (2D-3D δ-z discretization) to get the
      // concentration delta at each node per agent step.
      field.accumulateSource(agents.x[i], agents.y[i], h * inv_sigma * dt);
    }

    // 1b. Time-limited firing source (kicks the wave off at t̃ = 0).
    // Adds s_fire · dx² · dt to src[k] at every node inside r̃ < r_fire.
    // The solver consumes src with geometry-specific factors:
    //   2D-2D:  ΔL = s_fire · dt           (directSrc = src/dx²)
    //   2D-3D:  ΔL_z=0 = 2 s_fire · dt/h_0 (surfSrc = 2·src/(dx²·h_0))
    // matching parent_solver's `firing_strength · 2·δ(z̃)` in 2D-3D and a
    // flat bulk source in 2D-2D. After t ≥ t_fire the relay sustains the wave.
    if (t < t_fire && r_fire > 0 && s_fire > 0) {
      field.addFiringSource(r_fire, s_fire * field.dx * field.dx * dt);
    }

    // 2. Step field PDE (includes source injection and substepping).
    solver.step(dt, simParams);

    // 3. Step agents (SDE).
    stepAgents(agents, field, simParams, rng);
  }

  self.postMessage({ type: 'done' });
}
