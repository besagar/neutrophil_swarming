// setup4/worker.js  (Web Worker entry point, ES module)
// Batch loop for a single run + a multi-run c(σ̃) density sweep.
//
// The per-step physics lives in sim_core.js (createSim / advance /
// measureWaveSpeed) so the single run and the sweep share identical stepping.
// This file only owns worker plumbing: message dispatch, frame packing, and
// the sweep driver.
//
// Message protocol:
//   main → worker: { type: 'run',   params: SimParams }
//                  { type: 'sweep', params: SimParams, sigmas: number[] }
//                  { type: 'cancel' }
//   worker → main: { type: 'progress', pct, step, t }
//                  { type: 'frame', step, t, radialProfile, agentX, agentY,
//                    agentPx, agentPy, agentGx, agentGy, emitting, stuck, agentR,
//                    Lfield: Uint8Array, Lmax, LmaxNorm,
//                    Qfield?: Uint8Array, Qmax? }   (M6.2 only — second dish)
//                  { type: 'done' }
//                  { type: 'sweep-progress', i, n, sigma, c }
//                  { type: 'sweep-done', points: [{sigma, c}] }
//                  { type: 'error', message }

import { createSim, advance, measureWaveSpeed } from './sim_core.js';

let cancelled = false;

self.onmessage = function (e) {
  const msg = e.data;
  if (msg.type === 'cancel') { cancelled = true; return; }
  if (msg.type === 'run') {
    cancelled = false;
    try { runSimulation(msg.params); }
    catch (err) { self.postMessage({ type: 'error', message: String(err) }); }
  }
  if (msg.type === 'sweep') {
    cancelled = false;
    try { runSweep(msg.params, msg.sigmas); }
    catch (err) { self.postMessage({ type: 'error', message: String(err) }); }
  }
};

/**
 * Single batch run: steps the system and posts frames for time-scrubbing.
 * @param {SimParams} p
 */
function runSimulation(p) {
  const sim = createSim(p);
  const { field, solver, agents, model, dt, n_steps } = sim;

  if (n_steps > 1e6) {
    self.postMessage({ type: 'error', message: `Too many steps: ${n_steps}. Reduce t_max or increase dt.` });
    return;
  }

  // Target ~500 output frames; short runs emit every step.
  const n_frames      = Math.min(500, n_steps);
  const frameEvery    = Math.max(1, Math.floor(n_steps / n_frames));
  const progressEvery = Math.max(1, Math.floor(n_steps / 100));

  let frameCount = 0;

  for (let step = 0; step <= n_steps; step++) {
    if (cancelled) return;
    const t = step * dt;

    const isFrame = (step === 0) || (step % frameEvery === 0) || (step === n_steps);
    if (isFrame) {
      const radialProfile = field.getRadialProfile();
      const agentX   = new Float32Array(agents.x);
      const agentY   = new Float32Array(agents.y);
      const agentPx  = new Float32Array(agents.Px);
      const agentPy  = new Float32Array(agents.Py);
      const emitting = new Uint8Array(agents.emitting);
      // 1 = cell has engaged the target circle and no longer moves.
      const stuck    = new Uint8Array(agents.stuck);
      // Per-cell "R panel" readout: M2 → R̃_i, M6.1 → 𝓐_i, M6.2 → 𝓠_i,
      // M1 → zeros. (M6.2 also uses this channel for emission opacity.)
      const agentR = new Float32Array(agents.N);
      if ((model === 'M6.1' || model === 'M6.2') && solver.fieldAux) {
        for (let i = 0; i < agents.N; i++) {
          agentR[i] = Math.max(0, solver.fieldAux.sample(agents.x[i], agents.y[i]));
        }
      } else {
        agentR.set(agents.R);
      }
      // Cell-sampled ∇𝓛 for the bead-plot chemotactic tilt.
      const agentGx = new Float32Array(agents.N);
      const agentGy = new Float32Array(agents.N);
      for (let i = 0; i < agents.N; i++) {
        const g = field.sampleGrad(agents.x[i], agents.y[i]);
        agentGx[i] = g.gx; agentGy[i] = g.gy;
      }

      // Compress L to Uint8, normalized by the smooth radial-profile max (not
      // the grid max) so a single stochastically-hot cell can't dim the frame.
      const Lf = field.getLfield();
      let lmaxGrid = 0;
      for (let k = 0; k < Lf.length; k++) if (Lf[k] > lmaxGrid) lmaxGrid = Lf[k];
      let lmaxRadial = 0;
      for (let k = 0; k < radialProfile.length; k++) {
        if (radialProfile[k] > lmaxRadial) lmaxRadial = radialProfile[k];
      }
      const lmaxNorm = Math.max(lmaxRadial * 1.2, lmaxGrid * 1e-3, 1e-6);
      const inv = 255 / lmaxNorm;
      const Lfield = new Uint8Array(Lf.length);
      for (let k = 0; k < Lf.length; k++) {
        const v = Lf[k] * inv;
        Lfield[k] = v > 255 ? 255 : (v < 0 ? 0 : (v | 0));
      }

      // M6.2 only: also ship the 𝓠 grid so the page can draw a second petri
      // dish (cells + 𝓠). Normalized by its own grid max — 𝓠 is close to
      // uniform, so a radial-profile normalization would flatten it to a
      // featureless wash; the grid max keeps the per-cell hot spots and any
      // density bump visible. Not sent for other models (an extra N_grid²
      // bytes per frame is pure waste when nothing renders it).
      let Qfield = null, qmaxGrid = 0;
      if (model === 'M6.2' && solver.fieldAux) {
        const Qf = solver.fieldAux.getLfield();
        for (let k = 0; k < Qf.length; k++) if (Qf[k] > qmaxGrid) qmaxGrid = Qf[k];
        const qinv = 255 / Math.max(qmaxGrid, 1e-9);
        Qfield = new Uint8Array(Qf.length);
        for (let k = 0; k < Qf.length; k++) {
          const v = Qf[k] * qinv;
          Qfield[k] = v > 255 ? 255 : (v < 0 ? 0 : (v | 0));
        }
      }

      const outMsg = {
        type: 'frame', step, t, frameCount,
        radialProfile,
        agentX, agentY, agentPx, agentPy, agentGx, agentGy, emitting, stuck, agentR,
        Lfield, Lmax: lmaxGrid, LmaxNorm: lmaxNorm,
        Qfield, Qmax: qmaxGrid,
      };
      const transfer = [
        radialProfile.buffer, agentX.buffer, agentY.buffer,
        agentPx.buffer, agentPy.buffer, agentGx.buffer, agentGy.buffer,
        emitting.buffer, stuck.buffer, agentR.buffer, Lfield.buffer,
      ];
      if (Qfield) transfer.push(Qfield.buffer);
      self.postMessage(outMsg, transfer);
      frameCount++;
    }

    if (step > 0 && step % progressEvery === 0) {
      self.postMessage({ type: 'progress', pct: step / n_steps, step, t });
    }

    if (step === n_steps) break;
    advance(sim, t);
  }

  self.postMessage({ type: 'done' });
}

/**
 * Density sweep: for each σ̃, run a full (frame-less) sim and measure the
 * front speed c̃, posting a point as each finishes. N is recomputed per σ̃ so
 * the areal density is consistent (N = round(σ̃·π·R̃²)).
 * @param {SimParams} baseP - base params (already carries the sweep geometry:
 *                            capped R_dish / N_grid / t_max chosen by the UI)
 * @param {number[]} sigmas
 */
function runSweep(baseP, sigmas) {
  const isCancelled = () => cancelled;
  const points = [];
  for (let i = 0; i < sigmas.length; i++) {
    if (cancelled) return;
    const s = sigmas[i];
    // Mass-targeted ignition: keep the INJECTED 𝓛 mass (hence the ignition-halo
    // radius r̃_halo = r̃_fire√(s_fire σ̃ t̃_fire) = r̃_fire√K) the same at every
    // density, instead of letting it grow ∝ σ̃. Otherwise the halo alone floods
    // the dish at high σ̃ and measureWaveSpeed times diffusion, not relay.
    const s_fire = (baseP.fire_K != null)
      ? baseP.fire_K / Math.max(s * baseP.t_fire, 1e-9)
      : baseP.s_fire;
    const p = Object.assign({}, baseP, {
      sigma_tilde: s,
      s_fire,
      N: Math.max(1, Math.round(s * Math.PI * baseP.R_dish * baseP.R_dish)),
    });
    const res = measureWaveSpeed(p, isCancelled);
    if (res === null) return;   // cancelled mid-sim
    const pt = { sigma: s, c: res.c, status: res.status || (res.c > 0 ? 'ok' : 'noignite') };
    points.push(pt);
    self.postMessage({ type: 'sweep-progress', i: i + 1, n: sigmas.length, ...pt });
  }
  self.postMessage({ type: 'sweep-done', points });
}
