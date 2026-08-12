// setup4/sim_core.js
// Pure simulation core for Setup 4 — no DOM, no Web Worker `self` references,
// so it is unit-testable in node. Both the single-run worker loop and the
// multi-run c(σ̃) density sweep call the SAME stepping code here, so the two
// can never diverge physically.
//
// Public API:
//   createSim(p)            → sim context {field, solver, agents, rng, simParams, …}
//   advance(sim, t)         → advance the whole system one agent step at time t̃
//   measureWaveSpeed(p, cb) → run a full sim (no frames) and return the front speed c̃
//
// The per-step order (PIC consistency, matches docs/PLAN.md §4.4):
//   1. accumulate cell emissions (+ the basal auxiliary source for M6.1/M6.2),
//   2. add the time-limited firing source (M1/M6.1/M6.2),
//   3. step the field PDE(s),
//   4. step the agent SDE.

import { createField }  from './solvers/field.js';
import { createSolver } from './solvers/index.js';
import { initAgents, stepAgents, hillEmission, hillNeg, hillPos } from './agents.js';
import { makeRng }      from '../shared/rng.js';

/**
 * Build a simulation context from a nondim parameter object.
 * @param {Object} p - SimParams (same shape the worker receives)
 * @returns {Object} sim context (see fields assigned below)
 */
export function createSim(p) {
  const {
    geometry  = '2d2d',
    model     = 'M1',
    N         = 500,
    N_grid    = 128,
    R_dish    = 10,
    lam       = 1,
    nu        = 1,
    kap       = 0.075,
    mu        = 0.005,
    tht       = 1e-4,
    L_c_nd    = 1,
    n_L       = 10,
    gamma_L   = 0,
    sigma_tilde = 0.1,
    h_tilde   = 1,
    beta_R    = 0,
    gamma_R   = 0,
    L_r_nd    = 1,
    n_R       = 10,
    n_Lr      = 10,
    lam_A     = 0,
    D_A_nd    = 1,
    gamma_A   = 0.5,
    beta_Q    = 0,
    D_Q_nd    = 1,
    gamma_Q   = 0.5,
    m_Q       = 2,
    q_ic_ss   = true,
    stick_target = true,
    R_target     = 2.0,
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

  const n_steps = Math.ceil(t_max / dt);

  const field  = createField(N_grid, R_dish);
  const solver = createSolver(geometry, model, field, { N_z, h_0, alpha: alpha_z });
  solver.reset();

  const is2d2d = (geometry === '2d2d');

  // Auxiliary basal field (M6.1 𝓐 / M6.2 𝓠) — the solver is shared, so the
  // model-specific coefficients are mapped onto the generic aux keys here.
  //   aux_amp = per-cell source amplitude accumulated each step.
  // M6.1: amplitude 1 in BOTH geometries — its unit A_0 = b/(hD_L) is chosen
  //   per geometry to absorb the emission rate, so no h̃ factor appears.
  // M6.2: 𝓠 is measured in the PHYSICAL half-saturation constant Q_0, so the
  //   emission rate survives as β = b/(hD_LQ_0) and the h̃ bookkeeping is
  //   explicit: the 2D–3D surface-source coefficient is β·h̃ (same relation as
  //   the L source: α = 1/h̃ in 2D–2D vs 1 in 2D–3D). See catalog §7b.
  const isM62 = (model === 'M6.2');
  const D_aux     = isM62 ? D_Q_nd  : D_A_nd;
  const gamma_aux = isM62 ? gamma_Q : gamma_A;
  const aux_amp   = isM62 ? (beta_Q * (is2d2d ? 1 : h_tilde)) : 1;

  // M6.2: start 𝓠 at its mean-field steady tone so a run at small γ does not
  // silently measure unthrottled (M1) physics during the t̃ ~ 1/γ fill-up.
  //   2D–2D:  𝓠_ss = β σ̃ / γ                (uniform)
  //   2D–3D:  𝓠_ss(z̃) = (β h̃ σ̃ / √(Dγ)) e^{−z̃√(γ/D)}   (screened surface source)
  if (isM62 && q_ic_ss && solver.seedAux && gamma_aux > 0 && sigma_tilde > 0) {
    const tone = is2d2d
      ? (aux_amp * sigma_tilde) / gamma_aux
      : (aux_amp * sigma_tilde) / Math.sqrt(Math.max(D_aux * gamma_aux, 1e-12));
    solver.seedAux(tone, D_aux, gamma_aux);
  }

  // Sticking target: the target disk is excluded from the initial placement
  // (the pathogen occupies it), so no cell starts already engaged.
  const agents = initAgents(N, R_dish, seed, stick_target ? R_target : 0);
  const rng    = makeRng(seed + 1);   // separate stream from IC placement

  // M2: light up the forced-emit disk at t=0 so the first frame isn't blank.
  if (model === 'M2' && r_fire > 0) {
    const r2 = r_fire * r_fire;
    for (let i = 0; i < N; i++) {
      if (agents.x[i] * agents.x[i] + agents.y[i] * agents.y[i] < r2) agents.emitting[i] = 1;
    }
  }

  const simParams = {
    model,
    lam, nu, kap, mu, tht,
    L_c: L_c_nd, n_L, gamma_L, dt, R_dish,
    beta_R, gamma_R, L_r_nd, n_R, n_Lr, r_fire,
    stick_target, R_target,
    lam_A, D_A_nd, gamma_A,
    beta_Q, D_Q_nd, gamma_Q, m_Q,
    D_aux, gamma_aux,
  };

  return {
    field, solver, agents, rng, simParams,
    // cached loop scalars
    N, model, geometry, dt, n_steps, N_grid, R_dish,
    is2d2d,
    inv_h_tilde: 1 / h_tilde,
    n_L, n_R, lam_A, m_Q, aux_amp,
    r_fire, r_fire2: r_fire * r_fire, t_fire, s_fire, sigma_tilde,
  };
}

/**
 * Advance the whole system one agent step at (pre-step) time t̃.
 * Mutates sim.field / sim.solver / sim.agents in place.
 * @param {Object} sim - context from createSim()
 * @param {number} t   - current nondim time (for the firing-source window)
 */
export function advance(sim, t) {
  const { field, solver, agents, rng, simParams,
          N, model, is2d2d, inv_h_tilde, n_L, n_R, lam_A, m_Q, aux_amp,
          r_fire, r_fire2, t_fire, s_fire, sigma_tilde, dt } = sim;

  // 1. Cell emissions (+ basal auxiliary source for M6.1/M6.2) → grid.
  for (let i = 0; i < N; i++) {
    const xi = agents.x[i], yi = agents.y[i];
    const L_i = Math.max(0, field.sample(xi, yi));
    let w;
    if (model === 'M2') {
      const inFire = (xi * xi + yi * yi) < r_fire2;
      const wL = inFire ? 1 : hillEmission(L_i, n_L);
      const wR = hillNeg(agents.R[i], 1, n_R);
      w = wL * wR;
    } else if (model === 'M6.1') {
      // Relay gate with adenosine-shifted threshold 1 + λ_A·𝓐.
      const A_i = solver.fieldAux ? Math.max(0, solver.fieldAux.sample(xi, yi)) : 0;
      w = hillPos(L_i, 1 + lam_A * A_i, n_L);
      // Basal 𝓐 source: unconditional, amplitude 1 (no 1/h̃).
      solver.fieldAux.accumulateSource(xi, yi, aux_amp * dt);
    } else if (model === 'M6.2') {
      // Relay gate at the UNSHIFTED threshold, production throttled by the
      // quorum signal: w = H⁺(𝓛;1;n_L)·H⁻(𝓠;1;m).
      const Q_i = solver.fieldAux ? Math.max(0, solver.fieldAux.sample(xi, yi)) : 0;
      w = hillEmission(L_i, n_L) * hillNeg(Q_i, 1, m_Q);
      // Basal 𝓠 source: unconditional, amplitude β (×h̃ in 2D–3D; see createSim).
      solver.fieldAux.accumulateSource(xi, yi, aux_amp * dt);
    } else {
      w = hillEmission(L_i, n_L);
    }
    if (is2d2d) field.accumulateSource(xi, yi, w * dt * inv_h_tilde);
    else        field.accumulateSource(xi, yi, w * dt);
  }

  // 2. Time-limited firing source (M1/M6.1/M6.2 nucleation).
  // NOT throttled by 𝓠 in M6.2: this is an external stimulus (uncaging /
  // micropipette), not cell production — throttling it would make high-density
  // runs fail to ignite for an artefactual reason.
  if ((model === 'M1' || model === 'M6.1' || model === 'M6.2') &&
      t < t_fire && r_fire > 0 && s_fire > 0) {
    field.addFiringSource(r_fire, s_fire * sigma_tilde * field.dx * field.dx * dt);
  }

  // 3. Field PDE(s).
  solver.step(dt, simParams);

  // 4. Agent SDE (M6.1 needs 𝓐 for the shifted `emitting` flag).
  stepAgents(agents, field, simParams, rng, solver.fieldAux);
}

/**
 * Emission-front radius: the 90th-percentile radius of the cells that are
 * currently emitting (relaying). This is the RELAY wavefront — cells only emit
 * where 𝓛 exceeds their (adenosine-shifted) threshold, so a diffusing blob of
 * fired 𝓛 that no longer drives fresh cells over threshold does NOT extend it.
 * Contrast the 𝓛=1 contour, which a purely diffusive halo of the firing
 * source can inflate — the source of earlier false-positive wave speeds.
 * Returns 0 if fewer than 3 cells emit.
 */
function emissionFront(agents) {
  const rs = [];
  for (let i = 0; i < agents.N; i++) {
    if (agents.emitting[i]) rs.push(Math.hypot(agents.x[i], agents.y[i]));
  }
  if (rs.length < 3) return 0;
  rs.sort((a, b) => a - b);
  return rs[Math.min(rs.length - 1, Math.floor(0.9 * rs.length))];
}

/**
 * Run a full simulation (no frame output) and return the measured relay wave
 * speed. Returns c = 0 unless a genuine self-sustaining relay wave forms:
 *   (a) the emission front must spread well past the firing disk (spread), and
 *   (b) it must not have receded by the end (sustained — a diffusing/collapsing
 *       fired blob fails this).
 * Speed is a robust two-point average over the propagation phase (works for
 * both slow fronts and fast fills, unlike a windowed least-squares fit).
 * @param {Object} p - SimParams
 * @param {() => boolean} [isCancelled] - optional cancel poll; abort → returns null
 * @returns {{sigma_tilde:number, c:number, ignited:boolean}|null}
 */
export function measureWaveSpeed(p, isCancelled) {
  const sim = createSim(p);
  const { agents, dt, n_steps, R_dish, r_fire, t_fire, s_fire, sigma_tilde } = sim;
  const rMax = 0.9 * R_dish;
  const sampleEvery = Math.max(1, Math.floor(n_steps / 200));

  // Ignition-halo radius. The firing source injects a total nondim 𝓛 mass
  //   M = s_fire · σ̃ · π r̃_fire² · t̃_fire
  // (see field.addFiringSource). With Γ_L = 0 that mass is conserved, so a
  // purely diffusive halo — no relay at all — still carries ⟨𝓛⟩ > 1 out to
  //   r̃_halo = √(M/π) = r̃_fire √(s_fire σ̃ t̃_fire),
  // pushing cells over threshold and inflating the emission front. This is
  // invisible to the spread/sustained gates (the halo neither collapses nor
  // recedes), and it is what made a fully throttled M6.2 run (β → ∞, relay
  // impossible) report c ≈ 0.6. Timing therefore starts only OUTSIDE the halo.
  // In 2D–3D the same mass spreads into the half-space instead of staying in
  // the plane, so it dilutes as t^{-3/2} and the ⟨𝓛⟩>1 region is bounded by a
  // half-ball of the same mass: r̃_halo = (2 M/π)^{1/3}... in nondim terms
  // (π s_fire σ̃ r̃_fire² t̃_fire)^{1/3}. Using the 2D formula there would reject
  // every point.
  const fireMass = Math.max(0, s_fire * sigma_tilde * t_fire) * r_fire * r_fire;
  const r_halo  = sim.is2d2d
    ? Math.sqrt(fireMass)
    : Math.cbrt(Math.PI * fireMass);
  const r_start = Math.max(1.5 * r_fire, 1.15 * r_halo);

  const samples = [];   // {t, re}
  for (let step = 0; step <= n_steps; step++) {
    if (isCancelled && isCancelled()) return null;
    const t = step * dt;
    if (step % sampleEvery === 0) samples.push({ t, re: emissionFront(agents) });
    if (step === n_steps) break;
    advance(sim, t);
  }
  if (samples.length < 3) return { sigma_tilde, c: 0, ignited: false, status: 'noignite' };

  // Unmeasurable geometry: the halo alone would flood the dish, so nothing
  // outside it can be timed. Report this distinctly from "did not ignite" —
  // it is a diagnostic-setup failure, not a physical one.
  if (r_start > 0.75 * rMax) {
    return { sigma_tilde, c: 0, ignited: false, status: 'unmeasurable', r_halo };
  }

  // Self-sustainment gates: the relay must have carried emission out past the
  // ignition halo, and the front must not have collapsed back by the end.
  let maxRe = 0;
  for (const s of samples) if (s.re > maxRe) maxRe = s.re;
  const reEnd     = samples[samples.length - 1].re;
  const spread    = maxRe  > r_start;
  const sustained = reEnd  > 0.7 * maxRe;
  if (!spread || !sustained) {
    // Distinguish "no wave" from "a wave too slow to clear the halo within
    // t̃_max": a front that is still advancing at the end is real but
    // unresolved, and reporting it as 0 alongside genuine ignition failures
    // would draw a cliff where the physics has a slow tail (M6.2 at large β).
    const kLate    = Math.floor(0.6 * (samples.length - 1));
    const advancing = reEnd > samples[kLate].re + 0.05 * r_fire;
    return { sigma_tilde, c: 0, ignited: false, r_halo,
             status: advancing ? 'unresolved' : 'noignite' };
  }

  // Two-point average speed over the propagation phase: from when the front
  // clears the ignition halo (re > r_start — beyond that radius the stimulus
  // alone cannot hold 𝓛 above threshold, so further advance is genuine relay)
  // to when it nears the wall (or the last sample). No t_fire gate: fast fronts
  // fill during firing. Robust to both slow fronts and near-instant fills.
  let kIgnite = -1, kReach = -1;
  for (let k = 0; k < samples.length; k++) {
    if (kIgnite < 0 && samples[k].re > r_start) kIgnite = k;
    if (kIgnite >= 0 && samples[k].re > 0.8 * rMax) { kReach = k; break; }
  }
  if (kIgnite < 0) return { sigma_tilde, c: 0, ignited: false, r_halo, status: 'unresolved' };
  const kEnd = (kReach >= 0) ? kReach : samples.length - 1;
  const dtp = samples[kEnd].t - samples[kIgnite].t;
  // Reject too-short measurement arcs: a front that only clears the halo on the
  // very last samples gives a two-point average over a near-zero baseline,
  // which is pure noise and reads as a spuriously FAST wave (a slow M6.2 run
  // once reported 2.6× the M1 speed this way). Gate on the RADIAL span, not on
  // the sample count — a genuinely fast front crosses the whole window in few
  // samples and must still be measured.
  const span = samples[kEnd].re - samples[kIgnite].re;
  if (dtp < 1e-9 || span < Math.max(0.5 * r_fire, 0.1 * rMax)) {
    return { sigma_tilde, c: 0, ignited: false, r_halo, status: 'unresolved' };
  }
  const c = (samples[kEnd].re - samples[kIgnite].re) / dtp;
  return { sigma_tilde, c: c > 0 ? c : 0, ignited: c > 0, r_halo,
           status: c > 0 ? 'ok' : 'noignite' };
}
