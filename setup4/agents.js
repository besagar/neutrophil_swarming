// setup4/agents.js
// Per-cell state and Euler–Maruyama SDE step for Setup 4.
//
// GL polarization SDE (Euler–Maruyama, nondim):
//   dP_xi = [χ̃ ∂_x̃ 𝓛 + Λ((𝓛 − 𝓛_c) P_xi + λ(|P|²-|P|⁴) P_xi)] dt̃
//           + √(2 Λ ϑ 𝓛 dt̃) η_xi
//   dP_yi = [χ̃ ∂_ỹ 𝓛 + Λ((𝓛 − 𝓛_c) P_yi + λ(|P|²-|P|⁴) P_yi)] dt̃
//           + √(2 Λ ϑ 𝓛 dt̃) η_yi
//   dx̃_i/dt̃ = μ̃ P_xi
//   dỹ_i/dt̃ = μ̃ P_yi
//
// Note: noise prefactor is √(2 Λ ϑ 𝓛 dt̃) — Λ appears in both the activation
// term and noise because the cell time scale t_0 differs from Setup 3's 1/(r_0 L_c).
// See setup4_swarm3d.md §5.2.
//
// Boundary: reflective at r̃ = R̃_dish (hard outer boundary; no trap in Setup 4).
//
// Agents NEVER import the raw grid array. They access L only via:
//   field.sample(x, y)     → L value at position
//   field.sampleGrad(x, y) → {gx, gy} gradient
// Emissions are accumulated by the worker (not here) via field.accumulateSource().

import { makeRng } from '../shared/rng.js';

/**
 * @typedef {Object} AgentState
 * @property {Float32Array} x   - x positions [N]
 * @property {Float32Array} y   - y positions [N]
 * @property {Float32Array} Px  - x polarization [N]
 * @property {Float32Array} Py  - y polarization [N]
 * @property {Float32Array} R   - per-cell inhibitor R̃_i, nondim (M2 only; zeros in M1) [N]
 * @property {Uint8Array}  emitting - 1 if cell is currently emitting [N]
 * @property {number} N
 */

/**
 * @typedef {Object} SimParams
 * @property {number} N          - number of cells
 * @property {number} R_dish     - dish radius (nondim)
 * @property {number} Lambda     - Λ = wave time / cell activation time
 * @property {number} L_c        - 𝓛_c = cell threshold / relay threshold
 * @property {number} chi        - χ̃ chemotactic coupling
 * @property {number} mu         - μ̃ cell motility
 * @property {number} lam        - λ GL well depth
 * @property {number} tht        - ϑ noise strength
 * @property {number} n_L        - Hill exponent for H⁺(𝓛; 1; n_L)
 * @property {number} gamma_L    - Γ̃_L LTB4 degradation rate
 * @property {number} dt         - dt̃ agent SDE step
 * @property {number} seed       - RNG seed
 */

/**
 * Initialize agent state.
 * @param {number} N
 * @param {number} R_dish
 * @param {number} seed
 * @returns {AgentState}
 */
export function initAgents(N, R_dish, seed) {
  const rng = makeRng(seed);
  const x  = new Float32Array(N);
  const y  = new Float32Array(N);
  const Px = new Float32Array(N);
  const Py = new Float32Array(N);

  // Uniform distribution in the disk (rejection sampling).
  for (let i = 0; i < N; i++) {
    let rx, ry;
    do {
      rx = (rng.uniform() * 2 - 1) * R_dish;
      ry = (rng.uniform() * 2 - 1) * R_dish;
    } while (rx * rx + ry * ry > R_dish * R_dish);
    x[i] = rx;
    y[i] = ry;
    // Px[i] = Py[i] = 0 (default)
  }

  const emitting = new Uint8Array(N);
  // R̃_i starts at 0 for all cells (unsaturated). Used by M2 only; in M1 it
  // stays zero (no update path) and contributes no overhead.
  const R = new Float32Array(N);
  return { x, y, Px, Py, R, emitting, N };
}

/**
 * Step all agents one SDE step (Euler–Maruyama).
 * Worker must call field.accumulateSource() BEFORE calling this function
 * (for PIC consistency: emissions go onto the grid before the field advances,
 * then agents sense the updated field).
 *
 * @param {AgentState} agents
 * @param {Object} field - field API
 * @param {SimParams} params
 * @param {Object} rng - { gauss() }
 */
export function stepAgents(agents, field, params, rng) {
  const { x, y, Px, Py, R, emitting, N } = agents;
  const { Lambda, L_c, chi, mu, lam, tht, dt } = params;
  // M2 params (defaults make the R update a no-op when M1 is selected).
  const model   = params.model   || 'M1';
  const beta_R  = (model === 'M2') ? (params.beta_R  || 0) : 0;
  const gamma_R = (model === 'M2') ? (params.gamma_R || 0) : 0;
  const L_r_nd  = (model === 'M2') ? (params.L_r_nd  || 1) : 1;
  const n_Lr    = params.n_Lr    || 10;
  const n_R     = params.n_R     || 10;
  const n_L     = params.n_L     || 10;
  const r_fire  = params.r_fire  || 0;
  const r_fire2 = r_fire * r_fire;
  const R2 = params.R_dish * params.R_dish;

  // Explicit Euler–Maruyama stability for the GL drift requires
  //   dt_sub · Λ · |𝓛 − 𝓛_c| ≲ 0.3
  // In saturated regions where 𝓛 ≫ 𝓛_c (e.g. after wave nucleation in
  // a low-decay regime), Λ·(𝓛−𝓛_c) can reach 10³–10⁴ and one outer
  // dt = 0.01 step would blow P up to NaN within a handful of iterations.
  // We adaptively sub-step the SDE per cell so each cell stays stable.
  // L and ∇L are kept frozen across substeps (matches the worker's PIC
  // accumulation cadence — one L-sample per agent step).
  const SUB_CFL = 0.3;

  for (let i = 0; i < N; i++) {
    let xi = x[i], yi = y[i];

    const L_i = Math.max(0, field.sample(xi, yi));
    const { gx, gy } = field.sampleGrad(xi, yi);

    const drift_lin_mag = Lambda * Math.abs(L_i - L_c);
    // Cap n_sub defensively so a degenerate L (e.g. NaN propagated from
    // elsewhere) can't lock the worker in an unbounded inner loop.
    const n_sub  = Math.min(2000, Math.max(1, Math.ceil(drift_lin_mag * dt / SUB_CFL)));
    const dt_sub = dt / n_sub;
    // Noise: same total variance as one Euler step (variance scales with
    // dt_sub × n_sub = dt, so per-substep amplitude uses √dt_sub).
    const noise_sub = Math.sqrt(2 * Lambda * tht * L_i) * Math.sqrt(dt_sub);

    let px = Px[i], py = Py[i];
    for (let sub = 0; sub < n_sub; sub++) {
      const m2 = px * px + py * py;
      const m4 = m2 * m2;
      const drift_act = Lambda * ((L_i - L_c) + lam * (m2 - m4));
      px += (chi * gx + drift_act * px) * dt_sub + noise_sub * rng.gauss();
      py += (chi * gy + drift_act * py) * dt_sub + noise_sub * rng.gauss();
      // Guard against 0·∞ = NaN: with μ̃=0 we want positions frozen exactly,
      // but if the explicit Euler for P has overflowed to ±Infinity in a
      // saturated-L pocket (e.g. inside the firing disk with γ̃_L=0), `mu * px` would
      // be 0·∞ = NaN in JS — corrupting xi and making the cell invisible.
      if (mu !== 0) {
        xi += mu * px * dt_sub;
        yi += mu * py * dt_sub;
      }

      // Reflective BC per substep so a fast-moving cell can't shoot out.
      const rn2 = xi * xi + yi * yi;
      if (rn2 > R2) {
        const rn = Math.sqrt(rn2);
        const excess = rn - params.R_dish;
        const newR   = params.R_dish - excess;
        if (newR <= 0) {
          xi = (xi / rn) * params.R_dish;
          yi = (yi / rn) * params.R_dish;
        } else {
          xi = (xi / rn) * newR;
          yi = (yi / rn) * newR;
        }
      }
    }

    Px[i] = px;  Py[i] = py;
    x[i]  = xi;  y[i]  = yi;

    if (model === 'M2') {
      // M2 per-cell R̃_i ODE (explicit Euler, frozen L over the agent dt):
      //   dR̃_i/dt̃ = β̃ H⁺(𝓛(r̃_i); L̃_r; n_{Lr})  −  γ̃ R̃_i
      // L is sampled once per agent step (PIC cadence; consistent with the
      // SDE substeps using frozen-L). Explicit Euler is stable while
      // γ̃·dt < 2 — default γ̃ ≈ 1, dt = 0.01 → safe by 200×.
      const hLr = hillPos(L_i, L_r_nd, n_Lr);
      let r_new = R[i] + dt * (beta_R * hLr - gamma_R * R[i]);
      if (r_new < 0) r_new = 0;  // R is a concentration; floor at 0
      R[i] = r_new;

      // Cell is *emitting* iff its source weight is > 1/2. The source weight is
      //   w = H⁺(𝓛;1;n_L) · H⁻(R̃;1;n_R)            (outside r_fire)
      //   w =       1     · H⁻(R̃;1;n_R)            (inside  r_fire — forced)
      // Each Hill factor crosses 1/2 at its threshold, so w > 1/2 reduces to
      // every individual factor exceeding its 1/2 mark (matches the truth-table
      // convention in setup4_swarm3d.md §10).
      const inFire = (xi * xi + yi * yi) < r_fire2;
      const lOK = inFire ? true : (L_i > 1);
      const rOK = (R[i] < 1);
      emitting[i] = (lOK && rOK) ? 1 : 0;
    } else {
      // M1: cell is emitting iff sampled 𝓛 > 1.
      emitting[i] = (L_i > 1) ? 1 : 0;
    }
  }
}

/**
 * Compute Hill⁺ activation for source accumulation.
 * Worker calls this per-cell to get the emission weight before stepping field.
 * @param {number} L_val - current L at cell position
 * @param {number} n_L   - Hill exponent
 * @returns {number} H⁺(L_val; 1; n_L)
 */
export function hillEmission(L_val, n_L) {
  if (L_val <= 0) return 0;
  const xn  = Math.pow(L_val, n_L);
  const x0n = 1;  // x0 = 1 in nondim units
  return xn / (x0n + xn);
}

/**
 * General activating Hill: H⁺(x; x0; n) = x^n / (x0^n + x^n). → 1 as x → ∞.
 * Used by M2 for R-ODE production gate H⁺(𝓛; L̃_r; n_{Lr}).
 */
export function hillPos(x, x0, n) {
  if (x <= 0) return 0;
  const xn  = Math.pow(x, n);
  const x0n = Math.pow(x0, n);
  return xn / (x0n + xn);
}

/**
 * Inhibiting Hill: H⁻(x; x0; n) = x0^n / (x0^n + x^n). → 0 as x → ∞.
 * Used by M2 for R-inhibition gate H⁻(R̃; 1; n_R) on the L source.
 */
export function hillNeg(x, x0, n) {
  if (x <= 0) return 1;
  const xn  = Math.pow(x, n);
  const x0n = Math.pow(x0, n);
  return x0n / (x0n + xn);
}
