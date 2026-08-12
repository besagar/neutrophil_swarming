// setup4/agents.js
// Per-cell state and Euler–Maruyama SDE step for Setup 4.
//
// GL polarization SDE (Euler–Maruyama, nondim — intrinsic-units scheme):
//   dP_xi = [κ ∂_x̃ 𝓛 + λ(𝓛 − 𝓛_c) P_xi + ν(|P|²−|P|⁴) P_xi] dt̃
//           + √(2 ϑ 𝓛 dt̃) η_xi
//   dP_yi = [κ ∂_ỹ 𝓛 + λ(𝓛 − 𝓛_c) P_yi + ν(|P|²−|P|⁴) P_yi] dt̃
//           + √(2 ϑ 𝓛 dt̃) η_yi
//   dx̃_i/dt̃ = μ P_xi
//   dỹ_i/dt̃ = μ P_yi
//
// Note: noise prefactor is √(2 ϑ 𝓛 dt̃) — NO Λ factor. The three drift
// prefactors κ, λ, ν are independent (σ-free) nondim groups derived from
// intrinsic cell scales ℓ_0 = a/(L_0 D_L), t_0 = ℓ_0²/D_L.
// See docs/physics/setup4_intrinsic_units_implementation_plan.md §1 and §2 (A3).
//
// Boundary: reflective at r̃ = R̃_dish (hard outer boundary).
// Optional inner boundary ("stick to the target"): a circle of radius R̃_target
// at the origin stands for the pathogen (candida cluster / sterile injury).
// A cell that reaches it stops moving forever — it has engaged the target —
// while its polarization and its LTB4 emission keep running. See
// docs/physics/setup4_swarm3d.md §2c.
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
 * @property {Uint8Array}  stuck    - 1 if cell is stuck on the target circle [N]
 * @property {number} N
 */

/**
 * @typedef {Object} SimParams
 * @property {number} N      - number of cells
 * @property {number} R_dish - dish radius (nondim)
 * @property {number} lam    - λ = r_0·t_0·L_0  activation/deactivation rate (multiplies (𝓛−𝓛_c)·P)
 * @property {number} nu     - ν = u²·t_0/w      GL nonlinearity depth (multiplies (|P|²−|P|⁴)·P)
 * @property {number} kap    - κ = a·χ/(D_L²·p_0) chemotactic gradient coupling (multiplies ∇𝓛)
 * @property {number} mu     - μ = μ_dim·t_0/ℓ_0  cell motility (position-update coefficient)
 * @property {number} tht    - ϑ = (w·L_0·t_0/u)·θ noise strength (no Λ factor)
 * @property {number} L_c    - 𝓛_c = L_c/L_0     nondim polarization threshold
 * @property {number} n_L    - Hill exponent for H⁺(𝓛; 1; n_L)
 * @property {number} gamma_L - Γ_L LTB4 degradation rate (nondim)
 * @property {number} dt     - dt̃ agent SDE step
 * @property {number} seed   - RNG seed
 */

/**
 * Initialize agent state.
 * @param {number} N
 * @param {number} R_dish
 * @param {number} seed
 * @param {number} [R_target=0] - target radius; when > 0 the target disk is
 *   excluded from the initial placement (it is occupied by the pathogen), so
 *   cells start uniformly in the annulus R̃_target < r̃ < R̃_dish. N is NOT
 *   changed, so σ̃ = N/(πR̃_dish²) still names the same knob; the accessible
 *   area is smaller by (R̃_target/R̃_dish)², negligible at the defaults.
 * @returns {AgentState}
 */
export function initAgents(N, R_dish, seed, R_target = 0) {
  const rng = makeRng(seed);
  const x  = new Float32Array(N);
  const y  = new Float32Array(N);
  const Px = new Float32Array(N);
  const Py = new Float32Array(N);

  // Uniform distribution in the disk (rejection sampling), minus the target
  // disk when one is present. Guard R_in against a target that swallows the
  // dish, or the rejection loop would never terminate.
  const R_in  = Math.max(0, Math.min(R_target, 0.95 * R_dish));
  const R_in2 = R_in * R_in;
  for (let i = 0; i < N; i++) {
    let rx, ry, r2;
    do {
      rx = (rng.uniform() * 2 - 1) * R_dish;
      ry = (rng.uniform() * 2 - 1) * R_dish;
      r2 = rx * rx + ry * ry;
    } while (r2 > R_dish * R_dish || r2 < R_in2);
    x[i] = rx;
    y[i] = ry;
    // Px[i] = Py[i] = 0 (default)
  }

  const emitting = new Uint8Array(N);
  const stuck    = new Uint8Array(N);
  // R̃_i starts at 0 for all cells (unsaturated). Used by M2 only; in M1 it
  // stays zero (no update path) and contributes no overhead.
  const R = new Float32Array(N);
  return { x, y, Px, Py, R, emitting, stuck, N };
}

/**
 * Step all agents one SDE step (Euler–Maruyama).
 * Worker must call field.accumulateSource() BEFORE calling this function
 * (for PIC consistency: emissions go onto the grid before the field advances,
 * then agents sense the updated field).
 *
 * @param {AgentState} agents
 * @param {Object} field - field API (L field)
 * @param {SimParams} params
 * @param {Object} rng - { gauss() }
 * @param {Object} [fieldAux] - auxiliary basal field API (M6.1 𝓐 / M6.2 𝓠).
 *   Only M6.1 reads it here (threshold-shifted `emitting` flag); M6.2's
 *   throttle acts on the source weight, not on the flag. Ignored otherwise.
 */
export function stepAgents(agents, field, params, rng, fieldAux) {
  const { x, y, Px, Py, R, emitting, N } = agents;
  // `stuck` is optional so an AgentState built by older code still steps.
  const stuck = agents.stuck || (agents.stuck = new Uint8Array(N));
  const { lam, nu, kap, mu, tht, L_c, dt } = params;
  // M2 params (defaults make the R update a no-op when M1 is selected).
  const model   = params.model   || 'M1';
  // M6.1 threshold-shift coupling (0 makes M6.1 collapse to M1 emission logic).
  const lam_A   = (model === 'M6.1') ? (params.lam_A || 0) : 0;
  const beta_R  = (model === 'M2') ? (params.beta_R  || 0) : 0;
  const gamma_R = (model === 'M2') ? (params.gamma_R || 0) : 0;
  const L_r_nd  = (model === 'M2') ? (params.L_r_nd  || 1) : 1;
  const n_Lr    = params.n_Lr    || 10;
  const n_R     = params.n_R     || 10;
  const n_L     = params.n_L     || 10;
  const r_fire  = params.r_fire  || 0;
  const r_fire2 = r_fire * r_fire;
  const R2 = params.R_dish * params.R_dish;
  // Inner "stick to the target" boundary (absorbing, not reflecting).
  const stick  = !!params.stick_target && params.R_target > 0;
  const R_t    = stick ? Math.min(params.R_target, 0.95 * params.R_dish) : 0;
  const R_t2   = R_t * R_t;

  // Explicit Euler–Maruyama stability for the GL drift requires
  //   dt_sub · max(λ·|𝓛 − 𝓛_c|, ν) ≲ 0.3
  // Both lam and nu can independently drive blowup: lam when 𝓛 ≫ 𝓛_c
  // (e.g. after wave nucleation in a low-decay regime) and nu when the
  // GL nonlinearity dominates at large |P|. We take the max of both
  // contributions to set n_sub, ensuring each substep stays stable for
  // the full drift. L and ∇L are kept frozen across substeps (matches
  // the worker's PIC accumulation cadence — one L-sample per agent step).
  const SUB_CFL = 0.3;

  for (let i = 0; i < N; i++) {
    let xi = x[i], yi = y[i];

    const L_i = Math.max(0, field.sample(xi, yi));
    const { gx, gy } = field.sampleGrad(xi, yi);
    // M6.1: adenosine at this cell (pre-move position), for the shifted threshold.
    const A_i = (model === 'M6.1' && fieldAux) ? Math.max(0, fieldAux.sample(xi, yi)) : 0;

    const drift_lin_mag = lam * Math.abs(L_i - L_c);
    const drift_nl_mag  = nu;   // |P|²+|P|⁴ bound by ~2 at saturation; ν is the prefactor
    const drift_mag     = Math.max(drift_lin_mag, drift_nl_mag);
    // Cap n_sub defensively so a degenerate L (e.g. NaN propagated from
    // elsewhere) can't lock the worker in an unbounded inner loop.
    const n_sub  = Math.min(2000, Math.max(1, Math.ceil(drift_mag * dt / SUB_CFL)));
    const dt_sub = dt / n_sub;
    // Noise: same total variance as one Euler step (variance scales with
    // dt_sub × n_sub = dt, so per-substep amplitude uses √dt_sub).
    // Intrinsic-units scheme: NO Λ factor in the noise prefactor.
    const noise_sub = Math.sqrt(2 * tht * L_i) * Math.sqrt(dt_sub);

    // A cell already engaged with the target keeps polarizing and emitting,
    // but never moves again (adhesion to the pathogen, not a reflection).
    let frozen = stick && stuck[i] === 1;

    let px = Px[i], py = Py[i];
    for (let sub = 0; sub < n_sub; sub++) {
      const m2 = px * px + py * py;
      const m4 = m2 * m2;
      // drift_act is the scalar multiplying P: λ·(𝓛−𝓛_c) from activation +
      // ν·(|P|²−|P|⁴) from GL nonlinearity. Both terms then multiply px/py.
      const drift_act = lam * (L_i - L_c) + nu * (m2 - m4);
      px += (kap * gx + drift_act * px) * dt_sub + noise_sub * rng.gauss();
      py += (kap * gy + drift_act * py) * dt_sub + noise_sub * rng.gauss();
      // Guard against 0·∞ = NaN: with μ̃=0 we want positions frozen exactly,
      // but if the explicit Euler for P has overflowed to ±Infinity in a
      // saturated-L pocket (e.g. inside the firing disk with γ̃_L=0), `mu * px` would
      // be 0·∞ = NaN in JS — corrupting xi and making the cell invisible.
      if (mu !== 0 && !frozen) {
        xi += mu * px * dt_sub;
        yi += mu * py * dt_sub;
      }

      // Inner absorbing boundary: first crossing of r̃ = R̃_target pins the
      // cell onto the circle for the rest of the run.
      if (stick && !frozen) {
        const ri2 = xi * xi + yi * yi;
        if (ri2 < R_t2) {
          const ri = Math.sqrt(ri2);
          if (ri > 1e-12) { xi = (xi / ri) * R_t; yi = (yi / ri) * R_t; }
          else            { xi = R_t; yi = 0; }
          frozen = true;
        }
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
    stuck[i] = frozen ? 1 : 0;

    if (model === 'M2') {
      // M2 per-cell R̃_i ODE (explicit Euler, frozen L over the agent dt):
      //   dℛ_i/dt̃ = β H⁺(𝓛(r̃_i); 𝓛_r; n_{Lr})  −  Γ_R ℛ_i
      // where β = b·t_0/R_0 and Γ_R = γ_R·t_0 (intrinsic-units scheme;
      // t_0 = a²/(L_0² D_L³)). L is sampled once per agent step (PIC
      // cadence; consistent with the SDE substeps using frozen-L).
      // Explicit Euler stable while Γ_R·dt < 2; verify against current defaults.
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
    } else if (model === 'M6.1') {
      // M6.1: cell emits iff 𝓛 exceeds the adenosine-shifted relay threshold.
      // Gate H⁺(𝓛; 1 + λ_A 𝓐; n_L) crosses 1/2 at 𝓛 = 1 + λ_A 𝓐.
      emitting[i] = (L_i > 1 + lam_A * A_i) ? 1 : 0;
    } else if (model === 'M6.2') {
      // M6.2: DELIBERATE departure from the "source weight w > 1/2" truth-table
      // convention used by M2/M6.1. Here w = H⁺(𝓛;1;n_L)·H⁻(𝓠;1;m), so w > 1/2
      // would also require 𝓠 < 1 — but above the crossover σ̃★ = γ/β the steady
      // tone 𝓠_ss > 1 EVERYWHERE, which would flag every cell non-emitting and
      // make measureWaveSpeed() report c = 0 at exactly the densities this model
      // exists to describe. The relay front is defined by the L gate alone;
      // "throttled but still firing" is conveyed by emission opacity in the
      // renderer (frame.agentR carries 𝓠). See setup4_m6_2_implementation_plan §4.
      emitting[i] = (L_i > 1) ? 1 : 0;
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
