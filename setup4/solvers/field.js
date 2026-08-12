// setup4/solvers/field.js
// Grid data structure and PIC (particle-in-cell) operations for the L field.
//
// API (the ONLY interface agents and worker may use):
//   sample(x, y)              → bilinear interp of L at (x,y) in nondim coords
//   sampleGrad(x, y)          → {gx, gy} — finite-diff gradient of L at (x,y)
//   accumulateSource(x, y, w) → smear w units onto grid using same bilinear weights
//   step(dt, params)          → advance PDE one substep (called by solver)
//   getRadialProfile()        → Float32Array length N_grid: azimuthal avg of L vs r̃
//   reset(ic)                 → reinitialize L from IC object
//   getLfield()               → Float32Array of raw L values (for heatmap; read-only)
//
// PIC consistency note: sample() and accumulateSource() use identical bilinear
// weights so that the discrete emission operator is the adjoint of the
// interpolation operator. This prevents systematic flux imbalances.
//
// Geometry: the grid covers [-R_dish, R_dish]² in both x and y.
// Dirichlet BC: L = 0 at nodes outside r > R_dish (circular mask applied each step).

/**
 * @typedef {Object} FieldState
 * @property {Float32Array} L   - flat [N_grid × N_grid] nondim cue values
 * @property {Float32Array} src - flat accumulator for source smearing (cleared each step)
 * @property {number} N_grid
 * @property {number} dx       - grid spacing in nondim units (= 2 R_dish / (N_grid-1))
 * @property {number} R_dish
 */

export function createField(N_grid, R_dish) {
  if (!N_grid || N_grid < 4) throw new Error('field: N_grid must be >= 4');
  const N = N_grid;
  // Grid: N×N nodes covering [-R_dish, R_dish]²
  const dx = (2 * R_dish) / (N - 1);
  const L   = new Float32Array(N * N);
  const src = new Float32Array(N * N);  // source accumulator, zeroed before each PDE substep

  // Precompute circular mask: nodes outside r>R_dish are clamped to 0 (Dirichlet BC).
  const mask = new Uint8Array(N * N);
  for (let j = 0; j < N; j++) {
    for (let i = 0; i < N; i++) {
      const x = -R_dish + i * dx;
      const y = -R_dish + j * dx;
      mask[j * N + i] = (x * x + y * y <= R_dish * R_dish + 1e-9) ? 1 : 0;
    }
  }

  // Convert world coords to fractional grid indices.
  function toGrid(x, y) {
    const fi = (x + R_dish) / dx;
    const fj = (y + R_dish) / dx;
    return { fi, fj };
  }

  // Bilinear weights for fractional index (fi, fj).
  // Returns { i0, j0, w00, w10, w01, w11 } clamped to [0, N-1].
  function bilinear(fi, fj) {
    const i0 = Math.floor(fi), j0 = Math.floor(fj);
    const i1 = i0 + 1, j1 = j0 + 1;
    const tx = fi - i0, ty = fj - j0;
    // Clamp indices to valid range.
    const ci0 = Math.max(0, Math.min(N - 1, i0));
    const ci1 = Math.max(0, Math.min(N - 1, i1));
    const cj0 = Math.max(0, Math.min(N - 1, j0));
    const cj1 = Math.max(0, Math.min(N - 1, j1));
    return {
      ci0, ci1, cj0, cj1,
      w00: (1 - tx) * (1 - ty),
      w10: tx       * (1 - ty),
      w01: (1 - tx) * ty,
      w11: tx       * ty,
    };
  }

  /** Sample L at world position (x, y) using bilinear interpolation. */
  function sample(x, y) {
    const { fi, fj } = toGrid(x, y);
    const { ci0, ci1, cj0, cj1, w00, w10, w01, w11 } = bilinear(fi, fj);
    return (
      w00 * L[cj0 * N + ci0] +
      w10 * L[cj0 * N + ci1] +
      w01 * L[cj1 * N + ci0] +
      w11 * L[cj1 * N + ci1]
    );
  }

  /** Sample gradient of L at (x, y) using finite differences + bilinear interp. */
  function sampleGrad(x, y) {
    // Central difference with spacing dx; at boundaries, one-sided.
    const { fi, fj } = toGrid(x, y);
    const { ci0, ci1, cj0, cj1, w00, w10, w01, w11 } = bilinear(fi, fj);

    // Helper: get L at grid node (gi, gj) — clamped.
    function Lij(gi, gj) {
      const ii = Math.max(0, Math.min(N - 1, gi));
      const jj = Math.max(0, Math.min(N - 1, gj));
      return L[jj * N + ii];
    }

    // Gradient at each of the four surrounding nodes (central diff).
    function gx_at(gi, gj) {
      const lo = Math.max(0, gi - 1), hi = Math.min(N - 1, gi + 1);
      return (Lij(hi, gj) - Lij(lo, gj)) / ((hi - lo) * dx);
    }
    function gy_at(gi, gj) {
      const lo = Math.max(0, gj - 1), hi = Math.min(N - 1, gj + 1);
      return (Lij(gi, hi) - Lij(gi, lo)) / ((hi - lo) * dx);
    }

    const gx = w00 * gx_at(ci0, cj0) + w10 * gx_at(ci1, cj0) +
               w01 * gx_at(ci0, cj1) + w11 * gx_at(ci1, cj1);
    const gy = w00 * gy_at(ci0, cj0) + w10 * gy_at(ci1, cj0) +
               w01 * gy_at(ci0, cj1) + w11 * gy_at(ci1, cj1);
    return { gx, gy };
  }

  /**
   * Add a uniform contribution `amount` to src[k] for every grid node inside
   * a centred disk of radius `radius` (and inside the dish). Used by the
   * time-limited firing source IC.
   *
   * Under the intrinsic-units scheme the worker passes:
   *   amount = s_fire · σ̃ · dx² · dt
   * where s_fire scales the per-cell emission rate inside the firing disk
   * (s_fire = 1 means every point in the disk emits as strongly as one
   * saturated cell), and σ̃ = σ · ℓ_0² is the nondim cell density.
   *
   * The solver then recovers the per-step ΔL at each firing node as:
   *   2D-2D:  ΔL/dt = s_fire · σ̃         (effective bulk source of σ̃ saturated cells per area;
   *                                         solver_m1 computes directSrc = src[k]/dx²)
   *   2D-3D:  ΔL_{z=0}/dt = 2·s_fire·σ̃/h_0  (surface flux at z=0; solver_m1 computes
   *                                              surfSrc = 2·src[k]/(dx²·h_0))
   *
   * The geometric factor 2/h_0 in 2D-3D is the δ(z) discretization at the z=0 half-cell,
   * not an additional σ̃ scaling. The σ̃ factor comes entirely from the worker call-site.
   */
  function addFiringSource(radius, amount) {
    const r2 = radius * radius;
    for (let j = 0; j < N; j++) {
      for (let i = 0; i < N; i++) {
        if (!mask[j * N + i]) continue;
        const x = -R_dish + i * dx;
        const y = -R_dish + j * dx;
        if (x * x + y * y <= r2) {
          src[j * N + i] += amount;
        }
      }
    }
  }

  /**
   * Smear w units of emission from world position (x, y) onto the grid.
   * Uses the same bilinear weights as sample() — PIC consistency.
   * Accumulates into `src` (not L directly); caller calls applySource() next.
   */
  function accumulateSource(x, y, w) {
    const { fi, fj } = toGrid(x, y);
    const { ci0, ci1, cj0, cj1, w00, w10, w01, w11 } = bilinear(fi, fj);
    src[cj0 * N + ci0] += w00 * w;
    src[cj0 * N + ci1] += w10 * w;
    src[cj1 * N + ci0] += w01 * w;
    src[cj1 * N + ci1] += w11 * w;
  }

  /** Apply accumulated source to L and clear the accumulator. */
  function applySource() {
    for (let k = 0; k < N * N; k++) {
      L[k] += src[k];
      src[k] = 0;
    }
  }

  /** Apply Dirichlet BC: zero out nodes outside the dish. */
  function applyBC() {
    for (let k = 0; k < N * N; k++) {
      if (!mask[k]) L[k] = 0;
    }
  }

  /**
   * Azimuthal average of L vs r̃.
   * Returns Float32Array of length N_grid; index k → r̃ = k * dx.
   * Fills each radial bin with the average L over nodes in that annulus.
   */
  function getRadialProfile() {
    const nBins = N;
    const prof  = new Float64Array(nBins);
    const cnt   = new Float32Array(nBins);
    for (let j = 0; j < N; j++) {
      for (let i = 0; i < N; i++) {
        const x = -R_dish + i * dx;
        const y = -R_dish + j * dx;
        const r = Math.sqrt(x * x + y * y);
        const bin = Math.floor(r / dx);
        if (bin < nBins) {
          prof[bin] += L[j * N + i];
          cnt[bin]++;
        }
      }
    }
    const result = new Float32Array(nBins);
    for (let b = 0; b < nBins; b++) {
      result[b] = cnt[b] > 0 ? prof[b] / cnt[b] : 0;
    }
    return result;
  }

  /** Reset L and src to zero. The wave is launched by a time-limited
   * firing source in the worker time loop, not by a static IC. */
  function reset() {
    L.fill(0);
    src.fill(0);
  }

  /** Return the raw L array (read-only reference for heatmap rendering). */
  function getLfield() { return L; }

  /**
   * Set every in-dish node to `v` (out-of-dish nodes stay 0, Dirichlet).
   * Used to pre-seed an auxiliary field at its mean-field steady tone so a
   * run does not spend t̃ ~ 1/Γ filling it up (M6.x quorum fields; see
   * setup4_m6_2_implementation_plan.md §5).
   */
  function fillUniform(v) {
    for (let k = 0; k < N * N; k++) L[k] = mask[k] ? v : 0;
  }

  /** Internal accessor for solvers (src accumulator). Not for agents/worker. */
  function _getSrc() { return src; }

  return {
    N_grid: N,
    R_dish,
    dx,
    sample,
    sampleGrad,
    accumulateSource,
    addFiringSource,
    applySource,
    applyBC,
    getRadialProfile,
    reset,
    getLfield,
    fillUniform,
    _getSrc,
  };
}
