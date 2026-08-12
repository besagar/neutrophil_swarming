// setup4/solvers/solver_m1.js
// M1 — basic relay, no inhibitor.
//
// M1 nondim L equations (intrinsic-units scheme: ℓ_0 = a/(L_0 D_L), t_0 = ℓ_0²/D_L):
//   2D–2D:  ∂_t̃ 𝓛 = ∇̃²_{2D} 𝓛 + (1/h̃) Σ_i H⁺(𝓛_i;1;n_L) H⁻(ℛ_i;1;n_R) δ̃²(r̃−r̃_i) − Γ_L 𝓛
//           Per-cell prefactor 1/h̃ is baked into the worker's src[] accumulation
//           (weight = H⁺·H⁻·dt/h̃ per cell); the solver consumes src[] agnostically.
//   2D–3D:  ∂_t̃ 𝓛 = ∇̃²_{3D} 𝓛 + δ̃(z̃) Σ_i H⁺(𝓛_i;1;n_L) H⁻(ℛ_i;1;n_R) δ̃²(r̃−r̃_i) − Γ_L 𝓛
//           Per-cell prefactor is 1; the 2/h_0 in surfSrc comes purely from the
//           δ(z) control-volume discretization at the z=0 half-cell (see §1b), NOT
//           from any σ̃ scaling.
//
// Time integration:
//   2D–2D: explicit Euler with diffusion substepping.
//     CFL condition: dt̃_sub ≤ Δx̃²/4.  n_sub = ceil(dt̃_agent / (Δx̃²/4)).
//     Source accumulated into field._getSrc() before calling step(); applied on
//     first substep then cleared.
//
//   2D–3D: operator-split:
//     1. Explicit Euler for x/y diffusion (sub-stepped on Δx̃²/4 CFL).
//        Source at z=0 nodes with coefficient 2/h_0 (δ(z) discretization —
//        see setup4_cue_models.md §1b for 2/h_0 derivation).
//     2. Crank-Nicolson for z-diffusion (tridiagonal Thomas algorithm).
//        z-grid: non-uniform, Δz_j = h_0 · α^j. Neumann BCs at both ends.
//        System factorized inline per column call (N² columns per agent step;
//        Thomas is O(N_z) so total O(N² N_z) — fast enough for N=128, N_z=16).
//
// Hill function: H⁺(x; x0; n) = x^n / (x0^n + x^n).

// ─── 2D–2D stepper ─────────────────────────────────────────────────────────

/**
 * Create the 2D-2D M1 step function.
 * @param {Object} field - field API from createField()
 * @returns {Function} step(dt_agent, params) — advances L one agent step
 *
 * Assumes worker has already called field.accumulateSource() for all cells
 * before calling this function.  Source is in field._getSrc() and will be
 * consumed (zeroed) here.
 */
export function makeStepFn_2d2d(field) {
  // Pre-allocate scratch buffer (reused across calls; size N_grid²).
  let Lnew = new Float32Array(field.N_grid * field.N_grid);

  return function step2d2d(dt_agent, params) {
    const { gamma_L } = params;
    const N   = field.N_grid;
    const dx  = field.dx;
    const L   = field.getLfield();
    const src = field._getSrc();  // accumulated cell emissions (per nondim area, per nondim time)

    // Ensure scratch buffer size matches current N_grid (in case it changes between runs).
    if (Lnew.length !== N * N) Lnew = new Float32Array(N * N);

    // Number of diffusion substeps to satisfy CFL: dt_sub ≤ dx²/4.
    const dtCFL  = (dx * dx) / 4;
    const n_sub  = Math.ceil(dt_agent / dtCFL);
    const dt_sub = dt_agent / n_sub;

    for (let sub = 0; sub < n_sub; sub++) {

      for (let j = 0; j < N; j++) {
        for (let i = 0; i < N; i++) {
          const k  = j * N + i;
          const Lc = L[k];
          // 5-point Laplacian with Dirichlet zero at the boundary (out-of-range = 0).
          const Ln = (j > 0)     ? L[(j - 1) * N + i] : 0;
          const Ls = (j < N - 1) ? L[(j + 1) * N + i] : 0;
          const Lw = (i > 0)     ? L[j * N + i - 1]   : 0;
          const Le = (i < N - 1) ? L[j * N + i + 1]   : 0;
          const lap = (Ln + Ls + Lw + Le - 4 * Lc) / (dx * dx);

          // Source contribution:
          // src[k] = Σ_i w_ik * H⁺_i * H⁻_i * dt_agent / h̃  [worker-accumulated; PIC weights sum to 1 per cell].
          // To approximate δ²(r̃−r̃_i) on a 2D grid, the PIC density at node k is
          //   S_k = src[k] / dx²  (discrete 2D delta = PIC weight / dx²).
          // The 1/h̃ prefactor is already inside src[k] from the worker accumulation.
          // Apply all at once on sub=0 (total = src[k]/dx², regardless of n_sub).
          const directSrc = (sub === 0) ? src[k] / (dx * dx) : 0;
          Lnew[k] = Lc + dt_sub * (lap - gamma_L * Lc) + directSrc;
        }
      }

      // Clear source accumulator after first substep.
      if (sub === 0) src.fill(0);

      L.set(Lnew);
      field.applyBC();
    }
  };
}

// ─── 2D–3D z-grid and Crank-Nicolson ─────────────────────────────────────

/**
 * Build non-uniform z-grid.
 * Δz_j = h_0 · α^j,  j = 0..N_z-1.
 * Returns { dz: Float64Array[N_z] }
 */
function buildZgrid(N_z, h_0, alpha) {
  const dz = new Float64Array(N_z);
  for (let j = 0; j < N_z; j++) {
    dz[j] = h_0 * Math.pow(alpha, j);
  }
  return { dz };
}

/**
 * Build the finite-difference Laplacian coefficients for the non-uniform z-grid.
 * Second derivative at node j (interior):
 *   d²L/dz²|_j ≈ 2/(Δz_{j-1}+Δz_j) · [(L_{j+1}-L_j)/Δz_j - (L_j-L_{j-1})/Δz_{j-1}]
 * Neumann BCs (∂L/∂z = 0):
 *   j=0:      ghost node L[-1] = L[1]  →  d²L/dz²|_0 = 2*(L[1]-L[0])/Δz_0²
 *   j=N_z-1:  ghost node L[N_z] = L[N_z-2] → d²L/dz²|_{N_z-1} = 2*(L[N_z-2]-L[N_z-1])/Δz_{N_z-2}²
 * Returns { a, b, c } — sub/main/super diagonal arrays of the Laplacian operator.
 */
function buildZLaplacian(N_z, dz) {
  const a = new Float64Array(N_z);  // sub-diagonal (j-1 coefficient)
  const b = new Float64Array(N_z);  // main diagonal
  const c = new Float64Array(N_z);  // super-diagonal (j+1 coefficient)

  // j = 0 (Neumann ghost: L[-1] = L[1]).
  a[0] = 0;
  c[0] = 2 / (dz[0] * dz[0]);
  b[0] = -c[0];

  // Interior j = 1..N_z-2.
  for (let j = 1; j < N_z - 1; j++) {
    const dzL = dz[j - 1], dzR = dz[j];
    const denom = 0.5 * (dzL + dzR);
    a[j] =  1 / (dzL * denom);
    c[j] =  1 / (dzR * denom);
    b[j] = -(a[j] + c[j]);
  }

  // j = N_z-1 (Neumann ghost: L[N_z] = L[N_z-2]).
  {
    const j = N_z - 1;
    const dzL = dz[j - 1];
    a[j] = 2 / (dzL * dzL);
    c[j] = 0;
    b[j] = -a[j];
  }

  return { a, b, c };
}

// Thomas algorithm is inlined in makeStepFn_2d3d to avoid per-column allocations.

// ─── 2D–3D slab stepper (generalized reaction–diffusion) ────────────────────

/**
 * Create a 2D-3D thin-slab reaction–diffusion stepper for one field:
 *   ∂_t̃ F = D ∇̃²_{3D} F − decay·F + δ̃(z̃)·(surface source)
 * solved on a direct 3D grid (explicit x/y diffusion + Crank–Nicolson in z).
 * The z=0 slice is mirrored into `field2d` so agents can sample it, and the
 * surface source is read from `field2d._getSrc()` (worker-accumulated PIC
 * weights × dt, per-cell prefactor 1; the 2/h_0 in the flux is the δ(z)
 * control-volume discretization at the z=0 half-cell).
 *
 * Generalized over M1's original L stepper by two per-step accessors so the
 * SAME code drives both the LTB4 field (D=1, decay=Γ_L) and M6.1's adenosine
 * field (D=D_A/D_L, decay=Γ_A) — no duplicated 3D solver.
 *
 * @param {Object} field2d - 2D field API (z=0 slice, for agent I/O + src)
 * @param {Object} cfg - { N_z, h_0, alpha, getD?, getDecay? }
 *   getD(params)→diffusion coefficient (default 1), getDecay(params)→decay
 *   (default reads params.gamma_L). Both are constant within a run.
 * @returns {{ step, reset, getField }}
 */
export function makeSlabStepper(field2d, cfg) {
  const { N_z = 16, h_0 = 0.1, alpha = 1.4,
          getD = () => 1, getDecay = (p) => p.gamma_L || 0 } = cfg;
  const N      = field2d.N_grid;
  const dx     = field2d.dx;
  const R_dish = field2d.R_dish;
  const Fnew_xy = new Float32Array(N * N);

  const { dz } = buildZgrid(N_z, h_0, alpha);
  const { a: za, b: zb, c: zc } = buildZLaplacian(N_z, dz);

  // Full 3D array: F3d[kz * N * N + ky * N + kx], kz=0 is z=0 floor.
  const F3d = new Float32Array(N_z * N * N);

  const _col = new Float64Array(N_z);
  const _al  = new Float64Array(N_z);
  const _bl  = new Float64Array(N_z);
  const _cl  = new Float64Array(N_z);
  const _rhs = new Float64Array(N_z);
  const _c2  = new Float64Array(N_z);
  const _d2  = new Float64Array(N_z);
  const _sol = new Float64Array(N_z);
  // LHS depends on (dt, D); both constant within a run → built once.
  let _lhsDt = -1, _lhsD = -1;

  const mask = new Uint8Array(N * N);
  for (let j = 0; j < N; j++) {
    for (let i = 0; i < N; i++) {
      const x = -R_dish + i * dx;
      const y = -R_dish + j * dx;
      mask[j * N + i] = (x * x + y * y <= R_dish * R_dish + 1e-9) ? 1 : 0;
    }
  }

  function syncZ0toField2d() {
    const F2d = field2d.getLfield();
    for (let k = 0; k < N * N; k++) F2d[k] = F3d[k];  // kz=0 slice
  }

  function reset() {
    F3d.fill(0);
    syncZ0toField2d();
  }

  /**
   * Seed the slab with a prescribed z-profile: F(x,y,z) = profile(z̃) inside
   * the dish mask, 0 outside. Node j sits at z̃_j = Σ_{k<j} Δz_k (the same
   * cumulative non-uniform grid used by the Laplacian).
   * Used to start an auxiliary field at its screened steady tone instead of 0
   * (see setup4_m6_2_implementation_plan.md §5).
   */
  function seed(profile) {
    let z = 0;
    for (let kz = 0; kz < N_z; kz++) {
      const v = profile(z);
      const base = kz * N * N;
      for (let k = 0; k < N * N; k++) F3d[base + k] = mask[k] ? v : 0;
      z += dz[kz];
    }
    syncZ0toField2d();
  }

  function step(dt_agent, params) {
    const D     = Math.max(getD(params), 1e-12);
    const decay = getDecay(params);
    const src   = field2d._getSrc();

    // Step 1: explicit x/y diffusion + decay (sub-stepped on the D-scaled CFL).
    const dtCFL  = (dx * dx) / (4 * D);
    const n_sub  = Math.max(1, Math.ceil(dt_agent / dtCFL));
    const dt_sub = dt_agent / n_sub;

    for (let sub = 0; sub < n_sub; sub++) {
      for (let kz = 0; kz < N_z; kz++) {
        const base = kz * N * N;
        for (let ky = 0; ky < N; ky++) {
          for (let kx = 0; kx < N; kx++) {
            const k  = ky * N + kx;
            const Fc = F3d[base + k];
            const Fn = (ky > 0)     ? F3d[base + (ky - 1) * N + kx] : 0;
            const Fs = (ky < N - 1) ? F3d[base + (ky + 1) * N + kx] : 0;
            const Fw = (kx > 0)     ? F3d[base + ky * N + kx - 1]   : 0;
            const Fe = (kx < N - 1) ? F3d[base + ky * N + kx + 1]   : 0;
            const lap = (Fn + Fs + Fw + Fe - 4 * Fc) / (dx * dx);
            Fnew_xy[k] = Fc + dt_sub * (D * lap - decay * Fc);
          }
        }
        for (let k = 0; k < N * N; k++) {
          F3d[base + k] = mask[k] ? Fnew_xy[k] : 0;
        }
      }
    }

    // Step 2: Crank–Nicolson in z (operator scaled by D). Surface source at z=0.
    const col = _col, al = _al, bl = _bl, cl = _cl;
    const rhs = _rhs, c2 = _c2, d2 = _d2, sol = _sol;

    const h = dt_agent / 2;
    if (_lhsDt !== dt_agent || _lhsD !== D) {
      for (let j = 0; j < N_z; j++) {
        al[j] = -h * D * za[j];
        bl[j] =  1 - h * D * zb[j];
        cl[j] = -h * D * zc[j];
      }
      _lhsDt = dt_agent; _lhsD = D;
    }

    for (let ky = 0; ky < N; ky++) {
      for (let kx = 0; kx < N; kx++) {
        if (!mask[ky * N + kx]) continue;
        const k = ky * N + kx;
        for (let kz = 0; kz < N_z; kz++) col[kz] = F3d[kz * N * N + k];
        // surfSrc = 2·src[k]/(dx²·h_0): the 2/h_0 is the δ(z) discretization at
        // the z=0 half-cell (D-independent — it's a control-volume factor).
        const surfSrc = 2 * src[k] / (dx * dx * h_0);
        for (let j = 0; j < N_z; j++) {
          const jm = j > 0 ? j - 1 : 0;
          const jp = j < N_z - 1 ? j + 1 : N_z - 2;
          rhs[j] = col[j] + h * D * (za[j] * col[jm] + zb[j] * col[j] + zc[j] * col[jp]);
        }
        rhs[0] += surfSrc;
        c2[0] = cl[0] / bl[0];
        d2[0] = rhs[0] / bl[0];
        for (let i = 1; i < N_z; i++) {
          const denom = bl[i] - al[i] * c2[i - 1];
          c2[i] = cl[i] / denom;
          d2[i] = (rhs[i] - al[i] * d2[i - 1]) / denom;
        }
        sol[N_z - 1] = d2[N_z - 1];
        for (let i = N_z - 2; i >= 0; i--) sol[i] = d2[i] - c2[i] * sol[i + 1];
        for (let kz = 0; kz < N_z; kz++) F3d[kz * N * N + k] = Math.max(0, sol[kz]);
      }
    }

    src.fill(0);
    syncZ0toField2d();
  }

  return { step, reset, seed, getField: () => F3d };
}

/**
 * M1 2D-3D stepper — the LTB4 field slab (D=1, decay=Γ_L). Thin wrapper over
 * makeSlabStepper keeping the legacy {step, resetL3d, getL3d} interface.
 * @param {Object} field2d
 * @param {Object} params3d - { N_z, h_0, alpha }
 * @returns {{ step, resetL3d, getL3d }}
 */
export function makeStepFn_2d3d(field2d, params3d) {
  const s = makeSlabStepper(field2d, {
    ...params3d, getD: () => 1, getDecay: (p) => p.gamma_L || 0,
  });
  return { step: s.step, resetL3d: s.reset, getL3d: s.getField };
}
