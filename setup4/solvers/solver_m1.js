// setup4/solvers/solver_m1.js
// M1 — basic relay, no inhibitor.
//
// M1 nondim L equations (after ℓ_0 = D_L/c* choice):
//   2D–2D:  ∂_t̃ 𝓛 = ∇̃²_{2D} 𝓛 + S̃(x̃,ỹ) − Γ̃_L 𝓛
//           S̃ = Σ_i H⁺(𝓛_i; 1; n_L) δ(x̃-x̃_i)  (smeared onto grid via PIC)
//   2D–3D:  ∂_t̃ 𝓛 = ∇̃²_{3D} 𝓛 + S̃(x̃,ỹ) (2/h_0) δ(z̃) − Γ̃_L 𝓛
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
          // src[k] = Σ_i w_ik * H⁺(𝓛_i; 1; n_L) * dt_agent  [PIC weights sum to 1 per cell].
          // To approximate δ(x-x_i) on a 2D grid, the PIC density at node k is
          //   S_k = Σ_i w_ik * H⁺_i / dx²  (discrete delta = w/dx² in 2D).
          // Concentration change over dt_agent = S_k * dt_agent = src[k] / dx².
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

// ─── 2D–3D stepper ─────────────────────────────────────────────────────────

/**
 * Create the 2D-3D M1 step function.
 * @param {Object} field2d - 2D field API (z=0 slice, for agent I/O)
 * @param {Object} params3d - { N_z, h_0, alpha } fixed grid params (per-run)
 * @returns {{ step, resetL3d, getL3d }}
 */
export function makeStepFn_2d3d(field2d, params3d) {
  const { N_z = 16, h_0 = 0.1, alpha = 1.4 } = params3d;
  const N      = field2d.N_grid;
  const dx     = field2d.dx;
  const R_dish = field2d.R_dish;
  // Reusable scratch buffer for x/y explicit substeps.
  const Lnew_xy = new Float32Array(N * N);

  // Build z-grid.
  const { dz } = buildZgrid(N_z, h_0, alpha);
  const { a: za, b: zb, c: zc } = buildZLaplacian(N_z, dz);

  // Full 3D array: L3d[kz * N * N + ky * N + kx], kz=0 is z=0 floor.
  const L3d = new Float32Array(N_z * N * N);

  // Pre-allocate C-N scratch buffers at closure level (reused every step, no GC).
  const _col = new Float64Array(N_z);
  const _al  = new Float64Array(N_z);
  const _bl  = new Float64Array(N_z);
  const _cl  = new Float64Array(N_z);
  const _rhs = new Float64Array(N_z);
  const _c2  = new Float64Array(N_z);
  const _d2  = new Float64Array(N_z);
  const _sol = new Float64Array(N_z);
  // LHS coefficients depend on dt_agent (set on first call, cached for run).
  let _lhsDt = -1;  // sentinel: LHS not yet built

  // Circular mask for x/y Dirichlet BC.
  const mask = new Uint8Array(N * N);
  for (let j = 0; j < N; j++) {
    for (let i = 0; i < N; i++) {
      const x = -R_dish + i * dx;
      const y = -R_dish + j * dx;
      mask[j * N + i] = (x * x + y * y <= R_dish * R_dish + 1e-9) ? 1 : 0;
    }
  }

  /** Sync z=0 slice of L3d → field2d's L array (so agents can sample). */
  function syncZ0toField2d() {
    const L2d = field2d.getLfield();
    for (let k = 0; k < N * N; k++) L2d[k] = L3d[k];  // kz=0 slice
  }

  function resetL3d() {
    L3d.fill(0);
    syncZ0toField2d();
  }

  function step2d3d(dt_agent, params) {
    const { gamma_L } = params;
    const src = field2d._getSrc();  // cell-emission accumulator (accumulated by worker this step)

    // Step 1: explicit x/y diffusion (sub-stepped on CFL). NO source here.
    // Source enters in step 2 as a non-homogeneous Neumann flux at z=0.
    const dtCFL  = (dx * dx) / 4;
    const n_sub  = Math.ceil(dt_agent / dtCFL);
    const dt_sub = dt_agent / n_sub;

    for (let sub = 0; sub < n_sub; sub++) {
      for (let kz = 0; kz < N_z; kz++) {
        const base = kz * N * N;
        for (let ky = 0; ky < N; ky++) {
          for (let kx = 0; kx < N; kx++) {
            const k  = ky * N + kx;
            const Lc = L3d[base + k];
            const Ln = (ky > 0)     ? L3d[base + (ky - 1) * N + kx] : 0;
            const Ls = (ky < N - 1) ? L3d[base + (ky + 1) * N + kx] : 0;
            const Lw = (kx > 0)     ? L3d[base + ky * N + kx - 1]   : 0;
            const Le = (kx < N - 1) ? L3d[base + ky * N + kx + 1]   : 0;
            const lap = (Ln + Ls + Lw + Le - 4 * Lc) / (dx * dx);
            Lnew_xy[k] = Lc + dt_sub * (lap - gamma_L * Lc);
          }
        }
        // Apply Dirichlet BC (circular mask).
        for (let k = 0; k < N * N; k++) {
          L3d[base + k] = mask[k] ? Lnew_xy[k] : 0;
        }
      }
    }

    // Step 2: Crank-Nicolson in z for each (kx,ky) column.
    // The cell-emission source enters as a non-homogeneous Neumann flux BC at z=0:
    //   D_L ∂_z L|_{z=0} = S_2D(x,y)   where S_2D = src[k] / (dx² × dt_agent).
    // This is equivalent to the 2/h_0 δ(z) source term from the spec (§1b):
    //   The modified Laplacian at j=0 gains term +2 S_2D / h_0.
    // C-N RHS at j=0 gains: +h × 2 × S_2D / h_0  (CN θ=1/2, both sides).
    // We apply the full dt_agent source in this one C-N step (not split into substeps).
    //
    // Surface source per unit area per unit time: S_2D_k = src[k] / (dx² × dt_agent).
    // C-N contribution to rhs[0]: +2 × dt_agent/2 × S_2D_k / h_0
    //                            = src[k] / (dx² × h_0).

    // Use closure-level scratch buffers (no per-step allocation).
    const col = _col, al = _al, bl = _bl, cl = _cl;
    const rhs = _rhs, c2 = _c2, d2 = _d2, sol = _sol;

    const h = dt_agent / 2;
    // Rebuild LHS only if dt_agent changed (constant within a run).
    if (_lhsDt !== dt_agent) {
      for (let j = 0; j < N_z; j++) {
        al[j] = -h * za[j];
        bl[j] =  1 - h * zb[j];
        cl[j] = -h * zc[j];
      }
      _lhsDt = dt_agent;
    }

    for (let ky = 0; ky < N; ky++) {
      for (let kx = 0; kx < N; kx++) {
        if (!mask[ky * N + kx]) continue;
        const k = ky * N + kx;
        // Extract column.
        for (let kz = 0; kz < N_z; kz++) {
          col[kz] = L3d[kz * N * N + k];
        }
        // Build RHS: (I + h D_z²) L_old  +  surface-source correction at j=0.
        // Surface source per unit area per unit time: S_2D = src[k] / (dx² × dt_agent).
        // C-N correction to j=0: dt_agent × 2 S_2D / h_0 = 2 src[k] / (dx² h_0).
        // Derivation: control-volume balance at z=0 node (half-volume h_0/2):
        //   dL[0]/dt = (2/h_0²)(L[1]-L[0]) + 2 S_2D / h_0 (see spec §1b).
        // The C-N discretization adds dt × (2 S_2D / h_0) = 2 src/(dx² h_0) to rhs[0].
        const surfSrc = 2 * src[k] / (dx * dx * h_0);
        for (let j = 0; j < N_z; j++) {
          const jm = j > 0 ? j - 1 : 0;           // Neumann zero-flux ghost
          const jp = j < N_z - 1 ? j + 1 : N_z - 2;
          rhs[j] = col[j] + h * (za[j] * col[jm] + zb[j] * col[j] + zc[j] * col[jp]);
        }
        rhs[0] += surfSrc;
        // Thomas algorithm (in-place using pre-allocated scratch).
        c2[0] = cl[0] / bl[0];
        d2[0] = rhs[0] / bl[0];
        for (let i = 1; i < N_z; i++) {
          const denom = bl[i] - al[i] * c2[i - 1];
          c2[i] = cl[i] / denom;
          d2[i] = (rhs[i] - al[i] * d2[i - 1]) / denom;
        }
        sol[N_z - 1] = d2[N_z - 1];
        for (let i = N_z - 2; i >= 0; i--) sol[i] = d2[i] - c2[i] * sol[i + 1];

        for (let kz = 0; kz < N_z; kz++) {
          L3d[kz * N * N + k] = Math.max(0, sol[kz]);
        }
      }
    }

    // Clear source accumulator (consumed above).
    src.fill(0);

    // Sync z=0 back to 2D field for agent sampling.
    syncZ0toField2d();
  }

  return { step: step2d3d, resetL3d, getL3d: () => L3d };
}
