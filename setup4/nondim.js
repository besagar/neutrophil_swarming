// setup4/nondim.js
// Dimensional ↔ nondimensional linkage for Setup 4 (intrinsic-units scheme).
//
// Length and time scales: single-cell intrinsic — σ does NOT appear.
//   ℓ_0 = a / (L_0 D_L)                  [µm]
//   t_0 = a² / (L_0² D_L³)  (= ℓ_0²/D_L) [s]
//
// L PDE in nondim units (discrete-ABM form):
//   2D-3D:  ∂_t̃ 𝓛 = ∇̃²𝓛 + δ̃(z̃) Σ_i H⁺(𝓛_i;1;n_L) H⁻(ℛ_i;1;n_R) δ̃²(r̃−r̃_i) − Γ_L 𝓛
//   2D-2D:  ∂_t̃ 𝓛 = ∇̃²𝓛 + (1/h̃) Σ_i H⁺(𝓛_i;1;n_L) H⁻(ℛ_i;1;n_R) δ̃²(r̃−r̃_i) − Γ_L 𝓛
// per-cell prefactor = 1 in 2D-3D, 1/h̃ in 2D-2D. No 1/σ̃ anywhere.
//
// Cell SDE (no Λ wrapper):
//   Ṗ_i = κ ∇̃𝓛 + λ(𝓛 − 𝓛_c)P + ν(|P|² − |P|⁴)P + √(2ϑ𝓛) ξ
//   ṙ̃_i = μ P_i
//
// Nondim groups (all σ-independent except σ̃ itself):
//   λ   = r_0 · t_0 · L_0       (uses L_0, not L_c — see plan §1)
//   ν   = u² · t_0 / w
//   κ   = a χ / (D_L² · p_0)
//   μ   = μ_dim · t_0 / ℓ_0
//   ϑ   = (w L_0 t_0 / u) · θ
//   𝓛_c = L_c / L_0
//   σ̃   = σ · ℓ_0²              (TRUE cell density now — σ-independent units)
//   h̃   = h / ℓ_0               (only used in 2D-2D)
//   Γ_L = γ_L · t_0
//   Γ_R = γ_R · t_0
//   β   = b · t_0 / R_0         (R_0 = 1 implicit; β̃ = Beta_dim · t_0)
//   𝓛_r = L_r / L_0
//
// Reference: docs/physics/setup4_intrinsic_units_implementation_plan.md §1.

/**
 * Default dimensional parameter values.
 * Units: a [nM·µm³/s], σ [cells/µm²], D_L [µm²/s], L_0 [nM], L_c [nM],
 *        r_0 [1/(s·nM)], u [1/s], w [1/s], χ [µm²/(s·nM)], μ [µm/s],
 *        θ [nM/s], h [µm] (2D-2D only), Γ_L [1/s].
 *
 * Calibration is provisional (D1 of the plan) — these give:
 *   ℓ_0 = 10 µm, t_0 = 1 s
 *   λ = 1, ν = 1, κ = 0.075, μ = 0.005, ϑ = 1e-4
 *   σ̃ = 0.1 (discrete regime; raise σ or a for continuum)
 */
export const DIM_DEFAULTS = {
  a:    1000,   // nM·µm³/s per cell — bumped from 200 to give ℓ_0=10 µm,
                //   t_0=1 s with L_0=1, D_L=100. Lower → discrete regime.
  sigma: 0.01,  // cells/µm² (σ̃ ≈ 1 at default ℓ_0 — continuum/Dieterle regime)
  D_L:  100,    // µm²/s
  L_0:  1.0,    // nM
  L_c:  1.0,    // nM
  r_0:  1.0,    // 1/(s·nM)
  u:    1.0,    // 1/s
  w:    1.0,    // 1/s (p_0 = √(u/w) = 1)
  chi:  0.75,   // µm²/(s·nM) per polarization unit
  mu:   0.05,   // µm/s per polarization unit (μ_dim; m_dim ≡ 1/μ_dim notationally)
  theta: 1e-4,  // nM/s  (noise amplitude)
  h:    10,     // µm (2D-2D layer height)
  Gamma_L: 0,   // 1/s (LTB4 decay; 0 = no decay)
  // M2 per-cell inhibitor (intracellular R_i). R is normalized by R_0 ≡ 1
  // implicitly so β = Beta · t_0 directly.
  Beta:    0.075, // 1/s — per-cell R production rate (β̃ = 0.075 at t_0=1).
  Gamma_R: 1e-4,  // 1/s — per-cell R degradation rate.
  L_r:     0.01,  // nM — second activation threshold (𝓛_r = 0.01 since L_0 = 1).
};

/**
 * Compute nondim parameters from dimensional values (intrinsic-units scheme).
 *
 * Geometry only affects the L-PDE source prefactor (1 in 2D-3D, 1/h̃ in 2D-2D);
 * all cell-side groups (λ, ν, κ, μ, ϑ) are geometry-independent. This is the
 * whole point of the intrinsic-units rewrite — σ and the geometry are out of
 * the units.
 *
 * @param {Object} d - dimensional parameters (same keys as DIM_DEFAULTS)
 * @param {'2d2d'|'2d3d'} geometry  (geometry kept in signature for symmetry;
 *                                   used only to gate h̃ meaningfulness)
 * @returns {Object} nondim params
 */
export function dimToNondim(d, geometry) {
  const { a, sigma, D_L, L_0, L_c, r_0, u, w, chi, mu: mu_dim, theta, h, Gamma_L,
          Beta, Gamma_R, L_r } = d;

  // Guard against pathological dim params.
  if (!isFinite(a) || a <= 0 || !isFinite(L_0) || L_0 <= 0 ||
      !isFinite(D_L) || D_L <= 0) {
    return {
      lam: 1, nu: 1, kap: 0.075, mu: 0.005, tht: 1e-4,
      L_c_nd: 1, sigma_tilde: 0.1, h_tilde: 1,
      gamma_L: 0, beta_R: 0, gamma_R: 0, L_r_nd: 0.01,
      ell0: 1, t0: 1,
    };
  }

  const ell0 = a / (L_0 * D_L);
  const t0   = ell0 * ell0 / D_L;     // = a² / (L_0² D_L³)
  const p0   = Math.sqrt(u / w);

  return {
    // Cell-side nondim groups (σ-independent, geometry-independent).
    lam:    r_0 * t0 * L_0,
    nu:     (u * u) * t0 / w,
    kap:    a * chi / (D_L * D_L * p0),
    mu:     mu_dim * t0 / ell0,
    tht:    (w * L_0 * t0 / u) * theta,
    L_c_nd: L_c / L_0,
    // L-PDE side groups (σ̃ is a true cell density; h̃ only used in 2D-2D).
    sigma_tilde: sigma * ell0 * ell0,
    h_tilde:     h / ell0,
    gamma_L:     Gamma_L * t0,
    // M2 nondim groups (R_0 ≡ 1 implicit so β̃ = Beta · t_0).
    beta_R:  (Beta    !== undefined ? Beta    : 0) * t0,
    gamma_R: (Gamma_R !== undefined ? Gamma_R : 0) * t0,
    L_r_nd:  (L_r     !== undefined ? L_r     : L_0) / L_0,
    // Diagnostics.
    ell0, t0,
  };
}
