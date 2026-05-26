// setup4/nondim.js
// Dimensional ↔ nondimensional linkage for Setup 4.
//
// Geometry-dependent wave speed c* (nondim value = 1 by construction):
//   2D–2D:  c* = √(a σ D_L / (h L_0))      →  ℓ_0 = √(D_L h L_0 / (a σ))
//   2D–3D:  c* = a σ / L_0                  →  ℓ_0 = D_L L_0 / (a σ)
//
// Then: t_0 = D_L / c*² = ℓ_0 / c*   (= ℓ_0 in nondim since c*=1)
//
// Derived nondim groups:
//   Λ   = r_0 L_c D_L / c*²  = r_0 L_c t_0
//   𝓛_c = L_c / L_0
//   p_0 = √(u/w)
//   χ̃   = χ L_0 / (p_0 c*)
//   μ̃   = μ p_0 / c*
//   λ   = u² / (w r_0 L_c)
//   ϑ   = θ w / (u r_0)
//   Γ̃_L = Γ_L t_0   (optional LTB4 decay, default 0)
//   σ̃   = σ · ℓ_0²  (dim cell density in nondim area units; sets the
//                    per-cell emission prefactor 1/σ̃ — independent of the
//                    *simulated* cell count N; N only enters via Σ_i in the
//                    discrete source).
//
// The pattern is: dim object + recomputeFromDim() + recalibrate(), identical
// to the setup3/main.js idiom.

/**
 * Default dimensional parameter values.
 * Units: a [nM·µm³/s], σ [cells/µm²], D_L [µm²/s], L_0 [nM], L_c [nM],
 *        r_0 [1/(s·nM)], u [1/s], w [1/s], χ [µm²/(s·nM)], μ [µm/s],
 *        θ [nM/s], h [µm] (2D-2D only), Γ_L [1/s].
 */
export const DIM_DEFAULTS = {
  a:   200,    // nM·µm³/s per cell — tuned so ℓ₀ ≈ 15–20 µm at R_dim=200µm,
               //   N=3000 defaults → R̃ ≈ 10–14 in both geometries (visible wave
               //   in O(10) nondim time). Increase for sharper fronts / smaller
               //   ℓ₀; numerically requires more grid res or longer t_max.
  sigma: 1e-3, // cells/µm² (toy default; for experimental density use ~5e-3 = 5000 cells/mm²)
  D_L:  100,   // µm²/s
  L_0:  1.0,   // nM
  L_c:  1.0,   // nM
  r_0:  1.0,   // 1/(s·nM)
  u:    1.0,   // 1/s
  w:    1.0,   // 1/s (p_0 = √(u/w) = 1)
  chi:  0.75,  // µm²/(s·nM) per polarization unit
  mu:   0.05,  // µm/s per polarization unit
  theta: 1e-4, // nM/s  (noise amplitude)
  h:    10,    // µm (2D-2D layer height)
  Gamma_L: 0,  // 1/s (LTB4 decay; 0 = no decay)
  // M2 inhibitor (per-cell intracellular R_i). R is normalized by R_c, so the
  // *dim* knobs that survive into nondim are β, γ (R-ODE rates) and L_r
  // (second activation threshold for the R-ODE). Dim defaults below are tuned
  // so that under the 2D-3D geometry (t_0 ≈ 4.39 s at the default σ̃ etc.):
  //   β̃ ≈ 0.075, γ̃ ≈ 0 (effectively, set to a small ε for log-slider support),
  //   L̃_r = 0.01.
  // 2D-2D has t_0 ≈ 2.1 s so its nondim values are about half (β̃ ≈ 0.036).
  Beta:    0.01709, // 1/s — per-cell R production rate (β̃ = 0.075 in 2D-3D default).
  Gamma_R: 1e-4,    // 1/s — per-cell R degradation rate (γ̃ ≈ 4e-4 ≪ 1 — no decay).
  L_r:     0.01,    // nM — second activation threshold (L̃_r = 0.01 since L_0 = 1).
};

/**
 * Compute nondim parameters from dimensional values.
 *
 * No slaving: σ is the user-facing cell-density knob and all dim params
 * are independent. β = 1 is enforced at the ui.js level by deriving the
 * discrete cell count N = round(σ̃·πR̃²), so the simulated density always
 * matches the nondim density σ̃ that comes out of the dim → nondim map.
 *
 * @param {Object} d - dimensional parameters (same keys as DIM_DEFAULTS)
 * @param {'2d2d'|'2d3d'} geometry
 * @returns {Object} nondim params
 */
export function dimToNondim(d, geometry) {
  const { a, sigma, D_L, L_0, L_c, r_0, u, w, chi, mu: mu_dim, theta, h, Gamma_L,
          Beta, Gamma_R, L_r } = d;

  // Wave speed c* (nondim = 1 after rescaling).
  let cstar;
  if (geometry === '2d2d') {
    cstar = Math.sqrt(a * sigma * D_L / (h * L_0));
  } else {
    cstar = a * sigma / L_0;
  }

  // Guard against zero.
  if (!isFinite(cstar) || cstar <= 0) {
    return {
      Lambda: 1, L_c_nd: 1, chi_nd: 0.75, mu_nd: 0.05,
      lam: 1, tht: 1e-4, gamma_L: 0,
      sigma_tilde: 1,
      beta_R: 1, gamma_R: 0.5, L_r_nd: 1,
      cstar: 1, ell0: 1, t0: 1,
    };
  }

  const ell0 = D_L / cstar;
  const t0   = D_L / (cstar * cstar);  // = ell0 / cstar

  const p0 = Math.sqrt(u / w);

  return {
    Lambda:      r_0 * L_c * t0,
    L_c_nd:      L_c / L_0,
    chi_nd:      chi * L_0 / (p0 * cstar),
    mu_nd:       mu_dim * p0 / cstar,
    lam:         (u * u) / (w * r_0 * L_c),
    tht:         theta * w / (u * r_0),
    gamma_L:     Gamma_L * t0,
    // σ̃ = σ · ℓ_0² — fixes the per-cell nondim emission prefactor to 1/σ̃.
    // Holds for both 2D-2D and 2D-3D: in either geometry the discrete-ABM
    // source `(1/σ̃) Σ_i H⁺_i δ̃(x̃ − x̃_i) [δ̃(z̃)]` is the consistent
    // nondim form, with N appearing only in the sum.
    sigma_tilde: sigma * ell0 * ell0,
    // M2 nondim groups (β̃ = β·t_0, γ̃ = γ·t_0, L̃_r = L_r/L_0).
    // R̃ = R/R_c, so R_c ≡ 1 in nondim (no separate R_c knob).
    beta_R:   (Beta    !== undefined ? Beta    : 0) * t0,
    gamma_R:  (Gamma_R !== undefined ? Gamma_R : 0) * t0,
    L_r_nd:   (L_r     !== undefined ? L_r     : L_0) / L_0,
    cstar, ell0, t0,
  };
}
