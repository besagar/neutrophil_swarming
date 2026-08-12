// setup4/solvers/solver_auxfield.js
// M1 relay + ONE auxiliary basally-produced diffusing field.
//
// Two cue models share this solver — they differ only in how the auxiliary
// field enters the per-cell source WEIGHT, which is computed upstream in
// sim_core's emission loop, never here:
//
//   M6.1 (basal adenosine 𝓐, catalog §7)  — shifts the relay threshold:
//     ∂_t̃ 𝓛 = ∇̃²𝓛 + (1/h̃) Σ_i H⁺(𝓛_i; 1 + λ_A 𝓐_i; n_L) δ̃²(r̃−r̃_i) − Γ_L 𝓛
//     ∂_t̃ 𝓐 = D ∇̃²𝓐 +        Σ_i                        δ̃²(r̃−r̃_i) − Γ_A 𝓐
//
//   M6.2 (quorum signal 𝓠, catalog §7b) — throttles the production rate:
//     ∂_t̃ 𝓛 = ∇̃²𝓛 + (1/h̃) Σ_i H⁺(𝓛_i;1;n_L) H⁻(𝓠_i;1;m) δ̃²(r̃−r̃_i) − Γ_L 𝓛
//     ∂_t̃ 𝓠 = D ∇̃²𝓠 +    β Σ_i                          δ̃²(r̃−r̃_i) − γ 𝓠
//
// In both cases the auxiliary field obeys the SAME PDE — linear
// reaction–diffusion with a basal (unconditional) PIC source — so the solver
// reads its coefficients from two model-agnostic params keys:
//   params.D_aux      diffusivity ratio (D_A/D_L or D_Q/D_L)
//   params.gamma_aux  decay rate        (Γ_A or γ)
// The per-cell source amplitude (1 for M6.1, β for M6.2, times the geometry's
// h̃ factor) is applied by sim_core at accumulation time.
//
// The L PDE is structurally identical to M1, so the M1 steppers are reused
// verbatim. Geometry: both 2D–2D and 2D–3D.

import { makeStepFn_2d2d, makeSlabStepper } from './solver_m1.js';
import { createField }                      from './field.js';

/**
 * Auxiliary-field diffusion stepper (2D–2D).
 * ∂_t̃ F = D ∇̃² F + src − γ F, explicit Euler with CFL sub-stepping.
 *
 * The source `src[k]` is accumulated by sim_core (per-cell PIC weight ×
 * amplitude × dt). The discrete δ² on the grid is src[k]/dx², applied once
 * (matches solver_m1's L convention).
 *
 * BC: Dirichlet F = 0 outside the circular dish (same mask as L). With a
 * finite screening length ℓ = √(D/γ) ≪ R̃_dish, the interior steady tone is
 * unaffected by this edge choice.
 *
 * @param {Object} fieldAux - field API for the auxiliary field (from createField)
 * @returns {(dt_agent:number, D:number, gamma:number) => void}
 */
function makeAuxStep2d2d(fieldAux) {
  let Fnew = new Float32Array(fieldAux.N_grid * fieldAux.N_grid);

  return function stepAux(dt_agent, D, gamma) {
    const N   = fieldAux.N_grid;
    const dx  = fieldAux.dx;
    const F    = fieldAux.getLfield();  // fieldAux reuses the L-field storage
    const src  = fieldAux._getSrc();
    if (Fnew.length !== N * N) Fnew = new Float32Array(N * N);

    // CFL for the D-scaled diffusion: dt_sub ≤ dx²/(4D).
    const Dsafe  = Math.max(D, 1e-12);
    const dtCFL  = (dx * dx) / (4 * Dsafe);
    const n_sub  = Math.max(1, Math.ceil(dt_agent / dtCFL));
    const dt_sub = dt_agent / n_sub;

    for (let sub = 0; sub < n_sub; sub++) {
      for (let j = 0; j < N; j++) {
        for (let i = 0; i < N; i++) {
          const k  = j * N + i;
          const Fc = F[k];
          const Fn = (j > 0)     ? F[(j - 1) * N + i] : 0;
          const Fs = (j < N - 1) ? F[(j + 1) * N + i] : 0;
          const Fw = (i > 0)     ? F[j * N + i - 1]   : 0;
          const Fe = (i < N - 1) ? F[j * N + i + 1]   : 0;
          const lap = (Fn + Fs + Fw + Fe - 4 * Fc) / (dx * dx);
          // Source (already carries dt) applied once on sub 0, like solver_m1.
          const directSrc = (sub === 0) ? src[k] / (dx * dx) : 0;
          Fnew[k] = Fc + dt_sub * (D * lap - gamma * Fc) + directSrc;
        }
      }
      if (sub === 0) src.fill(0);
      F.set(Fnew);
      fieldAux.applyBC();
    }
  };
}

/**
 * Create the 2D–2D "M1 + one auxiliary basal field" solver.
 * Owns a second field alongside the shared L field.
 *
 * @param {Object} field2d - the L field API (from createField)
 * @returns {{ fieldAux, step, reset, seedAux }}
 */
export function makeStepFn_auxfield(field2d) {
  const fieldAux = createField(field2d.N_grid, field2d.R_dish);
  const stepL    = makeStepFn_2d2d(field2d);
  const stepAux  = makeAuxStep2d2d(fieldAux);

  return {
    fieldAux,
    /**
     * Advance the auxiliary field then 𝓛 by one agent step. Both source
     * accumulators are filled by sim_core BEFORE this call (PIC ordering).
     * The two fields do not read each other during the step: the auxiliary
     * field's effect on the L source (threshold shift or production throttle)
     * was already resolved upstream from the pre-step auxiliary values.
     */
    step(dt_agent, params) {
      stepAux(dt_agent, params.D_aux, params.gamma_aux);
      stepL(dt_agent, params);
    },
    reset() {
      field2d.reset();
      fieldAux.reset();
    },
    /** Pre-seed the auxiliary field at a uniform steady tone (2D–2D). */
    seedAux(tone) {
      fieldAux.fillUniform(tone);
    },
  };
}

/**
 * Create the 2D–3D "M1 + one auxiliary basal field" solver. Both 𝓛 and the
 * auxiliary field are full 3D thin-slab fields:
 *   ∂_t̃ 𝓛 = ∇̃²_{3D} 𝓛 + δ̃(z̃) Σ w_i        − Γ_L 𝓛
 *   ∂_t̃ F = D ∇̃²_{3D} F + δ̃(z̃) Σ (basal)  − γ F
 * driven by the same generalized slab stepper — 𝓛 with (D=1, decay=Γ_L), F
 * with (D_aux, γ_aux). Cells sit at z=0 and read F|_{z=0}; the z=0 slice of
 * each 3D field is mirrored into its 2D field for agent I/O.
 *
 * @param {Object} field2d  - the L 2D field (z=0 slice + src)
 * @param {Object} params3d - { N_z, h_0, alpha }
 * @returns {{ fieldAux, step, reset, seedAux, getL3d, getAux3d }}
 */
export function makeStepFn_auxfield_2d3d(field2d, params3d) {
  const fieldAux = createField(field2d.N_grid, field2d.R_dish);
  const stepL    = makeSlabStepper(field2d, {
    ...params3d, getD: () => 1, getDecay: (p) => p.gamma_L || 0,
  });
  const stepAux  = makeSlabStepper(fieldAux, {
    ...params3d, getD: (p) => p.D_aux || 1, getDecay: (p) => p.gamma_aux || 0,
  });

  return {
    fieldAux,
    step(dt_agent, params) {
      stepAux.step(dt_agent, params);
      stepL.step(dt_agent, params);
    },
    reset() {
      stepL.reset();
      stepAux.reset();
    },
    /**
     * Pre-seed the auxiliary field at its screened steady profile
     * F(z̃) = tone · exp(−z̃ √(γ/D)), the half-space solution of
     * D F'' = γ F with the basal surface flux fixing F(0) = tone.
     */
    seedAux(tone, D, gamma) {
      const k = Math.sqrt(Math.max(gamma, 0) / Math.max(D, 1e-12));
      stepAux.seed(z => tone * Math.exp(-k * z));
    },
    getL3d:   stepL.getField,
    getAux3d: stepAux.getField,
  };
}
