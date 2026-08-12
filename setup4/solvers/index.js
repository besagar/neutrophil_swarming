// setup4/solvers/index.js
// Factory: createSolver(geometry, model, field2d, extraParams)
//
// Returns a solver object with:
//   { step(dt_agent, params), reset(ic), field, fieldAux?, seedAux?, extras }
//
// Architecture note: model branching lives entirely here.
// Adding M3/M4 = import the new solver module + add a case below.
// worker.js/sim_core.js only call solver.step() and solver.field.*; they never
// import raw grid arrays.
//
// Supported now:
//   M1   × {2d2d, 2d3d}
//   M2   × {2d2d, 2d3d} — same PDE step as M1 (per-cell H⁻(R̃) gating is
//                         applied upstream in sim_core's source-accumulation
//                         loop; the L PDE is structurally identical).
//   M6.1 × {2d2d, 2d3d} — M1 relay + a basal adenosine field 𝓐 that SHIFTS the
//                         relay threshold (density-independent wave speed).
//   M6.2 × {2d2d, 2d3d} — M1 relay + a basal quorum field 𝓠 that THROTTLES the
//                         production rate by H⁻(𝓠;1;m).
// M6.1 and M6.2 share solver_auxfield.js: the auxiliary PDE is identical
// (linear reaction–diffusion with a basal PIC source; coefficients read from
// params.D_aux / params.gamma_aux). Only the per-cell L source weight differs,
// and that is computed in sim_core. Both expose solver.fieldAux.

import { makeStepFn_2d2d, makeStepFn_2d3d }   from './solver_m1.js';
import { makeStepFn_auxfield,
         makeStepFn_auxfield_2d3d }           from './solver_auxfield.js';

/** Models that carry one auxiliary basally-produced diffusing field. */
const AUX_MODELS = new Set(['M6.1', 'M6.2']);

/**
 * @param {'2d2d'|'2d3d'} geometry
 * @param {'M1'|'M2'|'M6.1'|'M6.2'} model
 * @param {Object} field2d - field API from createField()
 * @param {Object} [extraParams] - e.g. { N_z, h_0, alpha } for 2d3d
 * @returns {{ step(dt_agent, params): void, reset(ic): void, field: Object,
 *            fieldAux?: Object, seedAux?: Function, extras: Object }}
 */
export function createSolver(geometry, model, field2d, extraParams = {}) {
  if (model !== 'M1' && model !== 'M2' && !AUX_MODELS.has(model)) {
    throw new Error(`createSolver: model "${model}" not yet implemented (M3/M4/M5 are future work)`);
  }

  // M6.1 / M6.2: M1 relay + one basal auxiliary field.
  if (AUX_MODELS.has(model)) {
    if (geometry === '2d3d') {
      const { N_z = 16, h_0 = 0.1, alpha = 1.4 } = extraParams;
      const s = makeStepFn_auxfield_2d3d(field2d, { N_z, h_0, alpha });
      return {
        field:    field2d,
        fieldAux: s.fieldAux,
        step(dt_agent, params) { s.step(dt_agent, params); },
        reset() { s.reset(); },
        seedAux(tone, D, gamma) { s.seedAux(tone, D, gamma); },
        extras: { getL3d: s.getL3d, getAux3d: s.getAux3d },
      };
    }
    const s = makeStepFn_auxfield(field2d);
    return {
      field:    field2d,
      fieldAux: s.fieldAux,   // sim_core samples/accumulates it directly
      step(dt_agent, params) { s.step(dt_agent, params); },
      reset() { s.reset(); },
      seedAux(tone) { s.seedAux(tone); },
      extras: {},
    };
  }

  if (geometry === '2d2d') {
    const stepFn = makeStepFn_2d2d(field2d);
    return {
      field: field2d,
      step(dt_agent, params) {
        stepFn(dt_agent, params);
      },
      reset() {
        field2d.reset();
      },
      extras: {},
    };
  }

  if (geometry === '2d3d') {
    const { N_z = 16, h_0 = 0.1, alpha = 1.4 } = extraParams;
    const solver3d = makeStepFn_2d3d(field2d, { N_z, h_0, alpha });
    return {
      field: field2d,   // agents sample from field2d (z=0 slice is kept in sync)
      step(dt_agent, params) {
        solver3d.step(dt_agent, params);
      },
      reset() {
        solver3d.resetL3d();
      },
      extras: { getL3d: solver3d.getL3d },
    };
  }

  throw new Error(`createSolver: unknown geometry "${geometry}"`);
}
