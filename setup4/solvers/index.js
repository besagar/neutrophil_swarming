// setup4/solvers/index.js
// Factory: createSolver(geometry, model, field2d, extraParams)
//
// Returns a solver object with:
//   { step(dt_agent, params), reset(ic), field }
//
// Architecture note: model branching lives entirely here.
// Adding M3/M4 = import the new solver module + add a case below.
// worker.js only calls solver.step() and solver.field.*; it never
// imports raw grid arrays.
//
// Supported now:
//   M1 × {2d2d, 2d3d}
//   M2 × {2d2d, 2d3d} — same PDE step as M1 (per-cell H⁻(R̃) gating is
//                       applied upstream in the worker's source-accumulation
//                       loop; the L PDE is structurally identical).

import { makeStepFn_2d2d, makeStepFn_2d3d } from './solver_m1.js';

/**
 * @param {'2d2d'|'2d3d'} geometry
 * @param {'M1'|'M2'} model
 * @param {Object} field2d - field API from createField()
 * @param {Object} [extraParams] - e.g. { N_z, h_0, alpha } for 2d3d
 * @returns {{ step(dt_agent, params): void, reset(ic): void, field: Object, extras: Object }}
 */
export function createSolver(geometry, model, field2d, extraParams = {}) {
  if (model !== 'M1' && model !== 'M2') {
    throw new Error(`createSolver: model "${model}" not yet implemented (M3/M4 are future work)`);
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
