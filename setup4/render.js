// setup4/render.js
// Canvas2D rendering for Setup 4.
//
// Exports:
//   drawDish(canvasId, frame, params)       — L heatmap + cell dots
//   drawRadialProfile(canvasId, frame, params) — 1D azimuthal average
//   drawMeanRadius(canvasId, frames, params)   — ⟨r̃⟩_free vs t̃
//
// All plots are strictly nondimensional (no unit strings on axes).
// Axis labels are HTML overlays via decoratePlot() — called once at init.

import { autoFit, makeAxis, drawFrame, strokePath, dot, clipPlot } from '../shared/canvas.js';

// ─── magma colormap (25-stop sampled) ───────────────────────────────────────
// Matches the visual language of the Afonin et al. Ca²⁺-dye experimental
// panels (dark-purple background, magenta/orange mid, near-white bright
// spots). Used for the LTB4 (𝓛) heatmap so high-L regions read as "bright
// activity" just like the experimental fluorescent flashes — even though
// the experiment images Ca²⁺ rather than LTB4 directly.
const MAGMA = [
  [0,0,4],[6,5,21],[14,11,40],[26,16,60],[40,17,81],
  [56,17,99],[73,19,108],[89,25,113],[105,31,115],[121,37,114],
  [137,43,112],[153,49,108],[170,55,104],[186,62,99],[201,71,93],
  [216,82,86],[228,95,78],[239,110,72],[247,128,69],[251,148,72],
  [253,170,81],[253,191,97],[252,212,119],[252,232,146],[252,253,191],
];

// Pre-build a 256-entry magma LUT as Uint8Array[256*3] for fast heatmap rendering.
const MAGMA_LUT = new Uint8Array(256 * 3);
for (let i = 0; i < 256; i++) {
  const t = i / 255;
  const n = MAGMA.length - 1;
  const fi = t * n;
  const lo = Math.floor(fi), hi = Math.min(n, lo + 1);
  const f  = fi - lo;
  MAGMA_LUT[i * 3]     = MAGMA[lo][0] + f * (MAGMA[hi][0] - MAGMA[lo][0]);
  MAGMA_LUT[i * 3 + 1] = MAGMA[lo][1] + f * (MAGMA[hi][1] - MAGMA[lo][1]);
  MAGMA_LUT[i * 3 + 2] = MAGMA[lo][2] + f * (MAGMA[hi][2] - MAGMA[lo][2]);
}

// ── greens colormap (for the M6.2 quorum field 𝓠) ──
// Deliberately a different hue family from the magma 𝓛 dish so the two petri
// panels are never confused at a glance: dark blue-green background, mid
// emerald, pale yellow-green highlights. Monotone in lightness like magma, so
// it reads the same way (bright = more 𝓠).
const GREENS = [
  [2,10,8],[3,18,14],[4,27,20],[5,36,26],[6,45,32],
  [7,55,38],[8,65,44],[9,76,50],[10,87,56],[12,98,62],
  [16,110,68],[22,121,73],[30,132,78],[40,143,83],[52,154,88],
  [66,164,93],[82,175,99],[100,185,107],[121,195,118],[143,205,132],
  [166,215,149],[189,225,169],[210,234,192],[229,243,215],[244,251,236],
];

const GREENS_LUT = new Uint8Array(256 * 3);
for (let i = 0; i < 256; i++) {
  const t = i / 255;
  const n = GREENS.length - 1;
  const fi = t * n;
  const lo = Math.floor(fi), hi = Math.min(n, lo + 1);
  const f  = fi - lo;
  GREENS_LUT[i * 3]     = GREENS[lo][0] + f * (GREENS[hi][0] - GREENS[lo][0]);
  GREENS_LUT[i * 3 + 1] = GREENS[lo][1] + f * (GREENS[hi][1] - GREENS[lo][1]);
  GREENS_LUT[i * 3 + 2] = GREENS[lo][2] + f * (GREENS[hi][2] - GREENS[lo][2]);
}

// Outside-dish / no-data background. Same hue family as the magma low end so
// the dish blends seamlessly into the framing rather than sitting on a cold
// neutral. RGB ≈ darkest magma stop, slightly lifted.
const BG_R = 14, BG_G = 8, BG_B = 28;

// ─── dish / heatmap ──────────────────────────────────────────────────────────

/**
 * Draw the petri dish: L heatmap + circular dish boundary + cell dots.
 * @param {string} canvasId
 * @param {Object} frame - { Lfield?, Qfield?, agentX, agentY, emitting, agentPx?, agentPy? }
 * @param {Object} params - { N_grid, R_dish, t, L_max_display?,
 *   fieldKey?   — which frame array to render: 'Lfield' (default) or 'Qfield'
 *   palette?    — 'magma' (default, the 𝓛 dish) or 'greens' (the M6.2 𝓠 dish)
 *   haloScale?  — multiplier on the warm emission halo (0 disables it) }
 * The second (𝓠) dish on the M6.2 page is the same drawing routine with a
 * different field array and LUT, so cell positions/flags are pixel-identical
 * between the two panels and can be compared directly.
 */
export function drawDish(canvasId, frame, params) {
  const canvas = document.getElementById(canvasId);
  if (!canvas) return;
  const ctx = autoFit(canvas);
  const w = canvas.clientWidth, h = canvas.clientHeight;
  const R = params.R_dish;

  const ax = makeAxis({
    xMin: -R * 1.05, xMax: R * 1.05,
    yMin: -R * 1.05, yMax: R * 1.05,
    w, h, aspect: 1,
  });

  drawFrame(ctx, ax, { showGridX: false, showGridY: false, showTickLabelsX: false, showTickLabelsY: false });

  // No rectangular dark fill: areas outside the circular dish stay on the
  // canvas/page background. The dish itself is filled by the heatmap (or by
  // the faint no-data disc below if no field is available yet).

  // ── L heatmap ──
  // frame.Lfield is a Uint8Array already normalized 0..255 per-frame in the
  // worker. Render at native grid resolution (N×N) into an offscreen canvas
  // then drawImage-scale it to the plot area. Using drawImage (not
  // putImageData) is critical: putImageData IGNORES ctx.setTransform, so on
  // a Retina display the heatmap lands at half-scale in the upper-left of
  // the backing buffer. drawImage honors the dpr transform from autoFit.
  const fieldKey = params.fieldKey || 'Lfield';
  const LUT = (params.palette === 'greens') ? GREENS_LUT : MAGMA_LUT;
  if (frame[fieldKey]) {
    const N = params.N_grid;
    const Lf = frame[fieldKey];
    // Build a N×N ImageData (one byte per grid node mapped via viridis LUT).
    // Outside the dish: dark background colour. We use a tight 1.01·R² cutoff
    // here, in *grid* coordinates, to keep this consistent with the world-
    // coord cutoff used by drawImage destination clipping below.
    const imData = new ImageData(N, N);
    const buf = imData.data;
    const dx_grid = (2 * R) / (N - 1);
    // ImageData row 0 = canvas TOP (low canvas-y when drawImage'd). Our world
    // convention has world y = +R at canvas top (yToPx flips y). So map
    // image row j_img → field row j_field = (N-1) - j_img, giving world
    // y = +R for j_img = 0 (canvas top). Without this flip the heatmap
    // renders upside-down — the bug the user reported.
    for (let j_img = 0; j_img < N; j_img++) {
      const j_field = (N - 1) - j_img;
      const yw = -R + j_field * dx_grid;
      for (let i = 0; i < N; i++) {
        const xw = -R + i * dx_grid;
        const k_field = j_field * N + i;
        const idx = (j_img * N + i) * 4;
        if (xw*xw + yw*yw > R * R * 1.001) {
          // Outside the dish: transparent so the canvas/page bg shows through.
          buf[idx] = 0; buf[idx+1] = 0; buf[idx+2] = 0; buf[idx+3] = 0;
        } else {
          const li = Lf[k_field];
          buf[idx]   = LUT[li * 3];
          buf[idx+1] = LUT[li * 3 + 1];
          buf[idx+2] = LUT[li * 3 + 2];
          buf[idx+3] = 255;
        }
      }
    }
    // Stage in an offscreen canvas, then drawImage at world coords.
    const off = document.createElement('canvas');
    off.width = N; off.height = N;
    off.getContext('2d').putImageData(imData, 0, 0);
    // Source pixel CENTER alignment: drawImage treats source pixel (i,j) as
    // a cell occupying [i,i+1) — so its center is at i+0.5. Our grid nodes
    // live at world x = -R + i·dx_grid with dx_grid = 2R/(N-1). To make
    // the center of source pixel (0,0) land at canvas xToPx(-R) (i.e. on
    // grid node 0), the destination rect must extend a half-cell beyond
    // ±R in each direction.
    ctx.imageSmoothingEnabled = true;
    ctx.imageSmoothingQuality = 'medium';
    const halfGrid = R / (N - 1);    // = dx_grid / 2
    const dx_dest = ax.xToPx(-R - halfGrid);
    const dw_dest = ax.xToPx(R  + halfGrid) - dx_dest;
    const dy_dest = ax.yToPx(R  + halfGrid);
    const dh_dest = ax.yToPx(-R - halfGrid) - dy_dest;
    ctx.drawImage(off, 0, 0, N, N, dx_dest, dy_dest, dw_dest, dh_dest);
  } else {
    // No heatmap yet — draw faint dish fill.
    ctx.save();
    const cx = ax.xToPx(0), cy = ax.yToPx(0);
    const pr = ax.xToPx(R) - ax.xToPx(0);
    ctx.beginPath(); ctx.arc(cx, cy, Math.abs(pr), 0, 2*Math.PI);
    ctx.fillStyle = 'rgba(10,10,20,0.15)'; ctx.fill();
    ctx.restore();
  }

  // Clip to inner plot area for cell dots.
  clipPlot(ctx, ax);

  // ── dish boundary circle ──
  {
    ctx.save();
    const cx = ax.xToPx(0), cy = ax.yToPx(0);
    const pr = ax.xToPx(R) - ax.xToPx(0);
    ctx.beginPath(); ctx.arc(cx, cy, Math.abs(pr), 0, 2*Math.PI);
    ctx.strokeStyle = 'rgba(200,180,210,0.55)'; ctx.lineWidth = 1.5; ctx.stroke();
    ctx.restore();
  }

  // ── target circle (sticking boundary) ──
  // The pathogen the swarm converges on: cells that reach r̃ = R̃_target adhere
  // to it and stop. Drawn in cool teal so it never reads as part of the magma
  // 𝓛 field, with a faint fill to mark the region cells cannot enter.
  if (params.stick_target && params.R_target > 0) {
    ctx.save();
    const cx = ax.xToPx(0), cy = ax.yToPx(0);
    const pr = Math.abs(ax.xToPx(params.R_target) - ax.xToPx(0));
    ctx.beginPath(); ctx.arc(cx, cy, pr, 0, 2*Math.PI);
    ctx.fillStyle = 'rgba(70,190,190,0.16)'; ctx.fill();
    ctx.setLineDash([4, 3]);
    ctx.strokeStyle = 'rgba(120,235,225,0.85)'; ctx.lineWidth = 1.5; ctx.stroke();
    ctx.restore();
  }

  // ── cell dots ──
  if (frame.agentX && frame.agentY && frame.emitting) {
    const N = frame.agentX.length;
    // Slightly larger dots + dark stroke so cells stay visible against any
    // viridis value (especially the bright-yellow saturated regions where
    // a translucent fill would otherwise vanish).
    const dotR = Math.max(2.0, ax.plotW / (params.R_dish * 60));
    ctx.lineWidth = Math.max(0.5, dotR * 0.3);
    const tracked = (params.trackedCellIdx != null) ? (params.trackedCellIdx | 0) : -1;
    // Two visual classes, matched to the experimental Ca²⁺-dye look:
    //   emitting → hot near-white "flash", with a warm halo from the additive
    //              pass below;
    //   inactive → faint lavender ghost, almost the dish background colour.
    // Stroke colour is a deep plum so dots sit on the magma heatmap without
    // the harsh black ring we had against viridis.
    // M6.2: emission is throttled, not switched off — a cell over threshold
    // still fires, at rate H⁻(𝓠;1;m). frame.agentR carries 𝓠_i, so each
    // emitting dot is dimmed toward the inactive colour by its throttle factor
    // (see setup4_m6_2_implementation_plan.md §4: the `emitting` flag itself
    // must stay the pure L gate, or the front metric breaks above σ̃★).
    const throttled = (params.model === 'M6.2') && frame.agentR;
    const mQ = params.m_Q || 2;
    ctx.strokeStyle = 'rgba(30,15,40,0.85)';
    for (let i = 0; i < N; i++) {
      if (i === tracked) continue;  // draw tracked cell last (on top)
      const px = ax.xToPx(frame.agentX[i]);
      const py = ax.yToPx(frame.agentY[i]);
      const isEmit = frame.emitting[i] === 1;
      let fill = isEmit ? 'rgb(255,240,200)' : 'rgba(205,180,225,0.75)';
      if (isEmit && throttled) {
        const q = Math.max(0, frame.agentR[i]);
        const g = 1 / (1 + Math.pow(q, mQ));            // H⁻(𝓠;1;m)
        fill = `rgb(${(205 + 50 * g) | 0},${(180 + 60 * g) | 0},${(225 - 25 * g) | 0})`;
      }
      ctx.fillStyle = fill;
      ctx.beginPath();
      ctx.arc(px, py, dotR, 0, 2*Math.PI);
      ctx.fill();
      // Cells engaged with the target keep the same body but get a teal ring,
      // matching the target circle they are sitting on.
      if (frame.stuck && frame.stuck[i] === 1) {
        ctx.save();
        ctx.strokeStyle = 'rgba(120,235,225,0.95)';
        ctx.lineWidth = Math.max(0.8, dotR * 0.45);
        ctx.stroke();
        ctx.restore();
      } else {
        ctx.stroke();
      }
    }
    // Soft warm halo around emitting cells (additive blend), mimicking the
    // bloom of fluorescent emission in the experimental panels.
    const haloScale = (params.haloScale != null) ? params.haloScale : 1;
    ctx.save();
    ctx.globalCompositeOperation = 'lighter';
    for (let i = 0; i < N && haloScale > 0; i++) {
      if (i === tracked || frame.emitting[i] !== 1) continue;
      const px = ax.xToPx(frame.agentX[i]);
      const py = ax.yToPx(frame.agentY[i]);
      // Halo brightness ∝ actual emission rate (M6.2: dimmed by the throttle).
      const g = haloScale * (throttled
        ? 1 / (1 + Math.pow(Math.max(0, frame.agentR[i]), mQ)) : 1);
      const rg = ctx.createRadialGradient(px, py, dotR * 0.4, px, py, dotR * 2.6);
      rg.addColorStop(0,    `rgba(255,210,140,${(0.55 * g).toFixed(3)})`);
      rg.addColorStop(0.55, `rgba(220,90,80,${(0.18 * g).toFixed(3)})`);
      rg.addColorStop(1,    'rgba(120,30,80,0.0)');
      ctx.fillStyle = rg;
      ctx.beginPath();
      ctx.arc(px, py, dotR * 2.6, 0, 2*Math.PI);
      ctx.fill();
    }
    ctx.restore();
    // Tracked cell: same shape as the others, just bigger, always cyan, with
    // a soft white contour. Drawn last so it stays on top of every other dot
    // and halo.
    if (tracked >= 0 && tracked < N) {
      const px = ax.xToPx(frame.agentX[tracked]);
      const py = ax.yToPx(frame.agentY[tracked]);
      const Rt = dotR * 1.8;
      ctx.beginPath();
      ctx.arc(px, py, Rt, 0, 2*Math.PI);
      ctx.fillStyle = 'rgb(150,230,110)';
      ctx.fill();
      ctx.strokeStyle = 'rgba(30,15,40,0.9)';
      ctx.lineWidth = 1.2;
      ctx.stroke();
    }
  }

  ctx.restore();  // matches clipPlot

  // Time label (nondim).
  if (params.t != null) {
    ctx.save();
    ctx.fillStyle = '#ddd';
    ctx.font = '12px ui-monospace, monospace';
    ctx.fillText(`t̃ = ${params.t.toFixed(2)}`, ax.padL + 6, ax.padT + 14);
    ctx.restore();
  }
}

// ─── radial profile ──────────────────────────────────────────────────────────

/**
 * Dieterle planar-front analytical L̃(r̃) for the no-inhibitor (M1), 2D-3D
 * limit, in intrinsic-units nondim (ℓ_0 = a/(L_0 D_L), t_0 = ℓ_0²/D_L).
 *
 * In intrinsic units the continuum surface source is σ̃·H⁺, so the Dieterle
 * planar-front solution maps under ξ_old → σ̃·ξ_new (since ℓ_a/ℓ_old = σ̃).
 * The analytic profile in intrinsic-unit coordinate ξ ≡ r̃ − r̃_front is:
 *     𝓛(ξ) = √(−2 σ̃ ξ)                                  for ξ < 0  (behind)
 *     𝓛(ξ) = (π/2) · exp(−2 σ̃ ξ/π) / √(2 σ̃ ξ)          for ξ > 0  (ahead)
 *
 * Equivalently: z = σ̃·ξ, then the formulas are the canonical β=1 Dieterle
 * expressions in z. Wave speed: c̃_Dieterle = (2/π)·σ̃.
 * Front-layer width ~ 1/σ̃ (wide in discrete regime, narrow in continuum).
 *
 * @param {number} xi - signed distance from front in intrinsic-unit ℓ_0
 * @param {number} sigma_tilde - dimensionless cell density σ̃ = σ·ℓ_0²
 */
function dieterleProfile(xi, sigma_tilde = 1) {
  const z = sigma_tilde * xi;
  if (z < 0) return Math.sqrt(-2 * z);
  return (Math.PI / 2) * Math.exp(-2 * z / Math.PI) / Math.sqrt(2 * z);
}

/**
 * Draw 1D azimuthal average 𝓛(r̃).
 * @param {string} canvasId
 * @param {Object} frame - { radialProfile: Float32Array }
 * @param {Object} params - { N_grid, R_dish, model?, geometry?, L_r_nd? }
 *   When model === 'M1' && geometry === '2d3d', overlay the Dieterle
 *   analytical planar-front profile (dashed) anchored at the simulated
 *   r̃_front (where 𝓛 crosses 1). In M2 a second dashed line is drawn at
 *   L̃_r (the second activation threshold for the R-ODE), passed via
 *   params.L_r_nd.
 */
export function drawRadialProfile(canvasId, frame, params) {
  const canvas = document.getElementById(canvasId);
  if (!canvas) return;
  const ctx = autoFit(canvas);
  const w = canvas.clientWidth, h = canvas.clientHeight;

  if (!frame || !frame.radialProfile) {
    drawFrame(ctx, makeAxis({ xMin: 0, xMax: params.R_dish, yMin: 0, yMax: 5, w, h }));
    return;
  }

  const prof = frame.radialProfile;
  const N    = params.N_grid;
  const R    = params.R_dish;
  const dx   = (2 * R) / (N - 1);

  // Build xs, ys for r = 0..R_dish.
  const xs = [], ys = [];
  for (let k = 0; k < N; k++) {
    const r = k * dx;
    if (r > R * 1.01) break;
    xs.push(r);
    ys.push(prof[k]);
  }

  // Y range: 0 to max(prof)+pad.
  let ymax = 0.5;
  for (let k = 0; k < ys.length; k++) if (ys[k] > ymax) ymax = ys[k];
  ymax = ymax * 1.15;

  const ax = makeAxis({ xMin: 0, xMax: R, yMin: 0, yMax: ymax, w, h });
  drawFrame(ctx, ax);

  // Draw L=1 threshold line (relay threshold in nondim).
  ctx.save();
  ctx.strokeStyle = 'rgba(200,60,30,0.5)';
  ctx.lineWidth = 1;
  ctx.setLineDash([4, 4]);
  ctx.beginPath();
  ctx.moveTo(ax.xToPx(0), ax.yToPx(1));
  ctx.lineTo(ax.xToPx(R), ax.yToPx(1));
  ctx.stroke();
  ctx.restore();

  // Draw L̃_r threshold line (M2 only: second activation threshold for the
  // R-ODE). Skipped when L̃_r coincides with the relay threshold (= 1) or is
  // off-axis, to avoid drawing on top of the L=1 line.
  if (params.model === 'M2' && params.L_r_nd && Math.abs(params.L_r_nd - 1) > 1e-3
      && params.L_r_nd > 0 && params.L_r_nd <= ymax) {
    ctx.save();
    ctx.strokeStyle = 'rgba(120,80,170,0.6)';  // purple, distinct from the L=1 red
    ctx.lineWidth = 1;
    ctx.setLineDash([2, 3]);
    ctx.beginPath();
    ctx.moveTo(ax.xToPx(0), ax.yToPx(params.L_r_nd));
    ctx.lineTo(ax.xToPx(R), ax.yToPx(params.L_r_nd));
    ctx.stroke();
    ctx.fillStyle = 'rgba(120,80,170,0.95)';
    ctx.font = '10px ui-monospace, monospace';
    ctx.fillText(`L̃_r = ${params.L_r_nd.toFixed(2)}`,
                 ax.xToPx(R) - 70, ax.yToPx(params.L_r_nd) - 3);
    ctx.restore();
  }

  strokePath(ctx, ax, xs, ys, { color: '#2b6cb0', width: 1.8 });

  // ── Dieterle analytical overlay (M1 + 2D-3D only) ──
  // Planar-front solution in intrinsic units. Wave speed: c̃ = (2/π)·σ̃.
  // Front-layer width ~ 1/σ̃ — wide in the discrete regime, narrow in continuum.
  // β is no longer meaningful (β=1 by construction in intrinsic units).
  if (params.model === 'M1' && params.geometry === '2d3d') {
    const sigma_til = params.sigma_tilde || 0;
    // Find front position: largest r where 𝓛 crosses 1 going downward.
    let r_front = -1;
    for (let k = ys.length - 2; k >= 0; k--) {
      if (ys[k] >= 1 && ys[k + 1] < 1) {
        const r1 = xs[k], r2 = xs[k + 1];
        const y1 = ys[k], y2 = ys[k + 1];
        r_front = r1 + (1 - y1) / (y2 - y1) * (r2 - r1);
        break;
      }
    }
    if (r_front > 0) {
      // Buffer ±0.5/σ̃ around the front (both branches diverge there).
      // Scales with 1/σ̃: wide buffer at low density, narrow at high.
      const buf = 0.5 / Math.max(sigma_til, 1e-6);
      const NA = 400;
      const style = { color: 'rgba(20,20,20,0.85)', width: 1.6, dash: [5, 4] };
      // Behind the front: ξ ∈ [−r_front, −buf], i.e. r ∈ [0, r_front − buf].
      if (r_front > buf) {
        const xs_b = [], ys_b = [];
        for (let k = 0; k < NA; k++) {
          const xi = -r_front + (-buf - (-r_front)) * (k / (NA - 1));
          const r = r_front + xi;
          if (r < 0) continue;
          const Ly = dieterleProfile(xi, sigma_til);
          if (Ly <= ymax * 1.5 && isFinite(Ly)) { xs_b.push(r); ys_b.push(Ly); }
        }
        if (xs_b.length > 1) strokePath(ctx, ax, xs_b, ys_b, style);
      }
      // Ahead of the front: ξ ∈ [buf, R − r_front].
      const xiMax = R - r_front;
      if (xiMax > buf) {
        const xs_a = [], ys_a = [];
        for (let k = 0; k < NA; k++) {
          const xi = buf + (xiMax - buf) * (k / (NA - 1));
          const r = r_front + xi;
          const Ly = dieterleProfile(xi, sigma_til);
          if (Ly <= ymax * 1.5 && isFinite(Ly) && Ly > 0) { xs_a.push(r); ys_a.push(Ly); }
        }
        if (xs_a.length > 1) strokePath(ctx, ax, xs_a, ys_a, style);
      }

      // Annotate analytical wave speed c̃ = (2/π)·σ̃. Warn when the front
      // has reached the dish boundary — past that point the traveling-wave
      // description breaks down.
      const cTheory = (2 / Math.PI) * sigma_til;
      const saturated = (r_front >= R * 0.95);
      ctx.save();
      ctx.fillStyle = saturated ? '#a33' : '#222';
      ctx.font = '11px ui-monospace, monospace';
      ctx.fillText(`Dieterle (c̃ = (2/π)·σ̃ ≈ ${cTheory.toFixed(3)})`,
                   ax.padL + 6, ax.padT + 14);
      if (saturated) {
        ctx.fillText('⚠ wave hit boundary — traveling-wave profile invalid',
                     ax.padL + 6, ax.padT + 28);
      }
      // Regime warnings — overlay is the continuum-limit Dieterle prediction
      // and degrades at BOTH ends of σ̃:
      //   σ̃ ≪ 1 — discrete-cell regime (jagged hotspots, no smooth front).
      //   σ̃·dx̃ ≳ 0.5 — front layer 1/σ̃ is sub-grid; the √-cusp at the
      //     front renders as a linear ramp (NOT a physics bug; raise N_grid
      //     or lower σ̃ to resolve). See also Holmes et al. arXiv:2101.01181
      //     on density-paradox effects in 2D-3D LTB4 relays.
      const dx = (2 * R) / (N - 1);
      if (sigma_til < 0.3) {
        ctx.fillStyle = '#a33';
        ctx.fillText(`σ̃ = ${sigma_til.toFixed(2)} — discrete regime, overlay approximate`,
                     ax.padL + 6, ax.padT + 42);
      } else if (sigma_til * dx > 0.5) {
        ctx.fillStyle = '#a33';
        ctx.fillText(`σ̃·dx̃ = ${(sigma_til*dx).toFixed(2)} — front layer sub-grid; raise N_grid`,
                     ax.padL + 6, ax.padT + 42);
      }
      ctx.restore();
    }
  }
}

/**
 * Draw the per-cell inhibitor radial profile R̃(r̃) for the current frame.
 *
 * R̃ is a per-cell scalar (no grid), so the profile is built by binning cells
 * by their radial position and averaging R̃_i within each bin. Empty bins
 * (no cells) leave a gap in the line. A horizontal dashed line marks the
 * shut-off threshold R̃ = 1 (in nondim units, R_c ≡ 1 by construction).
 *
 * @param {string} canvasId
 * @param {Object} frame  - { agentX, agentY, agentR }
 * @param {Object} params - { R_dish, N_bins? } (default N_bins = 64)
 */
export function drawRadialR(canvasId, frame, params) {
  const canvas = document.getElementById(canvasId);
  if (!canvas) return;
  const ctx = autoFit(canvas);
  const w = canvas.clientWidth, h = canvas.clientHeight;
  const R = params.R_dish;
  const nBins = params.N_bins || 64;

  if (!frame || !frame.agentR || !frame.agentX) {
    drawFrame(ctx, makeAxis({ xMin: 0, xMax: R, yMin: 0, yMax: 2, w, h }));
    return;
  }

  const sums   = new Float64Array(nBins);
  const counts = new Uint32Array(nBins);
  const dr     = R / nBins;
  for (let i = 0; i < frame.agentR.length; i++) {
    const r = Math.hypot(frame.agentX[i], frame.agentY[i]);
    const b = Math.min(nBins - 1, Math.floor(r / dr));
    sums[b]   += frame.agentR[i];
    counts[b] += 1;
  }

  // M6.1 / M6.2: the panel shows the auxiliary basal field (𝓐 / 𝓠) sampled at
  // the cells; the guide line is its mean-field tone rather than the R̃ = 1
  // shutoff. 2D–3D uses the screened surface tone ∝ σ̃/√(Dγ), and M6.2's source
  // amplitude carries an extra h̃ there (catalog §7b).
  const isM61 = params.model === 'M6.1';
  const isM62 = params.model === 'M6.2';
  const is3dR = params.geometry === '2d3d';
  const sig   = params.sigma_tilde || 0;
  let guideY = 1;
  if (isM61) {
    const gA = Math.max(params.gamma_A || 0, 1e-12);
    guideY = is3dR ? sig / Math.sqrt(Math.max((params.D_A_nd || 1) * gA, 1e-12))
                   : sig / gA;
  } else if (isM62) {
    const gQ  = Math.max(params.gamma_Q || 0, 1e-12);
    const amp = (params.beta_Q || 0) * (is3dR ? (params.h_tilde || 1) : 1);
    guideY = is3dR ? amp * sig / Math.sqrt(Math.max((params.D_Q_nd || 1) * gQ, 1e-12))
                   : amp * sig / gQ;
  }

  // Y range: at least ~1.2× the guide line (so it sits within the frame even
  // when the profile stays small); expand if any bin exceeds it.
  let ymax = Math.max(1.2, guideY * 1.2, 1e-6);
  for (let b = 0; b < nBins; b++) {
    if (counts[b] > 0) {
      const r_avg = sums[b] / counts[b];
      if (r_avg * 1.15 > ymax) ymax = r_avg * 1.15;
    }
  }

  const ax = makeAxis({ xMin: 0, xMax: R, yMin: 0, yMax: ymax, w, h });
  drawFrame(ctx, ax);

  // Guide line: R̃=1 shutoff (M2), or the basal tone (M6.1 𝓐_ss / M6.2 𝓠_ss).
  if (guideY > 0) {
    ctx.save();
    ctx.strokeStyle = isM62 ? 'rgba(40,150,150,0.6)'
                    : isM61 ? 'rgba(60,120,200,0.55)' : 'rgba(200,60,30,0.5)';
    ctx.lineWidth = 1;
    ctx.setLineDash([4, 4]);
    ctx.beginPath();
    ctx.moveTo(ax.xToPx(0), ax.yToPx(guideY));
    ctx.lineTo(ax.xToPx(R), ax.yToPx(guideY));
    ctx.stroke();
    ctx.restore();
  }

  // Polyline through occupied bins; break the line at empty bins so a gap
  // doesn't get spanned by a misleading interpolated segment.
  ctx.save();
  ctx.strokeStyle = isM62 ? '#2a9d9d' : isM61 ? '#3a78c8' : '#a06030';
  ctx.lineWidth = 1.8;
  ctx.beginPath();
  let penDown = false;
  for (let b = 0; b < nBins; b++) {
    if (counts[b] === 0) { penDown = false; continue; }
    const r_avg = sums[b] / counts[b];
    const px = ax.xToPx((b + 0.5) * dr);
    const py = ax.yToPx(r_avg);
    if (penDown) ctx.lineTo(px, py);
    else { ctx.moveTo(px, py); penDown = true; }
  }
  ctx.stroke();
  ctx.restore();
}

/**
 * Angular (channelisation) spectrum Ψ_m vs m at one frame.
 *
 * Ψ_m is the shot-noise-corrected azimuthal Fourier amplitude of the cell
 * angular distribution (docs/physics/setup4_swarm3d.md §10):
 *   c_m = (1/N)Σ_i e^{imθ_i},   Ψ_m = √(max(0, (N|c_m|² − 1)/(N − 1))).
 * 0 = angularly uniform, 1 = all cells at one angle with m-fold periodicity.
 * The dashed line is the 5% shot-noise level √(2/(N−1)) — bars below it are
 * not significant. m = 1 is drawn in grey because it measures a bulk
 * off-centre drift of the whole swarm, not spokes; the dominant mode m* is
 * taken over m ≥ 2 and IS the channel count.
 *
 * @param {string} canvasId
 * @param {Object} frame - { chanPsi: Float32Array, chanMstar, chanNoise }
 */
export function drawAngularSpectrum(canvasId, frame) {
  const canvas = document.getElementById(canvasId);
  if (!canvas) return;
  const ctx = autoFit(canvas);
  const w = canvas.clientWidth, h = canvas.clientHeight;

  const psi = frame && frame.chanPsi;
  if (!psi || psi.length < 2) {
    drawFrame(ctx, makeAxis({ xMin: 0, xMax: 16, yMin: 0, yMax: 1, w, h }));
    return;
  }
  const M = psi.length - 1;
  const noise = frame.chanNoise || 0;

  let ymax = Math.max(noise * 1.6, 1e-3);
  for (let m = 1; m <= M; m++) if (psi[m] > ymax) ymax = psi[m];
  ymax *= 1.2;

  const ax = makeAxis({ xMin: 0.5, xMax: M + 0.5, yMin: 0, yMax: ymax, w, h });
  drawFrame(ctx, ax);

  clipPlot(ctx, ax);
  // Bars.
  const halfW = Math.abs(ax.xToPx(1) - ax.xToPx(0)) * 0.33;
  const y0 = ax.yToPx(0);
  for (let m = 1; m <= M; m++) {
    const xc = ax.xToPx(m);
    const yv = ax.yToPx(Math.min(psi[m], ymax));
    // m = 1 (bulk drift) muted; the dominant m ≥ 2 mode highlighted.
    ctx.fillStyle = (m === 1) ? 'rgba(150,150,160,0.55)'
                  : (m === frame.chanMstar) ? '#d9730d'
                  : 'rgba(43,108,176,0.75)';
    ctx.fillRect(xc - halfW, yv, 2 * halfW, y0 - yv);
  }
  // 5% shot-noise level.
  if (noise > 0) {
    strokePath(ctx, ax, [0.5, M + 0.5], [noise, noise],
               { color: 'rgba(120,120,120,0.7)', width: 1, dash: [4, 4] });
  }
  ctx.restore();

  // m* label — the channel count.
  if (frame.chanMstar >= 2 && psi[frame.chanMstar] > noise) {
    ctx.save();
    ctx.fillStyle = '#d9730d';
    ctx.font = '11px ui-monospace, monospace';
    ctx.textAlign = 'right'; ctx.textBaseline = 'top';
    ctx.fillText(`m* = ${frame.chanMstar}`, ax.padL + ax.plotW - 6, ax.padT + 4);
    ctx.restore();
  }
}

/**
 * Generic per-frame scalar vs t̃ time-series plot.
 *
 * Reads a precomputed scalar from each frame (worker frame handler attaches
 * `meanAbsVr`, `meanAbsP`, …) and plots it as a polyline. The vertical scale
 * auto-fits to the data; the horizontal scale runs to the last accumulated
 * frame time (NOT to t_max), so the trace grows as the worker streams.
 *
 * @param {string} canvasId
 * @param {Array} frames
 * @param {string} key   - frame[key] is the scalar to plot
 * @param {Object} [opts] - { color?, yMinAuto?: number, label?, currentT? }
 *   currentT: if provided and within the trace range, draw a vertical marker
 *             (current-frame indicator linked to the time-scrub slider).
 */
export function drawTimeSeries(canvasId, frames, key, opts = {}) {
  const canvas = document.getElementById(canvasId);
  if (!canvas) return;
  const ctx = autoFit(canvas);
  const w = canvas.clientWidth, h = canvas.clientHeight;

  if (!frames || frames.length === 0) {
    drawFrame(ctx, makeAxis({ xMin: 0, xMax: 1, yMin: 0, yMax: 1, w, h }));
    return;
  }

  const ts = [], ys = [];
  let ymax = opts.yMinAuto || 1e-6;
  for (const f of frames) {
    if (f[key] === undefined) continue;
    ts.push(f.t);
    ys.push(f[key]);
    if (f[key] > ymax) ymax = f[key];
  }
  if (ts.length === 0) {
    drawFrame(ctx, makeAxis({ xMin: 0, xMax: 1, yMin: 0, yMax: 1, w, h }));
    return;
  }

  const t_max = ts[ts.length - 1] || 1;
  const ax = makeAxis({ xMin: 0, xMax: t_max, yMin: 0, yMax: ymax * 1.15, w, h });
  drawFrame(ctx, ax);
  strokePath(ctx, ax, ts, ys, { color: opts.color || '#2b6cb0', width: 1.5 });

  // Current-time vertical indicator (time-scrub slider position).
  if (opts.currentT !== undefined && opts.currentT >= 0 && opts.currentT <= t_max) {
    ctx.save();
    ctx.strokeStyle = 'rgba(80,80,80,0.65)';
    ctx.lineWidth = 1;
    ctx.setLineDash([3, 3]);
    const px = ax.xToPx(opts.currentT);
    ctx.beginPath();
    ctx.moveTo(px, ax.padT);
    ctx.lineTo(px, ax.padT + ax.plotH);
    ctx.stroke();
    ctx.restore();
  }
}

// ─── bead-in-free-energy plots ────────────────────────────────────────────────
//
// Free energy (intrinsic-units nondim):
//   F(P) = −[ λ·(𝓛 − 𝓛_c)·|P|²/2  +  ν·(|P|⁴/4 − |P|⁶/6) ]  −  κ·∇𝓛·P
// Ring of minima at |P|² = (λ(𝓛-𝓛_c) + √(λ²(𝓛-𝓛_c)² + 4ν²·|P|⁰…)) / (2ν)
//   simplified for symmetric case: |P|_eq = √(λ(𝓛-𝓛_c)/ν)  when 𝓛 > 𝓛_c.
// Two views below: drawBead3D (wireframe surface + bead), drawBead2D
// (contour-ring phase plane + trajectory + bead).

// Effective free energy seen by one cell.
// Includes the GL Mexican-hat AND the chemotactic linear tilt −κ·(∇𝓛·P)
// which breaks the rotational symmetry of the well.
// Params: lam (λ), nu (ν), kap (κ) — intrinsic-units cell-side groups.
function freeEnergy(Px, Py, L, Lc, lam, nu, kap, gx, gy) {
  const p2 = Px*Px + Py*Py;
  const p4 = p2*p2;
  const p6 = p4*p2;
  const gl   = -(lam * (L - Lc) * p2 / 2 + nu * (p4/4 - p6/6));
  const tilt = -(kap || 0) * ((gx || 0) * Px + (gy || 0) * Py);
  return gl + tilt;
}

function sampleLAtCell(fNow, cellIdx, R_dish, N_grid) {
  if (!fNow.radialProfile) return 0;
  const cellR = Math.hypot(fNow.agentX[cellIdx], fNow.agentY[cellIdx]);
  const dxr = (2 * R_dish) / (N_grid - 1);
  const fbin = cellR / dxr;
  const b0 = Math.max(0, Math.min(fNow.radialProfile.length - 1, Math.floor(fbin)));
  const b1 = Math.min(fNow.radialProfile.length - 1, b0 + 1);
  const tb = fbin - b0;
  return (1 - tb) * fNow.radialProfile[b0] + tb * fNow.radialProfile[b1];
}

function autoPmaxForCell(frames, idx, cellIdx) {
  let Pmax = 0.5;
  for (let i = 0; i <= idx; i++) {
    if (!frames[i].Px) continue;
    const p = Math.hypot(frames[i].Px[cellIdx], frames[i].Py[cellIdx]);
    if (p > Pmax) Pmax = p;
  }
  return Pmax * 1.25;
}

// |P|_eq from the symmetric GL well (no tilt); 0 if sub-threshold.
// In new scheme: |P|_eq = √(λ(𝓛-𝓛_c)/ν). Guard for negative argument.
function gleEqRadius(L, Lc, lam, nu) {
  if (L <= Lc) return 0;
  const p2 = lam * (L - Lc) / Math.max(nu, 1e-12);
  if (p2 <= 0) return 0;
  return Math.sqrt(p2);
}

// Magnitude of the symmetric GL well, measured between F(P=0)=0 and F(|P|_eq).
// Used as a "natural scale" for the display rescaling of the chemotactic tilt.
function gleWellDepth(L, Lc, lam, nu, fallbackPmax) {
  if (L > Lc) {
    const p2 = lam * (L - Lc) / Math.max(nu, 1e-12);
    if (p2 > 0) {
      const p4 = p2 * p2;
      const p6 = p4 * p2;
      const F  = -(lam * (L - Lc) * p2 / 2 + nu * (p4/4 - p6/6));
      return Math.max(Math.abs(F), 1e-6);
    }
  }
  // Sub-threshold: F is a parabolic well at P = 0. Use its rise to fallback Pmax.
  const half = Math.abs(lam * (Lc - L) * fallbackPmax * fallbackPmax / 2);
  return Math.max(half, 1e-6);
}

// Compute a tilt-display rescaling factor: shrink the chemotactic-tilt
// contribution so its magnitude across the visible Pmax stays comparable to
// the GL well depth. Returns 1.0 when the tilt is naturally moderate (no
// rescaling needed). The 1D cross-section plot deliberately does NOT use this —
// it shows the true tilted potential along ∇𝓛.
// kap = κ (intrinsic-units chemotactic coupling).
function tiltDisplayScale(kap, gMag, Pmax, wellDepth, frac = 0.6) {
  const tiltMax = (kap || 0) * (gMag || 0) * Pmax;
  if (tiltMax <= 0) return 1;
  const target = frac * wellDepth;
  return tiltMax > target ? target / tiltMax : 1;
}

// Marching-squares contour extractor. Walks the (Nm-1)² cells of the F mesh
// over (P_x, P_y) ∈ [-Pmax, Pmax]², emitting line segments at F = level.
// Caller is responsible for ctx.beginPath() before / ctx.stroke() after.
function emitContourMS(ctx, ax, F, Nm, Pmax, level) {
  const dx = 2 * Pmax / (Nm - 1);
  function lerp(Fa, Fb, xa, ya, xb, yb) {
    const denom = Fb - Fa;
    if (Math.abs(denom) < 1e-30) return [xa, ya];
    const t = (level - Fa) / denom;
    return [xa + t*(xb - xa), ya + t*(yb - ya)];
  }
  function seg(p, q) {
    ctx.moveTo(ax.xToPx(p[0]), ax.yToPx(p[1]));
    ctx.lineTo(ax.xToPx(q[0]), ax.yToPx(q[1]));
  }
  for (let j = 0; j < Nm - 1; j++) {
    for (let i = 0; i < Nm - 1; i++) {
      const Fa = F[j*Nm + i];           // bottom-left
      const Fb = F[j*Nm + i + 1];       // bottom-right
      const Fc = F[(j+1)*Nm + i + 1];   // top-right
      const Fd = F[(j+1)*Nm + i];       // top-left
      const code = (Fa > level ? 1 : 0) | (Fb > level ? 2 : 0)
                 | (Fc > level ? 4 : 0) | (Fd > level ? 8 : 0);
      if (code === 0 || code === 15) continue;
      const xL = -Pmax + i*dx,     yB = -Pmax + j*dx;
      const xR = -Pmax + (i+1)*dx, yT = -Pmax + (j+1)*dx;
      const pAB = lerp(Fa, Fb, xL, yB, xR, yB);   // bottom
      const pBC = lerp(Fb, Fc, xR, yB, xR, yT);   // right
      const pDC = lerp(Fd, Fc, xL, yT, xR, yT);   // top
      const pAD = lerp(Fa, Fd, xL, yB, xL, yT);   // left
      switch (code) {
        case 1: case 14: seg(pAB, pAD); break;
        case 2: case 13: seg(pAB, pBC); break;
        case 3: case 12: seg(pBC, pAD); break;
        case 4: case 11: seg(pBC, pDC); break;
        case 6: case 9:  seg(pAB, pDC); break;
        case 7: case 8:  seg(pDC, pAD); break;
        case 5:  seg(pAB, pAD); seg(pBC, pDC); break;   // saddle (above-corners connected)
        case 10: seg(pAB, pBC); seg(pDC, pAD); break;   // saddle
      }
    }
  }
}

// Draws a 2D arrow from (x0, y0) to (x1, y1) in world coords, in current ctx
// stroke/fill style. Used for the ∇𝓛 indicator on bead plots.
function arrow2D(ctx, ax, x0, y0, x1, y1, color = 'rgb(255,200,40)') {
  const a = [ax.xToPx(x0), ax.yToPx(y0)];
  const b = [ax.xToPx(x1), ax.yToPx(y1)];
  ctx.save();
  ctx.strokeStyle = color; ctx.fillStyle = color; ctx.lineWidth = 2;
  ctx.beginPath(); ctx.moveTo(a[0], a[1]); ctx.lineTo(b[0], b[1]); ctx.stroke();
  const ang = Math.atan2(b[1] - a[1], b[0] - a[0]);
  const ah = 9;
  ctx.beginPath();
  ctx.moveTo(b[0], b[1]);
  ctx.lineTo(b[0] - ah*Math.cos(ang - Math.PI/7), b[1] - ah*Math.sin(ang - Math.PI/7));
  ctx.lineTo(b[0] - ah*Math.cos(ang + Math.PI/7), b[1] - ah*Math.sin(ang + Math.PI/7));
  ctx.closePath(); ctx.fill();
  ctx.restore();
}

/**
 * 3D wireframe surface F(P_x, P_y) with a bead at the tracked cell's current P.
 * Axonometric projection (no perspective, ~30° tilt), auto-scaled F axis.
 */
export function drawBead3D(canvasId, frames, idx, cellIdx, params) {
  const canvas = document.getElementById(canvasId);
  if (!canvas) return;
  const ctx = autoFit(canvas);
  const w = canvas.clientWidth, h = canvas.clientHeight;

  // Background.
  ctx.fillStyle = '#1a1a25';
  ctx.fillRect(0, 0, w, h);

  if (!frames || frames.length === 0 || !frames[0].Px) return;
  const N = frames[0].Px.length;
  cellIdx = Math.max(0, Math.min(N - 1, cellIdx | 0));

  const fNow = frames[idx];
  const L_cell = sampleLAtCell(fNow, cellIdx, params.R_dish, params.N_grid);
  const { lam, nu, kap, L_c } = params;
  const gx_cell = fNow.Gx ? fNow.Gx[cellIdx] : 0;
  const gy_cell = fNow.Gy ? fNow.Gy[cellIdx] : 0;
  const gMag = Math.hypot(gx_cell, gy_cell);

  // Pmax: always wide enough to contain the GL well minimum, even when the
  // bead is still near the origin. Without this, the surface auto-zooms to a
  // sub-|P|_eq region where the GL well isn't yet developed.
  const pEq = gleEqRadius(L_cell, L_c, lam, nu);
  const Pmax = Math.max(autoPmaxForCell(frames, idx, cellIdx), 1.5 * pEq, 0.5);

  // Chemotactic-tilt rescaling for display.
  // When κ·|∇𝓛|·Pmax ≫ well-depth (typical in saturated regimes), the
  // tilt term swamps the GL surface and the Mexican-hat geometry is invisible.
  // We shrink the displayed tilt so it remains comparable to the well depth;
  // the actual rescaling factor is reported in the label. The bead's height
  // and the 1D cross-section both use the un-rescaled tilt.
  const wellDepth  = gleWellDepth(L_cell, L_c, lam, nu, Pmax);
  const tiltScale  = tiltDisplayScale(kap, gMag, Pmax, wellDepth, 0.6);
  function F_display(px, py) {
    const gl   = freeEnergy(px, py, L_cell, L_c, lam, nu, 0, 0, 0); // GL only
    const tilt = -(kap || 0) * (gx_cell * px + gy_cell * py);
    return gl + tiltScale * tilt;
  }

  // Mesh of F on (Px, Py) using the display-rescaled tilt.
  const Nv = 28;
  const F = new Float32Array(Nv * Nv);
  let Fmin = Infinity, Fmax = -Infinity;
  for (let j = 0; j < Nv; j++) {
    const py = -Pmax + 2 * Pmax * j / (Nv - 1);
    for (let i = 0; i < Nv; i++) {
      const px = -Pmax + 2 * Pmax * i / (Nv - 1);
      const v = F_display(px, py);
      F[j*Nv + i] = v;
      if (v < Fmin) Fmin = v;
      if (v > Fmax) Fmax = v;
    }
  }
  const Frange = (Fmax - Fmin) || 1;

  // Axonometric projection. View angle 30°; F axis goes up the screen.
  const cosT = Math.cos(Math.PI / 6);   // ≈ 0.866
  const sinT = Math.sin(Math.PI / 6);   // = 0.5
  // Choose scales so the box fits inside the canvas with margins.
  const margin = 38;
  const screenW = w - 2 * margin;
  const screenH = h - 2 * margin;
  // Horizontal extent in screen: 2*Pmax * cosT * scaleXY * 2 (since x and y both contribute).
  const scaleXY = Math.min(screenW / (2 * Pmax * 2 * cosT),
                           screenH * 0.5 / (2 * Pmax * 2 * sinT * 0.5));
  const scaleZ = (screenH * 0.55) / Frange;
  const cx = w / 2;
  const cy = h * 0.78;  // bias down so the surface has room to rise.

  function project(px, py, f) {
    const sx = cx + (px - py) * cosT * scaleXY;
    const sy = cy + (px + py) * sinT * scaleXY * 0.5 - (f - Fmin) * scaleZ;
    return [sx, sy];
  }

  // ── Floor box edges + axes ──
  ctx.strokeStyle = 'rgba(180,180,200,0.4)';
  ctx.lineWidth = 1;
  // Four floor edges at z=Fmin (bottom of F range).
  function edge(p1, p2) {
    const a = project(...p1), b = project(...p2);
    ctx.beginPath(); ctx.moveTo(a[0], a[1]); ctx.lineTo(b[0], b[1]); ctx.stroke();
  }
  edge([-Pmax, -Pmax, Fmin], [ Pmax, -Pmax, Fmin]);
  edge([ Pmax, -Pmax, Fmin], [ Pmax,  Pmax, Fmin]);
  edge([ Pmax,  Pmax, Fmin], [-Pmax,  Pmax, Fmin]);
  edge([-Pmax,  Pmax, Fmin], [-Pmax, -Pmax, Fmin]);

  // ── Wireframe mesh on surface ──
  // Back-to-front (painter's): higher j+i = closer to viewer.
  ctx.lineWidth = 0.8;
  // Lines of constant Py (i varies).
  for (let j = 0; j < Nv; j++) {
    const py = -Pmax + 2 * Pmax * j / (Nv - 1);
    // Depth-based shading: rows farther back are dimmer.
    const depth = j / (Nv - 1);
    const alpha = 0.35 + 0.55 * depth;
    ctx.strokeStyle = `rgba(120,200,255,${alpha.toFixed(2)})`;
    ctx.beginPath();
    for (let i = 0; i < Nv; i++) {
      const px = -Pmax + 2 * Pmax * i / (Nv - 1);
      const [sx, sy] = project(px, py, F[j*Nv + i]);
      if (i === 0) ctx.moveTo(sx, sy); else ctx.lineTo(sx, sy);
    }
    ctx.stroke();
  }
  // Lines of constant Px (j varies).
  for (let i = 0; i < Nv; i++) {
    const px = -Pmax + 2 * Pmax * i / (Nv - 1);
    const depth = i / (Nv - 1);
    const alpha = 0.35 + 0.55 * depth;
    ctx.strokeStyle = `rgba(120,200,255,${alpha.toFixed(2)})`;
    ctx.beginPath();
    for (let j = 0; j < Nv; j++) {
      const py = -Pmax + 2 * Pmax * j / (Nv - 1);
      const [sx, sy] = project(px, py, F[j*Nv + i]);
      if (j === 0) ctx.moveTo(sx, sy); else ctx.lineTo(sx, sy);
    }
    ctx.stroke();
  }

  // ── ∇𝓛 vector on the floor ──
  // Drawn as an arrow from the origin in the (P_x, P_y) plane at z = Fmin.
  // Length scaled to ~40% of P-range for visibility, independent of |∇𝓛|.
  if (gMag > 0) {
    const arrowLen = Pmax * 0.7;
    const ghx = gx_cell / gMag * arrowLen;
    const ghy = gy_cell / gMag * arrowLen;
    const [oax, oay] = project(0, 0, Fmin);
    const [tipx, tipy] = project(ghx, ghy, Fmin);
    ctx.strokeStyle = 'rgb(255,200,40)';
    ctx.lineWidth = 2;
    ctx.beginPath(); ctx.moveTo(oax, oay); ctx.lineTo(tipx, tipy); ctx.stroke();
    // Arrowhead.
    const ang = Math.atan2(tipy - oay, tipx - oax);
    const ah = 8;
    ctx.fillStyle = 'rgb(255,200,40)';
    ctx.beginPath();
    ctx.moveTo(tipx, tipy);
    ctx.lineTo(tipx - ah*Math.cos(ang - Math.PI/7), tipy - ah*Math.sin(ang - Math.PI/7));
    ctx.lineTo(tipx - ah*Math.cos(ang + Math.PI/7), tipy - ah*Math.sin(ang + Math.PI/7));
    ctx.closePath(); ctx.fill();
    ctx.fillStyle = 'rgb(255,200,40)';
    ctx.font = '11px ui-monospace, monospace';
    ctx.fillText('∇𝓛', tipx + 6, tipy - 4);
  }

  // ── Bead ── 3D position (Px, Py, F at that point).
  // Use the SAME display-rescaled F as the surface so the bead sits on it.
  const Pxc = fNow.Px[cellIdx], Pyc = fNow.Py[cellIdx];
  const F_bead = F_display(Pxc, Pyc);
  const [bsx, bsy] = project(Pxc, Pyc, F_bead);
  // Connector from bead down to floor.
  const [bfx, bfy] = project(Pxc, Pyc, Fmin);
  ctx.strokeStyle = 'rgba(255,255,255,0.45)';
  ctx.setLineDash([3, 3]);
  ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(bsx, bsy); ctx.lineTo(bfx, bfy); ctx.stroke();
  ctx.setLineDash([]);
  // The bead itself.
  ctx.beginPath();
  ctx.arc(bsx, bsy, 7, 0, 2*Math.PI);
  ctx.fillStyle = 'rgb(255,80,30)';
  ctx.fill();
  ctx.strokeStyle = 'rgb(255,255,255)';
  ctx.lineWidth = 1.8;
  ctx.stroke();

  // Label.
  ctx.fillStyle = '#ddd';
  ctx.font = '12px ui-monospace, monospace';
  ctx.fillText(`F(P)  —  cell #${cellIdx}    𝓛 = ${L_cell.toFixed(2)}`, 8, 14);
  const tiltNote = tiltScale < 0.999
    ? `    tilt × ${tiltScale.toExponential(2)} for visibility`
    : '';
  ctx.fillText(`F_disp ∈ [${Fmin.toFixed(2)}, ${Fmax.toFixed(2)}]    |∇𝓛| = ${gMag.toFixed(2)}${tiltNote}`, 8, 28);
  ctx.fillText(`λ=${(lam||0).toFixed(2)}  ν=${(nu||0).toFixed(2)}  κ=${(kap||0).toFixed(2)}`, 8, 42);
}

/**
 * 2D phase plane (P_x, P_y) with:
 *   • polar grid + thin contour rings of F(|P|)
 *   • dashed white ring of minima at |P|_eq
 *   • trajectory polyline up to current frame
 *   • current bead
 * No background colour shading — structure is conveyed by lines only.
 */
export function drawBead2D(canvasId, frames, idx, cellIdx, params) {
  const canvas = document.getElementById(canvasId);
  if (!canvas) return;
  const ctx = autoFit(canvas);
  const w = canvas.clientWidth, h = canvas.clientHeight;

  if (!frames || frames.length === 0 || !frames[0].Px) {
    drawFrame(ctx, makeAxis({ xMin:-1, xMax:1, yMin:-1, yMax:1, w, h, aspect:1 }));
    return;
  }
  const N = frames[0].Px.length;
  cellIdx = Math.max(0, Math.min(N - 1, cellIdx | 0));

  const fNow = frames[idx];
  const L_cell = sampleLAtCell(fNow, cellIdx, params.R_dish, params.N_grid);
  const { lam, nu, kap, L_c } = params;
  const gx_cell = fNow.Gx ? fNow.Gx[cellIdx] : 0;
  const gy_cell = fNow.Gy ? fNow.Gy[cellIdx] : 0;
  const gMag = Math.hypot(gx_cell, gy_cell);

  // Same Pmax / tilt-rescaling logic as the 3D plot, so the two panels show
  // the same geometry.
  const pEq = gleEqRadius(L_cell, L_c, lam, nu);
  const Pmax = Math.max(autoPmaxForCell(frames, idx, cellIdx), 1.5 * pEq, 0.5);
  const wellDepth = gleWellDepth(L_cell, L_c, lam, nu, Pmax);
  const tiltScale = tiltDisplayScale(kap, gMag, Pmax, wellDepth, 0.6);
  const ax = makeAxis({ xMin:-Pmax, xMax:Pmax, yMin:-Pmax, yMax:Pmax, w, h, aspect:1 });

  // Dark canvas behind the plot for line contrast.
  ctx.save();
  ctx.fillStyle = '#15151c';
  ctx.fillRect(ax.padL, ax.padT, ax.plotW, ax.plotH);
  ctx.restore();
  drawFrame(ctx, ax);

  // Compute F (with DISPLAY-rescaled tilt) on a 2D mesh — contours are no
  // longer concentric circles, so we use marching squares to extract level sets.
  function F_display_2d(px, py) {
    const gl   = freeEnergy(px, py, L_cell, L_c, lam, nu, 0, 0, 0);
    const tilt = -(kap || 0) * (gx_cell * px + gy_cell * py);
    return gl + tiltScale * tilt;
  }
  const Nm = 90;
  const Fmesh = new Float32Array(Nm * Nm);
  let Fmin = Infinity, Fmax = -Infinity;
  for (let j = 0; j < Nm; j++) {
    const py = -Pmax + 2 * Pmax * j / (Nm - 1);
    for (let i = 0; i < Nm; i++) {
      const px = -Pmax + 2 * Pmax * i / (Nm - 1);
      const v = F_display_2d(px, py);
      Fmesh[j*Nm + i] = v;
      if (v < Fmin) Fmin = v;
      if (v > Fmax) Fmax = v;
    }
  }
  const cx_ax = ax.xToPx(0), cy_ax = ax.yToPx(0);
  const pxPerUnit = ax.xToPx(1) - ax.xToPx(0);

  // Marching-squares contours at 8 evenly-spaced levels.
  const nLevels = 8;
  for (let l = 0; l < nLevels; l++) {
    const t = (l + 0.5) / nLevels;
    const level = Fmin + t * (Fmax - Fmin);
    const colour = (t < 0.5)
      ? `rgba(255,170,90,${0.20 + 0.50*(0.5 - t)})`   // below mean — warm
      : `rgba(150,210,255,${0.20 + 0.50*(t - 0.5)})`; // above mean — cool
    ctx.save();
    ctx.strokeStyle = colour;
    ctx.lineWidth = 0.9;
    ctx.beginPath();
    emitContourMS(ctx, ax, Fmesh, Nm, Pmax, level);
    ctx.stroke();
    ctx.restore();
  }

  // Polar grid: faint radial spokes every 30°.
  ctx.save();
  ctx.strokeStyle = 'rgba(120,130,150,0.18)';
  ctx.lineWidth = 0.6;
  for (let k = 0; k < 12; k++) {
    const ang = (k / 12) * 2 * Math.PI;
    const r = Pmax * 0.98;
    ctx.beginPath();
    ctx.moveTo(cx_ax, cy_ax);
    ctx.lineTo(cx_ax + Math.cos(ang) * r * pxPerUnit,
               cy_ax - Math.sin(ang) * r * pxPerUnit);
    ctx.stroke();
  }
  ctx.restore();

  // Dashed white ring at the well minimum (current 𝓛, no-tilt case).
  // In intrinsic units: |P|_eq = √(λ(𝓛-𝓛_c)/ν). Guard for negative argument.
  if (L_cell > L_c) {
    const pEqSq = lam * (L_cell - L_c) / Math.max(nu, 1e-12);
    if (pEqSq > 0) {
      const pEqRing = Math.sqrt(pEqSq);
      if (pEqRing > 0 && pEqRing < Pmax) {
        ctx.save();
        ctx.beginPath();
        ctx.arc(cx_ax, cy_ax, Math.abs(pEqRing * pxPerUnit), 0, 2*Math.PI);
        ctx.strokeStyle = 'rgba(255,255,255,0.85)';
        ctx.lineWidth = 1.6;
        ctx.setLineDash([5, 4]);
        ctx.stroke();
        ctx.restore();
      }
    }
  }

  // ∇𝓛 arrow from origin in the gradient direction.
  if (gMag > 0) {
    const arrowLen = Pmax * 0.6;
    arrow2D(ctx, ax, 0, 0, gx_cell/gMag * arrowLen, gy_cell/gMag * arrowLen);
    // Label near the arrowhead.
    const tipx = ax.xToPx(gx_cell/gMag * arrowLen) + 5;
    const tipy = ax.yToPx(gy_cell/gMag * arrowLen) - 4;
    ctx.save();
    ctx.fillStyle = 'rgb(255,200,40)';
    ctx.font = '11px ui-monospace, monospace';
    ctx.fillText('∇𝓛', tipx, tipy);
    ctx.restore();
  }

  // Trajectory polyline of the tracked cell up to current frame.
  const xs = [], ys = [];
  for (let i = 0; i <= idx; i++) {
    if (frames[i].Px) {
      xs.push(frames[i].Px[cellIdx]);
      ys.push(frames[i].Py[cellIdx]);
    }
  }
  if (xs.length > 1) {
    strokePath(ctx, ax, xs, ys, { color: 'rgba(230,230,240,0.85)', width: 1.2 });
  }

  // Current bead.
  if (fNow.Px) {
    dot(ctx, ax, fNow.Px[cellIdx], fNow.Py[cellIdx], 5, 'rgb(255,80,30)');
  }

  // Label.
  ctx.save();
  ctx.fillStyle = '#ddd';
  ctx.font = '12px ui-monospace, monospace';
  const tiltNote = tiltScale < 0.999
    ? `    tilt × ${tiltScale.toExponential(2)} for visibility`
    : '';
  ctx.fillText(`(Pₓ, Pᵧ) — cell #${cellIdx}    𝓛 = ${L_cell.toFixed(2)}${tiltNote}`, ax.padL + 6, ax.padT + 14);
  ctx.fillText(`λ=${(lam||0).toFixed(2)}  ν=${(nu||0).toFixed(2)}  κ=${(kap||0).toFixed(2)}    |∇𝓛| = ${gMag.toFixed(2)}`,
               ax.padL + 6, ax.padT + 28);
  ctx.restore();
}

/**
 * 1D cross-section of F along the gradient direction.
 *   F(s) at P = s · ∇̂𝓛  =  −[λ(𝓛-𝓛_c) s²/2 + ν(s⁴/4 − s⁶/6)]  −  κ |∇𝓛| s
 * The chemotactic tilt term lowers one minimum and raises the other — that's
 * the broken-symmetry direction along which the cell preferentially polarizes.
 * Also shows the bead's projection onto this axis, s_cell = P · ∇̂𝓛.
 */
export function drawBead1D(canvasId, frames, idx, cellIdx, params) {
  const canvas = document.getElementById(canvasId);
  if (!canvas) return;
  const ctx = autoFit(canvas);
  const w = canvas.clientWidth, h = canvas.clientHeight;

  if (!frames || frames.length === 0 || !frames[0].Px) {
    drawFrame(ctx, makeAxis({ xMin:-1, xMax:1, yMin:-1, yMax:1, w, h }));
    return;
  }
  const N = frames[0].Px.length;
  cellIdx = Math.max(0, Math.min(N - 1, cellIdx | 0));

  const Pmax = autoPmaxForCell(frames, idx, cellIdx);
  const fNow = frames[idx];
  const L_cell = sampleLAtCell(fNow, cellIdx, params.R_dish, params.N_grid);
  const { lam, nu, kap, L_c } = params;
  const gx_cell = fNow.Gx ? fNow.Gx[cellIdx] : 0;
  const gy_cell = fNow.Gy ? fNow.Gy[cellIdx] : 0;
  const gMag = Math.hypot(gx_cell, gy_cell);

  // Sample F along s ∈ [-Pmax, Pmax] along ĝ.
  const ghx = gMag > 0 ? gx_cell / gMag : 1;
  const ghy = gMag > 0 ? gy_cell / gMag : 0;
  const Ns = 401;
  const ss = new Float64Array(Ns);
  const Fs = new Float64Array(Ns);
  let Fmin = Infinity, Fmax = -Infinity;
  for (let k = 0; k < Ns; k++) {
    const s = -Pmax + 2 * Pmax * k / (Ns - 1);
    const px = s * ghx, py = s * ghy;
    const v = freeEnergy(px, py, L_cell, L_c, lam, nu, kap, gx_cell, gy_cell);
    ss[k] = s; Fs[k] = v;
    if (v < Fmin) Fmin = v;
    if (v > Fmax) Fmax = v;
  }
  // Add a little vertical padding so the curve doesn't hug the axes.
  const Fpad = (Fmax - Fmin) * 0.08 || 1;
  const ax = makeAxis({
    xMin: -Pmax, xMax: Pmax,
    yMin: Fmin - Fpad, yMax: Fmax + Fpad,
    w, h,
  });
  drawFrame(ctx, ax);

  // Zero-axis horizontal & vertical for reference.
  ctx.save();
  ctx.strokeStyle = 'rgba(140,140,160,0.5)';
  ctx.lineWidth = 0.7;
  ctx.setLineDash([3, 3]);
  if (ax.yMin < 0 && ax.yMax > 0) {
    ctx.beginPath();
    ctx.moveTo(ax.xToPx(ax.xMin), ax.yToPx(0));
    ctx.lineTo(ax.xToPx(ax.xMax), ax.yToPx(0));
    ctx.stroke();
  }
  if (ax.xMin < 0 && ax.xMax > 0) {
    ctx.beginPath();
    ctx.moveTo(ax.xToPx(0), ax.yToPx(ax.yMin));
    ctx.lineTo(ax.xToPx(0), ax.yToPx(ax.yMax));
    ctx.stroke();
  }
  ctx.restore();

  // F(s) curve.
  strokePath(ctx, ax, Array.from(ss), Array.from(Fs), {
    color: 'rgb(120,200,255)', width: 1.8,
  });

  // Bead position s_cell, F(s_cell).
  if (fNow.Px) {
    const Pxc = fNow.Px[cellIdx], Pyc = fNow.Py[cellIdx];
    const s_cell = Pxc * ghx + Pyc * ghy;
    const F_at = freeEnergy(Pxc, Pyc, L_cell, L_c, lam, nu, kap, gx_cell, gy_cell);
    // Vertical drop line from bead down to the curve / down to F=0 baseline.
    ctx.save();
    ctx.strokeStyle = 'rgba(255,80,30,0.55)';
    ctx.lineWidth = 1;
    ctx.setLineDash([2, 3]);
    ctx.beginPath();
    ctx.moveTo(ax.xToPx(s_cell), ax.yToPx(ax.yMin));
    ctx.lineTo(ax.xToPx(s_cell), ax.yToPx(F_at));
    ctx.stroke();
    ctx.restore();
    dot(ctx, ax, s_cell, F_at, 5, 'rgb(255,80,30)');
  }

  // Label: which direction the slice is along.
  ctx.save();
  ctx.fillStyle = '#ddd';
  ctx.font = '12px ui-monospace, monospace';
  const dirTxt = gMag > 0
    ? `slice along ∇̂𝓛 = (${ghx.toFixed(2)}, ${ghy.toFixed(2)})    |∇𝓛| = ${gMag.toFixed(2)}`
    : `∇𝓛 ≈ 0 — symmetric well (no tilt)`;
  ctx.fillText(dirTxt, ax.padL + 6, ax.padT + 14);
  ctx.fillText(`cell #${cellIdx}    𝓛 = ${L_cell.toFixed(2)}    κ = ${(kap||0).toFixed(2)}`,
               ax.padL + 6, ax.padT + 30);
  ctx.restore();
}

// ─── mean radius plot ────────────────────────────────────────────────────────

/**
 * Draw ⟨r̃⟩_free(t̃) over accumulated frames.
 * @param {string} canvasId
 * @param {Array} frames - array of { t, agentX, agentY, emitting }
 * @param {Object} params - { R_dish }
 */
export function drawMeanRadius(canvasId, frames, params) {
  const canvas = document.getElementById(canvasId);
  if (!canvas || !frames || frames.length === 0) return;
  const ctx = autoFit(canvas);
  const w = canvas.clientWidth, h = canvas.clientHeight;

  const ts = [], rs = [];
  for (const f of frames) {
    if (!f.agentX) continue;
    let rsum = 0, n = 0;
    for (let i = 0; i < f.agentX.length; i++) {
      rsum += Math.hypot(f.agentX[i], f.agentY[i]);
      n++;
    }
    ts.push(f.t);
    rs.push(n > 0 ? rsum / n : 0);
  }

  const t_max = ts[ts.length - 1] || 1;
  let rmax = 1;
  for (const r of rs) if (r > rmax) rmax = r;

  const ax = makeAxis({ xMin: 0, xMax: t_max, yMin: 0, yMax: rmax * 1.1, w, h });
  drawFrame(ctx, ax);
  strokePath(ctx, ax, ts, rs, { color: '#2b6cb0', width: 1.5 });
}

/**
 * c(σ̃) density-sweep plot (M6.1 / M6.2): MEASURED front speeds only, on a
 * log-σ̃ axis.
 *
 * No analytical c(σ̃) curve is drawn. The mean-field pushed-front flux balance
 * (c²𝓛_c = α_eff σ̃ and its 2D–3D counterpart) does not describe this system
 * reliably, so plotting it invited reading the simulation as "off by X%" from
 * a law that was never right. The sweep reports what the ABM does; the density
 * dependence is read off the points themselves.
 *
 * All quantities nondimensional.
 *
 * @param {string} canvasId
 * @param {Object} d - { points, sigma_min, sigma_max, sigma_current }
 */
export function drawCsweep(canvasId, d) {
  const canvas = document.getElementById(canvasId);
  if (!canvas) return;
  const ctx = autoFit(canvas);
  const w = canvas.clientWidth, h = canvas.clientHeight;

  const smin = Math.max(1e-4, d.sigma_min || 0.1);
  const smax = Math.max(smin * 1.01, d.sigma_max || 10);

  const pts = (d.points || []).filter(p => isFinite(p.c));

  let ymax = 1e-6;
  for (const p of pts) if (p.c > ymax) ymax = p.c;
  ymax *= 1.25;

  const ax = makeAxis({ xMin: smin, xMax: smax, yMin: 0, yMax: ymax, w, h, logX: true });
  drawFrame(ctx, ax);

  // Current single-run σ̃: thin blue vertical marker (not a prediction — just
  // where the single run on this page sits along the swept axis).
  const sc = d.sigma_current;
  if (isFinite(sc) && sc > smin && sc < smax) {
    ctx.save();
    ctx.strokeStyle = 'rgba(43,108,176,0.4)';
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    ctx.moveTo(ax.xToPx(sc), ax.padT);
    ctx.lineTo(ax.xToPx(sc), ax.padT + ax.plotH);
    ctx.stroke();
    ctx.restore();
  }

  // Measured points, colour-coded by measurement status:
  //   ok           orange filled — a timed relay front
  //   noignite     red filled    — no self-sustaining wave (physical)
  //   unresolved   hollow amber  — a wave too slow to clear the ignition halo
  //                                within t̃_max (a window limit, not a zero)
  //   unmeasurable hollow gray   — the halo itself floods the dish (setup limit)
  for (const p of pts) {
    if (p.sigma < smin || p.sigma > smax) continue;
    const st = p.status || (p.c > 0 ? 'ok' : 'noignite');
    const y  = Math.min(p.c, ymax);
    if (st === 'ok')            { dot(ctx, ax, p.sigma, y, 3.5, '#d9730d'); continue; }
    if (st === 'noignite')      { dot(ctx, ax, p.sigma, y, 3.5, '#c0392b'); continue; }
    ctx.save();
    ctx.strokeStyle = (st === 'unresolved') ? '#d9a441' : '#8c8c8c';
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    ctx.arc(ax.xToPx(p.sigma), ax.yToPx(y), 3.5, 0, 2 * Math.PI);
    ctx.stroke();
    ctx.restore();
  }
}
