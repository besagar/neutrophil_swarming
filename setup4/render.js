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

// ─── viridis colormap (25-stop sampled) ─────────────────────────────────────
// Each entry: [r, g, b] in [0,255].
const VIRIDIS = [
  [68,1,84],[71,13,96],[72,28,109],[70,42,120],[67,56,130],
  [60,69,137],[53,82,142],[46,95,143],[40,107,142],[35,119,140],
  [30,131,137],[31,143,133],[40,155,126],[53,166,118],[70,177,106],
  [92,186,90],[116,195,72],[141,202,53],[166,209,35],[192,215,19],
  [217,220,7],[240,225,2],[253,231,37],[253,240,35],[253,231,37],
];

function viridis(t) {
  const v = Math.max(0, Math.min(1, t));
  const n = VIRIDIS.length - 1;
  const fi = v * n;
  const lo = Math.floor(fi), hi = Math.min(n, lo + 1);
  const f  = fi - lo;
  const r  = VIRIDIS[lo][0] + f * (VIRIDIS[hi][0] - VIRIDIS[lo][0]);
  const g  = VIRIDIS[lo][1] + f * (VIRIDIS[hi][1] - VIRIDIS[lo][1]);
  const b  = VIRIDIS[lo][2] + f * (VIRIDIS[hi][2] - VIRIDIS[lo][2]);
  return `rgb(${r|0},${g|0},${b|0})`;
}

// Pre-build a 256-entry viridis LUT as Uint8Array[256*3] for fast heatmap rendering.
const VIRIDIS_LUT = new Uint8Array(256 * 3);
for (let i = 0; i < 256; i++) {
  const t = i / 255;
  const n = VIRIDIS.length - 1;
  const fi = t * n;
  const lo = Math.floor(fi), hi = Math.min(n, lo + 1);
  const f  = fi - lo;
  VIRIDIS_LUT[i * 3]     = VIRIDIS[lo][0] + f * (VIRIDIS[hi][0] - VIRIDIS[lo][0]);
  VIRIDIS_LUT[i * 3 + 1] = VIRIDIS[lo][1] + f * (VIRIDIS[hi][1] - VIRIDIS[lo][1]);
  VIRIDIS_LUT[i * 3 + 2] = VIRIDIS[lo][2] + f * (VIRIDIS[hi][2] - VIRIDIS[lo][2]);
}

// ─── dish / heatmap ──────────────────────────────────────────────────────────

/**
 * Draw the petri dish: L heatmap + circular dish boundary + cell dots.
 * @param {string} canvasId
 * @param {Object} frame - { Lfield?, agentX, agentY, emitting, agentPx?, agentPy? }
 * @param {Object} params - { N_grid, R_dish, t, L_max_display? }
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

  // Fill the plot area with a dark background — areas outside the dish bbox
  // (or before a heatmap has been received) sit on this rather than canvas white.
  ctx.save();
  ctx.fillStyle = 'rgb(20,20,30)';
  ctx.fillRect(ax.padL, ax.padT, ax.plotW, ax.plotH);
  ctx.restore();

  // ── L heatmap ──
  // frame.Lfield is a Uint8Array already normalized 0..255 per-frame in the
  // worker. Render at native grid resolution (N×N) into an offscreen canvas
  // then drawImage-scale it to the plot area. Using drawImage (not
  // putImageData) is critical: putImageData IGNORES ctx.setTransform, so on
  // a Retina display the heatmap lands at half-scale in the upper-left of
  // the backing buffer. drawImage honors the dpr transform from autoFit.
  if (frame.Lfield) {
    const N = params.N_grid;
    const Lf = frame.Lfield;
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
          buf[idx] = 20; buf[idx+1] = 20; buf[idx+2] = 30; buf[idx+3] = 255;
        } else {
          const li = Lf[k_field];
          buf[idx]   = VIRIDIS_LUT[li * 3];
          buf[idx+1] = VIRIDIS_LUT[li * 3 + 1];
          buf[idx+2] = VIRIDIS_LUT[li * 3 + 2];
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
    ctx.strokeStyle = 'rgba(180,180,180,0.7)'; ctx.lineWidth = 1.5; ctx.stroke();
    ctx.restore();
  }

  // ── cell dots ──
  if (frame.agentX && frame.agentY && frame.emitting) {
    const N = frame.agentX.length;
    // Slightly larger dots + dark stroke so cells stay visible against any
    // viridis value (especially the bright-yellow saturated regions where
    // a translucent fill would otherwise vanish).
    const dotR = Math.max(2.0, ax.plotW / (params.R_dish * 60));
    ctx.lineWidth = Math.max(0.6, dotR * 0.35);
    ctx.strokeStyle = 'rgba(15,15,25,0.9)';
    const tracked = (params.trackedCellIdx != null) ? (params.trackedCellIdx | 0) : -1;
    for (let i = 0; i < N; i++) {
      if (i === tracked) continue;  // draw tracked cell last (on top)
      const px = ax.xToPx(frame.agentX[i]);
      const py = ax.yToPx(frame.agentY[i]);
      const isEmit = frame.emitting[i] === 1;
      // Emitting: saturated red-orange. Inactive: cool light grey-blue, fully
      // opaque so it doesn't blend into a saturated heatmap.
      ctx.fillStyle = isEmit ? 'rgb(245,70,30)' : 'rgb(210,220,230)';
      ctx.beginPath();
      ctx.arc(px, py, dotR, 0, 2*Math.PI);
      ctx.fill();
      ctx.stroke();
    }
    // Tracked cell: 2x size, cyan fill, thick white stroke + crosshair.
    if (tracked >= 0 && tracked < N) {
      const px = ax.xToPx(frame.agentX[tracked]);
      const py = ax.yToPx(frame.agentY[tracked]);
      const Rt = dotR * 2.4;
      ctx.beginPath();
      ctx.arc(px, py, Rt, 0, 2*Math.PI);
      ctx.fillStyle = 'rgb(80,220,255)';
      ctx.fill();
      ctx.lineWidth = 2;
      ctx.strokeStyle = 'rgb(255,255,255)';
      ctx.stroke();
      // Crosshair tick marks so it's easy to spot in a dense field.
      ctx.beginPath();
      ctx.moveTo(px - Rt*2.0, py); ctx.lineTo(px - Rt*1.0, py);
      ctx.moveTo(px + Rt*1.0, py); ctx.lineTo(px + Rt*2.0, py);
      ctx.moveTo(px, py - Rt*2.0); ctx.lineTo(px, py - Rt*1.0);
      ctx.moveTo(px, py + Rt*1.0); ctx.lineTo(px, py + Rt*2.0);
      ctx.strokeStyle = 'rgba(80,220,255,0.95)';
      ctx.lineWidth = 1.3;
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
 * limit. The nondim PDE  ∂_{t̃} 𝓛 = ∇̃²𝓛 + β·(2/h̃₀) δ(z̃) H⁺(𝓛;1;n)
 * with n → ∞ admits a traveling planar front at c̃ = β·(2/π); the z=0
 * profile in ξ ≡ r̃ − r̃_front is (β=1 reference Dieterle):
 *     𝓛(ξ) = √(−2ξ)                              for ξ < 0   (behind)
 *     𝓛(ξ) = (π/2) · exp(−2ξ/π) / √(2ξ)         for ξ > 0   (ahead)
 *
 * The β > 1 case maps to the β=1 problem under x→x/β, t→t/β² (z is mapped
 * the same as x; 𝓛 unchanged, threshold unchanged). So in original Setup 4
 * units we substitute ξ → β·ξ in both branches. β = ρ̃_sim/σ̃ — the ratio
 * of simulated cell density to the σ̃ implied by the dim parameters.
 */
function dieterleProfile(xi, beta = 1) {
  const z = beta * xi;
  if (z < 0) return Math.sqrt(-2 * z);
  return (Math.PI / 2) * Math.exp(-2 * z / Math.PI) / Math.sqrt(2 * z);
}

/**
 * Draw 1D azimuthal average 𝓛(r̃).
 * @param {string} canvasId
 * @param {Object} frame - { radialProfile: Float32Array }
 * @param {Object} params - { N_grid, R_dish, model?, geometry? }
 *   When model === 'M1' && geometry === '2d3d', overlay the Dieterle
 *   analytical planar-front profile (dashed) anchored at the simulated
 *   r̃_front (where 𝓛 crosses 1).
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

  strokePath(ctx, ax, xs, ys, { color: '#2b6cb0', width: 1.8 });

  // ── Dieterle analytical overlay (M1 + 2D-3D only) ──
  // Planar-front solution rescaled for the actual simulated cell density:
  // β = (N / π R̃²) / σ̃ — the per-area cell density in the simulation
  // divided by the σ̃ that defines the nondim units. β=1 is the canonical
  // continuum limit; β≠1 means the simulated swarm is denser/sparser than
  // what the dim parameters assume, which rescales the wave.
  if (params.model === 'M1' && params.geometry === '2d3d') {
    const N_cells   = (params.N_cells || 0) | 0;
    const sigma_til = params.sigma_tilde || 0;
    const beta = (N_cells > 0 && sigma_til > 0)
      ? N_cells / (Math.PI * R * R * sigma_til)
      : 1;
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
      // Sample the analytical profile, excluding a small buffer ±0.5/β around
      // the front (both branches diverge there). The buffer scales with 1/β
      // because the front layer is squeezed by the same factor.
      const buf = 0.5 / Math.max(beta, 1e-6);
      const NA = 400;
      const style = { color: 'rgba(20,20,20,0.85)', width: 1.6, dash: [5, 4] };
      // Behind the front: ξ ∈ [−r_front, −buf], i.e. r ∈ [0, r_front − buf].
      if (r_front > buf) {
        const xs_b = [], ys_b = [];
        for (let k = 0; k < NA; k++) {
          const xi = -r_front + (-buf - (-r_front)) * (k / (NA - 1));
          const r = r_front + xi;
          if (r < 0) continue;
          const Ly = dieterleProfile(xi, beta);
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
          const Ly = dieterleProfile(xi, beta);
          if (Ly <= ymax * 1.5 && isFinite(Ly) && Ly > 0) { xs_a.push(r); ys_a.push(Ly); }
        }
        if (xs_a.length > 1) strokePath(ctx, ax, xs_a, ys_a, style);
      }

      // Annotate analytical wave speed (rescaled by β). Warn when the front
      // has reached the dish boundary — past that point the traveling-wave
      // description breaks down and 𝓛 grows by half-space diffusion from a
      // saturated z=0 source (≈ 2β√(t̃/π) at the centre).
      const cTheory = beta * (2 / Math.PI);
      const saturated = (r_front >= R * 0.95);
      ctx.save();
      ctx.fillStyle = saturated ? '#a33' : '#222';
      ctx.font = '11px ui-monospace, monospace';
      const baseLbl = (Math.abs(beta - 1) < 1e-3)
        ? `Dieterle (c̃ = 2/π ≈ ${cTheory.toFixed(3)})`
        : `Dieterle, β=${beta.toFixed(2)}  (c̃ = ${cTheory.toFixed(3)})`;
      ctx.fillText(baseLbl, ax.padL + 6, ax.padT + 14);
      if (saturated) {
        ctx.fillText('⚠ wave hit boundary — traveling-wave profile invalid', ax.padL + 6, ax.padT + 28);
      }
      ctx.restore();
    }
  }
}

// ─── bead-in-free-energy plots ────────────────────────────────────────────────
//
// Free energy:
//   F(|P|) = -Λ [(𝓛 − 𝓛_c) |P|² / 2  +  λ (|P|⁴/4 − |P|⁶/6)]
// Ring of minima at |P|² = (1 + √(1 + 4(𝓛-𝓛_c)/λ)) / 2  when 𝓛 > 𝓛_c.
// Two views below: drawBead3D (wireframe surface + bead), drawBead2D
// (contour-ring phase plane + trajectory + bead).

// Effective free energy seen by one cell.
// Includes the GL Mexican-hat AND the chemotactic linear tilt −χ̃·(∇𝓛·P)
// which breaks the rotational symmetry of the well.
function freeEnergy(Px, Py, L, Lc, lam, Lambda, chi, gx, gy) {
  const p2 = Px*Px + Py*Py;
  const p4 = p2*p2;
  const p6 = p4*p2;
  const gl   = -Lambda * ((L - Lc) * p2 / 2 + lam * (p4/4 - p6/6));
  const tilt = -(chi || 0) * ((gx || 0) * Px + (gy || 0) * Py);
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

// |P|_eq from the symmetric GL well; 0 if sub-threshold.
function gleEqRadius(L, Lc, lam) {
  if (L <= Lc) return 0;
  const p2 = (1 + Math.sqrt(1 + 4 * (L - Lc) / lam)) / 2;
  return Math.sqrt(p2);
}

// Magnitude of the symmetric GL well, measured between F(P=0)=0 and F(|P|_eq).
// Used as a "natural scale" for the display rescaling of the chemotactic tilt.
function gleWellDepth(L, Lc, lam, Lambda, fallbackPmax) {
  if (L > Lc) {
    const p2 = (1 + Math.sqrt(1 + 4 * (L - Lc) / lam)) / 2;
    const p4 = p2 * p2;
    const p6 = p4 * p2;
    const F  = -Lambda * ((L - Lc) * p2 / 2 + lam * (p4/4 - p6/6));
    return Math.max(Math.abs(F), 1e-6);
  }
  // Sub-threshold: F is a parabolic well at P = 0. Use its rise to fallback Pmax.
  const half = Math.abs(Lambda * (Lc - L) * fallbackPmax * fallbackPmax / 2);
  return Math.max(half, 1e-6);
}

// Compute a tilt-display rescaling factor: shrink the chemotactic-tilt
// contribution so its magnitude across the visible Pmax stays comparable to
// the GL well depth. Returns 1.0 when the tilt is naturally moderate (no
// rescaling needed). The 1D cross-section plot deliberately does NOT use this —
// it shows the true tilted potential along ∇𝓛.
function tiltDisplayScale(chi, gMag, Pmax, wellDepth, frac = 0.6) {
  const tiltMax = (chi || 0) * (gMag || 0) * Pmax;
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
  const { Lambda, L_c, lam, chi } = params;
  const gx_cell = fNow.Gx ? fNow.Gx[cellIdx] : 0;
  const gy_cell = fNow.Gy ? fNow.Gy[cellIdx] : 0;
  const gMag = Math.hypot(gx_cell, gy_cell);

  // Pmax: always wide enough to contain the GL well minimum, even when the
  // bead is still near the origin. Without this, the surface auto-zooms to a
  // sub-|P|_eq region where the GL well isn't yet developed.
  const pEq = gleEqRadius(L_cell, L_c, lam);
  const Pmax = Math.max(autoPmaxForCell(frames, idx, cellIdx), 1.5 * pEq, 0.5);

  // Chemotactic-tilt rescaling for display.
  // When χ̃·|∇𝓛|·Pmax ≫ well-depth (typical in saturated regimes here), the
  // tilt term swamps the GL surface and the Mexican-hat geometry is invisible.
  // We shrink the displayed tilt so it remains comparable to the well depth;
  // the actual rescaling factor is reported in the label. The bead's height
  // and the 1D cross-section both use the un-rescaled tilt.
  const wellDepth  = gleWellDepth(L_cell, L_c, lam, Lambda, Pmax);
  const tiltScale  = tiltDisplayScale(chi, gMag, Pmax, wellDepth, 0.6);
  function F_display(px, py) {
    const gl   = freeEnergy(px, py, L_cell, L_c, lam, Lambda, 0, 0, 0); // GL only
    const tilt = -(chi || 0) * (gx_cell * px + gy_cell * py);
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
  const { Lambda, L_c, lam, chi } = params;
  const gx_cell = fNow.Gx ? fNow.Gx[cellIdx] : 0;
  const gy_cell = fNow.Gy ? fNow.Gy[cellIdx] : 0;
  const gMag = Math.hypot(gx_cell, gy_cell);

  // Same Pmax / tilt-rescaling logic as the 3D plot, so the two panels show
  // the same geometry.
  const pEq = gleEqRadius(L_cell, L_c, lam);
  const Pmax = Math.max(autoPmaxForCell(frames, idx, cellIdx), 1.5 * pEq, 0.5);
  const wellDepth = gleWellDepth(L_cell, L_c, lam, Lambda, Pmax);
  const tiltScale = tiltDisplayScale(chi, gMag, Pmax, wellDepth, 0.6);
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
    const gl   = freeEnergy(px, py, L_cell, L_c, lam, Lambda, 0, 0, 0);
    const tilt = -(chi || 0) * (gx_cell * px + gy_cell * py);
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
  if (L_cell > L_c) {
    const pEqSq = (1 + Math.sqrt(1 + 4 * (L_cell - L_c) / lam)) / 2;
    const pEq = Math.sqrt(pEqSq);
    if (pEq > 0 && pEq < Pmax) {
      ctx.save();
      ctx.beginPath();
      ctx.arc(cx_ax, cy_ax, Math.abs(pEq * pxPerUnit), 0, 2*Math.PI);
      ctx.strokeStyle = 'rgba(255,255,255,0.85)';
      ctx.lineWidth = 1.6;
      ctx.setLineDash([5, 4]);
      ctx.stroke();
      ctx.restore();
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
  ctx.restore();
}

/**
 * 1D cross-section of F along the gradient direction.
 *   F(s) at P = s · ∇̂𝓛  =  -Λ[(𝓛-𝓛_c) s²/2 + λ(s⁴/4 − s⁶/6)]  −  χ̃ |∇𝓛| s
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
  const { Lambda, L_c, lam, chi } = params;
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
    const v = freeEnergy(px, py, L_cell, L_c, lam, Lambda, chi, gx_cell, gy_cell);
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
    const F_at = freeEnergy(Pxc, Pyc, L_cell, L_c, lam, Lambda, chi, gx_cell, gy_cell);
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
  ctx.fillText(`cell #${cellIdx}    𝓛 = ${L_cell.toFixed(2)}    χ̃ = ${(chi||0).toFixed(2)}`,
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
