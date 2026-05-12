// Canvas2D plotting helpers — strictly nondimensional axes.

export function fitCanvas(canvas, cssWidth, cssHeight) {
  const dpr = window.devicePixelRatio || 1;
  canvas.style.width = cssWidth + 'px';
  canvas.style.height = cssHeight + 'px';
  canvas.width = Math.round(cssWidth * dpr);
  canvas.height = Math.round(cssHeight * dpr);
  const ctx = canvas.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return ctx;
}

const _fitCache = new WeakMap();
export function autoFit(canvas) {
  const r = canvas.getBoundingClientRect();
  const cssW = r.width, cssH = Math.max(120, r.height);
  const dpr = window.devicePixelRatio || 1;
  const bw = Math.round(cssW * dpr), bh = Math.round(cssH * dpr);
  const ctx = canvas.getContext('2d');
  const c = _fitCache.get(canvas);
  if (!c || c.bw !== bw || c.bh !== bh) {
    // Size changed: resize backing buffer (this clears the canvas and resets state).
    canvas.style.width = cssW + 'px'; canvas.style.height = cssH + 'px';
    canvas.width = bw; canvas.height = bh;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    _fitCache.set(canvas, { bw, bh, dpr });
  } else {
    // Same size: just reset transform and clear without touching the backing buffer.
    ctx.setTransform(c.dpr, 0, 0, c.dpr, 0, 0);
    ctx.clearRect(0, 0, cssW, cssH);
  }
  return ctx;
}

// Linear axis mapper. Padding leaves room for KaTeX label overlays
// (decoratePlot positions them in these pad regions).
// `aspect: 1` enforces equal x/y units-per-pixel by symmetrically expanding
// whichever data range would otherwise be over-sampled.
export function makeAxis({ xMin, xMax, yMin, yMax, w, h, padL = 44, padR = 10, padT = 18, padB = 30, aspect, logX = false }) {
  const plotW = w - padL - padR;
  const plotH = h - padT - padB;
  if (aspect === 1 && plotW > 0 && plotH > 0) {
    const sx = plotW / (xMax - xMin);
    const sy = plotH / (yMax - yMin);
    if (sx > sy) {
      const half = (plotW / sy) / 2;
      const c = (xMin + xMax) / 2;
      xMin = c - half; xMax = c + half;
    } else if (sy > sx) {
      const half = (plotH / sx) / 2;
      const c = (yMin + yMax) / 2;
      yMin = c - half; yMax = c + half;
    }
  }
  const xToPx = logX
    ? x => padL + (Math.log10(Math.max(x, 1e-300)) - Math.log10(xMin)) /
                  (Math.log10(xMax) - Math.log10(xMin)) * plotW
    : x => padL + ((x - xMin) / (xMax - xMin)) * plotW;
  const yToPx = y => padT + (1 - (y - yMin) / (yMax - yMin)) * plotH;
  return {
    xMin, xMax, yMin, yMax, w, h, padL, padR, padT, padB, plotW, plotH,
    logX, xToPx, yToPx,
  };
}

// Draws frame + ticks only. Axis labels and titles are LaTeX overlays
// rendered as DOM via decoratePlot() — never on canvas.
// Options (all default true): showGridX, showGridY, showTickLabelsX, showTickLabelsY.
export function drawFrame(ctx, ax, opts = {}) {
  const showGridX       = opts.showGridX       !== false;
  const showGridY       = opts.showGridY       !== false;
  const showTickLabelsX = opts.showTickLabelsX !== false;
  const showTickLabelsY = opts.showTickLabelsY !== false;

  ctx.save();
  ctx.strokeStyle = '#ececec';
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.rect(ax.padL, ax.padT, ax.plotW, ax.plotH);
  ctx.stroke();

  ctx.fillStyle = '#6b6b6b';
  ctx.font = '10px ui-monospace, monospace';

  const xTicks = ax.logX ? logTicks(ax.xMin, ax.xMax) : niceTicks(ax.xMin, ax.xMax, 5);
  ctx.textAlign = 'center'; ctx.textBaseline = 'top';
  for (const t of xTicks) {
    const x = ax.xToPx(t);
    if (showGridX) {
      ctx.strokeStyle = '#f4f4f4';
      ctx.beginPath(); ctx.moveTo(x, ax.padT); ctx.lineTo(x, ax.padT + ax.plotH); ctx.stroke();
    }
    ctx.strokeStyle = '#bbb';
    ctx.beginPath(); ctx.moveTo(x, ax.padT + ax.plotH); ctx.lineTo(x, ax.padT + ax.plotH + 3); ctx.stroke();
    if (showTickLabelsX) ctx.fillText(fmtTick(t), x, ax.padT + ax.plotH + 4);
  }
  const yTicks = niceTicks(ax.yMin, ax.yMax, 4);
  ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
  for (const t of yTicks) {
    const y = ax.yToPx(t);
    if (showGridY) {
      ctx.strokeStyle = '#f4f4f4';
      ctx.beginPath(); ctx.moveTo(ax.padL, y); ctx.lineTo(ax.padL + ax.plotW, y); ctx.stroke();
    }
    ctx.strokeStyle = '#bbb';
    ctx.beginPath(); ctx.moveTo(ax.padL, y); ctx.lineTo(ax.padL - 3, y); ctx.stroke();
    if (showTickLabelsY) ctx.fillText(fmtTick(t), ax.padL - 5, y);
  }
  ctx.restore();
}

// Clip subsequent drawing to the inner plot rectangle. Pair with ctx.restore().
export function clipPlot(ctx, ax) {
  ctx.save();
  ctx.beginPath();
  ctx.rect(ax.padL, ax.padT, ax.plotW, ax.plotH);
  ctx.clip();
}

export function strokePath(ctx, ax, xs, ys, style = {}) {
  if (xs.length === 0) return;
  ctx.save();
  ctx.strokeStyle = style.color || '#2b6cb0';
  ctx.lineWidth = style.width || 1.5;
  if (style.dash) ctx.setLineDash(style.dash);
  ctx.beginPath();
  ctx.moveTo(ax.xToPx(xs[0]), ax.yToPx(ys[0]));
  for (let i = 1; i < xs.length; i++) ctx.lineTo(ax.xToPx(xs[i]), ax.yToPx(ys[i]));
  ctx.stroke();
  ctx.restore();
}

export function dot(ctx, ax, x, y, r = 4, color = '#b34700') {
  ctx.save();
  ctx.fillStyle = color;
  ctx.beginPath();
  ctx.arc(ax.xToPx(x), ax.yToPx(y), r, 0, 2 * Math.PI);
  ctx.fill();
  ctx.restore();
}

function niceTicks(lo, hi, n = 5) {
  const span = hi - lo;
  if (!isFinite(span) || span <= 0) return [lo];
  const step0 = span / n;
  const mag = Math.pow(10, Math.floor(Math.log10(step0)));
  const norm = step0 / mag;
  let step;
  if (norm < 1.5) step = mag;
  else if (norm < 3) step = 2 * mag;
  else if (norm < 7) step = 5 * mag;
  else step = 10 * mag;
  const ticks = [];
  const start = Math.ceil(lo / step) * step;
  for (let v = start; v <= hi + 1e-9 * step; v += step) {
    ticks.push(Number(v.toFixed(12)));
  }
  return ticks;
}

// Log-axis ticks: decades within [lo, hi], plus 2..9 minor ticks per decade
// (minor ones get no label).
function logTicks(lo, hi) {
  const ticks = [];
  const dlo = Math.floor(Math.log10(lo));
  const dhi = Math.ceil(Math.log10(hi));
  for (let d = dlo; d <= dhi; d++) {
    const v = Math.pow(10, d);
    if (v >= lo - 1e-12 && v <= hi + 1e-12) ticks.push(v);
  }
  return ticks;
}

function fmtTick(v) {
  if (v === 0) return '0';
  const a = Math.abs(v);
  if (a >= 1e4 || a < 1e-3) return v.toExponential(1);
  return Number(v.toPrecision(3)).toString();
}
