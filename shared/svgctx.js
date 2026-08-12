// Canvas2D-compatible recorder that emits SVG instead of pixels.
//
// Why this exists: every plot on setups 1–2 is drawn through the small
// primitive set in canvas.js (lines, rects, arcs, text). None of it is
// inherently raster — no image data, no gradients, no per-pixel fields — so
// the same draw calls can be replayed into vector output for print.
//
// This implements only the subset the plotting code actually uses. It is NOT
// a general Canvas2D emulator; unimplemented members are absent on purpose so
// that a future draw call using them fails loudly rather than silently
// dropping marks from an exported figure.

const DEFAULT_STATE = {
  stroke: '#000000',
  fill: '#000000',
  lineWidth: 1,
  lineDash: [],
  lineCap: 'butt',
  lineJoin: 'miter',
  alpha: 1,
  font: '10px sans-serif',
  textAlign: 'start',
  textBaseline: 'alphabetic',
  clipId: null,
  matrix: null,          // null == identity
};

// 3 decimals is well below one device pixel at any sane figure size, and keeps
// the exported file an order of magnitude smaller than full float precision.
function fmt(v) {
  if (!isFinite(v)) return '0';
  return String(Math.round(v * 1000) / 1000);
}

// SVG (and Illustrator in particular) handles `rgb()` + a separate opacity
// attribute far more reliably than CSS `rgba()`, which the plotting code uses
// throughout for faint reference curves.
function parseColor(c) {
  if (typeof c !== 'string') return { color: '#000000', alpha: 1 };
  const m = c.match(/^rgba?\(([^)]+)\)$/i);
  if (!m) return { color: c, alpha: 1 };
  const p = m[1].split(',').map(s => parseFloat(s.trim()));
  const a = p.length > 3 && isFinite(p[3]) ? p[3] : 1;
  return { color: `rgb(${p[0] | 0},${p[1] | 0},${p[2] | 0})`, alpha: a };
}

function parseFont(font) {
  const m = String(font).match(/(\d+(?:\.\d+)?)px\s+(.*)$/);
  if (!m) return { size: 10, family: 'sans-serif' };
  return { size: parseFloat(m[1]), family: m[2] };
}

function escapeText(s) {
  return String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
}

const ANCHOR = { start: 'start', left: 'start', center: 'middle', middle: 'middle', right: 'end', end: 'end' };

// Canvas textBaseline resolved into an explicit y-shift rather than SVG's
// `dominant-baseline`, which Illustrator ignores. Factors are approximate
// ascent/x-height fractions — good to a fraction of a pixel at 10px type.
function baselineShift(baseline, size) {
  switch (baseline) {
    case 'top':        return 0.80 * size;
    case 'hanging':    return 0.72 * size;
    case 'middle':     return 0.35 * size;
    case 'bottom':     return -0.20 * size;
    default:           return 0;   // alphabetic
  }
}

export function createSVGContext(width, height, opts = {}) {
  const background = opts.background === undefined ? '#ffffff' : opts.background;
  const desc = opts.desc || null;

  const body = [];
  const defs = [];
  let clipSeq = 0;

  let st = { ...DEFAULT_STATE, lineDash: [] };
  const stack = [];

  let path = [];
  let hasCurrentPoint = false;

  // ─── attribute builders ────────────────────────────────────────────────
  function commonAttrs() {
    let s = '';
    if (st.clipId) s += ` clip-path="url(#${st.clipId})"`;
    if (st.matrix) s += ` transform="matrix(${st.matrix.map(fmt).join(' ')})"`;
    return s;
  }

  function strokeAttrs() {
    const { color, alpha } = parseColor(st.stroke);
    let s = ` fill="none" stroke="${color}" stroke-width="${fmt(st.lineWidth)}"`;
    const a = alpha * st.alpha;
    if (a < 1) s += ` stroke-opacity="${fmt(a)}"`;
    if (st.lineDash && st.lineDash.length) s += ` stroke-dasharray="${st.lineDash.map(fmt).join(' ')}"`;
    if (st.lineCap !== 'butt') s += ` stroke-linecap="${st.lineCap}"`;
    if (st.lineJoin !== 'miter') s += ` stroke-linejoin="${st.lineJoin}"`;
    return s + commonAttrs();
  }

  function fillAttrs() {
    const { color, alpha } = parseColor(st.fill);
    let s = ` fill="${color}" stroke="none"`;
    const a = alpha * st.alpha;
    if (a < 1) s += ` fill-opacity="${fmt(a)}"`;
    return s + commonAttrs();
  }

  // ─── the context object ────────────────────────────────────────────────
  const ctx = {
    // -- state properties (plain accessors onto `st`) --
    get strokeStyle() { return st.stroke; },  set strokeStyle(v) { st.stroke = v; },
    get fillStyle()   { return st.fill; },    set fillStyle(v)   { st.fill = v; },
    get lineWidth()   { return st.lineWidth; }, set lineWidth(v) { st.lineWidth = v; },
    get lineCap()     { return st.lineCap; },  set lineCap(v)    { st.lineCap = v; },
    get lineJoin()    { return st.lineJoin; }, set lineJoin(v)   { st.lineJoin = v; },
    get globalAlpha() { return st.alpha; },    set globalAlpha(v) { st.alpha = v; },
    get font()        { return st.font; },     set font(v)       { st.font = v; },
    get textAlign()   { return st.textAlign; }, set textAlign(v) { st.textAlign = v; },
    get textBaseline() { return st.textBaseline; }, set textBaseline(v) { st.textBaseline = v; },

    save() { stack.push({ ...st, lineDash: st.lineDash.slice() }); },
    restore() { if (stack.length) st = stack.pop(); },

    setLineDash(d) { st.lineDash = (d || []).slice(); },
    getLineDash() { return st.lineDash.slice(); },

    // -- transforms --
    setTransform(a, b, c, d, e, f) {
      st.matrix = (a === 1 && b === 0 && c === 0 && d === 1 && e === 0 && f === 0)
        ? null : [a, b, c, d, e, f];
    },
    translate(x, y) {
      const m = st.matrix || [1, 0, 0, 1, 0, 0];
      st.matrix = [m[0], m[1], m[2], m[3], m[4] + m[0] * x + m[2] * y, m[5] + m[1] * x + m[3] * y];
    },
    scale(sx, sy) {
      const m = st.matrix || [1, 0, 0, 1, 0, 0];
      st.matrix = [m[0] * sx, m[1] * sx, m[2] * sy, m[3] * sy, m[4], m[5]];
    },

    // -- path construction --
    beginPath() { path = []; hasCurrentPoint = false; },
    closePath() { if (path.length) path.push('Z'); },
    moveTo(x, y) { path.push(`M${fmt(x)} ${fmt(y)}`); hasCurrentPoint = true; },
    lineTo(x, y) {
      if (!hasCurrentPoint) return ctx.moveTo(x, y);
      path.push(`L${fmt(x)} ${fmt(y)}`);
    },
    rect(x, y, w, h) {
      path.push(`M${fmt(x)} ${fmt(y)}h${fmt(w)}v${fmt(h)}h${fmt(-w)}Z`);
      hasCurrentPoint = true;
    },
    arc(cx, cy, r, a0, a1, ccw = false) {
      const TAU = 2 * Math.PI;
      let delta = a1 - a0;
      if (!ccw) {
        if (delta < 0) delta = delta % TAU + TAU;
        if (delta > TAU) delta = TAU;
      } else {
        if (delta > 0) delta = delta % TAU - TAU;
        if (delta < -TAU) delta = -TAU;
      }
      const x0 = cx + r * Math.cos(a0), y0 = cy + r * Math.sin(a0);
      // Canvas draws an implicit line to the arc's start if a path is open.
      path.push(`${hasCurrentPoint ? 'L' : 'M'}${fmt(x0)} ${fmt(y0)}`);
      hasCurrentPoint = true;
      if (Math.abs(delta) >= TAU - 1e-9) {
        // A single elliptical-arc command cannot express a full circle
        // (start == end is degenerate); split into two half turns.
        const sweep = delta > 0 ? 1 : 0;
        const xm = cx - r * Math.cos(a0), ym = cy - r * Math.sin(a0);
        path.push(`A${fmt(r)} ${fmt(r)} 0 0 ${sweep} ${fmt(xm)} ${fmt(ym)}`);
        path.push(`A${fmt(r)} ${fmt(r)} 0 0 ${sweep} ${fmt(x0)} ${fmt(y0)}`);
        return;
      }
      const x1 = cx + r * Math.cos(a0 + delta), y1 = cy + r * Math.sin(a0 + delta);
      const large = Math.abs(delta) > Math.PI ? 1 : 0;
      const sweep = delta > 0 ? 1 : 0;
      path.push(`A${fmt(r)} ${fmt(r)} 0 ${large} ${sweep} ${fmt(x1)} ${fmt(y1)}`);
    },

    // -- painting --
    stroke() {
      if (!path.length) return;
      body.push(`<path d="${path.join('')}"${strokeAttrs()}/>`);
    },
    fill() {
      if (!path.length) return;
      body.push(`<path d="${path.join('')}"${fillAttrs()}/>`);
    },
    fillRect(x, y, w, h) {
      if (w === 0 || h === 0) return;
      const nx = w < 0 ? x + w : x, ny = h < 0 ? y + h : y;
      body.push(`<rect x="${fmt(nx)}" y="${fmt(ny)}" width="${fmt(Math.abs(w))}" height="${fmt(Math.abs(h))}"${fillAttrs()}/>`);
    },
    strokeRect(x, y, w, h) {
      body.push(`<rect x="${fmt(x)}" y="${fmt(y)}" width="${fmt(w)}" height="${fmt(h)}"${strokeAttrs()}/>`);
    },
    // No-op: the SVG starts from a blank (optionally white-filled) canvas, so
    // there is nothing to erase. Honouring it would mean masking earlier marks.
    clearRect() {},

    fillText(text, x, y) {
      const s = String(text);
      if (!s.length) return;
      const { size, family } = parseFont(st.font);
      const anchor = ANCHOR[st.textAlign] || 'start';
      const yy = y + baselineShift(st.textBaseline, size);
      body.push(
        `<text x="${fmt(x)}" y="${fmt(yy)}" font-family="${escapeText(family)}" ` +
        `font-size="${fmt(size)}" text-anchor="${anchor}"${fillAttrs()}>${escapeText(s)}</text>`
      );
    },
    measureText(text) {
      const { size } = parseFont(st.font);
      return { width: 0.6 * size * String(text).length };   // monospace approximation
    },

    clip() {
      if (!path.length) return;
      const id = `clip${++clipSeq}`;
      defs.push(`<clipPath id="${id}"><path d="${path.join('')}"/></clipPath>`);
      st.clipId = id;
    },

    // ─── output ──────────────────────────────────────────────────────────
    toSVG() {
      const parts = [
        `<svg xmlns="http://www.w3.org/2000/svg" width="${fmt(width)}" height="${fmt(height)}" ` +
        `viewBox="0 0 ${fmt(width)} ${fmt(height)}">`,
      ];
      if (desc) parts.push(`<desc>${escapeText(desc)}</desc>`);
      if (defs.length) parts.push(`<defs>${defs.join('')}</defs>`);
      if (background) parts.push(`<rect width="${fmt(width)}" height="${fmt(height)}" fill="${background}"/>`);
      parts.push(...body);
      parts.push('</svg>');
      return parts.join('\n');
    },
  };

  return ctx;
}
