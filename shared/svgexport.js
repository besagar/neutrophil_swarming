// Vector (SVG) export for the Canvas2D plots.
//
// The draw functions are reused verbatim: withContext() swaps the context that
// autoFit() hands out, so a redraw thunk that normally paints pixels instead
// records into an SVG. Nothing about the on-screen rendering path changes.
//
// Axis labels and titles are KaTeX DOM overlays (see CLAUDE.md), not canvas
// marks, so they are deliberately absent from the exported line art. Their
// LaTeX source is embedded in the SVG's <desc> for typesetting downstream.

import { withContext } from './canvas.js';
import { createSVGContext } from './svgctx.js';
import { plotLabels } from './dom.js';

function descFor(canvasId) {
  const l = plotLabels.get(canvasId);
  if (!l) return null;
  const bits = [];
  if (l.titleTex)  bits.push(`title: $${l.titleTex}$`);
  if (l.xLabelTex) bits.push(`x: $${l.xLabelTex}$`);
  if (l.yLabelTex) bits.push(`y: $${l.yLabelTex}$`);
  return bits.length ? `Axis labels (KaTeX source, not rendered in this file) — ${bits.join('  |  ')}` : null;
}

// Replays `redraw` into an SVG recorder and returns the SVG source.
// `redraw` must be a zero-arg thunk that performs the plot's normal drawing.
export function renderPlotSvg(canvasId, redraw, opts = {}) {
  const cv = document.getElementById(canvasId);
  if (!cv) return null;
  const w = opts.width  || cv.clientWidth;
  const h = opts.height || cv.clientHeight;
  if (!w || !h) return null;
  const ctx = createSVGContext(w, h, {
    background: opts.background === undefined ? '#ffffff' : opts.background,
    desc: descFor(canvasId),
  });
  withContext(ctx, redraw);
  return ctx.toSVG();
}

function download(filename, text) {
  const blob = new Blob([text], { type: 'image/svg+xml;charset=utf-8' });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  a.remove();
  // Revoke on the next tick: revoking synchronously can race the download in
  // WebKit.
  setTimeout(() => URL.revokeObjectURL(url), 0);
}

export function downloadPlotSvg(canvasId, redraw, opts = {}) {
  const svg = renderPlotSvg(canvasId, redraw, opts);
  if (!svg) return;
  download((opts.filename || canvasId) + '.svg', svg);
}

// Adds a small ⤓SVG button to the plot's .plot-wrap. Call after decoratePlot()
// so the wrapper exists.
export function attachSvgExport(canvasId, redraw, opts = {}) {
  const cv = document.getElementById(canvasId);
  if (!cv) return;
  const wrap = cv.closest('.plot-wrap') || cv.parentNode;
  if (!wrap) return;
  if (wrap.querySelector('.svg-export')) return;
  if (getComputedStyle(wrap).position === 'static') wrap.style.position = 'relative';

  const btn = document.createElement('button');
  btn.className = 'svg-export';
  btn.type = 'button';
  btn.textContent = '⤓ SVG';
  btn.title = 'Download this plot as vector SVG';
  btn.addEventListener('click', () => {
    downloadPlotSvg(canvasId, redraw, opts);
  });
  wrap.appendChild(btn);
}

// Convenience: wire a whole page at once from a {canvasId: redrawThunk} map.
export function attachSvgExports(map, opts = {}) {
  for (const [canvasId, redraw] of Object.entries(map)) {
    attachSvgExport(canvasId, redraw, {
      ...opts,
      filename: opts.prefix ? `${opts.prefix}-${canvasId.replace(/^cv-/, '')}` : canvasId,
    });
  }
}
