// Minimal widget helpers. Each `knob` is configured by:
//   { id, symbol, value, min, max, step, log?, units?, linkedLabel?, bind? }
// Returns { el, get, set, onChange(cb) }.
//
// Slider state and per-slider configuration (min, max, default value) are
// persisted to localStorage, scoped per page (setup1, setup2, ...). Each
// slider exposes a gear (⚙) button that pops up a small editor for those
// fields, plus reset buttons.
//
// `bind: [obj, key]` makes the slider the single source of truth for
// `obj[key]`: the built-in default is read from obj[key] (cfg.value is
// ignored when bind is present), and obj[key] is rewritten on every value
// change — including on construction. This makes it impossible for the
// displayed value and the simulated value to diverge.

const STORAGE_PREFIX = 'gl-motility';
function pageKey() {
  const m = (typeof window !== 'undefined' && window.location.pathname || '').match(/setup\d+/);
  return m ? m[0] : 'index';
}
function storageKey(id) { return `${STORAGE_PREFIX}:${pageKey()}:${id}`; }
function loadSliderState(id) {
  if (!id || typeof localStorage === 'undefined') return null;
  try { const r = localStorage.getItem(storageKey(id)); return r ? JSON.parse(r) : null; }
  catch { return null; }
}
function saveSliderState(id, state) {
  if (!id || typeof localStorage === 'undefined') return;
  try { localStorage.setItem(storageKey(id), JSON.stringify(state)); } catch {}
}
function clearSliderState(id) {
  if (!id || typeof localStorage === 'undefined') return;
  try { localStorage.removeItem(storageKey(id)); } catch {}
}

export function el(tag, attrs = {}, children = []) {
  const n = document.createElement(tag);
  for (const k in attrs) {
    if (k === 'class') n.className = attrs[k];
    else if (k === 'style') Object.assign(n.style, attrs[k]);
    else if (k.startsWith('on') && typeof attrs[k] === 'function') n.addEventListener(k.slice(2), attrs[k]);
    else n.setAttribute(k, attrs[k]);
  }
  for (const c of [].concat(children)) {
    if (c == null) continue;
    n.appendChild(typeof c === 'string' ? document.createTextNode(c) : c);
  }
  return n;
}

// Render `tex` into `node` with KaTeX, if it looks like LaTeX (contains a
// backslash, underscore, caret, or brace). Plain ASCII / unicode symbols are
// left as text. KaTeX may load after us, so retry on `katex-ready` and on a
// short polling schedule (older setup pages never dispatch that event).
function renderSymbol(node, tex) {
  if (!/[\\_^{}]/.test(tex)) return;
  const renderSym = () => {
    if (!window.katex) return false;
    try {
      window.katex.render(tex, node, { throwOnError: false });
      return true;
    } catch (_) { return false; }
  };
  if (renderSym()) return;
  window.addEventListener('katex-ready', renderSym, { once: true });
  let tries = 0;
  const tick = () => {
    if (renderSym() || ++tries > 20) return;
    setTimeout(tick, 150);
  };
  setTimeout(tick, 100);
}

// Linear or log slider. For log, the slider works in log-space internally
// but the public value is always the linear value.
//
// Persistence: if cfg.id is set, the current value and any user-customized
// {min, max, default} are saved to localStorage. On the next page load they
// override the cfg defaults. The gear button (⚙) opens a small popover
// where the user edits min/max/default explicitly.
export function makeSlider(cfg) {
  const log = !!cfg.log;
  const fmt = cfg.fmt || (v => Number(v).toPrecision(3));
  const id = cfg.id;
  const bind = cfg.bind || null;           // [obj, key]
  const transform = cfg.transform || null; // value -> canonicalized value (e.g. Math.round)

  // ─── load persisted state, if any ──────────────────────────────────────
  const saved = loadSliderState(id) || {};
  // Effective min/max/default (saved overrides cfg). The original cfg.value
  // is treated as the *built-in* default; the user can override it via gear.
  // When bind is set, the built-in default comes from obj[key] — so the
  // bound object's literal is the single source of truth and cfg.value is
  // ignored.
  const builtinDefault = bind ? bind[0][bind[1]] : cfg.value;
  let effMin = (saved.min  != null && isFinite(saved.min))  ? saved.min  : cfg.min;
  let effMax = (saved.max  != null && isFinite(saved.max))  ? saved.max  : cfg.max;
  let userDefault = (saved.default != null && isFinite(saved.default)) ? saved.default : builtinDefault;
  // On page load, always start at the (possibly user-customized) default
  // rather than the last value the slider held. Min/max/default config from
  // the gear popover still persists; only the live value is reset.
  let value = userDefault;
  // Sanity-clamp loaded value into current effMin..effMax
  value = Math.min(effMax, Math.max(effMin, value));

  // ─── DOM ───────────────────────────────────────────────────────────────
  const labelRow = el('div', { class: 'label-row' });
  const sym = el('span', { class: 'sym' }, cfg.symbol || cfg.id);
  renderSymbol(sym, cfg.symbol || '');
  const valSpan = el('span', { class: 'val' });
  const units = cfg.units ? el('span', { class: 'units' }, cfg.units) : null;
  const gearBtn = el('button', { class: 'knob-gear', type: 'button', title: 'configure slider' }, '⚙');
  const rightWrap = el('span', { class: 'right-wrap' });
  rightWrap.appendChild(valSpan);
  if (units) rightWrap.appendChild(units);
  rightWrap.appendChild(gearBtn);
  labelRow.appendChild(sym);
  labelRow.appendChild(rightWrap);

  const range = el('input', {
    type: 'range',
    step: cfg.step != null ? cfg.step : 'any',
  });
  // (min/max set after we know log/linear, see applyRange)

  const linked = cfg.linkedLabel ? el('div', { class: 'linked-readout' }, '') : null;

  // ─── popover (gear) ────────────────────────────────────────────────────
  const cfgPanel = el('div', { class: 'knob-config' });
  cfgPanel.style.display = 'none';
  const inMin = el('input', { type: 'number', class: 'cfg-min',  step: 'any' });
  const inMax = el('input', { type: 'number', class: 'cfg-max',  step: 'any' });
  const inDef = el('input', { type: 'number', class: 'cfg-default', step: 'any' });
  const lblMin = el('label', {}, 'min');
  const lblMax = el('label', {}, 'max');
  const lblDef = el('label', {}, 'default');
  const btnResetVal = el('button', { type: 'button', class: 'cfg-action' }, 'reset value');
  const btnResetCfg = el('button', { type: 'button', class: 'cfg-action ghost' }, 'reset config');
  const actions = el('div', { class: 'actions' });
  actions.appendChild(btnResetVal);
  actions.appendChild(btnResetCfg);
  cfgPanel.appendChild(lblMin); cfgPanel.appendChild(inMin);
  cfgPanel.appendChild(lblMax); cfgPanel.appendChild(inMax);
  cfgPanel.appendChild(lblDef); cfgPanel.appendChild(inDef);
  cfgPanel.appendChild(actions);

  function syncCfgInputs() {
    inMin.value = String(effMin);
    inMax.value = String(effMax);
    inDef.value = String(userDefault);
  }

  // ─── range/value mechanics ─────────────────────────────────────────────
  function applyRange() {
    const lo = log ? Math.log(effMin) : effMin;
    const hi = log ? Math.log(effMax) : effMax;
    range.min = String(lo);
    range.max = String(hi);
    range.value = String(log ? Math.log(value) : value);
  }
  applyRange();

  const subs = [];

  function refresh() {
    // transform canonicalizes the raw slider value (e.g. Math.round for integer
    // knobs) so display, bound params, and downstream consumers see exactly the
    // same number — no rounding-via-onChange surprises on initial construction.
    if (transform) value = transform(value);
    // Write-through to bound params object on every value change (including
    // construction). This is what guarantees displayed-value === simulated-value.
    if (bind) bind[0][bind[1]] = value;
    valSpan.textContent = fmt(value);
    if (linked && cfg.linkedLabel) linked.textContent = cfg.linkedLabel(value);
  }
  refresh();

  function persist() {
    saveSliderState(id, {
      value,
      min: effMin === cfg.min ? null : effMin,
      max: effMax === cfg.max ? null : effMax,
      default: userDefault === builtinDefault ? null : userDefault,
    });
  }

  const relSubs = [];
  range.addEventListener('input', () => {
    const raw = parseFloat(range.value);
    value = log ? Math.exp(raw) : raw;
    refresh();
    persist();
    for (const s of subs) s(value);
  });
  // 'change' fires once when the user releases the slider (mouseup / touchend).
  range.addEventListener('change', () => { for (const s of relSubs) s(value); });

  // ─── gear popover open/close ──────────────────────────────────────────
  gearBtn.addEventListener('click', e => {
    e.stopPropagation();
    const isHidden = cfgPanel.style.display === 'none';
    if (isHidden) syncCfgInputs();
    cfgPanel.style.display = isHidden ? '' : 'none';
  });

  function commitCfgInputs() {
    const newMin = parseFloat(inMin.value);
    const newMax = parseFloat(inMax.value);
    const newDef = parseFloat(inDef.value);
    if (!(isFinite(newMin) && isFinite(newMax) && newMin < newMax)) {
      // Reject: revert inputs to current
      syncCfgInputs();
      return;
    }
    effMin = newMin; effMax = newMax;
    if (isFinite(newDef)) userDefault = Math.min(effMax, Math.max(effMin, newDef));
    applyRange();
    // Re-clamp current value into the new range
    const clamped = Math.min(effMax, Math.max(effMin, value));
    if (clamped !== value) {
      value = clamped;
      range.value = String(log ? Math.log(value) : value);
    }
    refresh();
    persist();
    for (const s of subs) s(value);
    syncCfgInputs();
  }
  inMin.addEventListener('change', commitCfgInputs);
  inMax.addEventListener('change', commitCfgInputs);
  inDef.addEventListener('change', commitCfgInputs);

  btnResetVal.addEventListener('click', () => {
    value = Math.min(effMax, Math.max(effMin, userDefault));
    range.value = String(log ? Math.log(value) : value);
    refresh();
    persist();
    for (const s of subs) s(value);
  });
  btnResetCfg.addEventListener('click', () => {
    effMin = cfg.min; effMax = cfg.max; userDefault = builtinDefault;
    value = Math.min(effMax, Math.max(effMin, value));
    applyRange();
    refresh();
    clearSliderState(id);
    syncCfgInputs();
    for (const s of subs) s(value);
  });

  const wrap = el('div', { class: 'knob' }, [labelRow, range]);
  if (linked) wrap.appendChild(linked);
  wrap.appendChild(cfgPanel);

  const api = {
    el: wrap,
    get value() { return value; },
    set(v) {
      // clamp into slider range so internal value stays in sync with DOM
      const clamped = Math.min(effMax, Math.max(effMin, v));
      value = clamped;
      range.value = String(log ? Math.log(clamped) : clamped);
      refresh();
      persist();
      for (const s of subs) s(value);
    },
    onChange(cb) { subs.push(cb); },
    onRelease(cb) { relSubs.push(cb); },
    setLinkedText(text) { if (linked) linked.textContent = text; },
    setMinMax(newMin, newMax) {
      // Programmatic min/max change (e.g. from another slider's linkage).
      // Counts as a *config* change just like the gear popover.
      effMin = newMin; effMax = newMax;
      cfg.min = newMin; cfg.max = newMax;
      applyRange();
      api.set(value); // re-clamp + persist + notify
    },
  };
  return api;
}

export function makeToggle(cfg) {
  // cfg: { label, options: [{id, label}], value, onChange }
  const wrap = el('div');
  if (cfg.label) wrap.appendChild(el('div', { class: 'label-row' }, [el('span', { class: 'sym' }, cfg.label)]));
  const row = el('div', { class: 'toggle-row' });
  const buttons = {};
  let current = cfg.value;
  function refresh() {
    for (const id in buttons) buttons[id].classList.toggle('on', id === current);
  }
  for (const opt of cfg.options) {
    const b = el('button', { type: 'button' }, opt.label);
    b.addEventListener('click', () => {
      if (current === opt.id) return;
      current = opt.id;
      refresh();
      cfg.onChange && cfg.onChange(current);
    });
    buttons[opt.id] = b;
    row.appendChild(b);
  }
  refresh();
  wrap.appendChild(row);
  return {
    el: wrap,
    get value() { return current; },
    set(v) { current = v; refresh(); },
  };
}

// Typed numeric field. Same `bind: [obj, key]` contract as makeSlider — the
// bound value is the single source of truth and is rewritten on every commit,
// so the displayed number and the simulated number can never diverge.
//
// Unlike a slider there is no drag, so nothing is committed while typing: the
// value lands on Enter or blur. Out-of-range or unparseable input is clamped /
// reverted and the field is redrawn from the committed value.
//
// cfg: { id, symbol, bind:[obj,key], min, max, fmt?, hint?, onCommit? }
export function makeNumberField(cfg) {
  const fmt = cfg.fmt || (v => String(v));
  const bind = cfg.bind || null;
  const lo = cfg.min != null ? cfg.min : -Infinity;
  const hi = cfg.max != null ? cfg.max : Infinity;

  const saved = loadSliderState(cfg.id) || {};
  let value = (saved.value != null && isFinite(saved.value))
    ? Math.min(hi, Math.max(lo, saved.value))
    : (bind ? bind[0][bind[1]] : cfg.value);

  const sym = el('span', { class: 'sym' });
  sym.textContent = cfg.symbol || cfg.id;
  renderSymbol(sym, cfg.symbol || '');
  const input = el('input', { type: 'text', inputmode: 'decimal', class: 'numfield-input' });
  const labelRow = el('div', { class: 'label-row' }, [sym, input]);
  const hint = cfg.hint ? el('div', { class: 'linked-readout' }, '') : null;
  const wrap = el('div', { class: 'knob numfield' }, hint ? [labelRow, hint] : [labelRow]);

  const subs = [];
  function write(v) {
    value = v;
    if (bind) bind[0][bind[1]] = v;
    input.value = fmt(v);
    if (cfg.id) saveSliderState(cfg.id, { value: v });
  }
  function commit() {
    const parsed = parseFloat(input.value);
    const next = isFinite(parsed) ? Math.min(hi, Math.max(lo, parsed)) : value;
    const changed = next !== value;
    write(next);
    if (changed) for (const s of subs) s(value);
  }
  input.addEventListener('keydown', e => { if (e.key === 'Enter') { commit(); input.blur(); } });
  input.addEventListener('blur', commit);
  write(value);

  return {
    el: wrap,
    get value() { return value; },
    set(v) {
      const clamped = Math.min(hi, Math.max(lo, v));
      const changed = clamped !== value;
      write(clamped);
      if (changed) for (const s of subs) s(value);
    },
    onChange(cb) { subs.push(cb); },
    setHintText(text) { if (hint) hint.textContent = text; },
  };
}

export function makeButtonRow(buttons) {
  // buttons: [{label, onClick, ghost?}]
  const row = el('div', { class: 'btnrow' });
  const refs = {};
  for (const b of buttons) {
    const btn = el('button', { type: 'button', class: b.ghost ? 'ghost' : '' }, b.label);
    btn.addEventListener('click', b.onClick);
    refs[b.label] = btn;
    row.appendChild(btn);
  }
  return { el: row, refs };
}

export function makeKpis(items) {
  // items: [{id, label, init?}]
  const wrap = el('div', { class: 'kpis' });
  const refs = {};
  for (const it of items) {
    const v = el('span', { class: 'v' }, it.init || '–');
    const node = el('div', { class: 'kpi' }, [
      el('span', { class: 'k' }, it.label),
      v,
    ]);
    refs[it.id] = v;
    wrap.appendChild(node);
  }
  return {
    el: wrap,
    set(id, text) { if (refs[id]) refs[id].textContent = text; },
  };
}

// Wrap an existing canvas in a div and overlay KaTeX-rendered axis labels.
// `opts`: { titleTex, xLabelTex, yLabelTex } — TeX source strings.
// Call once at startup, after KaTeX is loaded.
// Records the LaTeX passed to decoratePlot, keyed by canvas id. The labels are
// DOM overlays, not canvas marks, so they cannot appear in a vector export of
// the canvas itself — svgexport.js writes them into the SVG's <desc> so the
// source strings travel with the figure into the poster/figure tool.
export const plotLabels = new Map();

export function decoratePlot(canvasId, opts = {}) {
  const canvas = document.getElementById(canvasId);
  if (!canvas) return;
  plotLabels.set(canvasId, opts);
  if (canvas.parentNode.classList.contains('plot-wrap')) return;
  const wrap = el('div', { class: 'plot-wrap' });
  canvas.parentNode.insertBefore(wrap, canvas);
  wrap.appendChild(canvas);
  const tex = s => (window.katex
    ? window.katex.renderToString(s, { throwOnError: false, displayMode: false })
    : s);
  if (opts.titleTex) {
    const t = el('div', { class: 'ax-title' });
    t.innerHTML = tex(opts.titleTex);
    wrap.appendChild(t);
  }
  if (opts.xLabelTex) {
    const x = el('div', { class: 'ax-x' });
    x.innerHTML = tex(opts.xLabelTex);
    wrap.appendChild(x);
  }
  if (opts.yLabelTex) {
    const y = el('div', { class: 'ax-y' });
    y.innerHTML = tex(opts.yLabelTex);
    wrap.appendChild(y);
  }
}

export function section(title, children) {
  return el('div', { class: 'section' }, [
    el('h4', {}, title),
    ...[].concat(children),
  ]);
}

// Like section() but collapsible. Hidden by default unless open=true.
export function detailsSection(title, children, open = false) {
  const det = el('details', { class: 'section section-details' });
  if (open) det.setAttribute('open', '');
  det.appendChild(el('summary', {}, title));
  for (const c of [].concat(children)) det.appendChild(c);
  return det;
}
