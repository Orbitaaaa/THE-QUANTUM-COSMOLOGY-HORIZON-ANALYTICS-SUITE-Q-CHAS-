/**
 * DOM helpers for panel authors: element builder, panel chrome, stat tiles, controls, tables,
 * provenance badges, collapsible paper notes and the "bundle not built" state.
 *
 * Strings are always inserted as text (never as HTML), so bundle content such as claim
 * statements or dataset descriptions cannot inject markup.
 */

import { BUILD_COMMAND, type LoadFailure } from './data';

// --- element builder ---------------------------------------------------------------------

export type Child = Node | string | number | null | undefined | false;
type AttrValue = string | number | boolean | null | undefined | EventListener | Record<string, string>;
export type Attrs = Record<string, AttrValue>;

/**
 * `h('div', { class: 'x', onclick: fn, 'aria-label': 'y' }, 'text', child)`.
 * `on*` keys attach listeners, `style` takes a string or object, `dataset` an object.
 */
export function h<K extends keyof HTMLElementTagNameMap>(
  tag: K,
  attrs?: Attrs | null,
  ...children: (Child | Child[])[]
): HTMLElementTagNameMap[K] {
  const el = document.createElement(tag);
  if (attrs) {
    for (const [k, v] of Object.entries(attrs)) {
      if (v === null || v === undefined || v === false) continue;
      if (k.startsWith('on') && typeof v === 'function') el.addEventListener(k.slice(2), v as EventListener);
      else if (k === 'style' && typeof v === 'object') Object.assign(el.style, v);
      else if (k === 'dataset' && typeof v === 'object') Object.assign(el.dataset, v);
      else if (v === true) el.setAttribute(k, '');
      else el.setAttribute(k, String(v));
    }
  }
  append(el, ...children);
  return el;
}

export function append(el: Element, ...children: (Child | Child[])[]): void {
  for (const c of children.flat()) {
    if (c === null || c === undefined || c === false) continue;
    el.append(c instanceof Node ? c : document.createTextNode(String(c)));
  }
}

export const clear = (el: Element): void => el.replaceChildren();

// --- icons -------------------------------------------------------------------------------

const ICONS = {
  database: 'M12 3c4.4 0 8 1.3 8 3v12c0 1.7-3.6 3-8 3s-8-1.3-8-3V6c0-1.7 3.6-3 8-3zm8 3c0 1.7-3.6 3-8 3S4 7.7 4 6m16 6c0 1.7-3.6 3-8 3s-8-1.3-8-3',
  warning: 'M12 3l10 18H2L12 3zm0 6v6m0 3v.5',
  check: 'M20 6L9 17l-5-5',
  info: 'M12 3a9 9 0 110 18 9 9 0 010-18zm0 8v6m0-9v.5',
  critical: 'M12 3a9 9 0 110 18 9 9 0 010-18zm-4 5l8 8m0-8l-8 8',
  sun: 'M12 7a5 5 0 110 10 5 5 0 010-10zm0-5v2m0 16v2M4.2 4.2l1.4 1.4m12.8 12.8l1.4 1.4M2 12h2m16 0h2M4.2 19.8l1.4-1.4M18.4 5.6l1.4-1.4',
  moon: 'M20 14.5A8.5 8.5 0 019.5 4a8.5 8.5 0 1010.5 10.5z',
  menu: 'M4 6h16M4 12h16M4 18h16',
  external: 'M14 4h6v6m0-6L10 14M18 14v5a1 1 0 01-1 1H5a1 1 0 01-1-1V7a1 1 0 011-1h5',
  wrench: 'M14.7 6.3a4 4 0 00-5.4 5.4L3 18l3 3 6.3-6.3a4 4 0 005.4-5.4l-2.5 2.5-2.4-.6-.6-2.4 2.5-2.5z',
} as const;

export type IconName = keyof typeof ICONS;

export function icon(name: IconName, label?: string): SVGSVGElement {
  const ns = 'http://www.w3.org/2000/svg';
  const svg = document.createElementNS(ns, 'svg');
  svg.setAttribute('viewBox', '0 0 24 24');
  svg.setAttribute('class', 'icon');
  svg.setAttribute('fill', 'none');
  svg.setAttribute('stroke', 'currentColor');
  svg.setAttribute('stroke-width', '1.8');
  svg.setAttribute('stroke-linecap', 'round');
  svg.setAttribute('stroke-linejoin', 'round');
  if (label) {
    svg.setAttribute('role', 'img');
    svg.setAttribute('aria-label', label);
  } else svg.setAttribute('aria-hidden', 'true');
  const path = document.createElementNS(ns, 'path');
  path.setAttribute('d', ICONS[name]);
  svg.append(path);
  return svg;
}

// --- formatting --------------------------------------------------------------------------

const isNum = (x: unknown): x is number => typeof x === 'number' && Number.isFinite(x);

/** Significant-figure formatting with thousands separators; '—' for null/NaN. */
export function fmtNum(x: number | null | undefined, sig = 4): string {
  if (!isNum(x)) return '—';
  if (x === 0) return '0';
  const a = Math.abs(x);
  if (a >= 1e6 || a < 1e-3) return fmtSci(x, sig);
  const digits = Math.max(0, sig - 1 - Math.floor(Math.log10(a)));
  return x.toLocaleString('en-US', { maximumFractionDigits: digits, minimumFractionDigits: 0 });
}

export const fmtFixed = (x: number | null | undefined, d = 2): string => (isNum(x) ? x.toFixed(d) : '—');

export function fmtSci(x: number | null | undefined, sig = 3): string {
  if (!isNum(x)) return '—';
  const [m, e] = x.toExponential(sig - 1).split('e');
  const exp = Number(e);
  const sup = String(exp).replace(/[-0-9]/g, (c) => '⁻⁰¹²³⁴⁵⁶⁷⁸⁹'['-0123456789'.indexOf(c)]!);
  return `${m}×10${sup}`;
}

export const fmtInt = (x: number | null | undefined): string => (isNum(x) ? Math.round(x).toLocaleString('en-US') : '—');

export function fmtBytes(n: number | null | undefined): string {
  if (!isNum(n)) return '—';
  const units = ['B', 'KB', 'MB', 'GB', 'TB'];
  let v = n;
  let i = 0;
  while (v >= 1024 && i < units.length - 1) {
    v /= 1024;
    i++;
  }
  return `${i === 0 ? v : v.toFixed(v < 10 ? 2 : 1)} ${units[i]}`;
}

/** ISO timestamp -> "2026-10-09 12:00 UTC" (input returned unchanged if unparseable). */
export function fmtUtc(iso: string | null | undefined, seconds = false): string {
  if (!iso) return '—';
  const d = new Date(iso);
  if (Number.isNaN(d.getTime())) return iso;
  return d.toISOString().slice(0, seconds ? 19 : 16).replace('T', ' ') + ' UTC';
}

// --- provenance ----------------------------------------------------------------------------

export type ProvenanceKind = 'proved' | 'computed' | 'measured' | 'assumed' | 'open' | 'quoted' | 'live';

const KIND_LABEL: Record<ProvenanceKind, string> = {
  proved: 'Proved',
  computed: 'Computed',
  measured: 'Measured',
  assumed: 'Assumed',
  open: 'Open',
  quoted: 'Quoted',
  live: 'Live',
};

export const provenanceLabel = (kind: ProvenanceKind): string => KIND_LABEL[kind];

/** Provenance badge, e.g. `badge('quoted', '[TOP] Sec. 9')` renders QUOTED [TOP] Sec. 9. */
export function badge(kind: ProvenanceKind, source?: string, title?: string): HTMLSpanElement {
  return h(
    'span',
    { class: `badge badge-${kind}`, title: title ?? (source ? `${KIND_LABEL[kind]}: ${source}` : KIND_LABEL[kind]) },
    KIND_LABEL[kind],
    source ? h('span', { class: 'badge-src' }, source) : null,
  );
}

/**
 * Map a status string from a paper's status table ("Theorem, unconditional", "Numerical, ...",
 * "Pre-registered data analysis", "Conditional: ...", "Must be measured; ...") to a provenance
 * kind. The leading clause decides; qualifiers after ';' or ',' refine only when it is silent.
 * A statement "derived by the gravitational path integral" ([BH] Sec. 16, H* for an evaporating
 * hole) is not a theorem of the paper, so it counts as assumed, never as proved.
 */
export function classifyStatus(status: string): ProvenanceKind {
  const s = status.toLowerCase().trim();
  const lead = s.split(/[;,:(]/)[0] ?? s;
  if (/^(must be|unknown|open|undetermined)/.test(lead)) return 'open';
  if (/^(conditional|hypothes|assum|not used|conjectur|estimate|derived|postulat)/.test(lead)) return 'assumed';
  if (/theorem|proved|proof|lemma|corollary|proposition/.test(lead)) return 'proved';
  if (/data|measur|observ/.test(lead)) return 'measured';
  if (/comput|numer|simulat|calculat/.test(lead)) return 'computed';
  if (/hypothes|assum|conditional|estimate/.test(s)) return 'assumed';
  if (/theorem/.test(s)) return 'proved';
  return 'quoted';
}

export type StatusLevel = 'good' | 'warning' | 'serious' | 'critical';

/** Status with icon + label (status colour never carries meaning alone). */
export function statusPill(level: StatusLevel, label: string): HTMLSpanElement {
  const name: IconName = level === 'good' ? 'check' : level === 'critical' ? 'critical' : 'warning';
  return h('span', { class: `status-pill ${level}` }, icon(name), label);
}

// --- panel chrome ------------------------------------------------------------------------------

export type Span = 3 | 4 | 5 | 6 | 7 | 8 | 9 | 12;

export interface PanelOptions {
  title: string;
  subtitle?: string;
  /** provenance badges in the header, e.g. `badge('quoted', '[TOP] Sec. 9')` */
  badges?: Node[];
  /** extra header controls (buttons, selects) */
  actions?: Node[];
  span?: Span;
  /** stretch to the row height */
  stretch?: boolean;
  /** body without padding (tables, sky views) */
  flush?: boolean;
  foot?: Child;
  id?: string;
}

export interface PanelHandle {
  root: HTMLElement;
  head: HTMLElement;
  body: HTMLElement;
  foot: HTMLElement | null;
}

/** A titled panel; add content to `.body`. */
export function panel(opts: PanelOptions): PanelHandle {
  const body = h('div', { class: `panel-body${opts.flush ? ' flush' : ''}` });
  const head = h(
    'div',
    { class: 'panel-head' },
    h(
      'div',
      { class: 'panel-titles' },
      h('h2', { class: 'panel-title' }, opts.title),
      opts.subtitle ? h('div', { class: 'panel-sub' }, opts.subtitle) : null,
    ),
    opts.badges?.length || opts.actions?.length
      ? h('div', { class: 'panel-actions' }, ...(opts.badges ?? []), ...(opts.actions ?? []))
      : null,
  );
  const foot = opts.foot ? h('div', { class: 'panel-foot' }, opts.foot) : null;
  const root = h(
    'section',
    { class: `panel span-${opts.span ?? 12}${opts.stretch ? ' stretch' : ''}`, id: opts.id },
    head,
    body,
    foot,
  );
  return { root, head, body, foot };
}

/** A 12-column grid row holding panels. */
export const grid = (...children: (Child | Child[])[]): HTMLDivElement => h('div', { class: 'grid' }, ...children);

// --- stat tiles ----------------------------------------------------------------------------------

export interface TileOptions {
  label: string;
  value: string | number | null | undefined;
  unit?: string;
  sub?: Child | Child[];
  badge?: Node;
  /** monospace value (identifiers, timestamps) */
  mono?: boolean;
}

export function statTile(o: TileOptions): HTMLDivElement {
  const absent = o.value === null || o.value === undefined || o.value === '—';
  const value = typeof o.value === 'number' ? fmtNum(o.value) : (o.value ?? '—');
  return h(
    'div',
    { class: 'tile' },
    h('div', { class: 'tile-label' }, o.label),
    h(
      'div',
      { class: `tile-value${o.mono ? ' mono' : ''}${absent ? ' absent' : ''}` },
      value,
      o.unit && !absent ? h('span', { class: 'unit' }, o.unit) : null,
    ),
    o.sub || o.badge ? h('div', { class: 'tile-sub' }, o.badge ?? null, ...[o.sub ?? null].flat()) : null,
  );
}

export const tiles = (...items: HTMLElement[]): HTMLDivElement => h('div', { class: 'tiles' }, ...items);

// --- controls --------------------------------------------------------------------------------------

export interface SliderOptions {
  label: string;
  min: number;
  max: number;
  step: number;
  value: number;
  format?: (v: number) => string;
  onInput: (v: number) => void;
}

export interface SliderHandle {
  root: HTMLElement;
  input: HTMLInputElement;
  get value(): number;
  set(v: number, notify?: boolean): void;
}

/** Range slider with a live value readout. */
export function slider(o: SliderOptions): SliderHandle {
  const fmt = o.format ?? ((v: number) => fmtNum(v));
  const out = h('output', { class: 'control-value' }, fmt(o.value));
  const input = h('input', { type: 'range', min: o.min, max: o.max, step: o.step, value: o.value, 'aria-label': o.label });
  const paint = (): void => {
    const v = Number(input.value);
    out.textContent = fmt(v);
    input.style.setProperty('--fill', `${((v - o.min) / (o.max - o.min)) * 100}%`);
  };
  input.addEventListener('input', () => {
    paint();
    o.onInput(Number(input.value));
  });
  paint();
  const root = h('label', { class: 'control' }, h('div', { class: 'control-head' }, h('span', { class: 'control-label' }, o.label), out), input);
  return {
    root,
    input,
    get value() {
      return Number(input.value);
    },
    set(v: number, notify = false) {
      input.value = String(v);
      paint();
      if (notify) o.onInput(Number(input.value));
    },
  };
}

export interface SelectOption {
  value: string;
  label: string;
}

export function select(o: {
  label?: string;
  options: (SelectOption | string)[];
  value?: string;
  onChange: (v: string) => void;
}): { root: HTMLElement; select: HTMLSelectElement } {
  const sel = h('select', { 'aria-label': o.label ?? 'select' });
  for (const opt of o.options) {
    const { value, label } = typeof opt === 'string' ? { value: opt, label: opt } : opt;
    sel.append(h('option', { value, selected: value === o.value }, label));
  }
  sel.addEventListener('change', () => o.onChange(sel.value));
  const root = o.label
    ? h('label', { class: 'control' }, h('span', { class: 'control-label' }, o.label), sel)
    : h('span', null, sel);
  return { root, select: sel };
}

export function toggle(o: { label: string; checked?: boolean; onChange: (v: boolean) => void }): {
  root: HTMLElement;
  input: HTMLInputElement;
} {
  const input = h('input', { type: 'checkbox', role: 'switch', checked: !!o.checked });
  input.addEventListener('change', () => o.onChange(input.checked));
  return { root: h('label', { class: 'switch' }, input, o.label), input };
}

/** Segmented button group (single choice). */
export function segmented(o: {
  options: SelectOption[];
  value: string;
  label: string;
  onChange: (v: string) => void;
}): { root: HTMLElement; set(v: string): void } {
  const buttons = o.options.map((opt) =>
    h('button', { type: 'button', class: 'btn', 'aria-pressed': String(opt.value === o.value), dataset: { value: opt.value } }, opt.label),
  );
  const set = (v: string): void => buttons.forEach((b) => b.setAttribute('aria-pressed', String(b.dataset.value === v)));
  buttons.forEach((b) =>
    b.addEventListener('click', () => {
      set(b.dataset.value!);
      o.onChange(b.dataset.value!);
    }),
  );
  return { root: h('div', { class: 'btn-group', role: 'group', 'aria-label': o.label }, ...buttons), set };
}

// --- tables -----------------------------------------------------------------------------------------

export interface Column<T> {
  key: string;
  label: string;
  /** 'num' right-aligns with tabular figures */
  align?: 'num' | 'text';
  format?: (row: T) => Child | Child[];
  className?: string;
  title?: string;
}

/** Render rows into a table; values come from `format(row)` or `row[key]` as text. */
export function table<T extends object>(columns: Column<T>[], rows: T[], empty = 'No rows'): HTMLDivElement {
  const thead = h(
    'thead',
    null,
    h('tr', null, ...columns.map((c) => h('th', { class: c.align === 'num' ? 'num' : null, title: c.title, scope: 'col' }, c.label))),
  );
  const tbody = h('tbody');
  if (!rows.length) tbody.append(h('tr', null, h('td', { class: 'empty', colspan: columns.length }, empty)));
  for (const row of rows) {
    tbody.append(
      h(
        'tr',
        null,
        ...columns.map((c) => {
          const raw = c.format ? c.format(row) : ((row as Record<string, unknown>)[c.key] as Child);
          const cls = [c.align === 'num' ? 'num' : '', c.className ?? ''].join(' ').trim() || null;
          return h('td', { class: cls }, ...[raw].flat());
        }),
      ),
    );
  }
  return h('div', { class: 'table-wrap' }, h('table', { class: 'data' }, thead, tbody));
}

// --- notes and states ----------------------------------------------------------------------------------

/** Collapsible note, by default titled "What the paper proves". */
export function note(content: Child | Child[], title = 'What the paper proves', open = false): HTMLDetailsElement {
  const items = Array.isArray(content) ? content : [content];
  const body =
    items.length > 1 && items.every((x) => typeof x === 'string')
      ? h('ul', null, ...items.map((x) => h('li', null, x)))
      : h('div', null, ...items);
  return h('details', { class: 'note', open }, h('summary', null, title), h('div', { class: 'note-body' }, body));
}

/** The explicit "Data bundle not built yet" state for one or more missing bundles. */
export function missingBundle(names: string | string[], detail?: string): HTMLDivElement {
  const list = [names].flat();
  return h(
    'div',
    { class: 'empty-state', role: 'status' },
    icon('database'),
    h(
      'div',
      null,
      h('div', { class: 'empty-title' }, `Data bundle not built yet: run ${BUILD_COMMAND}`),
      h(
        'div',
        { class: 'empty-text' },
        list.length === 1 ? 'Missing ' : 'Missing bundles ',
        ...list.flatMap((n, i) => [i ? ', ' : '', h('code', null, `data/${n}`)]),
        '. ',
        detail ?? 'This view shows only real pipeline output; it stays empty until the Python build has produced it.',
      ),
      h('div', { class: 'cmd' }, BUILD_COMMAND),
    ),
  );
}

/** State for any failed load: missing bundles get `missingBundle`, others an error message. */
export function loadFailure(f: LoadFailure): HTMLDivElement {
  if (f.reason === 'missing') return missingBundle(f.name);
  return h(
    'div',
    { class: 'empty-state', role: 'alert' },
    icon('warning'),
    h(
      'div',
      null,
      h('div', { class: 'empty-title' }, f.reason === 'invalid' ? 'Data bundle is malformed' : 'Data bundle could not be loaded'),
      h('div', { class: 'empty-text' }, f.message),
    ),
  );
}

export function copyButton(text: string, label = 'Copy'): HTMLButtonElement {
  const btn = h('button', { type: 'button', class: 'copy-btn', title: `Copy ${text}` }, label);
  btn.addEventListener('click', () => {
    void navigator.clipboard?.writeText(text).then(
      () => {
        btn.textContent = 'Copied';
        setTimeout(() => (btn.textContent = label), 1200);
      },
      () => undefined,
    );
  });
  return btn;
}
