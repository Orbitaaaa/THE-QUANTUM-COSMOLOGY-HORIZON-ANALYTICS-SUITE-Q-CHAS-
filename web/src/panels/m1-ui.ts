/**
 * Shared DOM pieces of the Module 1 panels (manifolds, lattice, bounce, ghost-twin): the
 * provenance-led statement list, the ours-vs-paper validation summary, tile rows, a raw data
 * table with row classes, and in-place line-chart updates for slider-driven live physics.
 * Styles live in m1.css.
 */

import { type ChartHandle, type LineSeries, loadPlotly } from '../lib/charts';
import type { LoadResult, ValidationBundle, ValidationRow } from '../lib/data';
import type { PanelContext } from '../lib/panel';
import { type Child, type ProvenanceKind, badge, fmtNum, h, statusPill, tiles } from '../lib/ui';

/**
 * Fixed palette slots of the lattice shapes across every Module 1 chart and sky view, so a
 * shape keeps its colour everywhere: SC/cubic orange, FCC green, BCC violet. Slot 1 (blue) is
 * the main computed curve.
 */
export const SHAPE_SLOT: Record<string, number> = { SC: 2, cubic: 2, FCC: 3, BCC: 7 };

export interface ProofItem {
  kind: ProvenanceKind;
  src: string;
  text: Child | Child[];
}

/** A list of statements, each led by its provenance badge. */
export function proofList(items: ProofItem[], stacked = false): HTMLUListElement {
  return h(
    'ul',
    { class: `m1-proofs${stacked ? ' stacked' : ''}` },
    ...items.map((it) => h('li', null, badge(it.kind, it.src), h('span', null, ...[it.text].flat()))),
  );
}

/** Register a disposable for unmount and return it. */
export function keep<T extends { destroy(): void }>(ctx: PanelContext, x: T): T {
  ctx.onCleanup(() => x.destroy());
  return x;
}

/** Stat tiles in a fixed number of columns (two below 1100 px). */
export function tilesN(n: 2 | 3 | 4 | 5, ...items: HTMLElement[]): HTMLDivElement {
  const t = tiles(...items);
  t.classList.add(`m1-cols-${n}`);
  return t;
}

export const sectionLabel = (text: string): HTMLDivElement => h('div', { class: 'm1-section-label' }, text);

/** Number in a table cell with the paper's value beneath it (muted). */
export function oursPaper(ours: string, paper: string | null, title?: string): HTMLElement {
  return h(
    'span',
    { class: 'm1-op', title },
    h('span', { class: 'm1-op-ours' }, ours),
    paper !== null ? h('span', { class: 'm1-op-paper' }, `paper ${paper}`) : null,
  );
}

/** Rows of validation.json in the given groups. */
export function validationRows(res: LoadResult<ValidationBundle>, groups: string[], filter?: (r: ValidationRow) => boolean): ValidationRow[] {
  return res.ok ? res.data.rows.filter((r) => groups.includes(r.group) && (filter?.(r) ?? true)) : [];
}

/** "N of M match the printed values" pill with the validation groups named. */
export function validationSummary(rows: ValidationRow[], groups: string[]): HTMLElement | null {
  if (!rows.length) return null;
  const nMatch = rows.filter((r) => r.verdict === 'match').length;
  const nDiff = rows.length - nMatch;
  return h(
    'div',
    { class: 'row m1-validation' },
    statusPill(nDiff ? 'warning' : 'good', `${nMatch} of ${rows.length} match the printed values`),
    nDiff ? h('span', { class: 'small muted' }, `${nDiff} differ, each with a stated reason`) : null,
    h('span', { class: 'small muted' }, `validation.json · ${groups.join(' · ')}`),
  );
}

export interface RawColumn<T> {
  label: string;
  num?: boolean;
  title?: string;
  cell: (row: T) => Child | Child[];
  className?: string;
}

/** A data table whose rows can carry classes (survivor / excluded highlighting). */
export function rawTable<T>(columns: RawColumn<T>[], rows: T[], rowClass?: (r: T) => string | null): { root: HTMLDivElement; tbody: HTMLTableSectionElement; render(rows: T[]): void } {
  const tbody = h('tbody');
  const render = (rs: T[]): void => {
    tbody.replaceChildren(
      ...rs.map((r) =>
        h(
          'tr',
          { class: rowClass?.(r) ?? null },
          ...columns.map((c) => h('td', { class: [c.num ? 'num' : '', c.className ?? ''].join(' ').trim() || null }, ...[c.cell(r)].flat())),
        ),
      ),
    );
  };
  render(rows);
  const root = h(
    'div',
    { class: 'table-wrap' },
    h(
      'table',
      { class: 'data' },
      h('thead', null, h('tr', null, ...columns.map((c) => h('th', { scope: 'col', class: c.num ? 'num' : null, title: c.title }, c.label)))),
      tbody,
    ),
  );
  return { root, tbody, render };
}

/** Signed difference with a fixed number of significant digits, e.g. "+7.6×10⁻⁴". */
export const fmtSigned = (x: number, sig = 2): string => (x > 0 ? '+' : x < 0 ? '−' : '') + fmtNum(Math.abs(x), sig);

/**
 * Update a mounted line chart in place: mutate the series objects (so the Table twin and theme
 * re-renders show the same numbers) and restyle the given traces.
 */
export async function restyleSeries(chart: ChartHandle, series: LineSeries[], indices: number[]): Promise<void> {
  await chart.ready;
  if (!chart.plot.isConnected || chart.plot.hidden) return;
  const P = await loadPlotly();
  const upd = { x: indices.map((i) => series[i]!.x), y: indices.map((i) => series[i]!.y) };
  await P.restyle(chart.plot, upd as Parameters<typeof P.restyle>[1], indices);
}

/** A small "label: value" definition list. */
export function kv(pairs: [string, Child | Child[]][]): HTMLDListElement {
  return h('dl', { class: 'kv' }, ...pairs.flatMap(([k, v]) => [h('dt', null, k), h('dd', null, ...[v].flat())]));
}
