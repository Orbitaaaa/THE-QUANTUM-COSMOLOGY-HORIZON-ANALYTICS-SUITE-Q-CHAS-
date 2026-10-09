/**
 * Shared DOM pieces of the Module 2 panels (capacity, page-curve, qes, hayden-preskill,
 * gw-horizons): the "What the paper proves" list with provenance badges, the ours-vs-paper
 * validation block filtered by group, and the diary mutual-information chart of [BH] Fig. 2.
 */

import { type ChartHandle, type ChartOptions, type LineSeries, lineChart, loadPlotly } from '../lib/charts';
import type { HorizonBundle, LoadResult, ValidationBundle, ValidationRow } from '../lib/data';
import type { PanelContext } from '../lib/panel';
import { type Child, type ProvenanceKind, badge, fmtNum, h, note, statusPill, table, tiles } from '../lib/ui';

export interface ProofItem {
  kind: ProvenanceKind;
  src: string;
  text: Child | Child[];
}

/** A list of statements, each led by its provenance badge (PROVED / COMPUTED / ...). */
export function proofList(items: ProofItem[], stacked = false): HTMLUListElement {
  return h(
    'ul',
    { class: `m2-proofs${stacked ? ' stacked' : ''}` },
    ...items.map((it) => h('li', null, badge(it.kind, it.src), h('span', null, ...[it.text].flat()))),
  );
}

/** The collapsible "What the paper proves" box. */
export const provesNote = (items: ProofItem[], open = true, title = 'What the paper proves'): HTMLDetailsElement =>
  note(proofList(items), title, open);

const fmtVal = (v: number): string => (Math.abs(v) >= 1e6 || (Math.abs(v) < 1e-3 && v !== 0) ? fmtNum(v, 4) : fmtNum(v, 6));

/**
 * Ours-vs-paper rows of validation.json for the given groups: a summary pill and a table.
 * Returns null when validation.json is not built (the panel then shows its own tables only).
 */
export function validationBlock(
  res: LoadResult<ValidationBundle>,
  groups: string[],
  opts: { maxRows?: number; filter?: (r: ValidationRow) => boolean; compact?: boolean; summaryOnly?: boolean } = {},
): HTMLElement | null {
  if (!res.ok) return null;
  const rows = res.data.rows.filter((r) => groups.includes(r.group) && (opts.filter?.(r) ?? true));
  if (!rows.length) return null;
  const nMatch = rows.filter((r) => r.verdict === 'match').length;
  const nDiff = rows.length - nMatch;
  const shown = opts.maxRows ? rows.slice(0, opts.maxRows) : rows;
  if (opts.summaryOnly)
    return h(
      'div',
      { class: 'row m2-validation-head' },
      statusPill(nDiff ? 'warning' : 'good', `${nMatch} of ${rows.length} match the printed values`),
      h('span', { class: 'small muted' }, 'validation.json · ', groups.join(' · ')),
    );
  return h(
    'div',
    { class: 'm2-validation' },
    h(
      'div',
      { class: 'row m2-validation-head' },
      statusPill(nDiff ? 'warning' : 'good', `${nMatch} of ${rows.length} match the printed values`),
      nDiff ? h('span', { class: 'small muted' }, `${nDiff} differ; the note on each row says why`) : null,
      h('span', { class: 'small muted' }, 'validation.json · ', groups.join(' · ')),
    ),
    table(
      opts.compact
        ? [
            { key: 'quantity', label: 'Quantity', className: 'wrap' },
            { key: 'ours', label: 'Q-CHAS', align: 'num', format: (r) => fmtVal(r.ours) },
            { key: 'paper', label: 'Paper', align: 'num', format: (r) => fmtVal(r.paper) },
            { key: 'verdict', label: '', format: (r) => h('span', { title: `relative difference ${fmtNum(r.rel_diff, 2)}${r.note ? `; ${r.note}` : ''}` }, statusPill(r.verdict === 'match' ? 'good' : 'warning', r.verdict)) },
          ]
        : [
            { key: 'quantity', label: 'Quantity', className: 'wrap' },
            { key: 'ours', label: 'Q-CHAS', align: 'num', format: (r) => fmtVal(r.ours) },
            { key: 'paper', label: 'Paper', align: 'num', format: (r) => fmtVal(r.paper) },
            { key: 'rel_diff', label: 'Rel. diff.', align: 'num', format: (r) => fmtNum(r.rel_diff, 2) },
            {
              key: 'verdict',
              label: 'Verdict',
              format: (r) => [statusPill(r.verdict === 'match' ? 'good' : 'warning', r.verdict), r.note ? h('div', { class: 'small muted m2-note' }, r.note) : null],
              className: 'wrap',
            },
          ],
      shown,
    ),
    opts.maxRows && rows.length > opts.maxRows ? h('div', { class: 'chart-caption' }, `First ${opts.maxRows} of ${rows.length} rows`) : null,
  );
}

/** Stat tiles in a fixed number of columns (2 below 1100 px). */
export function tilesN(n: 2 | 3 | 4 | 5, ...items: HTMLElement[]): HTMLDivElement {
  const t = tiles(...items);
  t.classList.add(`m2-cols-${n}`);
  return t;
}

/** A small "label: value" definition list. */
export function kv(pairs: [string, Child | Child[]][]): HTMLDListElement {
  return h('dl', { class: 'kv' }, ...pairs.flatMap(([k, v]) => [h('dt', null, k), h('dd', null, ...[v].flat())]));
}

/**
 * [BH] Fig. 2 right: diary mutual information I(Ref : R) in bits for the early and late diary,
 * with the throw-in steps, the Page time and the paper's recovery steps as reference lines.
 */
export function diaryChart(host: HTMLElement, ev: HorizonBundle['evaporation'], height = 300): ChartHandle {
  const early = ev.early_diary;
  const late = ev.late_diary;
  const p = ev.paper;
  const pageStep = p.max_step ?? early.max_step;
  return lineChart(
    host,
    [
      { name: `diary thrown in at step ${early.diary_step ?? '—'}`, x: early.steps, y: early.diary_mi, color: 3, markers: true },
      { name: `diary thrown in at step ${late.diary_step ?? '—'}`, x: late.steps, y: late.diary_mi, color: 7, markers: true },
    ],
    {
      xTitle: 'emitted qubits',
      yTitle: 'diary information in radiation (bits)',
      yRange: [-0.05, 2.1],
      height,
      refLines: [
        ...(early.diary_step ? [{ axis: 'x' as const, value: early.diary_step, label: 'early in' }] : []),
        ...(pageStep ? [{ axis: 'x' as const, value: pageStep, label: 'Page time' }] : []),
        ...(late.diary_step ? [{ axis: 'x' as const, value: late.diary_step, label: 'late in' }] : []),
        // recovery threshold, unlabelled (the caption names it) so it does not cover the markers
        { axis: 'y' as const, value: p.recovery_threshold_bits ?? 1.9 },
      ],
      format: (v) => fmtNum(v, 4),
      ariaLabel: 'Diary mutual information against emitted qubits for an early and a late diary',
    },
  );
}

/** Register a chart for disposal and return it. */
export function keep<T extends { destroy(): void }>(ctx: PanelContext, x: T): T {
  ctx.onCleanup(() => x.destroy());
  return x;
}

/**
 * Update a mounted line chart in place (slider-driven live physics) instead of rebuilding it:
 * the series objects and reference lines passed to lineChart() are mutated, so a later theme
 * re-render or the Table twin show the same numbers, and Plotly restyles the traces now.
 * `refLines` must list the chart's reference lines in their original order.
 */
export async function updateLineChart(
  chart: ChartHandle,
  series: LineSeries[],
  next: { y?: (number | null)[][]; x?: number[][]; refLines?: ChartOptions['refLines']; yRange?: [number, number] },
  opts?: ChartOptions,
): Promise<void> {
  next.y?.forEach((y, i) => series[i] && (series[i].y = y));
  next.x?.forEach((x, i) => series[i] && (series[i].x = x));
  if (opts && next.refLines) opts.refLines = next.refLines;
  if (opts && next.yRange) opts.yRange = next.yRange;
  await chart.ready;
  if (!chart.plot.isConnected) return;
  const P = await loadPlotly();
  const upd: Record<string, unknown> = { y: series.map((s) => s.y) };
  if (next.x) upd.x = series.map((s) => s.x);
  await P.restyle(chart.plot, upd as Parameters<typeof P.restyle>[1]);
  if (next.refLines || next.yRange) {
    const rel: Record<string, unknown> = {};
    if (next.yRange) rel['yaxis.range'] = next.yRange;
    let a = 0;
    next.refLines?.forEach((r, i) => {
      const key = r.axis === 'x' ? 'x' : 'y';
      rel[`shapes[${i}].${key}0`] = r.value;
      rel[`shapes[${i}].${key}1`] = r.value;
      if (r.label) {
        rel[`annotations[${a}].${key}`] = r.value;
        rel[`annotations[${a}].text`] = r.label;
        a++;
      }
    });
    await P.relayout(chart.plot, rel as Parameters<typeof P.relayout>[1]);
  }
}
