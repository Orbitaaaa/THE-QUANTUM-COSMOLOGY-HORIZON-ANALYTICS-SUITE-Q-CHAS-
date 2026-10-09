/**
 * Plotly wrappers with the shared palette and chart anatomy: line, scatter with error bars,
 * bar and heatmap. Every chart
 *
 * * takes series colours from the validated categorical palette in fixed slot order
 *   (`color: 1..8`, a CSS colour, or 'context' for de-emphasized gray),
 * * uses transparent backgrounds, hairline solid grid, 2px lines, >= 8px markers with a 2px
 *   surface ring, a unified crosshair tooltip on line charts,
 * * re-renders itself on theme change and resizes with its container,
 * * carries a "Table" toggle with the same numbers (the accessible twin of the chart),
 * * disposes itself when its element leaves the document (or on `destroy()`).
 *
 * Plotly (cartesian bundle) is loaded lazily on the first chart.
 */

import type * as PlotlyNS from 'plotly.js';
import { divergingLut, plotlyScale, sequentialLut } from './colormap';
import { cssVar, currentTheme, onThemeChange, seriesColors } from './theme';
import { fmtNum, h, table } from './ui';

type PlotlyModule = typeof PlotlyNS;
type Layout = Partial<PlotlyNS.Layout>;
type Trace = Partial<PlotlyNS.PlotData> & Record<string, unknown>;

let plotlyPromise: Promise<PlotlyModule> | null = null;

/** Load the Plotly cartesian bundle once. */
export function loadPlotly(): Promise<PlotlyModule> {
  plotlyPromise ??= import('plotly.js-cartesian-dist-min').then(
    (m) => ((m as unknown as { default?: PlotlyModule }).default ?? m) as PlotlyModule,
  );
  return plotlyPromise;
}

export type SeriesColor = number | string | 'context';

export interface LineSeries {
  name: string;
  x: number[];
  y: (number | null)[];
  /** palette slot 1..8 (default: series index + 1), CSS colour, or 'context' */
  color?: SeriesColor;
  dash?: 'solid' | 'dot' | 'dash' | 'dashdot';
  /** soft ~10% area fill under the line */
  fill?: boolean;
  width?: number;
  /** marker on every point (accessibility markers) */
  markers?: boolean;
  hidden?: boolean;
}

export interface ScatterSeries {
  name: string;
  x: number[];
  y: (number | null)[];
  /** symmetric y error bars */
  error?: number[];
  /** asymmetric y error bars: [minus, plus] */
  errorAsym?: [number[], number[]];
  errorX?: number[];
  color?: SeriesColor;
  mode?: 'markers' | 'lines+markers';
  size?: number;
  symbol?: string;
  /** per-point hover labels */
  text?: string[];
}

export interface BarSeries {
  name: string;
  x: (number | string)[];
  y: number[];
  color?: SeriesColor;
  error?: number[];
}

export interface HeatmapData {
  z: (number | null)[][];
  x?: (number | string)[];
  y?: (number | string)[];
  /** 'sequential' for magnitude, 'diverging' (symmetric about zmid) for signed values */
  scale?: 'sequential' | 'diverging';
  zmid?: number;
  zmin?: number;
  zmax?: number;
  colorbarTitle?: string;
  /** hover label for z */
  zLabel?: string;
}

/** A vertical or horizontal reference line with a short label (e.g. Page time). */
export interface RefLine {
  axis: 'x' | 'y';
  value: number;
  label?: string;
}

export interface ChartOptions {
  xTitle?: string;
  yTitle?: string;
  xType?: 'linear' | 'log';
  yType?: 'linear' | 'log';
  xRange?: [number, number];
  yRange?: [number, number];
  /** plot height in px including the axis band (default 300) */
  height?: number;
  /** accessible name of the figure */
  ariaLabel?: string;
  caption?: string;
  refLines?: RefLine[];
  /** add a Table toggle (default true) */
  tableToggle?: boolean;
  /** formatter for table cells */
  format?: (v: number) => string;
  /** escape hatch merged into the generated layout */
  layout?: Layout;
}

export interface ChartHandle {
  /** the <figure> holding plot, table and toolbar */
  el: HTMLElement;
  plot: HTMLElement;
  ready: Promise<void>;
  destroy(): void;
}

interface Palette {
  text1: string;
  text2: string;
  text3: string;
  grid: string;
  axis: string;
  surface: string;
  surface2: string;
  line2: string;
  series: string[];
  other: string;
  sans: string;
  mono: string;
}

function palette(): Palette {
  return {
    text1: cssVar('--text-1'),
    text2: cssVar('--text-2'),
    text3: cssVar('--text-3'),
    grid: cssVar('--chart-grid'),
    axis: cssVar('--chart-axis'),
    surface: cssVar('--surface'),
    surface2: cssVar('--surface-2'),
    line2: cssVar('--line-2'),
    series: seriesColors(),
    other: cssVar('--series-other'),
    sans: cssVar('--font-sans'),
    mono: cssVar('--font-mono'),
  };
}

/** Resolve a series colour: slot number, 'context', explicit CSS colour, or default slot. */
export function resolveColor(c: SeriesColor | undefined, index: number, p = palette()): string {
  if (c === 'context') return p.other;
  if (typeof c === 'string') return c.startsWith('--') ? cssVar(c) : c;
  const slot = (c ?? index + 1) - 1;
  // never generate a 9th hue: past eight series fold into the "other" gray
  return p.series[slot] ?? p.other;
}

function axis(p: Palette, title: string | undefined, type: 'linear' | 'log' | undefined, range?: [number, number]) {
  return {
    title: title ? { text: title, font: { size: 12, color: p.text2 }, standoff: 8 } : undefined,
    type: type ?? 'linear',
    range: range && type === 'log' ? range.map(Math.log10) : range,
    gridcolor: p.grid,
    gridwidth: 1,
    zeroline: false,
    showline: true,
    linecolor: p.axis,
    linewidth: 1,
    ticks: 'outside' as const,
    ticklen: 4,
    tickcolor: p.axis,
    tickfont: { family: p.mono, size: 11, color: p.text3 },
    automargin: true,
    separatethousands: true,
    exponentformat: 'power' as const,
    showspikes: false,
  };
}

function baseLayout(p: Palette, o: ChartOptions, legend: boolean): Layout {
  const layout: Layout = {
    paper_bgcolor: 'rgba(0,0,0,0)',
    plot_bgcolor: 'rgba(0,0,0,0)',
    font: { family: p.sans, size: 12, color: p.text2 },
    margin: { l: 56, r: 18, t: legend ? 36 : 12, b: 46, pad: 0 },
    height: o.height ?? 300,
    xaxis: axis(p, o.xTitle, o.xType, o.xRange) as Layout['xaxis'],
    yaxis: axis(p, o.yTitle, o.yType, o.yRange) as Layout['yaxis'],
    showlegend: legend,
    legend: {
      orientation: 'h',
      x: 0,
      xanchor: 'left',
      y: 1.02,
      yanchor: 'bottom',
      bgcolor: 'rgba(0,0,0,0)',
      font: { color: p.text2, size: 12 },
      itemclick: 'toggleothers',
      itemdoubleclick: 'toggle',
    },
    hoverlabel: {
      bgcolor: p.surface2,
      bordercolor: p.line2,
      font: { family: p.mono, size: 12, color: p.text1 },
      align: 'left',
    },
    dragmode: 'zoom',
  };
  const refs = o.refLines ?? [];
  layout.shapes = refs.map((r) => ({
    type: 'line',
    xref: r.axis === 'x' ? 'x' : 'paper',
    yref: r.axis === 'y' ? 'y' : 'paper',
    x0: r.axis === 'x' ? r.value : 0,
    x1: r.axis === 'x' ? r.value : 1,
    y0: r.axis === 'y' ? r.value : 0,
    y1: r.axis === 'y' ? r.value : 1,
    line: { color: p.text3, width: 1 },
    layer: 'below',
  })) as Layout['shapes'];
  layout.annotations = refs
    .filter((r) => r.label)
    .map((r) => ({
      text: r.label,
      showarrow: false,
      xref: r.axis === 'x' ? 'x' : 'paper',
      yref: r.axis === 'y' ? 'y' : 'paper',
      x: r.axis === 'x' ? r.value : 1,
      y: r.axis === 'y' ? r.value : 1,
      xanchor: r.axis === 'x' ? 'left' : 'right',
      yanchor: r.axis === 'x' ? 'top' : 'bottom',
      xshift: r.axis === 'x' ? 4 : 0,
      font: { family: p.mono, size: 10.5, color: p.text3 },
    })) as Layout['annotations'];
  return layout;
}

function mergeLayout(a: Layout, b?: Layout): Layout {
  if (!b) return a;
  const out: Record<string, unknown> = { ...a };
  for (const [k, v] of Object.entries(b)) {
    const cur = out[k];
    out[k] = cur && typeof cur === 'object' && v && typeof v === 'object' && !Array.isArray(v) ? { ...cur, ...v } : v;
  }
  return out as Layout;
}

const CONFIG: Partial<PlotlyNS.Config> = {
  displaylogo: false,
  responsive: false,
  displayModeBar: 'hover',
  modeBarButtonsToRemove: ['lasso2d', 'select2d', 'autoScale2d', 'toggleSpikelines', 'hoverClosestCartesian', 'hoverCompareCartesian'],
  toImageButtonOptions: { format: 'png', scale: 2 },
};

type Builder = (p: Palette) => { data: Trace[]; layout: Layout };

/** Generic mount: figure with toolbar, plot, optional table twin; theme + resize aware. */
function mountChart(
  el: HTMLElement,
  build: Builder,
  o: ChartOptions,
  tableFor: () => HTMLElement,
): ChartHandle {
  const plot = h('div', { class: 'chart-plot', role: 'img', 'aria-label': o.ariaLabel ?? o.yTitle ?? 'chart' });
  plot.style.minHeight = `${o.height ?? 300}px`;
  const tableBox = h('div', { class: 'chart-table', hidden: true });
  const toggleBtn = h('button', { type: 'button', class: 'btn', 'aria-pressed': 'false', title: 'Show the numbers as a table' }, 'Table');
  const figure = h(
    'figure',
    { class: 'chart-figure' },
    o.tableToggle === false ? null : h('div', { class: 'chart-toolbar' }, toggleBtn),
    plot,
    tableBox,
    o.caption ? h('figcaption', { class: 'chart-caption' }, o.caption) : null,
  );
  el.append(figure);

  let disposed = false;
  let Plotly: PlotlyModule | null = null;
  const render = async (): Promise<void> => {
    if (disposed) return;
    if (!figure.isConnected && Plotly) return destroy();
    const { data, layout } = build(palette());
    Plotly ??= await loadPlotly();
    if (disposed) return;
    await Plotly.react(plot, data as PlotlyNS.Data[], layout, CONFIG);
  };

  toggleBtn.addEventListener('click', () => {
    const showTable = tableBox.hidden;
    tableBox.hidden = !showTable;
    plot.hidden = showTable;
    toggleBtn.setAttribute('aria-pressed', String(showTable));
    if (showTable) tableBox.replaceChildren(tableFor());
    else void render();
  });

  let raf = 0;
  const ro = new ResizeObserver(() => {
    cancelAnimationFrame(raf);
    raf = requestAnimationFrame(() => {
      if (!figure.isConnected) return destroy();
      if (Plotly && !plot.hidden && plot.offsetWidth > 0) void Plotly.Plots.resize(plot);
    });
  });
  ro.observe(figure);
  const offTheme = onThemeChange(() => void render());

  function destroy(): void {
    if (disposed) return;
    disposed = true;
    ro.disconnect();
    offTheme();
    cancelAnimationFrame(raf);
    if (Plotly) Plotly.purge(plot);
  }

  const ready = render().catch((err: unknown) => {
    plot.replaceChildren(h('div', { class: 'empty-state' }, `Chart failed to render: ${String(err)}`));
  });
  return { el: figure, plot, ready, destroy };
}

const fmtDefault = (v: number): string => fmtNum(v, 5);

function seriesTable(
  series: { name: string; x: (number | string)[]; y: (number | null)[]; error?: number[] }[],
  xLabel: string,
  fmt: (v: number) => string,
): HTMLElement {
  const fmtCell = (v: number | string | null | undefined): string =>
    v === null || v === undefined ? '—' : typeof v === 'number' ? fmt(v) : v;
  const first = series[0];
  const shared = !!first && series.every((s) => s.x.length === first.x.length && s.x.every((v, i) => v === first.x[i]));
  if (shared && first) {
    const rows = first.x.map((x, i) => ({ x, i }));
    return table(
      [
        { key: 'x', label: xLabel, align: 'num', format: (r) => fmtCell(r.x) },
        ...series.map((s, k) => ({
          key: `s${k}`,
          label: s.name,
          align: 'num' as const,
          format: (r: { i: number }) => fmtCell(s.y[r.i]) + (s.error ? ` ± ${fmtCell(s.error[r.i])}` : ''),
        })),
      ],
      rows,
    );
  }
  const rows = series.flatMap((s) => s.x.map((x, i) => ({ series: s.name, x, y: s.y[i] ?? null, e: s.error?.[i] })));
  return table(
    [
      { key: 'series', label: 'Series' },
      { key: 'x', label: xLabel, align: 'num', format: (r) => fmtCell(r.x) },
      { key: 'y', label: 'Value', align: 'num', format: (r) => fmtCell(r.y) + (r.e !== undefined ? ` ± ${fmtCell(r.e)}` : '') },
    ],
    rows,
  );
}

/** Multi-series line chart with a crosshair readout of every series at the hovered x. */
export function lineChart(el: HTMLElement, series: LineSeries[], o: ChartOptions = {}): ChartHandle {
  const build: Builder = (p) => {
    const data: Trace[] = series.map((s, i) => {
      const color = resolveColor(s.color, i, p);
      return {
        type: 'scatter',
        mode: s.markers ? 'lines+markers' : 'lines',
        name: s.name,
        x: s.x,
        y: s.y,
        visible: s.hidden ? 'legendonly' : true,
        line: { color, width: s.width ?? 2, dash: s.dash ?? 'solid', shape: 'linear' },
        marker: { color, size: 8, line: { color: p.surface, width: 2 } },
        fill: s.fill ? 'tozeroy' : 'none',
        fillcolor: s.fill ? withAlpha(color, 0.1) : undefined,
        hovertemplate: `%{y:.4~g}<extra>${escapeHover(s.name)}</extra>`,
        connectgaps: false,
      };
    });
    const layout = baseLayout(p, o, series.length >= 2);
    layout.hovermode = 'x unified';
    (layout.xaxis as Record<string, unknown>).showspikes = true;
    Object.assign(layout.xaxis as object, { spikemode: 'across', spikethickness: 1, spikecolor: p.text3, spikedash: 'solid', spikesnap: 'data' });
    return { data, layout: mergeLayout(layout, o.layout) };
  };
  return mountChart(el, build, o, () => seriesTable(series, o.xTitle ?? 'x', o.format ?? fmtDefault));
}

/** Scatter plot with optional (symmetric or asymmetric) error bars. */
export function scatterChart(el: HTMLElement, series: ScatterSeries[], o: ChartOptions = {}): ChartHandle {
  const build: Builder = (p) => {
    const data: Trace[] = series.map((s, i) => {
      const color = resolveColor(s.color, i, p);
      const errY: Record<string, unknown> | undefined = s.errorAsym
        ? { type: 'data', symmetric: false, array: s.errorAsym[1], arrayminus: s.errorAsym[0] }
        : s.error
          ? { type: 'data', array: s.error }
          : undefined;
      return {
        type: 'scatter',
        mode: s.mode ?? 'markers',
        name: s.name,
        x: s.x,
        y: s.y,
        text: s.text,
        marker: { color, size: s.size ?? 8, symbol: s.symbol ?? 'circle', line: { color: p.surface, width: 2 } },
        line: { color, width: 2 },
        error_y: errY ? ({ ...errY, color, thickness: 1.5, width: 0 } as Trace['error_y']) : undefined,
        error_x: s.errorX ? ({ type: 'data', array: s.errorX, color, thickness: 1.5, width: 0 } as Trace['error_x']) : undefined,
        hovertemplate: `${s.text ? '%{text}<br>' : ''}x %{x:.4~g}<br>y %{y:.4~g}<extra>${escapeHover(s.name)}</extra>`,
      };
    });
    const layout = baseLayout(p, o, series.length >= 2);
    layout.hovermode = 'closest';
    layout.hoverdistance = 24;
    return { data, layout: mergeLayout(layout, o.layout) };
  };
  return mountChart(el, build, o, () => seriesTable(series, o.xTitle ?? 'x', o.format ?? fmtDefault));
}

/** Vertical bar chart (grouped when several series), thin bars with rounded data ends. */
export function barChart(el: HTMLElement, series: BarSeries[], o: ChartOptions = {}): ChartHandle {
  const build: Builder = (p) => {
    const data: Trace[] = series.map((s, i) => {
      const color = resolveColor(s.color, i, p);
      return {
        type: 'bar',
        name: s.name,
        x: s.x,
        y: s.y,
        marker: { color, cornerradius: 4, line: { width: 0 } },
        error_y: s.error ? { type: 'data', array: s.error, color: p.text3, thickness: 1.5, width: 0 } : undefined,
        hovertemplate: `%{x}: %{y:.4~g}<extra>${escapeHover(s.name)}</extra>`,
      };
    });
    const layout = baseLayout(p, o, series.length >= 2);
    layout.hovermode = 'closest';
    layout.barmode = 'group';
    const n = Math.max(1, ...series.map((s) => s.x.length));
    layout.bargap = n <= 4 ? 0.6 : n <= 10 ? 0.45 : 0.25;
    layout.bargroupgap = 0.12;
    return { data, layout: mergeLayout(layout, o.layout) };
  };
  return mountChart(el, build, o, () => seriesTable(series, o.xTitle ?? 'Category', o.format ?? fmtDefault));
}

/** Heatmap with a one-hue sequential or blue-gray-red diverging scale for the current theme. */
export function heatmapChart(el: HTMLElement, d: HeatmapData, o: ChartOptions = {}): ChartHandle {
  const build: Builder = (p) => {
    const theme = currentTheme();
    const diverging = d.scale === 'diverging';
    const lut = diverging ? divergingLut(theme) : sequentialLut(theme);
    const cols = d.z[0]?.length ?? 0;
    const gap = cols <= 40 && d.z.length <= 40 ? 1 : 0;
    let zmin = d.zmin;
    let zmax = d.zmax;
    if (diverging && (zmin === undefined || zmax === undefined)) {
      const mid = d.zmid ?? 0;
      let lim = 0;
      for (const row of d.z) for (const v of row) if (v !== null && Number.isFinite(v)) lim = Math.max(lim, Math.abs(v - mid));
      zmin = mid - lim;
      zmax = mid + lim;
    }
    const trace: Trace = {
      type: 'heatmap',
      z: d.z,
      x: d.x,
      y: d.y,
      zmin,
      zmax,
      colorscale: plotlyScale(lut, 17),
      xgap: gap,
      ygap: gap,
      hoverongaps: false,
      hovertemplate: `x %{x}<br>y %{y}<br>${escapeHover(d.zLabel ?? 'value')} %{z:.4~g}<extra></extra>`,
      colorbar: {
        thickness: 10,
        outlinewidth: 0,
        tickfont: { family: p.mono, size: 10.5, color: p.text3 },
        title: d.colorbarTitle ? { text: d.colorbarTitle, font: { size: 11, color: p.text2 }, side: 'right' } : undefined,
        len: 0.9,
      },
    };
    const layout = baseLayout(p, o, false);
    layout.hovermode = 'closest';
    (layout.xaxis as Record<string, unknown>).showgrid = false;
    (layout.yaxis as Record<string, unknown>).showgrid = false;
    return { data: [trace], layout: mergeLayout(layout, o.layout) };
  };
  const tableFor = (): HTMLElement => {
    const maxRows = 200;
    const rows = d.z.slice(0, maxRows).map((row, i) => ({ y: d.y?.[i] ?? i, row }));
    const xs = d.x ?? (d.z[0] ?? []).map((_, j) => j);
    const fmt = o.format ?? fmtDefault;
    return h(
      'div',
      null,
      table(
        [
          { key: 'y', label: o.yTitle ?? 'y', align: 'num', format: (r) => String(r.y) },
          ...xs.map((x, j) => ({
            key: `c${j}`,
            label: String(x),
            align: 'num' as const,
            format: (r: { row: (number | null)[] }) => (r.row[j] === null || r.row[j] === undefined ? '—' : fmt(r.row[j]!)),
          })),
        ],
        rows,
      ),
      d.z.length > maxRows ? h('div', { class: 'chart-caption' }, `First ${maxRows} of ${d.z.length} rows`) : null,
    );
  };
  return mountChart(el, build, o, tableFor);
}

const escapeHover = (s: string): string => s.replace(/[<>&]/g, (c) => ({ '<': '&lt;', '>': '&gt;', '&': '&amp;' })[c]!);

/** `#rrggbb` + alpha -> rgba(); other colour syntaxes are returned unchanged. */
export function withAlpha(color: string, alpha: number): string {
  const m = /^#([0-9a-f]{6})$/i.exec(color.trim());
  if (!m) return color;
  const v = parseInt(m[1]!, 16);
  return `rgba(${(v >> 16) & 255},${(v >> 8) & 255},${v & 255},${alpha})`;
}
