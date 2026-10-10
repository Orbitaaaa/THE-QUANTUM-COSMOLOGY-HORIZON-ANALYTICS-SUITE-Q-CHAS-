/**
 * 1.4 Resonance Shell Inspector ([TOP] Secs. 10-11.8).
 *
 * The law of limitation (KL of the l <= 30 CMB against d/D with the Pinsker bound, and the
 * Sec. 11.2 information ceiling of every probe), topological resonance on dual-lattice shells
 * (KL with a primordial line at known orientation, ours against the printed table), a shell
 * explorer with the shell directions computed live, the rank-N Woodbury benchmark, the
 * pre-registered resonance search re-run on the build map (WMAP 9-yr ILC) beside the paper's
 * Planck values, and the Sec. 11.8 exits.
 *
 * Bundle fields: resonance.{kl_table, shells, woodbury_benchmark, search, search_paper,
 * single_field, ceiling, exits_3d}; cmb.kl.{d_over_D, cubic, BCC, FCC, paper, pinsker_bound,
 * pinsker_bound_paper}; validation rows of the KL and resonance groups.
 */

import './m1.css';
import './resonance.css';
import { type LineSeries, type ScatterSeries, barChart, lineChart, scatterChart } from '../lib/charts';
import type { CmbBundle, PaperValue, ResonanceBundle, ResonanceKlRow, ResonanceShell, SingleFieldRow } from '../lib/data';
import type { PanelContext } from '../lib/panel';
import { cssVar } from '../lib/theme';
import { badge, clear, fmtFixed, fmtNum, grid, h, loadFailure, note, panel, segmented, select, statTile, statusPill } from '../lib/ui';
import { type RawColumn, SHAPE_SLOT, keep, kv, oursPaper, proofList, rawTable, sectionLabel, tilesN, validationRows, validationSummary } from './m1-ui';
import { azEl, directionFamily, familyCounts, histogram, median, shellKChi, shellVectors } from './resonance-logic';

export const title = 'Resonance Shell Inspector';
export const subtitle =
  'The law of limitation, and its one physical loophole: a primordial line on a dual-lattice shell. KL gain at known orientation, the rank-N Woodbury test, and the pre-registered search on a real sky.';

const SHAPES = ['cubic', 'BCC', 'FCC'] as const;
const SHAPE_LABEL: Record<string, string> = { cubic: 'Cubic', BCC: 'BCC', FCC: 'FCC' };
const VAL_GROUPS = ['Detectability KL ([TOP] Theorem P)', 'Topological resonance ([TOP] Sec. 11.6)'];

/** Paper cell as text ("—" when the paper prints nothing). */
const paperText = (v: PaperValue | undefined, d = 2): string | null =>
  v === null || v === undefined ? null : typeof v === 'number' ? fmtFixed(v, d) : v;
const sign = (v: number, d = 2): string => (v > 0 ? '+' : v < 0 ? '−' : '') + Math.abs(v).toFixed(d);
const pretty = (s: string): string =>
  s
    .replace(/\|m\|\^2/g, '|m|²')
    .replace(/\bpi\b/g, 'π')
    .replace(/<=/g, '≤');
/** "BCC_1.05" -> "BCC · 1.05 D" */
const keyLabel = (k: string): string => {
  const [s, d] = k.split('_');
  return `${SHAPE_LABEL[s ?? ''] ?? s} · ${d} D`;
};

export async function mount(el: HTMLElement, ctx: PanelContext): Promise<void> {
  ctx.setHeaderExtras([badge('proved', '[TOP] Theorems P, V, L'), badge('computed', 'qchas.cmb.resonance'), badge('quoted', '[TOP] Secs. 11.2–11.8')]);
  const [res, cmb, val] = await Promise.all([ctx.load('resonance.json'), ctx.load('cmb.json'), ctx.load('validation.json')]);
  if (ctx.unmounted) return;
  if (!res.ok) {
    const p = panel({ title: 'Resonance Shell Inspector', span: 12 });
    p.body.append(loadFailure(res));
    el.append(grid(p.root));
    return;
  }
  const R = res.data;
  const C = cmb.ok ? cmb.data : null;

  mountHeadline(el, R, validationSummary(validationRows(val, VAL_GROUPS), VAL_GROUPS));
  mountLimitation(el, ctx, R, C, cmb.ok ? null : loadFailure(cmb));
  mountResonanceTable(el, ctx, R);
  mountShellExplorer(el, ctx, R);
  mountSearch(el, ctx, R);
  mountWoodbury(el, ctx, R);
  mountExits(el, R);
}

// --- headline -------------------------------------------------------------------------------------

function mountHeadline(el: HTMLElement, R: ResonanceBundle, validation: HTMLElement | null): void {
  const s = R.search;
  const sp = s?.paper ?? R.search_paper ?? null;
  const gainRow = [...R.kl_table].sort((a, b) => b.gain - a.gain)[0];
  const sf = R.single_field;
  const p = panel({
    title: 'Verdict',
    subtitle: 'What the pre-registered search finds, on the paper’s Planck sky and on our independent WMAP re-analysis',
    badges: [badge('measured', 'Planck 2018 SMICA · [TOP] Sec. 11.6'), badge('computed', s?.map ?? 'build map')],
    span: 12,
  });
  const verdictTitle = s
    ? `No resonance: not in Planck (paper, S_max = ${fmtFixed(sp?.S_max, 2)}, p = ${fmtFixed(sp?.p_max, 2)}), not in our ${s.map} re-analysis (S_max = ${fmtFixed(s.S_max, 2)}, p = ${fmtFixed(s.p_max, 2)}).`
    : `No resonance in Planck (paper, S_max = ${fmtFixed(sp?.S_max, 2)}, p = ${fmtFixed(sp?.p_max, 2)}). Our re-analysis is not in this build.`;
  p.body.append(
    h(
      'div',
      { class: 'res-headline' },
      h(
        'div',
        { class: 'res-verdict' },
        h('h3', { class: 'res-verdict-title' }, verdictTitle),
        proofList([
          {
            kind: 'measured',
            src: '[TOP] Sec. 11.6',
            text: s
              ? `The Planck sky is typical of skies without topology: S_max = ${fmtFixed(sp?.S_max, 2)} lies below the null median ${fmtFixed(sp?.null_median_S_max, 2)} and S_avg = ${sign(sp?.S_avg ?? NaN)}. On ${s.map} S_max = ${fmtFixed(s.S_max, 2)} ${s.S_max < (s.null_median_S_max ?? median(s.null.S_max)) ? 'also lies below' : 'lies above'} our null median ${fmtFixed(s.null_median_S_max ?? median(s.null.S_max), 2)}, S_avg = ${sign(s.S_avg)}. No resonance is seen; it is not excluded either.`
              : `The Planck sky is typical of skies without topology: S_max = ${fmtFixed(sp?.S_max, 2)} lies below the null median ${fmtFixed(sp?.null_median_S_max, 2)}. No resonance is seen; it is not excluded either.`,
          },
          {
            kind: 'computed',
            src: '[TOP] Sec. 11.6',
            text: [h('strong', null, 'Up to 64× more information only at known orientation. '), 'With the orientation unknown the Bayes factor averages over all of SO(3); with the 160-model grid this consumes most of the gain, so a resonant torus at 1.2 D would be found in only 30–35 % of skies.'],
          },
          {
            kind: 'computed',
            src: '[TOP] Sec. 11.8',
            text: [h('strong', null, 'A torus as large as 3D needs a multi-field line mapped in 3D. '), 'The CMB alone cannot decide it; single-field features cannot either. It is conditional on nature: current data show no line.'],
          },
          {
            kind: 'assumed',
            src: 'this build',
            text: /^planck/i.test(s?.map ?? '')
              ? `Our numbers are an independent recomputation on ${s?.map}, run with our own code and null skies; the paper's values are quoted beside them.`
              : `Our numbers are an independent re-analysis of ${s?.map ?? 'the build map'} (Planck PR3 was not reachable from the build host), not a reproduction of the Planck values.`,
          },
        ]),
      ),
      h(
        'div',
        null,
        tilesN(
          2,
          statTile({
            label: 'S_max · Planck (paper)',
            value: fmtFixed(sp?.S_max, 2),
            sub: [badge('quoted', '[TOP] Sec. 11.6'), `p = ${fmtFixed(sp?.p_max, 2)} · null median ${fmtFixed(sp?.null_median_S_max, 2)}`],
          }),
          statTile({
            label: `S_max · ${s?.map ?? 'build map'} (ours)`,
            value: s ? fmtFixed(s.S_max, 2) : null,
            sub: s ? [badge('computed'), `p = ${fmtFixed(s.p_max, 2)} · null median ${fmtFixed(s.null_median_S_max ?? median(s.null.S_max), 2)}`] : 'not in this build',
          }),
          statTile({
            label: 'Largest KL gain, known orientation',
            value: gainRow ? `${fmtNum(gainRow.gain, 3)}×` : null,
            sub: gainRow
              ? [badge('computed'), `cubic ${gainRow.d_over_D} D · paper ${paperText(gainRow.paper?.gain, 0) ?? '—'}×`]
              : undefined,
          }),
          statTile({
            label: 'Shells where single-field inflation caps the bump at f < 0.05',
            value: sf ? `${sf['n_f_below_0.05']} of ${sf.n_shells}` : null,
            sub: sf
              ? [badge('computed', '[TOP] Sec. 11.8'), `paper ${sf.paper?.['shells_with_f_below_0.05']?.count ?? '—'} of ${sf.paper?.['shells_with_f_below_0.05']?.of ?? '—'}: no gain there`]
              : undefined,
          }),
        ),
        validation,
      ),
    ),
  );
  el.append(grid(p.root));
}

// --- (a) law of limitation and the ceiling -----------------------------------------------------------

/** Label of a horizontal line at y on a log axis (Plotly annotations take log10 there). */
function logLineLabel(y: number, text: string): Record<string, unknown> {
  return {
    text,
    showarrow: false,
    xref: 'paper',
    yref: 'y',
    x: 0,
    y: Math.log10(y),
    xanchor: 'left',
    yanchor: 'bottom',
    xshift: 4,
    // no colour: inherits the theme-aware layout font colour on every re-render
    font: { family: cssVar('--font-mono'), size: 10.5 },
  };
}

function mountLimitation(el: HTMLElement, ctx: PanelContext, R: ResonanceBundle, C: CmbBundle | null, cmbFailure: HTMLElement | null): void {
  const a = panel({
    title: 'Law of limitation · KL of the ℓ ≤ 30 CMB',
    subtitle: 'Exact Gaussian KL divergence, torus against infinite ΛCDM space, smooth spectrum; ours (lines) against the Theorem P table (diamonds)',
    badges: [badge('proved', '[TOP] Theorem P'), badge('computed')],
    span: 7,
    stretch: true,
  });
  const b = panel({
    title: 'Information ceiling · every probe',
    subtitle: 'Rigorous lower bounds on the KL of a perfect map (cubic torus, smooth spectrum)',
    badges: [badge('proved', '[TOP] Theorems V, L'), badge('quoted', '[TOP] Sec. 11.2')],
    span: 5,
    stretch: true,
  });
  el.append(grid(a.root, b.root));

  // (a) KL / Pinsker
  if (!C) a.body.append(cmbFailure ?? h('div'));
  else {
    const K = C.kl;
    const x = K.d_over_D;
    const host = h('div');
    let view = 'kl';
    const seg = segmented({
      label: 'Quantity',
      value: view,
      options: [
        { value: 'kl', label: 'KL divergence' },
        { value: 'pinsker', label: 'Pinsker bound on power − size' },
      ],
      onChange: (v) => {
        view = v;
        draw();
      },
    });
    let chart: { destroy(): void } | null = null;
    const draw = (): void => {
      chart?.destroy();
      clear(host);
      const series: ScatterSeries[] = [];
      for (const s of SHAPES) {
        const ours = view === 'kl' ? K[s] : (K.pinsker_bound[s] ?? []);
        series.push({ name: `${SHAPE_LABEL[s]} · ours`, x, y: ours, color: SHAPE_SLOT[s], mode: 'lines+markers', size: 7 });
      }
      for (const s of SHAPES) {
        const paper = view === 'kl' ? x.map((d) => K.paper[s]?.[String(d)] ?? K.paper[s]?.[d.toFixed(1)] ?? null) : (K.pinsker_bound_paper?.[s] ?? []);
        series.push({ name: `${SHAPE_LABEL[s]} · paper`, x, y: paper, color: SHAPE_SLOT[s], symbol: 'diamond', size: 10 });
      }
      chart = keep(
        ctx,
        scatterChart(host, series, {
          xTitle: 'torus size d / D (D = 27.74 Gpc)',
          yTitle: view === 'kl' ? 'KL divergence (nats)' : 'bound on power − size',
          yType: view === 'kl' ? 'log' : 'linear',
          yRange: view === 'kl' ? undefined : [0, 1.05],
          height: 330,
          refLines: view === 'kl' ? [{ axis: 'y', value: 1 }] : [],
          layout: view === 'kl' ? { annotations: [logLineLabel(1, 'KL = 1 nat')] as never, yaxis: { dtick: 1 } as never } : undefined,
          ariaLabel: view === 'kl' ? 'KL divergence of the CMB against torus size' : 'Pinsker bound against torus size',
          format: (v) => fmtNum(v, 4),
        }),
      );
    };
    a.body.append(h('div', { class: 'm1-toolbar' }, seg.root), host);
    draw();
    a.body.append(
      h(
        'p',
        { class: 'm1-takeaway' },
        'Power minus false-alarm rate of any test is at most √(KL/2) (Pinsker). With a smooth spectrum the CMB information collapses from 38 nats (FCC at 1.0 D) to below 1 nat by 1.15–1.2 D: past that size no CMB test can decide the topology. ',
        badge('proved', '[TOP] Theorem P'),
      ),
    );
  }

  // (b) ceiling
  const ceil = R.ceiling ?? C?.ceiling ?? null;
  if (!ceil) b.body.append(h('p', { class: 'muted small' }, 'The ceiling table is not in this build.'));
  else {
    const series: LineSeries[] = [
      { name: 'Causal ball, perfect 3D map (every probe)', x: ceil.d_over_D, y: ceil.causal_ball, color: 1, markers: true },
      { name: '… its boundary sphere alone', x: ceil.d_over_D, y: ceil.boundary_sphere, color: 4, markers: true, dash: 'dash' },
      { name: 'Last-scattering ball (galaxy, 21-cm volume)', x: ceil.d_over_D, y: ceil.lss_ball, color: 5, markers: true },
      { name: 'CMB ℓ ≤ 30 (Theorem P)', x: ceil.d_over_D, y: ceil.cmb_l30, color: SHAPE_SLOT.cubic, markers: true },
    ];
    keep(
      ctx,
      lineChart(b.body, series, {
        xTitle: 'torus size d / D',
        yTitle: 'KL divergence (nats)',
        yType: 'log',
        height: 330,
        refLines: [{ axis: 'y', value: 1 }],
        layout: { annotations: [logLineLabel(1, 'KL = 1 nat')] as never, yaxis: { dtick: 1 } as never },
        ariaLabel: 'Information ceiling of every probe against torus size',
        format: (v) => fmtNum(v, 3),
      }),
    );
    b.body.append(
      h(
        'p',
        { class: 'm1-takeaway' },
        'Theorem L: with a smooth spectrum and Einstein causality no observation of any kind exceeds these values. Strong evidence is possible only for d ≲ 1.1 D, a sliver above the present exclusion limits. ',
        badge('quoted', ceil.source ?? '[TOP] Sec. 11.2'),
      ),
    );
  }
}

// --- (b) resonance KL table ----------------------------------------------------------------------------

function mountResonanceTable(el: HTMLElement, ctx: PanelContext, R: ResonanceBundle): void {
  const rows = R.kl_table;
  const sizes = [...new Set(rows.map((r) => r.d_over_D))].sort((a, b) => a - b);
  const combos = [...new Map(rows.map((r) => [`${r.shell}_${r.f}`, { shell: r.shell, n: r.n_vectors, f: r.f }])).values()].sort(
    (a, b) => a.shell - b.shell || a.f - b.f,
  );
  const find = (d: number, shell: number, f: number): ResonanceKlRow | undefined =>
    rows.find((r) => r.d_over_D === d && r.shell === shell && r.f === f);
  const byD = (d: number): ResonanceKlRow | undefined => rows.find((r) => r.d_over_D === d);

  const c = panel({
    title: 'Topological resonance · KL with a primordial line',
    subtitle: 'Cubic torus, known orientation, ℓ ≤ 30: a narrow line of fractional bump f on one dual-lattice shell; same C_ℓ under both hypotheses',
    badges: [badge('computed'), badge('quoted', '[TOP] Sec. 11.6')],
    span: 5,
    stretch: true,
  });
  const t = panel({
    title: 'Ours against the printed table',
    subtitle: `KL (nats). Columns “|m|² · f”: a line on shell |m|² (${combos
      .filter((k, i, a) => a.findIndex((x) => x.shell === k.shell) === i)
      .map((k) => `${k.shell}: ${k.n} vectors`)
      .join(', ')}) with bump f. Our value, the paper’s beneath; gain = best KL with a line ÷ KL without`,
    badges: [badge('computed'), badge('quoted', '[TOP] Sec. 11.6')],
    span: 7,
    stretch: true,
  });
  el.append(grid(c.root, t.root));

  const series: LineSeries[] = [
    { name: 'no line', x: sizes, y: sizes.map((d) => byD(d)?.no_line ?? null), color: 'context', markers: true },
    ...combos.map(
      (k): LineSeries => ({
        name: `|m|² = ${k.shell} (${k.n} vectors), f = ${k.f}`,
        x: sizes,
        y: sizes.map((d) => find(d, k.shell, k.f)?.kl ?? null),
        color: k.shell === 9 ? 1 : 4,
        dash: k.f < 0.4 ? 'dash' : 'solid',
        markers: true,
      }),
    ),
  ];
  keep(
    ctx,
    lineChart(c.body, series, {
      xTitle: 'torus size d / D',
      yTitle: 'KL divergence (nats)',
      yType: 'log',
      height: 320,
      refLines: [{ axis: 'y', value: 1 }],
      layout: { annotations: [logLineLabel(1, 'KL = 1 nat')] as never, yaxis: { dtick: 1 } as never },
      ariaLabel: 'KL divergence with and without a resonant line',
      format: (v) => fmtNum(v, 4),
    }),
  );
  const top = [...rows].sort((a, b) => b.gain - a.gain)[0];
  c.body.append(
    h(
      'p',
      { class: 'm1-takeaway' },
      top
        ? `At ${top.d_over_D} D a torus invisible with a smooth spectrum (KL = ${fmtFixed(top.no_line, 2)}) becomes clearly visible (KL = ${fmtFixed(
            Math.max(...rows.filter((r) => r.d_over_D === top.d_over_D).map((r) => r.kl)),
            2,
          )}) if a line sits on the sparse 6-vector shell: a ${fmtNum(top.gain, 3)}-fold gain, at known orientation only. `
        : '',
      badge('computed'),
    ),
  );

  const cols: RawColumn<number>[] = [
    { label: 'd / D', num: true, cell: (d) => h('strong', null, d.toFixed(2)) },
    { label: 'no line', num: true, cell: (d) => oursPaper(fmtFixed(byD(d)?.no_line, 2), paperText(byD(d)?.paper?.no_line)) },
    ...combos.map(
      (k): RawColumn<number> => ({
        label: `${k.shell} · ${k.f}`,
        title: `shell |m|² = ${k.shell} (${k.n} vectors), f = ${k.f}`,
        num: true,
        cell: (d) => {
          const r = find(d, k.shell, k.f);
          return r ? oursPaper(fmtFixed(r.kl, 2), paperText(r.paper?.kl)) : '—';
        },
      }),
    ),
    {
      label: 'gain',
      num: true,
      cell: (d) => {
        const r = byD(d);
        return r ? oursPaper(`${fmtNum(r.gain, 3)}×`, r.paper?.gain !== undefined ? `${paperText(r.paper.gain, 0)}×` : null) : '—';
      },
    },
  ];
  t.body.append(h('div', { class: 'res-pivot' }, rawTable(cols, sizes).root));
  t.body.append(
    h(
      'p',
      { class: 'm1-takeaway' },
      'Every KL agrees with the printed value to the printed digits (largest relative difference ≈ 0.5 %). The paper prints the gain rounded down to an integer (3×, 16×, 28×, 64×); the unrounded ratios are 3.35, 16.5, 28.1 and 65.0. ',
      badge('computed', 'validation.json'),
    ),
    note(
      [
        'Theorem V(c): a spectral line breaks the dichotomy. On a torus, power at |k| = k₀ exists only on the dual-lattice vectors of that shell: a finite set of plane waves with fixed directions, coherent across all multipoles. In infinite space the same power is spread over all directions.',
        'The addition theorem makes the C_ℓ of the two hypotheses identical, so the divergence measures purely topological structure. ([TOP] Sec. 11.6)',
        'The Planck spectrum itself allows such lines: its C_ℓ power caps an isotropic line at f ≈ 0.3–0.9 (95 %), depending on k₀. ([TOP] Sec. 11.6)',
      ],
      'Why a line helps',
    ),
  );
}

// --- shell explorer -----------------------------------------------------------------------------------

function mountShellExplorer(el: HTMLElement, ctx: PanelContext, R: ResonanceBundle): void {
  const keys = Object.keys(R.shells);
  if (!keys.length) return;
  const shapes = SHAPES.filter((s) => keys.some((k) => k.startsWith(`${s}_`)));
  const sizesOf = (s: string): string[] =>
    keys
      .filter((k) => k.startsWith(`${s}_`))
      .map((k) => k.split('_')[1] ?? '')
      .sort((a, b) => Number(a) - Number(b));
  const state: { shape: string; size: string; idx: number } = { shape: shapes.includes('cubic') ? 'cubic' : (shapes[0] ?? 'cubic'), size: '1.2', idx: 0 };
  if (!sizesOf(state.shape).includes(state.size)) state.size = sizesOf(state.shape)[0] ?? '';

  const p = panel({
    title: 'Shell explorer · the dual-lattice shells of the search grid',
    subtitle: `Every shell with at most 12 vectors and 2 ≤ k₀χ* ≤ 14 (${R.n_shells ?? '—'} shells in all). Directions are recomputed live from the lattice; click a row.`,
    badges: [badge('computed', 'resonance.json shells'), badge('live', 'resonance-logic.ts')],
    span: 12,
  });
  el.append(grid(p.root));

  const shapeSeg = segmented({
    label: 'Torus shape',
    value: state.shape,
    options: shapes.map((s) => ({ value: s, label: SHAPE_LABEL[s] ?? s })),
    onChange: (v) => {
      state.shape = v;
      if (!sizesOf(v).includes(state.size)) state.size = sizesOf(v)[0] ?? '';
      state.idx = 0;
      renderSizes();
      render();
    },
  });
  const sizeHost = h('span');
  const renderSizes = (): void => {
    clear(sizeHost);
    sizeHost.append(
      segmented({
        label: 'Torus size',
        value: state.size,
        options: sizesOf(state.shape).map((d) => ({ value: d, label: `${d} D` })),
        onChange: (v) => {
          state.size = v;
          state.idx = 0;
          render();
        },
      }).root,
    );
  };
  renderSizes();

  const left = h('div');
  const right = h('div');
  p.body.append(
    h('div', { class: 'm1-toolbar' }, h('span', { class: 'control-label' }, 'Torus'), shapeSeg.root, sizeHost),
    h('div', { class: 'res-shell-layout' }, left, right),
  );

  let chart: { destroy(): void } | null = null;
  const models = R.search?.models ?? [];
  const fGrid = R.single_field?.f_max_grid ?? [];

  function render(): void {
    const key = `${state.shape}_${state.size}`;
    const shells = R.shells[key] ?? [];
    const d = Number(state.size);
    const lnB = (sh: ResonanceShell, f: number): number | null =>
      models.find((m) => m.shape === state.shape && Math.abs(m.d_over_D - d) < 1e-9 && Math.abs(m.shell_m2 - sh.m2) < 1e-9 && Math.abs(m.f - f) < 1e-9)?.lnB ??
      null;
    const fmax = (sh: ResonanceShell): number | null =>
      fGrid.find((g) => g.shape === state.shape && Math.abs(g.d_over_D - d) < 1e-9 && Math.abs(g.m2 - sh.m2) < 1e-9)?.f_max ?? null;
    const fs = [...new Set(models.map((m) => m.f))].sort((a, b) => a - b);
    const vecs = shells.map((sh) => shellVectors(state.shape, sh.m2));

    const cols: RawColumn<number>[] = [
      { label: '#', num: true, title: 'shell index in the dual lattice (the paper counts all shells)', cell: (i) => String(shells[i]?.index ?? i + 1) },
      { label: '|m|²', num: true, cell: (i) => fmtNum(shells[i]!.m2, 4) },
      { label: 'k₀χ*', num: true, title: 'k₀χ* = π|m| / (d/D)', cell: (i) => fmtFixed(shells[i]!.k_chi, 2) },
      { label: 'N', num: true, title: 'number of dual-lattice vectors (rank of the line term)', cell: (i) => String(shells[i]!.n_vectors) },
      {
        label: 'directions',
        cell: (i) =>
          h(
            'span',
            { class: 'res-family' },
            familyCounts(vecs[i] ?? [])
              .map(([f, n]) => `${n}⟨${f}⟩`)
              .join(' + '),
          ),
      },
      ...fs.map(
        (f): RawColumn<number> => ({
          label: `ln B · f ${f}`,
          num: true,
          title: `orientation-averaged ln B on ${R.search?.map ?? 'the build map'} (ours)`,
          cell: (i) => {
            const v = lnB(shells[i]!, f);
            return v === null ? '—' : sign(v);
          },
        }),
      ),
      {
        label: 'f_max 1-field',
        num: true,
        title: 'largest bump single-field inflation allows on this shell ([TOP] Sec. 11.8 growth bound), computed',
        cell: (i) => {
          const v = fmax(shells[i]!);
          return v === null ? '—' : h('span', { class: v < 0.05 ? 'res-dim' : null }, fmtNum(v, 2));
        },
      },
    ];
    const tbl = rawTable(cols, shells.map((_, i) => i), (i) => (i === state.idx ? 'res-row-selected' : null));
    tbl.root.classList.add('res-clickable');
    tbl.tbody.addEventListener('click', (e) => {
      const tr = (e.target as Element).closest('tr');
      const i = tr ? [...tbl.tbody.children].indexOf(tr) : -1;
      if (i >= 0 && i !== state.idx) {
        state.idx = i;
        render();
      }
    });
    const shellSel = select({
      label: 'Shell',
      value: String(state.idx),
      options: shells.map((sh, i) => ({ value: String(i), label: `|m|² = ${fmtNum(sh.m2, 4)} · ${sh.n_vectors} vectors` })),
      onChange: (v) => {
        state.idx = Number(v);
        render();
      },
    });
    clear(left);
    left.append(
      tbl.root,
      h(
        'p',
        { class: 'm1-takeaway' },
        `ln B: the orientation-averaged Bayes factor of torus + line against infinite space + the same isotropic line, from our search on ${R.search?.map ?? 'the build map'}; none is significant. Grey f_max < 0.05: single-field inflation leaves no room for a line there. `,
        badge('computed'),
      ),
    );

    const sh = shells[state.idx];
    const vs = vecs[state.idx] ?? [];
    chart?.destroy();
    clear(right);
    if (!sh) return;
    const pts = vs.map((v) => azEl(v));
    right.append(h('div', { class: 'm1-toolbar' }, shellSel.root));
    chart = keep(
      ctx,
      scatterChart(
        right,
        [
          {
            name: `${keyLabel(key)} · |m|² = ${fmtNum(sh.m2, 4)}`,
            x: pts.map((q) => q[0]),
            y: pts.map((q) => q[1]),
            color: SHAPE_SLOT[state.shape],
            size: 10,
            text: vs.map((v) => `⟨${directionFamily(v)}⟩ m = (${v.map((c) => fmtFixed(c, 2)).join(', ')})`),
          },
        ],
        {
          xTitle: 'azimuth in the cube frame (deg)',
          yTitle: 'elevation (deg)',
          xRange: [-190, 190],
          yRange: [-95, 95],
          height: 270,
          layout: { xaxis: { dtick: 45 } as never, yaxis: { dtick: 45 } as never },
          ariaLabel: 'Directions of the shell vectors',
          format: (v) => fmtNum(v, 4),
        },
      ),
    );
    right.append(
      kv([
        ['Wave number', [`k₀χ* = ${fmtFixed(sh.k_chi, 3)} `, badge('computed'), h('span', { class: 'small muted' }, ` live π|m|/(d/D) = ${fmtFixed(shellKChi(sh.m2, d), 3)}`)]],
        ['Plane waves', `${vs.length} (bundle: ${sh.n_vectors}) — ${familyCounts(vs).map(([f, n]) => `${n} along ⟨${f}⟩`).join(', ')}`],
        ['Rank of the line term', `N = ${sh.n_vectors}: C = C₀ + U W Uᵀ with U of size 957 × ${sh.n_vectors}`],
      ]),
    );
  }
  render();
}

// --- (d) the pre-registered search -----------------------------------------------------------------------

function mountSearch(el: HTMLElement, ctx: PanelContext, R: ResonanceBundle): void {
  const s = R.search;
  const sp = s?.paper ?? R.search_paper ?? null;
  const head = panel({
    title: 'Pre-registered resonance search',
    subtitle: sp
      ? `FCC, BCC and cubic tori at 1.05–1.3 D; every shell with ≤ 12 vectors and 2 ≤ k₀χ* ≤ 14; f = 0.2 and 0.4: ${s?.n_models ?? 160} models. S_max = max ln B over the models, S_avg = ln of their mean B; calibrated on ${s?.n_null ?? sp.n_null ?? 100} isotropic skies.`
      : undefined,
    badges: [badge('measured', sp?.map ?? 'Planck'), badge('computed', s?.map ?? 'build map')],
    span: 12,
  });
  el.append(grid(head.root));
  if (!s) {
    head.body.append(h('p', { class: 'muted' }, 'Our re-run of the search is not in this build (resonance.json search = null); only the paper’s Planck values are shown below.'));
  }
  const pb = sp?.best_model;
  const ob = s?.best;
  head.body.append(
    h(
      'div',
      { class: 'res-tile-group' },
      sectionLabel(`S_max and S_avg · ${s?.map ?? 'build map'} (ours) against Planck (paper)`),
      tilesN(
        4,
        statTile({ label: `S_max · ${s?.map ?? '—'}`, value: s ? fmtFixed(s.S_max, 2) : null, sub: s ? [badge('computed'), `p = ${fmtFixed(s.p_max, 2)}`] : undefined }),
        statTile({ label: 'S_max · Planck SMICA', value: fmtFixed(sp?.S_max, 2), sub: [badge('quoted', '[TOP] Sec. 11.6'), `p = ${fmtFixed(sp?.p_max, 2)}`] }),
        statTile({ label: `S_avg · ${s?.map ?? '—'}`, value: s ? sign(s.S_avg, 2) : null, sub: s ? [badge('computed'), `p = ${fmtFixed(s.p_avg, 2)}`] : undefined }),
        statTile({ label: 'S_avg · Planck SMICA', value: sp ? sign(sp.S_avg, 2) : null, sub: [badge('quoted', '[TOP] Sec. 11.6'), `p = ${fmtFixed(sp?.p_avg, 2)}`] }),
      ),
    ),
    h(
      'div',
      { class: 'res-tile-group' },
      sectionLabel('Best single model (before the look-elsewhere correction)'),
      tilesN(
        2,
        statTile({
          label: `${s?.map ?? '—'} (ours)`,
          value: ob ? `${ob.shape} ${ob.d_over_D} D · ln B = ${sign(ob.lnB)}` : null,
          mono: true,
          sub: ob ? [badge('computed'), `k₀χ* = ${fmtFixed(ob.k_chi, 2)}, N = ${ob.n_vectors}, f = ${ob.f}, local p = ${fmtFixed(ob.p_local, 2)}`, statusPill('warning', 'not significant')] : undefined,
        }),
        statTile({
          label: 'Planck SMICA (paper)',
          value: pb ? `${pb.shape} ${pb.d_over_D} D · ln B = ${sign(pb.lnB)}` : null,
          mono: true,
          sub: pb ? [badge('quoted', '[TOP] Sec. 11.6'), `k₀χ* = ${fmtFixed(pb.k0_chi, 2)}, N = ${pb.n_vectors}, local p = ${fmtFixed(pb.p_local, 2)}`, statusPill('warning', 'not significant')] : undefined,
        }),
      ),
    ),
  );
  if (!s) return;

  // null histograms with the observed values
  const two = h('div', { class: 'm1-two' });
  const hA = h('div');
  const hB = h('div');
  two.append(hA, hB);
  head.body.append(two);
  const hist = (vals: number[], step: number): ReturnType<typeof histogram> => {
    const lo = Math.floor(Math.min(...vals) / step) * step;
    const hi = Math.ceil(Math.max(...vals) / step) * step;
    return histogram(vals, Math.max(1, Math.round((hi - lo) / step)), lo, hi);
  };
  const hmax = hist(s.null.S_max, 0.25);
  const havg = hist(s.null.S_avg, 0.1);
  keep(
    ctx,
    barChart(hA, [{ name: 'isotropic ΛCDM skies', x: hmax.x, y: hmax.counts, color: 1 }], {
      xTitle: 'S_max = max ln B over 160 models',
      yTitle: `skies (of ${s.null.S_max.length})`,
      height: 260,
      refLines: [{ axis: 'x', value: s.S_max, label: `${s.map} ${fmtFixed(s.S_max, 2)} · p = ${fmtFixed(s.p_max, 2)}` }],
      layout: { bargap: 0.08 },
      ariaLabel: 'Null distribution of S_max with the observed value',
      format: (v) => fmtNum(v, 3),
    }),
  );
  keep(
    ctx,
    barChart(hB, [{ name: 'isotropic ΛCDM skies', x: havg.x, y: havg.counts, color: 1 }], {
      xTitle: 'S_avg = ln (mean B over 160 models)',
      yTitle: `skies (of ${s.null.S_avg.length})`,
      height: 260,
      refLines: [{ axis: 'x', value: s.S_avg, label: `${s.map} ${sign(s.S_avg)} · p = ${fmtFixed(s.p_avg, 2)}` }],
      layout: { bargap: 0.08 },
      ariaLabel: 'Null distribution of S_avg with the observed value',
      format: (v) => fmtNum(v, 3),
    }),
  );
  head.body.append(
    h(
      'p',
      { class: 'm1-takeaway' },
      `Bars: our ${s.null.S_max.length} isotropic skies processed by the identical code and grid (bin centres). The vertical line is our ${s.map} value. The paper’s Planck values (S_max = ${fmtFixed(sp?.S_max, 2)}, S_avg = ${sign(sp?.S_avg ?? 0)}) were calibrated on the paper’s own null skies, with null median ${fmtFixed(sp?.null_median_S_max, 2)}; ours is ${fmtFixed(s.null_median_S_max ?? median(s.null.S_max), 2)}. `,
      badge('computed'),
    ),
  );

  // every model's ln B
  const mHost = h('div');
  head.body.append(sectionLabel(`ln B of every model on ${s.map}`), mHost);
  const series: ScatterSeries[] = SHAPES.map((shape) => {
    const ms = s.models.filter((m) => m.shape === shape);
    return {
      name: SHAPE_LABEL[shape] ?? shape,
      x: ms.map((m) => m.k_chi ?? shellKChi(m.shell_m2, m.d_over_D)),
      y: ms.map((m) => m.lnB),
      color: SHAPE_SLOT[shape],
      symbol: 'circle',
      size: 8,
      text: ms.map((m) => `${shape} ${m.d_over_D} D · |m|² ${fmtNum(m.shell_m2, 4)} (N ${m.n_vectors ?? '—'}) · f ${m.f}`),
    };
  });
  keep(
    ctx,
    scatterChart(mHost, series, {
      xTitle: 'shell wave number k₀χ*',
      yTitle: 'ln B (torus + line vs isotropic line)',
      height: 280,
      refLines: [{ axis: 'y', value: 0 }],
      ariaLabel: 'Bayes factor of every model of the search',
      format: (v) => fmtNum(v, 3),
    }),
  );
  const gc = s.grid_check;
  head.body.append(
    note(
      [
        h(
          'ul',
          { class: 'res-notes' },
          h('li', null, `Statistic: ${s.statistic ?? 'ln B per model, S_max and S_avg over the models'}.`),
          s.grid ? h('li', null, `SO(3) grid ${s.grid.n_alpha} × ${s.grid.n_beta} × ${s.grid.n_gamma} (α, β, γ); the paper’s search grid is 36 × 19 × 72.`) : null,
          gc
            ? h(
                'li',
                null,
                `Grid check on the same data: S_max = ${gc.S_max.map((v, i) => `${fmtFixed(v, 2)} (${(gc.grids[i] ?? []).join('×')})`).join(', ')}. The coarse 36 × 19 × 72 grid under-resolves the narrow orientation peaks, so we use the finer grid for every number and calibration.`,
              )
            : null,
          sp?.power
            ? h(
                'li',
                null,
                `Power of the test (paper, ${sp.power.skies ?? 20} skies, 5 % level): cubic 1.2 D, f = 0.4 on shell 7: ${Math.round(100 * (sp.power['cubic_1.2D_f0.4_shell7'] ?? NaN))} %; FCC 1.2 D, f = 0.4 on shell 5: ${Math.round(100 * (sp.power['FCC_1.2D_f0.4_shell5'] ?? NaN))} %. The search is not blind to resonance, but far from decisive. `,
                badge('quoted', '[TOP] Sec. 11.6'),
              )
            : null,
          ...(s.notes ?? []).map((n) => h('li', null, n)),
        ),
      ],
      'How the search is computed',
    ),
  );
}

// --- (c) Woodbury -------------------------------------------------------------------------------------------

function mountWoodbury(el: HTMLElement, ctx: PanelContext, R: ResonanceBundle): void {
  const w = R.woodbury_benchmark;
  const p = panel({
    title: 'Rank-N Woodbury update',
    subtitle: `ln N(a; 0, C₀ + U W Uᵀ) for n = ${w.n} multipoles (2 ≤ ℓ ≤ 30), ${w.shape ?? 'cubic'} torus at ${w.d_over_D ?? '—'} D, mean over ${w.n_orientations ?? '—'} random orientations`,
    badges: [badge('computed', 'qchas.cmb.resonance.woodbury_benchmark')],
    span: 5,
    stretch: true,
  });
  const q = panel({
    title: 'Single-field limits · exit 1 examined',
    subtitle: 'The largest bump single-field inflation allows on the lowest shells, and what it buys at known orientation',
    badges: [badge('computed'), badge('quoted', '[TOP] Sec. 11.8')],
    span: 7,
    stretch: true,
  });
  el.append(grid(p.root, q.root));

  const labels = w.ranks.map((r) => `N = ${r}`);
  const speed = w.ranks.map((_, i) => (w.direct_ms[i] ?? NaN) / (w.woodbury_ms[i] ?? NaN));
  p.body.append(
    tilesN(
      2,
      statTile({
        label: 'Speed-up per orientation',
        value: `${fmtNum(Math.min(...speed), 3)}–${fmtNum(Math.max(...speed), 3)}×`,
        sub: [badge('computed'), `${fmtNum(w.direct_ms[0], 3)} → ${fmtNum(w.woodbury_ms[0], 3)} ms at N = ${w.ranks[0]}`],
      }),
      statTile({ label: 'max |Δ ln L| vs dense Cholesky', value: fmtNum(w.max_abs_diff, 2), unit: 'nats', sub: [badge('computed'), 'machine precision'] }),
    ),
  );
  keep(
    ctx,
    barChart(
      p.body,
      [
        { name: 'dense Cholesky of C₀ + UWUᵀ', x: labels, y: w.direct_ms, color: 2 },
        { name: 'Woodbury + determinant lemma', x: labels, y: w.woodbury_ms, color: 1 },
      ],
      {
        xTitle: 'rank of the line term (shell vectors)',
        yTitle: 'time per orientation (ms)',
        height: 250,
        layout: { xaxis: { type: 'category' } as never },
        ariaLabel: 'Direct against Woodbury time',
        format: (v) => fmtNum(v, 3),
      },
    ),
  );
  p.body.append(
    h(
      'p',
      { class: 'res-explain' },
      `A line on a shell of N vectors adds only a rank-N term U W Uᵀ to the ${w.n} × ${w.n} torus covariance C₀, so the Woodbury identity and the matrix-determinant lemma give the likelihood from one factorization of C₀ plus an N × N solve, O(n²N) instead of a fresh O(n³) Cholesky for every orientation. `,
      'That is what makes the orientation-averaged search over SO(3) and 160 models feasible, and it reproduces the dense computation to machine precision, as the paper reports for its own run. ',
      badge('computed'),
      badge('quoted', '[TOP] Sec. 11.6'),
    ),
  );

  // single-field rows
  const sf = R.single_field;
  if (!sf) {
    q.body.append(h('p', { class: 'muted small' }, 'Not in this build.'));
    return;
  }
  const gb = sf.paper?.growth_bound;
  q.body.append(
    rawTable<SingleFieldRow>(
      [
        {
          label: 'Torus · shell',
          cell: (r) => h('span', { class: 'res-nowrap' }, h('strong', null, `${r.torus} ${r.d_over_D} D`), h('br'), h('span', { class: 'small muted' }, `k₀χ* ${fmtFixed(r.k0_chi, 2)} · N ${r.n_vectors}`)),
        },
        { label: 'max f, 1-field', num: true, cell: (r) => oursPaper(fmtFixed(r.single_field_max_f, 2), paperText(r.paper?.single_field_max_f)) },
        { label: 'Planck cap', num: true, title: 'quoted Planck cap on f', cell: (r) => fmtFixed(r.planck_cap, 2) },
        { label: 'KL no line', num: true, cell: (r) => oursPaper(fmtFixed(r.kl_no_line, 2), paperText(r.paper?.kl_no_line)) },
        { label: 'KL with line', num: true, cell: (r) => oursPaper(fmtFixed(r.kl_allowed_line, 2), paperText(r.paper?.kl_allowed_line)) },
        { label: 'gain', num: true, cell: (r) => oursPaper(`${fmtNum(r.gain, 3)}×`, r.paper?.gain !== undefined ? `${paperText(r.paper.gain, 1)}×` : null) },
      ],
      sf.rows,
    ).root,
    h(
      'p',
      { class: 'm1-takeaway' },
      gb
        ? `Growth bound d ln P / d ln k ≤ ${gb.dlnP_dlnk_max} (${gb.refined} in refined analyses): a line reaching its peak within a half-width has strength a ≤ ${gb.a_max} and bump f = ${gb.f}. `
        : '',
      `On ${sf['n_f_below_0.05']} of ${sf.n_shells} shells this leaves f < 0.05 and the resonance gains nothing; on the lowest shells a single-field bump gains at most ${fmtNum(Math.max(...sf.rows.map((r) => r.gain)), 3)}×. The 64-fold gains need narrow lines on sparse high shells, which only multi-field dynamics can produce. `,
      badge('computed'),
      badge('quoted', '[TOP] Sec. 11.8'),
    ),
  );
}

// --- (e) the 3D exits ---------------------------------------------------------------------------------------

function mountExits(el: HTMLElement, R: ResonanceBundle): void {
  const p = panel({
    title: 'Exit 1 in three dimensions · a perfect map of the causal ball',
    subtitle: 'KL divergence (nats), cubic torus, lines on six-vector shells, known orientation; five spheres out to the causal horizon',
    badges: [badge('quoted', '[TOP] Sec. 11.8')],
    span: 12,
  });
  el.append(grid(p.root));
  const cols = ['d1.5', 'd2', 'd3'].filter((k) => R.exits_3d.some((r) => k in r));
  const cell = (v: number | string | null | undefined): HTMLElement => {
    if (typeof v !== 'number') return h('span', { class: 'res-dim', title: 'the paper prints a dash' }, '—');
    // quoted values, shown exactly as printed
    return h('span', { class: v >= 1 ? 'res-strong' : null }, String(v));
  };
  p.body.append(
    rawTable<ResonanceBundle['exits_3d'][number]>(
      [
        { label: 'Primordial spectrum', cell: (r) => pretty(r.spectrum), className: 'wrap' },
        ...cols.map((k): RawColumn<ResonanceBundle['exits_3d'][number]> => ({ label: `d = ${k.slice(1)} D`, num: true, cell: (r) => cell(r[k]) })),
      ],
      R.exits_3d,
      (r) => (/multi-field/.test(r.spectrum) ? 'm1-row-survivor' : null),
    ).root,
    h(
      'p',
      { class: 'm1-takeaway' },
      'Bold: KL ≥ 1 nat. Only a multi-field narrow line, mapped in 3D (for example by dark-ages 21-cm tomography), carries enough information to decide a torus three times larger than the observable universe, where a smooth spectrum gives 0.0003 nats. This is conditional on nature: it needs a narrow primordial line on a lattice shell, and current data show none. ',
      badge('quoted', '[TOP] Sec. 11.8'),
    ),
  );
}
