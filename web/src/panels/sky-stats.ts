/**
 * 1.5 Sky Statistics ([TOP] Secs. 8-9, 14).
 *
 * The large-angle sky against LambdaCDM and torus models: D_l for 2 <= l <= 30 (published WMAP
 * spectrum with errors, Planck when the build has it, our own a_lm extraction, CAMB LambdaCDM
 * and the Sec. 8 torus models), the anomaly statistics D_2 and S_1/2 with LambdaCDM p-values
 * beside the paper's Planck numbers, C(theta), the quadrupole / octopole axes on a Mollweide,
 * and the orientation-scanned torus likelihood (our independent re-analysis of the build map)
 * next to [TOP] Tables 9a/9b.
 *
 * Bundle fields: cmb.{ell, observed, map_D, lcdm_D, torus_D, stats.{D2, S12, p_C2, p_S12,
 * labels, C_theta, S12_lcdm_hist, S12_lcdm_quantiles, axes, paper, paper_check}, sec8,
 * likelihood, likelihood_meta}; validation rows of the Sec. 8 groups.
 */

import './m1.css';
import './resonance.css';
import { type LineSeries, type ScatterSeries, barChart, lineChart, scatterChart } from '../lib/charts';
import type { CmbBundle, LBPair, LikelihoodRow, Sec8Row } from '../lib/data';
import type { PanelContext } from '../lib/panel';
import type { SkyView } from '../lib/sky';
import { append, badge, clear, fmtFixed, fmtInt, fmtNum, grid, h, loadFailure, note, panel, segmented, statTile, statusPill } from '../lib/ui';
import { type RawColumn, SHAPE_SLOT, keep, oursPaper, proofList, rawTable, sectionLabel, tilesN, validationRows, validationSummary } from './m1-ui';
import { antipode, axisAngleDeg, isotropicAlignmentP, lbVec } from './resonance-logic';

export const title = 'Sky Statistics';
export const subtitle =
  'The large-angle CMB sky against ΛCDM and torus models: low-ℓ spectra, the quadrupole and S₁/₂ anomalies, aligned axes, and the orientation-scanned torus likelihood on a real map.';

const VAL_GROUPS = ['Large-angle CMB ([TOP] Sec. 8)', 'Low-l anomalies ([TOP] Sec. 8)'];
const SHAPE_LABEL: Record<string, string> = { cubic: 'Cubic', BCC: 'BCC', FCC: 'FCC' };
/** [TOP] Table 9b: 95 % lower bounds on the shortest loop d / D (quoted). */
const LOWER_BOUND: Record<string, number> = { FCC: 1.04, BCC: 1.0, cubic: 0.985 };

const sign = (v: number | null | undefined, d = 2): string =>
  v === null || v === undefined || !Number.isFinite(v) ? '—' : (v > 0 ? '+' : v < 0 ? '−' : '') + Math.abs(v).toFixed(d);
const fmtP = (v: number | null | undefined): string => (v === null || v === undefined ? '—' : v < 0.01 ? fmtNum(v, 2) : fmtFixed(v, 3));
const lbText = (lb: readonly number[]): string => `(${fmtFixed(lb[0], 1)}°, ${sign(lb[1] ?? NaN, 1)}°)`;
/** "FCC_0.985" -> ['FCC', 0.985] */
const parseKey = (k: string): [string, number] => {
  const [s, d] = k.split('_');
  return [s ?? k, Number(d)];
};

interface Verdict {
  excluded: boolean;
  label: string;
  src: string;
}

/** Table 9b verdict of a torus model, or the quoted 95 % lower bound when the table has no row. */
function verdictOf(C: CmbBundle, shape: string, d: number): Verdict {
  const row = C.likelihood?.find((r) => r.shape === shape && Math.abs(r.d_over_D - d) < 1e-6);
  if (row?.paper?.verdict) return { excluded: row.paper.verdict.startsWith('excluded'), label: row.paper.verdict, src: '[TOP] Table 9b' };
  const lb = LOWER_BOUND[shape];
  if (lb === undefined) return { excluded: false, label: 'no bound quoted', src: '[TOP] Table 9b' };
  return d < lb
    ? { excluded: true, label: `excluded: below the 95 % bound d ≥ ${lb} D`, src: '[TOP] Table 9b' }
    : { excluded: false, label: `allowed: above the 95 % bound ${lb} D`, src: '[TOP] Table 9b' };
}

const verdictPill = (v: Verdict, short = false): HTMLElement => {
  const pill = statusPill(v.excluded ? 'critical' : 'good', short ? (v.excluded ? 'excluded' : 'allowed') : v.label);
  pill.title = `${v.label} (${v.src})`;
  return pill;
};

/** Section label whose units / formula keep their case (the label style is uppercase). */
const unitLabel = (text: string, tail: string): HTMLElement => {
  const el = sectionLabel(text);
  el.append(' ', h('span', { style: 'text-transform: none; letter-spacing: 0.02em' }, tail));
  return el;
};

export async function mount(el: HTMLElement, ctx: PanelContext): Promise<void> {
  ctx.setHeaderExtras([badge('measured', 'WMAP 9-yr'), badge('computed', 'qchas.cmb'), badge('quoted', '[TOP] Secs. 8–9')]);
  const [cmb, val] = await Promise.all([ctx.load('cmb.json'), ctx.load('validation.json')]);
  if (ctx.unmounted) return;
  if (!cmb.ok) {
    const p = panel({ title: 'Sky Statistics', span: 12 });
    p.body.append(loadFailure(cmb));
    el.append(grid(p.root));
    return;
  }
  const C = cmb.data;
  mountAnomalies(el, C, validationSummary(validationRows(val, VAL_GROUPS), VAL_GROUPS));
  mountSpectrum(el, ctx, C);
  mountCorrelation(el, ctx, C);
  await mountAxes(el, ctx, C);
  if (ctx.unmounted) return;
  mountLikelihood(el, C);
}

// --- anomaly tiles ----------------------------------------------------------------------------------------

function mountAnomalies(el: HTMLElement, C: CmbBundle, validation: HTMLElement | null): void {
  const S = C.stats;
  const P = S.paper ?? {};
  const chk = S.paper_check ?? {};
  const ax = S.axes;
  const n = S.n_realizations;
  const p = panel({
    title: 'Large-angle anomalies against ΛCDM',
    subtitle: `One-sided p = P(statistic ≤ observed) under ΛCDM (CAMB, Planck 2018), ${n ? fmtInt(n) : '—'} Gaussian realizations; full sky, 2 ≤ ℓ ≤ 30`,
    badges: [badge('measured', 'WMAP 9-yr'), badge('computed'), badge('quoted', '[TOP] Sec. 8')],
    span: 12,
  });
  const planckTile = (label: string, value: string, sub: (Node | string)[]): HTMLElement =>
    statTile({ label, value, sub: [badge('quoted', P.source?.split(' and ')[0] ?? '[TOP] Sec. 8'), ...sub] });
  const mapLabel = C.source ?? 'map';

  const gD2 = h(
    'div',
    { class: 'res-tile-group' },
    unitLabel('Quadrupole', 'D₂ (μK²)'),
    tilesN(
      4,
      statTile({ label: 'WMAP 9-yr published spectrum', value: fmtNum(S.D2?.wmap, 4), sub: [badge('measured', 'Hinshaw et al. 2013'), `p = ${fmtP(S.p_C2.wmap)}`] }),
      statTile({ label: `Our extraction from ${mapLabel}`, value: fmtNum(S.D2?.map, 4), sub: [badge('computed'), `p = ${fmtP(S.p_C2.map)}`] }),
      planckTile('Planck 2018 (paper)', fmtNum(P.D2_planck, 4), [`p = ${fmtP(P.p_C2_planck)}`, h('span', null, `· our sims at this value: ${fmtP(chk.p_C2_planck_value)}`)]),
      statTile({ label: 'ΛCDM expectation', value: fmtNum(S.D2_lcdm ?? C.lcdm_D[0], 5), sub: [badge('computed', 'CAMB'), `paper ${fmtNum(P.D2_lcdm, 4)}`] }),
    ),
  );
  const gS = h(
    'div',
    { class: 'res-tile-group' },
    unitLabel('Missing large-angle correlation', 'S₁/₂ = ∫ C(θ)² d cos θ over θ ≥ 60° (μK⁴)'),
    tilesN(
      4,
      statTile({ label: 'WMAP 9-yr published spectrum', value: fmtInt(S.S12?.wmap), sub: [badge('computed', 'from the published D_ℓ'), `p = ${fmtP(S.p_S12.wmap)}`] }),
      statTile({ label: `Our extraction from ${mapLabel}`, value: fmtInt(S.S12?.map), sub: [badge('computed'), `p = ${fmtP(S.p_S12.map)}`] }),
      planckTile('Planck 2018 (paper)', fmtInt(P.S12_planck), [`p = ${fmtP(P.p_S12_planck)}`, h('span', null, `· our sims at this value: ${fmtP(chk.p_S12_planck_value)}`)]),
      statTile({ label: 'ΛCDM median', value: fmtInt(S.S12_lcdm_median), sub: [badge('computed'), `paper ${fmtInt(P.S12_lcdm_median)}`] }),
    ),
  );
  const live = ax?.paper ? axisAngleDeg(ax.paper.quadrupole, ax.paper.octopole) : null;
  const gA = ax
    ? h(
        'div',
        { class: 'res-tile-group' },
        sectionLabel('Quadrupole–octopole alignment'),
        tilesN(
          3,
          statTile({ label: `Angle between the axes · ${ax.source}`, value: fmtFixed(ax.angle_deg, 1), unit: '°', sub: [badge('computed'), `isotropic p = ${fmtP(ax.p_isotropic)}`] }),
          statTile({
            label: `Angle · ${ax.paper?.map ?? 'Planck'} (paper)`,
            value: fmtFixed(ax.paper?.angle_deg, 1),
            unit: '°',
            sub: [badge('quoted', ax.paper?.source ?? '[TOP] Sec. 9.1'), `isotropic p = ${fmtP(ax.paper?.p_isotropic)}, ΛCDM sims ${fmtP(ax.paper?.p_lcdm_sims)}`],
          }),
          statTile({
            label: 'Recomputed from the printed axes',
            value: live === null ? null : fmtFixed(live, 1),
            unit: '°',
            sub: [badge('live', 'resonance-logic.ts'), live === null ? '' : `1 − cos θ = ${fmtP(isotropicAlignmentP(live))}`],
          }),
        ),
      )
    : null;
  const ps = [S.p_C2.wmap, S.p_C2.map, S.p_C2.planck, S.p_S12.wmap, S.p_S12.map, S.p_S12.planck, P.p_C2_planck, P.p_S12_planck].filter(
    (v): v is number => typeof v === 'number',
  );
  append(
    p.body,
    gD2,
    gS,
    gA,
    h(
      'p',
      { class: 'm1-takeaway' },
      ps.length
        ? `Across the published spectra, our extraction and the paper’s Planck numbers, the low quadrupole and the missing large-angle correlation have ΛCDM p-values between ${fmtP(Math.min(...ps))} and ${fmtP(Math.max(...ps))}: unusual, not decisive. `
        : '',
      'These statistics motivated small tori; the orientation-scanned likelihood below finds no evidence for any torus. ',
      badge('quoted', '[TOP] Secs. 8–9'),
    ),
    validation,
  );
  el.append(grid(p.root));
}

// --- D_l ---------------------------------------------------------------------------------------------------

function mountSpectrum(el: HTMLElement, ctx: PanelContext, C: CmbBundle): void {
  const a = panel({
    title: 'Low-ℓ power spectrum',
    subtitle: 'D_ℓ = ℓ(ℓ+1)C_ℓ / 2π; published spectra with their errors, our extraction from the map, ΛCDM and one torus model at a time',
    badges: [badge('measured', 'WMAP 9-yr TT'), badge('computed', 'CAMB, transfer functions')],
    span: 6,
    stretch: true,
  });
  const b = panel({
    title: 'Torus models · [TOP] Sec. 8 table',
    subtitle: 'Our recomputation, the paper’s value beneath; D₂ in μK², KL of the ℓ ≤ 30 CMB in nats; status from Table 9b',
    badges: [badge('computed'), badge('quoted', '[TOP] Sec. 8, Table 9b')],
    span: 6,
    stretch: true,
  });
  el.append(grid(a.root, b.root));

  const keys = Object.keys(C.torus_D);
  let model = keys.includes('FCC_0.70') ? 'FCC_0.70' : (keys[0] ?? 'none');
  const host = h('div');
  const statusHost = h('div', { class: 'res-model-note' });
  const seg = segmented({
    label: 'Torus model',
    value: model,
    options: [
      ...keys.map((k) => {
        const [s, d] = parseKey(k);
        return { value: k, label: `${s} ${d} D` };
      }),
      { value: 'none', label: 'none' },
    ],
    onChange: (v) => {
      model = v;
      draw();
    },
  });
  a.body.append(h('div', { class: 'm1-toolbar' }, h('span', { class: 'control-label' }, 'Torus model'), seg.root), host, statusHost);

  let chart: { destroy(): void } | null = null;
  const ell = C.ell;
  function draw(): void {
    chart?.destroy();
    clear(host);
    clear(statusHost);
    const series: ScatterSeries[] = [];
    const W = C.observed.wmap;
    if (W)
      series.push({
        name: W.label ?? 'WMAP 9-yr published',
        x: W.ell,
        y: W.D,
        ...(W.err_minus && W.err_plus ? { errorAsym: [W.err_minus, W.err_plus] as [number[], number[]] } : { error: W.err }),
        color: 1,
      });
    const Pl = C.observed.planck;
    if (Pl)
      series.push({
        name: Pl.label ?? 'Planck 2018 published',
        x: Pl.ell,
        y: Pl.D,
        ...(Pl.err_minus && Pl.err_plus ? { errorAsym: [Pl.err_minus, Pl.err_plus] as [number[], number[]] } : { error: Pl.err }),
        color: 5,
      });
    if (C.map_D) series.push({ name: `our extraction (${C.source ?? 'map'})`, x: ell, y: C.map_D, color: 4, symbol: 'diamond', size: 9 });
    series.push({ name: 'ΛCDM (CAMB, Planck 2018)', x: ell, y: C.lcdm_D, color: 'context', mode: 'lines+markers', size: 5 });
    if (model !== 'none' && C.torus_D[model]) {
      const [s, d] = parseKey(model);
      const v = verdictOf(C, s, d);
      series.push({
        name: `${SHAPE_LABEL[s] ?? s} torus ${d} D${v.excluded ? ' (excluded)' : ''}`,
        x: ell,
        y: C.torus_D[model]!,
        color: SHAPE_SLOT[s] ?? 2,
        mode: 'lines+markers',
        size: 5,
      });
      statusHost.append(h('span', null, `${SHAPE_LABEL[s] ?? s} ${d} D:`), verdictPill(v), badge('quoted', v.src));
    }
    chart = keep(
      ctx,
      scatterChart(host, series, {
        xTitle: 'multipole ℓ',
        yTitle: 'D_ℓ (μK²)',
        xType: 'log',
        xRange: [1.8, 32],
        height: 340,
        layout: { xaxis: { tickvals: [2, 3, 4, 5, 7, 10, 15, 20, 30] } as never },
        ariaLabel: 'Low-ell power spectrum',
        format: (v) => fmtNum(v, 4),
      }),
    );
  }
  draw();
  const ex = (C.stats.extraction_check as { paper?: { D_extracted?: number[]; D_published?: number[] } } | undefined)?.paper;
  a.body.append(
    h(
      'p',
      { class: 'm1-takeaway' },
      `Our full-sky extraction gives D₂ = ${fmtNum(C.map_D?.[0], 4)} μK² against the published maximum-likelihood ${fmtNum(C.observed.wmap?.D[0], 4)} μK² (the published estimator marginalizes over foregrounds). `,
      ex?.D_extracted && ex.D_published ? `The paper’s own Planck check: ${fmtNum(ex.D_extracted[0], 3)} against ${fmtNum(ex.D_published[0], 3)}. ` : '',
      badge('computed'),
      ex ? badge('quoted', '[TOP] Sec. 9') : null,
    ),
  );

  // Sec. 8 model table
  const rows = C.sec8?.rows ?? [];
  if (!rows.length) {
    b.body.append(h('p', { class: 'muted small' }, 'The Sec. 8 table is not in this build.'));
    return;
  }
  b.body.append(
    rawTable<Sec8Row>(
      [
        { label: 'Torus', cell: (r) => h('strong', { class: 'res-nowrap' }, `${r.shape} ${r.d_over_D} D`) },
        { label: 'V / D³', num: true, cell: (r) => oursPaper(fmtFixed(r.V_over_D3, 3), r.paper?.V_over_D3 !== undefined ? fmtFixed(r.paper.V_over_D3, 3) : null) },
        { label: 'D₂', num: true, cell: (r) => oursPaper(fmtFixed(r.D2, 0), r.paper?.D2 !== undefined ? fmtFixed(r.paper.D2, 0) : null) },
        { label: 'KL', num: true, cell: (r) => oursPaper(fmtNum(r.kl, 3), r.paper?.kl !== undefined ? fmtNum(r.paper.kl, 3) : null) },
        { label: 'Status', cell: (r) => verdictPill(verdictOf(C, r.shape, r.d_over_D), true) },
      ],
      rows,
      (r) => (verdictOf(C, r.shape, r.d_over_D).excluded ? 'm1-row-excl' : null),
    ).root,
    h(
      'p',
      { class: 'm1-takeaway' },
      `ΛCDM D₂ = ${fmtNum(C.sec8?.lcdm_D[0], 5)} μK² (paper 1016). A torus suppresses the quadrupole only when it is small, and the small ones are excluded: FCC needs d ≥ 1.04 D, BCC 1.00 D, cubic 0.985 D at 95 %. `,
      badge('quoted', '[TOP] Table 9b'),
    ),
  );
}

// --- C(theta) and the S1/2 distribution ------------------------------------------------------------------

function mountCorrelation(el: HTMLElement, ctx: PanelContext, C: CmbBundle): void {
  const S = C.stats;
  const ct = S.C_theta;
  const a = panel({
    title: 'Two-point correlation C(θ)',
    subtitle: 'C(θ) = Σ (2ℓ+1)/4π C_ℓ P_ℓ(cos θ), 2 ≤ ℓ ≤ 30',
    badges: [badge('computed')],
    span: 6,
    stretch: true,
  });
  const b = panel({
    title: 'S₁/₂ under ΛCDM',
    subtitle: `Distribution over ${S.n_realizations ? fmtInt(S.n_realizations) : '—'} Gaussian realizations of the CAMB spectrum`,
    badges: [badge('computed')],
    span: 6,
    stretch: true,
  });
  el.append(grid(a.root, b.root));
  if (!ct) a.body.append(h('p', { class: 'muted small' }, 'C(θ) is not in this build.'));
  else {
    const series: LineSeries[] = [{ name: 'ΛCDM', x: ct.theta_deg, y: ct.lcdm, color: 'context', dash: 'dash' }];
    if (ct.wmap) series.push({ name: 'WMAP 9-yr published spectrum', x: ct.theta_deg, y: ct.wmap, color: 1 });
    if (ct.planck) series.push({ name: 'Planck 2018 published spectrum', x: ct.theta_deg, y: ct.planck, color: 5 });
    if (ct.map) series.push({ name: `our extraction (${C.source ?? 'map'})`, x: ct.theta_deg, y: ct.map, color: 4 });
    keep(
      ctx,
      lineChart(a.body, series, {
        xTitle: 'angular separation θ (deg)',
        yTitle: 'C(θ) (μK²)',
        xRange: [0, 180],
        height: 300,
        refLines: [
          { axis: 'x', value: 60, label: 'S₁/₂ integrates θ ≥ 60°' },
          { axis: 'y', value: 0 },
        ],
        layout: { xaxis: { dtick: 30 } as never },
        ariaLabel: 'Two-point correlation function',
        format: (v) => fmtNum(v, 4),
      }),
    );
    a.body.append(
      h(
        'p',
        { class: 'm1-takeaway' },
        'Beyond 60° the observed C(θ) stays close to zero, while ΛCDM predicts correlations of a few hundred μK² there; S₁/₂ measures exactly this lack of large-angle correlation. ',
        badge('computed'),
      ),
    );
  }

  const hist = S.S12_lcdm_hist;
  if (!hist) {
    b.body.append(h('p', { class: 'muted small' }, 'The ΛCDM S₁/₂ distribution is not in this build.'));
    return;
  }
  const cut = 160000;
  const total = hist.counts.reduce((x, y) => x + y, 0) + hist.overflow;
  const xs: number[] = [];
  const ys: number[] = [];
  hist.counts.forEach((c, i) => {
    const lo = hist.edges[i] ?? 0;
    const hi = hist.edges[i + 1] ?? lo;
    if (lo < cut) {
      xs.push(0.5 * (lo + hi));
      ys.push(c / total);
    }
  });
  const obs = S.S12?.wmap ?? null;
  keep(
    ctx,
    barChart(b.body, [{ name: 'ΛCDM realizations', x: xs, y: ys, color: 1 }], {
      xTitle: 'S₁/₂ (μK⁴)',
      yTitle: 'fraction of realizations per bin',
      height: 300,
      refLines: [
        ...(obs !== null ? [{ axis: 'x' as const, value: obs, label: `WMAP ${fmtInt(obs)}` }] : []),
        { axis: 'x', value: S.S12_lcdm_median, label: `ΛCDM median ${fmtInt(S.S12_lcdm_median)}` },
      ],
      layout: { bargap: 0.06, xaxis: { exponentformat: 'none' } as never },
      ariaLabel: 'Distribution of S1/2 under LambdaCDM',
      format: (v) => fmtNum(v, 4),
    }),
  );
  const q = S.S12_lcdm_quantiles ?? {};
  b.body.append(
    h(
      'p',
      { class: 'm1-takeaway' },
      `Bins of ${fmtInt((hist.edges[1] ?? 0) - (hist.edges[0] ?? 0))} μK⁴, shown up to ${fmtInt(cut)} (the tail holds the rest). ΛCDM 5 % / 16 % / 50 % quantiles: ${fmtInt(q['0.05'])} / ${fmtInt(q['0.16'])} / ${fmtInt(q['0.5'])}. Our map gives ${fmtInt(S.S12?.map)}; the paper’s Planck value is ${fmtInt(S.paper?.S12_planck)} (Copi et al. range ${(S.paper?.S12_copi_range ?? []).map((v) => fmtInt(v)).join('–')}). `,
      badge('computed'),
      badge('quoted', '[TOP] Sec. 8'),
    ),
  );
}

// --- axes on a Mollweide --------------------------------------------------------------------------------------

async function mountAxes(el: HTMLElement, ctx: PanelContext, C: CmbBundle): Promise<void> {
  const ax = C.stats.axes;
  const meta = C.likelihood_meta;
  const best = meta?.best_model;
  const paperFcc = meta?.paper_fcc_105_axes_gal;
  if (!ax && !best) return;
  const a = panel({
    title: 'Axes on the sky',
    subtitle: 'Galactic Mollweide, l = 0 at the centre and increasing to the left; each axis is drawn with its antipode',
    badges: [badge('computed', C.source ?? 'map'), badge('quoted', '[TOP] Sec. 9')],
    span: 7,
    stretch: true,
  });
  const b = panel({
    title: 'Ours against the paper',
    subtitle: 'Directions in Galactic (l, b); separations between headless axes',
    badges: [badge('computed'), badge('quoted', '[TOP] Sec. 9')],
    span: 5,
    stretch: true,
  });
  b.root.classList.add('dense-table');
  el.append(grid(a.root, b.root));

  const views = [
    ...(ax ? [{ value: 'multipoles', label: 'Quadrupole & octopole' }] : []),
    ...(best ? [{ value: 'fcc', label: `${best.shape} ${best.d_over_D} D best fit` }] : []),
  ];
  let view = views[0]?.value ?? 'multipoles';
  const seg = segmented({ label: 'Axes', value: view, options: views, onChange: (v) => ((view = v), draw()) });
  const skyHost = h('div', { class: 'm1-sky-host' });
  const legend = h('div', { class: 'res-axes-legend' });
  const explain = h('p', { class: 'm1-takeaway' });
  a.body.append(h('div', { class: 'm1-toolbar' }, seg.root), skyHost, legend, explain);

  const { createMollweide } = await import('../lib/sky');
  if (ctx.unmounted) return;
  const sky: SkyView = createMollweide(skyHost, { hud: true, colorbar: false, gridStepDeg: 30 });
  ctx.onCleanup(() => sky.dispose());

  const dot = (color: string, text: string): HTMLElement => h('span', null, h('span', { class: 'res-dot', style: `background: var(${color})` }), text);
  const ring = (color: string, text: string): HTMLElement => h('span', { style: `color: var(${color})` }, h('span', { class: 'res-ring' }), h('span', { style: 'color: var(--text-2)' }, text));

  function draw(): void {
    sky.clearOverlays();
    clear(legend);
    clear(explain);
    if (view === 'multipoles' && ax) {
      const q: LBPair = [ax.quadrupole.l, ax.quadrupole.b];
      const o: LBPair = [ax.octopole.l, ax.octopole.b];
      for (const [lb, col] of [
        [q, '--series-1'],
        [o, '--series-4'],
      ] as const) {
        sky.addPoint(lb, col);
        sky.addPoint(antipode(lb[0], lb[1]), col);
      }
      if (ax.paper) {
        for (const [lb, col] of [
          [ax.paper.quadrupole, '--series-1'],
          [ax.paper.octopole, '--series-4'],
        ] as const) {
          sky.addSmallCircle(lb, 3, col, { width: 2 });
          sky.addSmallCircle(antipode(lb[0], lb[1]), 3, col, { width: 2 });
        }
      }
      legend.append(
        dot('--series-1', `quadrupole · ${ax.source}`),
        dot('--series-4', `octopole · ${ax.source}`),
        ring('--series-1', `quadrupole · ${ax.paper?.map ?? 'Planck'} (paper)`),
        ring('--series-4', `octopole · ${ax.paper?.map ?? 'Planck'} (paper)`),
      );
      explain.append(
        `Each axis n maximizes Σ m²|a_ℓm(n)|² (${ax.method?.replace(/^.*\((.*)\)$/, '$1') ?? 'angular-momentum dispersion'}), searched on an Nside ${ax.quadrupole.nside_search ?? 64} grid and refined. The two ${ax.source} axes lie ${fmtFixed(ax.angle_deg, 1)}° apart; ${ax.paper?.map ?? 'Planck'}’s lie ${fmtFixed(ax.paper?.angle_deg, 1)}° apart (filled: ours, rings: paper). `,
        badge('computed'),
      );
    } else if (best) {
      best.axes_gal.forEach((lb, i) => {
        sky.addPoint(lb, '--series-3', String(i + 1));
        sky.addPoint(antipode(lb[0], lb[1]), '--series-3');
      });
      for (const lb of paperFcc?.axes ?? []) {
        sky.addSmallCircle(lb, 3, '--series-5', { width: 2 });
        sky.addSmallCircle(antipode(lb[0], lb[1]), 3, '--series-5', { width: 2 });
      }
      legend.append(
        dot('--series-3', `best-fit ${best.shape} ${best.d_over_D} D lattice directions · ${C.source ?? 'map'} (ours)`),
        ring('--series-5', `${paperFcc?.map ?? 'Planck'} best fit (paper)`),
        statusPill('warning', 'not significant'),
      );
      explain.append(
        `The orientation of the ${best.shape} ${best.d_over_D} D torus that maximizes Δln L on ${C.source ?? 'the map'} (filled, numbered) and on Planck in the paper (rings). Two maps of the same sky give nearly the same best fit, but neither is significant: these are not coordinates of detected structures. `,
        badge('computed'),
        badge('quoted', paperFcc?.source ?? '[TOP] Sec. 9'),
      );
    }
  }
  draw();

  // right: tables
  if (ax) {
    const rows = [
      { name: 'Quadrupole ℓ = 2', ours: [ax.quadrupole.l, ax.quadrupole.b] as LBPair, paper: ax.paper?.quadrupole ?? null },
      { name: 'Octopole ℓ = 3', ours: [ax.octopole.l, ax.octopole.b] as LBPair, paper: ax.paper?.octopole ?? null },
    ];
    b.body.append(
      unitLabel('Multipole axes', 'max over n of Σ m²|a_ℓm(n)|²'),
      rawTable<(typeof rows)[number]>(
        [
          { label: 'Axis', cell: (r) => r.name },
          { label: `${ax.source} (ours)`, num: true, cell: (r) => lbText(r.ours) },
          { label: `${ax.paper?.map ?? 'Planck'} (paper)`, num: true, cell: (r) => (r.paper ? lbText(r.paper) : '—') },
          { label: 'separation', num: true, cell: (r) => (r.paper ? `${fmtFixed(axisAngleDeg(r.ours, r.paper), 1)}°` : '—') },
        ],
        rows,
      ).root,
      h(
        'p',
        { class: 'm1-takeaway' },
        `Alignment ${fmtFixed(ax.angle_deg, 1)}° on ${ax.source} (isotropic p = ${fmtP(ax.p_isotropic)}) against ${fmtFixed(ax.paper?.angle_deg, 1)}° on ${ax.paper?.map ?? 'Planck'} (p = ${fmtP(ax.paper?.p_isotropic)}). Isotropic p = 1 − cos θ for two random axes. `,
        badge('computed'),
        badge('quoted', ax.paper?.source ?? '[TOP] Sec. 9.1'),
      ),
    );
  }
  if (best) {
    const pAx = meta?.fcc_105_axes_vs_paper?.paper_axes ?? paperFcc?.axes ?? [];
    // Pair each printed axis with our nearest headless axis, shown in the paper axis' hemisphere,
    // so every row compares the same direction (the two lists come in different orders).
    const pairs = pAx.length
      ? pAx.map((pa, i) => {
          let best_j = 0;
          let bestSep = Infinity;
          best.axes_gal.forEach((o, j) => {
            const d = axisAngleDeg(o, pa);
            if (d < bestSep) [best_j, bestSep] = [j, d];
          });
          const o = best.axes_gal[best_j]!;
          const u = lbVec(o[0], o[1]);
          const v = lbVec(pa[0], pa[1]);
          const direct = u[0] * v[0] + u[1] * v[1] + u[2] * v[2] >= 0;
          return { n: i + 1, paper: pa as [number, number] | null, ours: (direct ? o : antipode(o[0], o[1])) as [number, number], sep: bestSep };
        })
      : best.axes_gal.map((o, i) => ({ n: i + 1, paper: null, ours: o as [number, number], sep: NaN }));
    b.body.append(
      unitLabel(`${best.shape} ${best.d_over_D} D best-fit lattice directions`, '(not significant)'),
      rawTable<(typeof pairs)[number]>(
        [
          { label: '#', num: true, cell: (r) => String(r.n) },
          { label: 'paper (Planck)', num: true, cell: (r) => (r.paper ? lbText(r.paper) : '—') },
          { label: 'nearest of ours', num: true, cell: (r) => lbText(r.ours) },
          { label: 'separation', num: true, cell: (r) => (Number.isFinite(r.sep) ? `${fmtFixed(r.sep, 2)}°` : '—') },
        ],
        pairs,
      ).root,
      h(
        'p',
        { class: 'm1-takeaway' },
        `max Δln L = ${fmtFixed(best.max_dlnL, 2)} with p_max = ${fmtFixed(best.p_max, 2)} on isotropic skies: `,
        statusPill('warning', 'not significant'),
        ` ${meta?.fcc_105_axes_vs_paper?.note ?? best.note ?? ''} `,
        badge('computed'),
      ),
    );
  }
}

// --- torus likelihood --------------------------------------------------------------------------------------------

function mountLikelihood(el: HTMLElement, C: CmbBundle): void {
  const L = C.likelihood;
  const M = C.likelihood_meta;
  const p = panel({
    title: 'Orientation-scanned torus likelihood · independent re-analysis',
    subtitle: M
      ? `${M.map}, full-sky a_ℓm, 2 ≤ ℓ ≤ 30; SO(3) grid ${M.grid.n_alpha} × ${M.grid.n_beta} × ${M.grid.n_gamma}; calibrated on ${M.n_null} isotropic skies and ${M.n_true ?? '—'} torus-true skies per model (${M.n_true_marginal ?? '—'} where marginal). Paper: Planck 2018 Commander, Tables 9a/9b.`
      : undefined,
    badges: [badge('computed', C.source ?? 'map'), badge('quoted', '[TOP] Tables 9a, 9b')],
    span: 12,
  });
  el.append(grid(p.root));
  if (!L?.length) {
    p.body.append(h('p', { class: 'muted' }, 'Our orientation-scanned likelihood is not in this build (cmb.json likelihood = null).'));
    return;
  }
  const g = M?.global;
  const f = M?.fcc_family;
  const best = M?.best_model;
  p.body.append(
    tilesN(
      3,
      statTile({
        label: `Average over all ${g?.n_models ?? '—'} models of this table`,
        value: g ? `ln B = ${sign(g.lnB_avg)}` : null,
        sub: g ? [badge('computed'), `p = ${fmtFixed(g.p, 2)} · paper ${sign(g.paper?.lnB_avg)} (p = ${fmtFixed(g.paper?.p, 2)}, ${g.paper?.n_models ?? '—'} models)`] : undefined,
      }),
      statTile({
        label: 'Pre-registered FCC family (d ≥ 1 D)',
        value: f ? `ln B = ${sign(f.lnB)}` : null,
        sub: f ? [badge('computed'), `p = ${fmtFixed(f.p, 2)} · paper ${sign(f.paper?.lnB)} (p = ${fmtFixed(f.paper?.p, 2)})`] : undefined,
      }),
      statTile({
        label: 'Largest single-model max Δln L',
        value: best ? `${best.shape} ${best.d_over_D} D · ${fmtFixed(best.max_dlnL, 2)}` : null,
        mono: true,
        sub: best ? [badge('computed'), `p_max = ${fmtFixed(best.p_max, 2)}`, statusPill('warning', 'not significant')] : undefined,
      }),
    ),
  );

  let shape = 'all';
  const shapes = [...new Set(L.map((r) => r.shape))];
  const host = h('div');
  const seg = segmented({
    label: 'Shape',
    value: shape,
    options: [{ value: 'all', label: 'All' }, ...shapes.map((s) => ({ value: s, label: SHAPE_LABEL[s] ?? s }))],
    onChange: (v) => ((shape = v), render()),
  });
  p.body.append(h('div', { class: 'm1-toolbar' }, h('span', { class: 'control-label' }, 'Shape'), seg.root, h('span', { class: 'small muted m1-grow' }, 'Each cell: ours on the build map, the paper’s Planck value beneath (— where Table 9a or 9b has no entry).')), host);

  const num = (v: unknown, d = 2): string | null => (typeof v === 'number' ? sign(v, d) : null);
  const cols: RawColumn<LikelihoodRow>[] = [
    { label: 'Torus', cell: (r) => h('strong', null, `${r.shape} ${r.d_over_D} D`) },
    { label: 'max Δln L', num: true, cell: (r) => oursPaper(sign(r.max_dlnL), num(r.paper?.max_dlnL)) },
    { label: 'ln B (optimal)', num: true, cell: (r) => oursPaper(sign(r.lnB), num(r.paper?.lnB ?? r.paper?.lnB_planck)) },
    { label: '95 % of ln B, no topology', num: true, cell: (r) => oursPaper(sign(r.null_q95), num(r.paper?.null95)) },
    { label: 'p (no topology)', num: true, title: 'fraction of isotropic skies with ln B ≥ observed', cell: (r) => oursPaper(fmtFixed(r.p, 2), typeof r.paper?.p === 'number' ? fmtFixed(r.paper.p, 2) : null) },
    { label: 'p_max', num: true, title: 'fraction of isotropic skies with max Δln L ≥ observed (ours)', cell: (r) => fmtFixed(r.p_max, 2) },
    {
      label: 'P(ln B ≤ obs | torus)',
      num: true,
      title: 'torus-true calibration (Table 9b)',
      cell: (r) => oursPaper(r.exclusion ? r.exclusion.p_le_label.replace('<=', '≤ ') : '—', typeof r.paper?.p_le_planck === 'string' ? r.paper.p_le_planck.replace('<=', '≤ ') : null),
    },
    {
      label: 'Verdict (95 %)',
      cell: (r) => {
        const ours = r.exclusion?.verdict;
        const paper = r.paper?.verdict;
        return h(
          'span',
          { class: 'm1-op', style: 'align-items: flex-start' },
          ours ? statusPill(ours.startsWith('excluded') ? 'critical' : 'good', ours) : '—',
          paper ? h('span', { class: 'm1-op-paper' }, `paper: ${paper}`) : null,
        );
      },
      className: 'wrap',
    },
  ];
  const tbl = rawTable(cols, L, (r) => (r.exclusion?.verdict.startsWith('excluded') ? 'm1-row-excl' : null));
  host.append(tbl.root);
  const render = (): void => tbl.render(L.filter((r) => shape === 'all' || r.shape === shape));

  const diffs = M?.verdict_differences ?? [];
  const maxRow = L.reduce((a, r) => (r.lnB > a.lnB ? r : a));
  const minRow = L.reduce((a, r) => (r.p < a.p ? r : a));
  p.body.append(
    proofList(
      [
        {
          kind: 'computed',
          src: C.source ?? 'map',
          text: `No torus shows evidence on this map either: the largest ln B is ${sign(maxRow.lnB)} (${maxRow.shape} ${maxRow.d_over_D} D), far below the Kass–Raftery “positive” level ln B ≥ 1; the smallest single-model p is ${fmtFixed(minRow.p, 2)} (${minRow.shape} ${minRow.d_over_D} D, ln B = ${sign(minRow.lnB)}) before the look-elsewhere correction${g ? `, and averaged over all models p = ${fmtFixed(g.p, 2)}` : ''}. The best-fit orientations are reported for the record only; they are not significant and are not coordinates of detected structures.`,
        },
        ...diffs.map((d) => ({
          kind: 'computed' as const,
          src: 'verdict differs',
          text: `${d.shape} ${d.d_over_D} D: our ${M?.map ?? 'map'} analysis gives P(ln B ≤ observed | torus true) = ${fmtFixed(d.p_le_obs, 3)} over ${d.n_true} torus-true skies, hence “${d.ours}”; the paper’s Planck analysis gives ${d.paper_p_le_planck}, “${d.paper}”. Different maps (foreground residuals near the Galactic plane enter them differently), so both values are shown side by side; neither is a correction of the other.`,
        })),
      ],
      true,
    ),
    note(h('ul', { class: 'res-notes' }, ...(M?.notes ?? []).map((n) => h('li', null, n))), 'How this re-analysis differs from the paper'),
  );
}
