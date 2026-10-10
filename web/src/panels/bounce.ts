/**
 * 1.6 Quantum Big Bang ([QG] Theorems Q, Q', Lambda, N; Sec. 19.3).
 *
 * The closed radiation universe quantized as a Wheeler-DeWitt oscillator with Psi(0) = 0 evolves
 * unitarily through a = 0. We show <a>(T) and <a^-2>(T) for a selectable a_max against the
 * classical a_max |cos T| that hits a = 0, |Psi(a, T)|^2 before, at and after the bounce, the
 * Theorem Lambda table (exact solution, the paper's own finite-difference grid, printed values)
 * with the grid bias explained, the spectrum check E_j = 2j + 3/2, the finite graviton
 * amplitudes, and the quoted mixmaster and Planck no-feature results.
 *
 * Bundle fields: bounce.{runs, snapshots, snapshots_paper_grid, snapshots_a_max_2, spectrum_check,
 * spectrum_check_evolution_grid, graviton, limits, mixmaster, planck_bounce_test, units, notes, claims}.
 */

import './m1.css';
import { type ChartHandle, type LineSeries, type ScatterSeries, lineChart, scatterChart } from '../lib/charts';
import type { BounceBundle, BounceRun, BounceSnapshots, BounceSpectrumCheck } from '../lib/data';
import type { PanelContext } from '../lib/panel';
import { append, badge, classifyStatus, clear, fmtNum, fmtSci, grid, h, loadFailure, note, panel, segmented, statTile } from '../lib/ui';
import { classicalScaleFactor } from './m1-logic';
import { type RawColumn, fmtSigned, keep, kv, oursPaper, proofList, rawTable, sectionLabel, tilesN } from './m1-ui';

export const title = 'Quantum Big Bang';
export const subtitle =
  'The quantized closed radiation universe evolves unitarily through a = 0: the bounce scale, its duration, finite graviton amplitudes, and what Planck says.';

const fmtAmax = (a: number): string => (Number.isInteger(a) ? String(a) : a.toFixed(1));
const num = (x: unknown): number | null => (typeof x === 'number' && Number.isFinite(x) ? x : null);

export async function mount(el: HTMLElement, ctx: PanelContext): Promise<void> {
  ctx.setHeaderExtras([badge('proved', '[QG] Theorem Q'), badge('computed', 'bounce.json')]);
  const res = await ctx.load('bounce.json');
  if (ctx.unmounted) return;
  if (!res.ok) {
    const p = panel({ title: 'Quantum Big Bang', span: 12 });
    p.body.append(loadFailure(res));
    el.append(grid(p.root));
    return;
  }
  const B = res.data;
  if (!B.runs.length) {
    const p = panel({ title: 'Quantum Big Bang', span: 12 });
    p.body.append(h('div', { class: 'empty-state' }, 'bounce.json has no runs.'));
    el.append(grid(p.root));
    return;
  }
  mountEvolution(el, ctx, B);
  mountSnapshots(el, ctx, B);
  mountLambda(el, ctx, B);
  mountChecks(el, ctx, B);
  mountQuoted(el, B);
}

// --- <a>(T), <a^-2>(T) for a selectable a_max ------------------------------------------------------

function mountEvolution(el: HTMLElement, ctx: PanelContext, B: BounceBundle): void {
  const runs = [...B.runs].sort((a, b) => a.a_max - b.a_max);
  const ell = B.units.ell_over_lP;
  let run = runs.find((r) => r.a_max === 2) ?? runs[0]!;

  const p = panel({
    title: 'Expectation values through the bounce · [QG] Theorem Q',
    subtitle: 'Exact solution for a packet released at rest at a_max (odd coherent-state combination). Classically a = a_max |cos T| reaches a = 0 at T = π/2; quantum mechanically ⟨a⟩ turns around at a finite ⟨a⟩_min and ⟨a⁻²⟩ stays finite.',
    badges: [badge('proved', '[QG] Theorem Q'), badge('computed', 'bounce.runs')],
    span: 8,
  });
  const side = panel({
    title: 'What is proved, what is computed',
    subtitle: '[QG] status of each bounce statement, with the units used here.',
    span: 4,
  });
  el.append(grid(p.root, side.root));

  const sel = segmented({
    label: 'a_max',
    options: runs.map((r) => ({ value: String(r.a_max), label: fmtAmax(r.a_max) })),
    value: String(run.a_max),
    onChange: (v) => {
      run = runs.find((r) => String(r.a_max) === v) ?? run;
      render();
    },
  });
  const tileHost = h('div');
  const chartsHost = h('div', { class: 'm1-two' });
  p.body.append(
    h('div', { class: 'm1-toolbar' }, h('span', { class: 'control-label' }, 'Release radius a_max (units of ℓ)'), sel.root),
    tileHost,
    chartsHost,
  );
  let charts: ChartHandle[] = [];
  ctx.onCleanup(() => charts.forEach((c) => c.destroy()));

  function render(): void {
    charts.forEach((c) => c.destroy());
    clear(chartsHost);
    const T = run.T;
    const classical = T.map((t) => classicalScaleFactor(run.a_max, t));
    const aHost = h('div');
    const iHost = h('div');
    chartsHost.append(aHost, iHost);
    const bounce = run.t_bounce ?? Math.PI / 2;
    charts = [
      lineChart(
        aHost,
        [
          { name: '⟨a⟩ quantum (exact)', x: T, y: run.a_mean, color: 1 },
          { name: 'a_max |cos T| classical', x: T, y: classical, color: 'context', dash: 'dash' },
        ],
        {
          xTitle: 'conformal time T (rad)',
          yTitle: 'scale factor a (units of ℓ)',
          height: 280,
          refLines: [{ axis: 'x', value: bounce, label: 'T = π/2' }],
          ariaLabel: 'Mean scale factor through the bounce',
        },
      ),
      lineChart(
        iHost,
        [
          { name: '⟨a⁻²⟩ quantum (exact)', x: T, y: run.a_inv2, color: 1 },
          { name: '1/(a_max cos T)² classical', x: T, y: classical.map((a) => (a > 1e-3 ? 1 / (a * a) : null)), color: 'context', dash: 'dash' },
        ],
        {
          xTitle: 'conformal time T (rad)',
          yTitle: '⟨a⁻²⟩ (units of ℓ⁻²)',
          yType: 'log',
          yRange: [0.01, Math.max(10, 20 * run.peak_a_inv2)],
          height: 280,
          refLines: [{ axis: 'x', value: bounce }],
          ariaLabel: 'Inverse-square scale factor through the bounce',
          caption: 'Classically 1/a² diverges at T = π/2 (vertical line); the quantum value peaks and stays finite.',
        },
      ),
    ];

    const pp = run.paper;
    const printedMin = num(pp.a_min) ?? num(pp.a_min_theorem_Q) ?? num(pp.a_min_theorem_N);
    const printedMinSrc = num(pp.a_min) !== null ? 'Theorem Λ table' : num(pp.a_min_theorem_Q) !== null ? 'Theorem Q(iii)' : 'Theorem N';
    const g = run.paper_grid;
    tileHost.replaceChildren(
      tilesN(
        4,
        statTile({
          label: '⟨a⟩_min',
          value: fmtNum(run.a_min, 5),
          unit: 'ℓ',
          sub: [badge('computed', 'exact'), `${fmtNum(run.a_min * ell, 4)} l_P${printedMin !== null ? `; printed ${fmtNum(printedMin, 4)} (${printedMinSrc})` : ''}`],
        }),
        statTile({
          label: 'Peak ⟨a⁻²⟩',
          value: fmtNum(run.peak_a_inv2, 5),
          unit: 'ℓ⁻²',
          sub: [badge('computed', 'exact'), num(pp.peak_a_inv2) !== null ? `printed ${fmtNum(num(pp.peak_a_inv2), 4)}, paper's grid ${fmtNum(g?.peak_a_inv2, 4)}` : 'finite: no singularity'],
        }),
        statTile({
          label: 'Bounce duration ΔT (FWHM of ⟨a⁻²⟩)',
          value: run.fwhm_T === null ? '—' : fmtNum(run.fwhm_T, 4),
          sub: [badge('computed'), `ΔT · a_max = ${fmtNum(run.fwhm_times_amax ?? null, 4)}${num(pp.fwhm_times_amax) !== null ? `; printed ${fmtNum(num(pp.fwhm_times_amax), 4)}` : ''}`],
        }),
        statTile({
          label: 'Unitarity',
          value: run.norm_max_dev === undefined ? '—' : fmtSci(run.norm_max_dev, 1),
          sub: [badge('computed'), `max |‖Ψ‖ − 1|; period ${fmtNum(run.period ?? Math.PI, 5)} = π (Theorem Q(ii))`],
        }),
      ),
    );
  }
  render();

  const claims = B.claims ?? [];
  append(
    side.body,
    claims.length
      ? proofList(
          claims.map((c) => ({ kind: classifyStatus(c.status), src: c.source, text: [c.statement, h('span', { class: 'small muted' }, ` (${c.status})`)] })),
          true,
        )
      : null,
    sectionLabel('Units'),
    kv([
      ['ℓ', [`${fmtNum(ell, 4)} l_P `, badge('quoted', B.units.source ?? '[QG] Sec. 4'), h('span', { class: 'small muted' }, ` ${B.units.definition ?? ''}`)]],
      ['⟨a⟩_min limit', [`${fmtNum(B.limits?.a_min ?? null, 4)} ℓ = ${fmtNum(num(B.units.a_min_lP) ?? (B.limits?.a_min ?? NaN) * ell, 3)} l_P `, badge('computed', 'a_max → ∞')]],
      ['T', 'conformal time; the bounce is at T = π/2, the period is π'],
    ]),
  );
}

// --- |Psi(a, T)|^2 snapshots --------------------------------------------------------------------

function mountSnapshots(el: HTMLElement, ctx: PanelContext, B: BounceBundle): void {
  const sets: { key: string; label: string; s: BounceSnapshots; what: string }[] = [
    { key: 'exact6', label: 'a_max = 6 exact', s: B.snapshots, what: 'exact solution, a_max = 6 ([QG] Fig. 1 right)' },
  ];
  if (B.snapshots_paper_grid)
    sets.push({
      key: 'grid6',
      label: "a_max = 6 paper's grid",
      s: B.snapshots_paper_grid,
      what: `the paper's finite differences on (0, ${B.snapshots_paper_grid.grid?.L ?? 40}) with ${B.snapshots_paper_grid.grid?.n ?? 8000} points, as used for Fig. 1 (right)`,
    });
  if (B.snapshots_a_max_2) sets.push({ key: 'exact2', label: 'a_max = 2 exact', s: B.snapshots_a_max_2, what: 'exact solution, a_max = 2' });
  let cur = sets[0]!;

  const p = panel({
    title: 'The wave packet through a = 0 · [QG] Fig. 1 (right)',
    subtitle: '|Ψ(a, T)|² shortly before, at and after the bounce. Ψ(0) = 0 at all times: the packet is reflected at a = 0, not lost.',
    badges: [badge('computed', 'bounce.snapshots')],
    span: 12,
  });
  const sel = segmented({
    label: 'Snapshot set',
    options: sets.map((s) => ({ value: s.key, label: s.label })),
    value: cur.key,
    onChange: (v) => {
      cur = sets.find((s) => s.key === v) ?? cur;
      render();
    },
  });
  const host = h('div');
  const cap = h('p', { class: 'm1-takeaway' });
  p.body.append(h('div', { class: 'm1-toolbar' }, sel.root), host, cap);
  el.append(grid(p.root));
  let chart: ChartHandle | null = null;
  ctx.onCleanup(() => chart?.destroy());

  function render(): void {
    chart?.destroy();
    clear(host);
    const s = cur.s;
    const half = Math.PI / 2;
    const name = (t: number): string => {
      const d = t - half;
      return Math.abs(d) < 1e-3 ? 'T = π/2 (bounce)' : `T = π/2 ${d < 0 ? '−' : '+'} ${Math.abs(d).toFixed(2)} (${d < 0 ? 'before' : 'after'})`;
    };
    const series: LineSeries[] = s.T.map((t, i) => ({
      name: name(t),
      x: s.a,
      y: s.psi2[i] ?? [],
      color: Math.abs(t - half) < 1e-3 ? 1 : 'context',
      dash: t > half + 1e-3 ? 'dash' : 'solid',
    }));
    chart = lineChart(host, series, {
      xTitle: 'scale factor a (units of ℓ)',
      yTitle: '|Ψ(a, T)|² (per ℓ)',
      height: 300,
      ariaLabel: 'Probability density of the scale factor',
    });
    const before = s.psi2[0] ?? [];
    const after = s.psi2[2] ?? [];
    const maxDiff = before.reduce((m, v, i) => Math.max(m, Math.abs(v - (after[i] ?? 0))), 0);
    cap.textContent =
      `Source: ${cur.what}. ` +
      (maxDiff < 1e-9
        ? 'Before and after coincide exactly: the exact evolution is time-symmetric about the bounce.'
        : `Before and after differ by up to ${fmtNum(maxDiff, 2)} per ℓ: the finite-difference grid disperses the packet, which also delays the bounce slightly.`);
  }
  render();
}

// --- Theorem Lambda table ----------------------------------------------------------------------

interface LambdaRow {
  run: BounceRun;
}

function mountLambda(el: HTMLElement, ctx: PanelContext, B: BounceBundle): void {
  const runs = [...B.runs].sort((a, b) => a.a_max - b.a_max);
  const tab = runs.filter((r) => num(r.paper.peak_a_inv2) !== null || r.paper_grid);
  const lim = B.limits;
  const ell = B.units.ell_over_lP;

  const p = panel({
    title: 'Theorem Λ table: exact solution, paper\'s grid, printed',
    subtitle: 'The bounce scale and duration for the a_max the paper tabulates. Our finite-difference rerun on the paper\'s grid reproduces every printed number; the exact solution differs in the direction a grid error predicts.',
    badges: [badge('computed', 'runs[].paper_grid'), badge('quoted', '[QG] Theorem Λ')],
    span: 12,
  });
  p.root.classList.add('dense-table');
  const cols: RawColumn<LambdaRow>[] = [
    { label: 'a_max (ℓ)', num: true, cell: (r) => fmtAmax(r.run.a_max) },
    { label: '⟨a⟩_min exact', num: true, cell: (r) => fmtNum(r.run.a_min, 5) },
    { label: '⟨a⟩_min grid / printed', num: true, cell: (r) => oursPaper(fmtNum(r.run.paper_grid?.a_min, 5), fmtNum(num(r.run.paper.a_min), 4)) },
    { label: 'peak ⟨a⁻²⟩ exact', num: true, cell: (r) => fmtNum(r.run.peak_a_inv2, 5) },
    { label: 'peak ⟨a⁻²⟩ grid / printed', num: true, cell: (r) => oursPaper(fmtNum(r.run.paper_grid?.peak_a_inv2, 5), fmtNum(num(r.run.paper.peak_a_inv2), 4)) },
    { label: 'exact / printed', num: true, title: 'Ratio of the exact peak to the printed peak', cell: (r) => (num(r.run.paper.peak_a_inv2) ? fmtNum(r.run.peak_a_inv2 / num(r.run.paper.peak_a_inv2)!, 4) : '—') },
    { label: 'ΔT·a_max exact', num: true, cell: (r) => fmtNum(r.run.fwhm_times_amax ?? null, 4) },
    { label: 'ΔT·a_max grid / printed', num: true, cell: (r) => oursPaper(fmtNum(r.run.paper_grid?.fwhm_times_amax ?? null, 4), fmtNum(num(r.run.paper.fwhm_times_amax), 4)) },
    { label: 'bounce delay on grid', num: true, title: 'T_bounce on the grid minus π/2', cell: (r) => (r.run.paper_grid ? fmtSigned(r.run.paper_grid.bounce_delay, 2) : '—') },
  ];
  const tbl = rawTable(cols, tab.map((run) => ({ run })));

  const chartHost = h('div');
  const exactX = runs.map((r) => r.a_max);
  const printed = tab.filter((r) => num(r.paper.a_min) !== null);
  const gridRuns = tab.filter((r) => r.paper_grid);
  keep(
    ctx,
    scatterChart(
      chartHost,
      [
        { name: 'exact solution', x: exactX, y: runs.map((r) => r.a_min), color: 1, mode: 'lines+markers' },
        { name: "paper's grid (our rerun)", x: gridRuns.map((r) => r.a_max), y: gridRuns.map((r) => r.paper_grid!.a_min), color: 4, symbol: 'square' },
        { name: 'printed', x: printed.map((r) => r.a_max), y: printed.map((r) => num(r.paper.a_min)!), color: 'context', symbol: 'diamond', size: 10 },
      ],
      {
        xTitle: 'release radius a_max (units of ℓ)',
        yTitle: '⟨a⟩_min (units of ℓ)',
        xType: 'log',
        height: 290,
        refLines: lim ? [{ axis: 'y', value: lim.a_min }] : [],
        caption: lim ? `The horizontal line is the exact large-a_max limit ⟨a⟩_min → ${fmtNum(lim.a_min, 4)} ℓ.` : undefined,
      },
    ),
  );

  const notes = B.notes ?? [];
  p.body.append(
    tbl.root,
    h('div', { class: 'chart-caption' }, "Grid: [QG] Sec. 21 finite differences on (0, 24) with 6000 points (our rerun). Under each grid value, the paper's printed value."),
    h(
      'div',
      { class: 'm1-two' },
      chartHost,
      h(
        'div',
        { class: 'callout m1-callout' },
        h('strong', null, 'Why the printed table drifts with a_max. '),
        h(
          'p',
          null,
          'The grid sum misses an O(h) boundary term at a = 0, exactly where |Ψ|²/a² is largest at the bounce, so the printed peaks (262, 533, 1035, 1917) fall increasingly below the exact values; second-order dispersion also slows the packet, delaying the bounce on the grid by up to 0.063 at a_max = 16.',
        ),
        h(
          'p',
          null,
          `Exactly, ⟨a⟩_min tends to a universal ${fmtNum(lim?.a_min ?? null, 4)} ℓ = ${fmtNum((lim?.a_min ?? NaN) * ell, 3)} l_P (√(12/π)/12), so the printed decrease to 0.1518 at a_max = 16 is the same discretisation effect, not physics. ΔT · a_max tends to ${fmtNum(lim?.fwhm_times_amax ?? null, 3)} (printed ≈ ${fmtNum(lim?.paper_fwhm_times_amax ?? null, 2)}).`,
        ),
        h('p', null, 'The paper\'s qualitative claims, a universal Planck-scale bounce of duration ∝ 1/a_max, hold with either set of numbers.'),
        h('div', { class: 'row small' }, badge('computed', 'exact and grid'), badge('quoted', '[QG] Theorem Λ table')),
      ),
    ),
    notes.length ? note(notes, 'Pipeline notes (bounce.json)') : '',
  );
  el.append(grid(p.root));
}

// --- spectrum check and gravitons -----------------------------------------------------------------

function relErr(s: BounceSpectrumCheck): { x: number[]; y: (number | null)[] } {
  const x: number[] = [];
  const y: (number | null)[] = [];
  s.j.forEach((j, i) => {
    if (j < 1) return;
    const e = s.E_exact[i]!;
    x.push(j);
    y.push(Math.abs((s.E_numeric[i]! - e) / e));
  });
  return { x, y };
}

function mountChecks(el: HTMLElement, ctx: PanelContext, B: BounceBundle): void {
  const S = B.spectrum_check;
  const S2 = B.spectrum_check_evolution_grid;
  const p = panel({
    title: 'Spectrum check: E_j = 2j + 3/2 · [QG] Theorem Q(ii)',
    subtitle: 'Relative error of the finite-difference eigenvalues against the exact spectrum. Integer spacing 2 makes the evolution exactly periodic with period π.',
    badges: [badge('proved', '[QG] Theorem Q(ii)'), badge('computed', 'spectrum_check')],
    span: 6,
  });
  const series: ScatterSeries[] = [{ name: `grid L = ${S.grid?.L ?? '?'}, n = ${S.grid?.n ?? '?'}`, ...relErr(S), color: 1, mode: 'lines+markers' }];
  if (S2) series.push({ name: `grid L = ${S2.grid?.L ?? '?'}, n = ${S2.grid?.n ?? '?'}`, ...relErr(S2), color: 4, mode: 'lines+markers', symbol: 'square' });
  keep(
    ctx,
    scatterChart(p.body, series, {
      xTitle: 'level j',
      yTitle: '|E_numeric − E_exact| / E_exact',
      xType: 'log',
      yType: 'log',
      height: 280,
    }),
  );
  const e0 = S.E_numeric[0];
  p.body.append(
    h(
      'p',
      { class: 'm1-takeaway' },
      `Ground level ${fmtNum(e0, 8)} against 3/2. The error grows like (h E)², from the second-order differences, not from the theory; the packets above use at most ${Math.max(...B.runs.map((r) => r.paper_grid?.levels_used ?? 0)) || '—'} levels on the evolution grid, and the exact runs use the closed form.`,
    ),
  );

  const G = B.graviton;
  const q = panel({
    title: "Gravitons through the bounce · [QG] Theorem Q'",
    subtitle: `Every graviton mode is an exact oscillator of integer frequency n; its amplitude stays finite: ${G.formula ?? '⟨h_n²⟩ = (N + 1/2)/n ⟨a⁻²⟩'} at the peak, vacuum N = ${G.occupation ?? 0}.`,
    badges: [badge('proved', "[QG] Theorem Q'"), badge('computed', 'graviton.h2_peak_by_a_max')],
    span: 6,
  });
  const by = G.h2_peak_by_a_max ?? {};
  const amaxes = Object.keys(by).map(Number).sort((a, b) => a - b);
  if (amaxes.length) {
    const slots = [1, 4, 5, 6];
    keep(
      ctx,
      lineChart(
        q.body,
        G.n.map((n, k) => ({ name: `n = ${n}`, x: amaxes, y: amaxes.map((a) => by[String(a)]?.[k] ?? null), color: slots[k % slots.length], markers: true })),
        { xTitle: 'release radius a_max (units of ℓ)', yTitle: 'peak ⟨h_n²⟩ (dimensionless)', xType: 'log', yType: 'log', height: 280 },
      ),
    );
  } else {
    keep(
      ctx,
      scatterChart(q.body, [{ name: `a_max = ${G.a_max ?? 2}`, x: G.n, y: G.h2_peak, color: 1 }], { xTitle: 'graviton mode n', yTitle: 'peak ⟨h_n²⟩ (dimensionless)', height: 280 }),
    );
  }
  const pk = G.paper_a_inv2_peak ?? {};
  q.body.append(
    h(
      'p',
      { class: 'm1-takeaway' },
      'Finite at every a_max: classically h_n would blow up with 1/a² at a = 0. The amplitudes scale with the peak ⟨a⁻²⟩, which we take from the exact solution',
      Object.keys(pk).length ? ` (the paper's printed peaks ${Object.entries(pk).map(([a, v]) => `${v} at a_max = ${a}`).join(', ')} are its grid values; see the table above).` : '.',
    ),
  );
  el.append(grid(p.root, q.root));
}

// --- quoted: mixmaster and Planck ------------------------------------------------------------------

function mountQuoted(el: HTMLElement, B: BounceBundle): void {
  const M = B.mixmaster;
  const P = B.planck_bounce_test;
  const p = panel({
    title: 'Beyond isotropy: the quantum mixmaster bounce · [QG] Theorem N',
    subtitle: 'Quoted from the paper (not recomputed here): with anisotropy channels the universe still bounces at a Planck-scale ⟨a⟩_min and re-expands.',
    badges: [badge('quoted', M.source)],
    span: 6,
  });
  p.body.append(
    tilesN(
      3,
      statTile({ label: '⟨a⟩_min', value: fmtNum(M.a_min, 3), unit: 'ℓ', sub: [badge('quoted', '[QG] Theorem N'), M.a_min_range ? `range ${M.a_min_range.join('–')} over runs` : ''] }),
      statTile({ label: '⟨a⟩ at T = π', value: fmtNum(M.a_mean_at_pi ?? null, 4), unit: 'ℓ', sub: [badge('quoted'), `released at a_max = ${M.a_max ?? '—'}`] }),
      statTile({ label: 'Norm drift', value: fmtSci(M.norm_drift ?? null, 1), sub: [badge('quoted'), 'unitary evolution'] }),
    ),
  );
  if (M.runs?.length)
    p.body.append(
      rawTable(
        [
          { label: 'Run', cell: (r) => r.run },
          { label: '⟨a⟩_min', num: true, cell: (r) => fmtNum(r.a_min, 4) },
          { label: '⟨a⟩(π)', num: true, cell: (r) => fmtNum(r.a_at_pi, 4) },
          { label: 'max ⟨a⁻²⟩', num: true, cell: (r) => fmtNum(r.max_a_inv2, 4) },
          { label: '|β|² at π', num: true, cell: (r) => fmtNum(r.beta2_at_pi, 3) },
        ],
        M.runs,
      ).root,
      h('div', { class: 'chart-caption' }, `${M.source}. Our isotropic exact value for comparison: ⟨a⟩_min = ${fmtNum(B.runs.find((r) => r.a_max === 1.5)?.a_min ?? null, 4)} at a_max = 1.5.`),
    );

  const q = panel({
    title: 'Does the CMB show a bounce? · [QG] Sec. 19.3',
    subtitle: P.data ?? 'Planck 2018 CMB power spectra',
    badges: [badge('quoted', P.source)],
    span: 6,
  });
  const tmpl = P.templates ?? [];
  q.body.append(
    tmpl.length
      ? rawTable(
          [
            { label: 'Template', cell: (t) => h('strong', null, t.name) },
            { label: 'Primordial modulation', cell: (t) => h('span', { class: 'm1-mono' }, t.form), className: 'wrap' },
            { label: 'Extra params', num: true, cell: (t) => String(t.extra_params) },
            { label: 'Δχ² max', num: true, cell: (t) => fmtNum(t.dchi2_max, 2) },
            { label: 'p (chance)', num: true, cell: (t) => fmtNum(t.p_chance, 2) },
          ],
          tmpl,
        ).root
      : h('div', { class: 'empty-state' }, 'No template table in the bundle.'),
    tilesN(
      2,
      statTile({ label: 'Cutoff scale excluded above', value: fmtSci(P.kc_excluded_above, 2), unit: 'Mpc⁻¹', sub: [badge('quoted', P.source), `Δχ² > ${P.kc_exclusion_dchi2 ?? 3.84}`] }),
      statTile({ label: 'Null skies', value: P.n_null_skies ?? '—', sub: [badge('quoted'), 'calibrate the chance probabilities'] }),
    ),
    h(
      'div',
      { class: 'callout m1-callout' },
      h('strong', null, 'No bounce feature in Planck 2018. '),
      h(
        'p',
        null,
        `The best template improves χ² by at most ${fmtNum(Math.max(...tmpl.map((t) => t.dchi2_max), P.dchi2_cutoff_max), 2)} for up to three extra parameters, which chance does ${tmpl.length ? `${Math.round(100 * Math.min(...tmpl.map((t) => t.p_chance)))}–${Math.round(100 * Math.max(...tmpl.map((t) => t.p_chance)))} %` : 'often'} of the time. This is consistent with the framework, which predicts no Bogoliubov mixing at any scale; it is not a detection of anything.`,
      ),
    ),
  );
  el.append(grid(p.root, q.root));
}
