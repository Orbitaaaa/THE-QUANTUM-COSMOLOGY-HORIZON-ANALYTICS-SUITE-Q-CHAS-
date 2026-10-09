/**
 * 2.2 Dynamic Page Curve Plotter ([BH] Theorems B-F, Sec. 11, Figs. 1-2).
 *
 * Hawking's curve, the capacity A/4G, the Page curve min{S_H, A/4G} and the island's minimal
 * entropy, recomputed live (physics/page.ts) for any beta; the interior volume of Theorem F; the
 * exactly unitary 20-qubit evaporation with diaries thrown in early and late; and an absolute
 * time axis for a real LIGO remnant. Bundle fields: horizon.page.{page_time, page_mass_fraction,
 * beta}, horizon.evaporation.{no_diary, early_diary, late_diary, paper, paper_comparable},
 * gw.events[].{name, mf, t_evap_log10_yr, page_time_log10_yr, T_H_K, lifetime_model,
 * mass_frame}; validation rows "Page transfer law ([BH] Fig. 1)", "Unitary evaporation ([BH] Sec. 11)".
 */

import './m2.css';
import { type ChartHandle, type ChartOptions, type LineSeries, lineChart } from '../lib/charts';
import type { EvaporationRun, GwEvent } from '../lib/data';
import type { PanelContext } from '../lib/panel';
import { BETA_RADIATION, N_QUANTA_COEFF, PAGE_MASS_FRACTION_PAPER, PAGE_TIME_PAPER } from '../lib/physics/constants';
import { curves, interiorVolume } from '../lib/physics/page';
import { badge, clear, fmtNum, grid, h, loadFailure, panel, segmented, select, slider, statTile, statusPill } from '../lib/ui';
import { fmtPow10, log10InteriorVolumeEndM3, MSUN_PLANCK } from './horizon-calc';
import { diaryChart, keep, proofList, tilesN, updateLineChart, validationBlock } from './m2-ui';

export const title = 'Dynamic Page Curve Plotter';
export const subtitle =
  "Field theory alone forces Hawking's ever-rising curve; finite capacity A/4G forces the island and the Page curve. Recomputed live, and checked against an exactly unitary evaporation of qubits.";

const RUNS = [
  { value: 'early_diary', label: 'Diary in early (step 5)' },
  { value: 'late_diary', label: 'Diary in late (step 14)' },
  { value: 'no_diary', label: 'No diary' },
] as const;
type RunKey = (typeof RUNS)[number]['value'];

export async function mount(el: HTMLElement, ctx: PanelContext): Promise<void> {
  ctx.setHeaderExtras([badge('proved', '[BH] Theorems B–F'), badge('computed', '[BH] Sec. 11'), badge('live')]);
  const [hor, gw, val] = await Promise.all([ctx.load('horizon.json'), ctx.load('gw.json'), ctx.load('validation.json')]);
  if (ctx.unmounted) return;
  const events: GwEvent[] = gw.ok ? gw.data.events.filter((e) => Number.isFinite(e.t_evap_log10_yr)) : [];
  const bundlePage = hor.ok ? hor.data.page : null;

  // --- Page transfer law -----------------------------------------------------------------------------------
  const lawP = panel({
    title: 'Page transfer law · [BH] Fig. 1',
    subtitle: 'Entropy of the radiation in units of the initial Bekenstein–Hawking entropy S_BH(M₀), massless scalar',
    badges: [badge('proved', '[BH] Theorem E'), badge('live', 'physics/page.ts')],
    span: 8,
  });
  const state = { beta: bundlePage?.beta ?? BETA_RADIATION.value, axis: 'fraction' };
  const betaSlider = slider({
    label: 'β = S_rad / S_BH (Hawking endpoint)',
    min: 0.5,
    max: 4,
    step: 0.005,
    value: state.beta,
    format: (v) => v.toFixed(3),
    onInput: (v) => {
      state.beta = v;
      schedule();
    },
  });
  const resetBtn = h('button', { type: 'button', class: 'btn', title: 'Massless scalar, [BH] Theorem C' }, `β = ${BETA_RADIATION.value} (massless scalar)`);
  resetBtn.addEventListener('click', () => {
    betaSlider.set(BETA_RADIATION.value);
    state.beta = BETA_RADIATION.value;
    schedule();
  });
  const axisSel = select({
    label: 'Time axis',
    options: [{ value: 'fraction', label: 't / t_evap (any hole)' }, ...events.map((e) => ({ value: e.name, label: `years · ${e.name} remnant` }))],
    value: 'fraction',
    onChange: (v) => {
      state.axis = v;
      buildLaw();
      buildVolume();
    },
  });
  lawP.body.append(
    h('div', { class: 'm2-controls' }, betaSlider.root, h('div', { class: 'control' }, h('span', { class: 'control-label' }, 'Reset'), resetBtn), axisSel.root),
  );
  const lawTiles = h('div');
  const lawHost = h('div');
  const lawCaption = h('div', { class: 'chart-caption' });
  lawP.body.append(lawTiles, lawHost, lawCaption);

  const provesP = panel({ title: 'What the paper proves', subtitle: '[BH] Secs. 4–8 and the status table of Sec. 16', span: 4, stretch: true });
  provesP.body.append(
    proofList(
      [
        { kind: 'proved', src: '[BH] Theorem B', text: 'Field theory alone forces Hawking’s curve: on a fixed background the radiation never learns what fell in, its entanglement with the interior only grows, and the final state is mixed. No field-theory construction gives Page’s curve.' },
        { kind: 'proved', src: '[BH] Theorem C', text: '4D radiation is a genuine density matrix with entropy β S_BH(M₀), β = 1.855 for a massless scalar.' },
        { kind: 'proved', src: '[BH] Theorem D', text: 'Finite capacity forces the island: once S_H > A/4G, at least S_H − A/4G of the partners’ entropy must sit in the radiation. Given capacity and H*.' },
        { kind: 'proved', src: '[BH] Theorem E', text: 'S(R) ≤ min{S_H, A/4G}, equality within e^(−|S_H − A/4G|) for 2-design scrambling; after the Page time dS/du = (1/4G) dA/du. Given capacity and H*.' },
        { kind: 'proved', src: '[BH] Theorem F', text: 'The interior volume grows to ~M₀⁵ while the capacity goes to zero.' },
        { kind: 'assumed', src: '[BH] Sec. 16', text: 'H* for an evaporating hole is derived by the gravitational path integral, not proved (a theorem for eternal holes, Theorem K).' },
      ],
      true,
    ),
  );
  el.append(grid(lawP.root, provesP.root));

  let lawChart: ChartHandle | null = null;
  let lawSeries: LineSeries[] = [];
  let lawOpts: ChartOptions = {};
  const event = (): GwEvent | undefined => events.find((e) => e.name === state.axis);

  function lawTilesRender(c: ReturnType<typeof curves>): void {
    clear(lawTiles);
    const e = event();
    const atDefault = Math.abs(state.beta - BETA_RADIATION.value) < 1e-9;
    lawTiles.append(
      tilesN(
        4,
        statTile({
          label: 'Page time t_P / t_evap',
          value: fmtNum(c.page_time, 5),
          badge: badge('live'),
          sub: atDefault
            ? [badge('quoted', PAGE_TIME_PAPER.source), ` paper ${PAGE_TIME_PAPER.value}`, bundlePage ? ` · pipeline ${fmtNum(bundlePage.page_time, 5)}` : '']
            : `paper ${PAGE_TIME_PAPER.value} at β = 1.855`,
        }),
        statTile({
          label: 'Mass left at the Page time',
          value: fmtNum(100 * c.page_mass_fraction, 4),
          unit: '%',
          badge: badge('live'),
          sub: atDefault ? [badge('quoted', PAGE_MASS_FRACTION_PAPER.source), ` paper ${fmtNum(100 * PAGE_MASS_FRACTION_PAPER.value, 3)} %`] : `paper 80.6 % at β = 1.855`,
        }),
        statTile({ label: 'Hawking entropy at the end', value: fmtNum(state.beta, 4), unit: '× S_BH(M₀)', badge: badge('proved', '[BH] Thm C'), sub: 'more than the hole ever held' }),
        e
          ? statTile({
              label: `${e.name}: Page time`,
              value: fmtPow10(e.t_evap_log10_yr + Math.log10(c.page_time)),
              unit: 'yr',
              badge: badge('assumed', 'scalar emission'),
              sub: `t_evap = ${fmtPow10(e.t_evap_log10_yr)} yr; gw.json Page time ${fmtPow10(e.page_time_log10_yr)} yr at β = 1.855`,
            })
          : statTile({ label: 'Quanta emitted in all', value: fmtNum(N_QUANTA_COEFF.value, 4), unit: 'M₀²/m_P²', badge: badge('quoted', N_QUANTA_COEFF.source), sub: 'finite: a density matrix on Fock space' }),
      ),
    );
  }

  function buildLaw(): void {
    const c = curves(401, state.beta);
    const e = event();
    const scale = e ? 10 ** e.t_evap_log10_yr : 1;
    const x = c.t.map((t) => t * scale);
    lawSeries = [
      { name: "Hawking's curve S_H (field theory, Theorems B, C)", x, y: c.hawking, color: 2 },
      { name: 'black-hole entropy A/4G (capacity)', x, y: c.bh_entropy, color: 'context', dash: 'dash' },
      { name: 'Page curve min{S_H, A/4G} (Theorem E)', x, y: c.page, color: 1, width: 2.5 },
      { name: 'island minimal entropy S_H − A/4G (Theorem D)', x, y: c.island_min, color: 5, fill: true },
    ];
    lawOpts = {
      xTitle: e ? `time since formation (years) · ${e.name} remnant` : 'time / lifetime t_evap',
      yTitle: 'entropy / S_BH(M₀)',
      height: 360,
      yRange: [0, Math.max(2, state.beta * 1.05)],
      refLines: [{ axis: 'x', value: c.page_time * scale, label: `Page time ${fmtNum(c.page_time, 3)}` }],
      ariaLabel: "Hawking's curve, black-hole entropy, Page curve and island entropy against time",
      format: (v) => fmtNum(v, 5),
    };
    lawChart?.destroy();
    clear(lawHost);
    lawChart = keep(ctx, lineChart(lawHost, lawSeries, lawOpts));
    lawTilesRender(c);
    lawCaption.textContent = e
      ? `Absolute time for the ${e.name} remnant: t_evap = ${fmtPow10(e.t_evap_log10_yr)} yr from gw.json — massless-scalar Hawking emission, isolated hole (${e.lifetime_model ?? 'Schwarzschild, α = 7.44×10⁻⁵/M²'}), ${e.mass_frame ?? 'detector'}-frame mass ${fmtNum(e.mf, 4)} M☉, spin ignored. With T_H = ${fmtNum(e.T_H_K, 3)} K, far below today’s 2.7 K microwave background, the real remnant currently absorbs more than it emits; the clock describes an isolated hole.`
      : 'Field theory on the fixed background (Theorems B, C) gives Hawking’s curve, which ends at β S_BH(M₀). Finite capacity gives the Page curve of Theorem E: it turns over at the Page time and then falls with the area to zero. The shaded island entropy is the part of the interior that must be encoded in the radiation (Theorem D).';
  }

  let raf = 0;
  function schedule(): void {
    cancelAnimationFrame(raf);
    raf = requestAnimationFrame(() => {
      if (!lawChart) return buildLaw();
      const c = curves(401, state.beta);
      const e = event();
      const scale = e ? 10 ** e.t_evap_log10_yr : 1;
      lawTilesRender(c);
      void updateLineChart(
        lawChart,
        lawSeries,
        {
          y: [c.hawking, c.bh_entropy, c.page, c.island_min],
          refLines: [{ axis: 'x', value: c.page_time * scale, label: `Page time ${fmtNum(c.page_time, 3)}` }],
          // keep the y range wide enough for the Hawking endpoint
          yRange: [0, Math.max(2, state.beta * 1.05)],
        },
        lawOpts,
      );
    });
  }
  ctx.onCleanup(() => cancelAnimationFrame(raf));

  // --- Theorem F: area, not volume ------------------------------------------------------------------------------
  const volP = panel({
    title: 'Area, not volume · Theorem F',
    subtitle: 'Christodoulou–Rovelli interior volume against the capacity A/4G during evaporation',
    badges: [badge('proved', '[BH] Theorem F'), badge('live')],
    span: 5,
    stretch: true,
  });
  const volTiles = h('div');
  const volHost = h('div');
  volP.body.append(volTiles, volHost);
  let volChart: ChartHandle | null = null;
  function buildVolume(): void {
    const n = 201;
    const t = Array.from({ length: n }, (_, i) => i / (n - 1));
    const vEnd = interiorVolume(1);
    const e = event();
    clear(volTiles);
    volTiles.append(
      tilesN(
        2,
        statTile({ label: 'V_end = (9√3π/5) M₀² t_evap', value: fmtNum(vEnd, 4), unit: 'M₀² t_evap', badge: badge('proved', '[BH] Theorem F') }),
        e
          ? statTile({
              label: `${e.name}: V_end = 3√3π M₀⁵/(5α)`,
              value: fmtPow10(log10InteriorVolumeEndM3(e.mf * MSUN_PLANCK)),
              unit: 'm³',
              badge: badge('assumed', 'massless scalar'),
              sub: 'while its capacity A/4G → 0',
            })
          : statTile({ label: 'Capacity at the endpoint', value: '0', unit: 'nats', badge: badge('proved', '[BH] Theorem E(iii)'), sub: 'S(R) → 0 as A → 0' }),
      ),
    );
    volChart?.destroy();
    clear(volHost);
    volChart = keep(
      ctx,
      lineChart(
        volHost,
        [
          { name: 'interior volume V / V_end', x: t, y: t.map((x) => interiorVolume(x) / vEnd), color: 3 },
          { name: 'capacity A/A₀ = S_BH / S_BH(M₀)', x: t, y: t.map((x) => (1 - x) ** (2 / 3)), color: 'context', dash: 'dash' },
        ],
        {
          xTitle: 'time / lifetime t_evap',
          yTitle: 'fraction of V_end / of A₀',
          height: 280,
          yRange: [0, 1.02],
          ariaLabel: 'Interior volume and horizon area against time',
          caption: 'Both indexed to 1 on one axis: the volume grows monotonically to V_end while the area, and with it the information the hole can hold, falls to zero.',
        },
      ),
    );
  }

  // --- exact unitary evaporation -------------------------------------------------------------------------------
  const evP = panel({
    title: 'Exactly unitary evaporation · [BH] Fig. 2',
    subtitle: '20 qubits, Haar-scrambled before every emission (computed by qchas, seeds in horizon.json)',
    badges: [badge('computed', '[BH] Sec. 11'), badge('proved', '[BH] Theorem E(ii)')],
    span: 7,
    stretch: true,
  });
  el.append(grid(volP.root, evP.root));
  buildLaw();
  buildVolume();

  if (!hor.ok) {
    evP.body.append(loadFailure(hor));
  } else {
    const ev = hor.data.evaporation;
    const paper = ev.paper;
    let run: RunKey = 'early_diary';
    const runSeg = segmented({ label: 'Run', options: RUNS.map((r) => ({ ...r })), value: run, onChange: (v) => ((run = v as RunKey), drawRun()) });
    const evTiles = h('div');
    const evHost = h('div');
    const evCaption = h('div', { class: 'chart-caption' });
    evP.body.append(h('div', { class: 'm2-pill-row' }, runSeg.root, ev.paper_comparable === false ? statusPill('warning', 'fast build: run size differs from the paper') : null), evTiles, evHost, evCaption);
    let evChart: ChartHandle | null = null;
    const drawRun = (): void => {
      const r: EvaporationRun = ev[run];
      const isPaperRun = run === 'early_diary';
      clear(evTiles);
      evTiles.append(
        tilesN(
          3,
          statTile({
            label: 'Peak radiation entropy',
            value: fmtNum(r.max_bits, 4),
            unit: `bits at step ${r.max_step ?? '—'}`,
            badge: badge('computed'),
            sub: isPaperRun ? [badge('quoted', '[BH] Sec. 11'), ` paper ${paper.max_bits} at step ${paper.max_step ?? 11}`] : 'paper’s figure: the early-diary run (21 emissions)',
          }),
          statTile({
            label: 'Largest deviation from min{S_H, S_BH}',
            value: fmtNum(r.max_dev, 3),
            unit: 'bits',
            badge: badge('computed'),
            sub: [badge('quoted', '[BH] Sec. 11'), ` paper ${paper.max_dev}, at the Page time`],
          }),
          statTile({
            label: 'Radiation entropy at the end',
            value: fmtNum(Math.abs(r.final_bits ?? NaN) < 1e-9 ? 0 : r.final_bits, 4),
            unit: 'bits',
            badge: badge('computed'),
            sub: r.diary_step ? `= the diary’s entanglement with its reference (paper ${paper.final_bits_with_diary ?? 1})` : 'pure final state: all information out',
          }),
        ),
      );
      evChart?.destroy();
      clear(evHost);
      evChart = keep(
        ctx,
        lineChart(
          evHost,
          [
            { name: 'exact unitary evaporation S(R)', x: r.steps, y: r.s_rad, color: 1, markers: true, width: 0 },
            { name: 'envelope min{S_H, S_BH} (Theorem E)', x: r.steps, y: r.envelope, color: 'context' },
            { name: "Page's exact average", x: r.steps, y: r.page_average, color: 2, dash: 'dash' },
            { name: 'Hawking S_H', x: r.steps, y: r.hawking, color: 2, dash: 'dot', hidden: true },
          ],
          {
            xTitle: 'emitted qubits',
            yTitle: 'radiation entropy (bits)',
            height: 300,
            refLines: r.max_step ? [{ axis: 'x', value: r.max_step, label: 'Page time' }] : [],
            ariaLabel: 'Radiation entropy of the exact unitary evaporation against emitted qubits',
            format: (v) => fmtNum(v, 5),
          },
        ),
      );
      evCaption.textContent =
        `Run ${run.replace('_', ' ')}: ${r.n_qubits} qubits${r.diary_step ? ` plus a diary qubit maximally entangled with an outside reference, thrown in at step ${r.diary_step}` : ''}; seed ${r.seed}. ` +
        'Dots are the exact entropies; they follow S_H up, turn over at the Page time and fall along S_BH, within Page’s average. Hawking’s rising curve is in the legend (click to show).';
    };
    drawRun();

    // diary information and validation
    const diaryP = panel({
      title: 'Where the diary goes · Hayden–Preskill',
      subtitle: 'Bits of the diary recoverable from the radiation, I(Ref : R), for a diary thrown in early and after the Page time',
      badges: [badge('computed', '[BH] Sec. 11'), badge('proved', '[BH] Theorem E(iv)')],
      span: 7,
    });
    const er = ev.early_diary;
    const lr = ev.late_diary;
    diaryP.body.append(
      tilesN(
        2,
        statTile({
          label: `Early diary (in at ${er.diary_step}) recovered at`,
          value: er.recovered_step ?? '—',
          unit: 'emitted qubits',
          badge: badge('computed'),
          sub: [badge('quoted', '[BH] Sec. 11'), ` paper ${paper.early_recovered_step}: only after the Page time`],
        }),
        statTile({
          label: `Late diary (in at ${lr.diary_step}) recovered at`,
          value: lr.recovered_step ?? '—',
          unit: 'emitted qubits',
          badge: badge('computed'),
          sub: [badge('quoted', '[BH] Sec. 11'), ` paper ${paper.late_recovered_step}: ${(lr.recovered_step ?? 0) - (lr.diary_step ?? 0)} emissions later, the Hayden–Preskill mirror`],
        }),
      ),
    );
    const dHost = h('div', { style: { marginTop: '12px' } });
    diaryP.body.append(dHost);
    keep(ctx, diaryChart(dHost, ev, 280));
    diaryP.body.append(
      h(
        'div',
        { class: 'chart-caption' },
        `Recovered means at least ${paper.recovery_threshold_bits ?? 1.9} of the diary’s 2 bits of mutual information with its reference are in the radiation. Thrown in early, the diary waits for the Page time; thrown in late, it comes out after about k + log₂(1/ε) further emissions (Theorem E(iv)).`,
      ),
    );
    const valP = panel({ title: 'Ours vs paper', subtitle: 'validation.json rows for this tool', span: 5, stretch: true, badges: [badge('computed')] });
    const vb = validationBlock(val, ['Page transfer law ([BH] Fig. 1)', 'Unitary evaporation ([BH] Sec. 11)'], { compact: true });
    valP.body.append(vb ?? (val.ok ? h('div', { class: 'muted small' }, 'No validation rows for this tool.') : loadFailure(val)));
    el.append(grid(diaryP.root, valP.root));
  }
}
