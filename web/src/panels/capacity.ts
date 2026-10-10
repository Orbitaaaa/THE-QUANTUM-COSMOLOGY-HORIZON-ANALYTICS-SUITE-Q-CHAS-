/**
 * 2.1 Crossed-Product Capacity Monitor ([BH] Theorem G, [QG] Theorems G_inf, C, U).
 *
 * The algebra ladder III_1 -> II_inf -> II_1 as a stepper, a live horizon calculator (custom mass
 * or a LIGO remnant from gw.json) giving area, A/4G, Hawking temperature and the number of area
 * quanta, the energy-window capacity ln tr Pi against beta DeltaE, and the de Sitter entropy of
 * our universe. Bundle fields: horizon.capacity.{ladder, de_sitter, min_area_lp2,
 * window_examples, window_mass_planck}; gw.events[].{name, mf, af, S_BH_nats, T_H_K,
 * t_evap_log10_yr, page_time_log10_yr, mass_frame}; validation rows of group
 * "Crossed product ([QG] Theorem U)".
 */

import './m2.css';
import './capacity.css';
import { type ChartHandle, type LineSeries, lineChart } from '../lib/charts';
import type { GwEvent, HorizonBundle } from '../lib/data';
import type { Layout } from 'plotly.js';
import type { PanelContext } from '../lib/panel';
import { ALGEBRA_LADDER, deSitterEntropy, minimalArea } from '../lib/physics/capacity';
import { BARBERO_IMMIRZI_SU2, MIN_AREA_LP2, PLANCK2018 } from '../lib/physics/constants';
import { pageTimeFraction } from '../lib/physics/page';
import {
  type Child,
  append,
  badge,
  clear,
  fmtNum,
  fmtSci,
  grid,
  h,
  loadFailure,
  missingBundle,
  panel,
  segmented,
  select,
  slider,
  statTile,
  statusPill,
} from '../lib/ui';
import { blackHole, fmtPow10, logspace, windowCurve } from './horizon-calc';
import { keep, kv, proofList, tilesN, validationBlock } from './m2-ui';

export const title = 'Crossed-Product Capacity Monitor';
export const subtitle =
  'How gravitational dressing turns the traceless type III₁ field algebra into a type II₁ factor whose maximal entropy is A/4G, and what that capacity is for real black holes and for our universe.';

type Rung = HorizonBundle['capacity']['ladder'][number];

/** "III_1" -> "III₁", "II_inf" -> "II∞"; ASCII formulas of the bundle -> symbols. */
const prettyType = (t: string): string => t.replace(/_1$/, '₁').replace(/_inf$/, '∞');
const prettyText = (s: string): string =>
  s
    .replace(/e\^\{([^}]*)\}/g, 'exp($1)')
    .replace(/x_sigma R/g, '⋊σ ℝ')
    .replace(/\(A x R\)/g, '(A ⋊ ℝ)')
    .replace(/\bbeta\b ?/g, 'β')
    .replace(/\bPi\b/g, 'Π')
    .replace(/\bint\b/g, '∫')
    .replace(/<=/g, '≤')
    .replace(/ ~ /g, ' ≈ ')
    .replace(/II_inf/g, 'II∞')
    .replace(/II_1/g, 'II₁')
    .replace(/III_1/g, 'III₁');

/** What changes on each rung, paraphrasing [BH] Sec. 12 and [QG] Sec. 10 (Theorem G_inf). */
const RUNG_DETAIL: { how: string; trace: string; entropy: string; refs: string[] }[] = [
  {
    how: 'The quantum fields outside the horizon form a type III₁ factor (Buchholz–D’Antoni–Fredenhagen; Araki for free fields). The Hartle–Hawking state is KMS at β = 2π/κ, so its modular operator is exp(−βH) with H the Killing energy.',
    trace: 'There is no trace at all, hence no density matrices.',
    entropy: 'The von Neumann entropy of the exterior is undefined (the entanglement across the horizon is UV-infinite). [QG] Theorem G∞(i): the same holds for every region with a non-empty causal complement in every globally hyperbolic spacetime.',
    refs: ['[BH] Sec. 12', '[QG] Theorem G∞(i)'],
  },
  {
    how: 'The Hamiltonian constraint ties H to the ADM energy H + X measured by a clock at infinity. Dressing every exterior observable to the clock time, a ↦ a(q̂) = exp(iq̂H) a exp(−iq̂H), gives the crossed product A ⋊σ ℝ by the modular group; Takesaki’s structure theorem makes it a type II∞ factor (Theorem G(i)–(ii)).',
    trace: 'tr â = ⟨Ψ ⊗ exp(βX/2), â (Ψ ⊗ exp(βX/2))⟩ is faithful, normal and semifinite (G(iii)); the identity has infinite trace. The weight exp(βX) = exp(S_BH(M) − S_BH(M*)) is the Bekenstein–Hawking density of states.',
    entropy: 'Entropies are defined up to one state-independent constant (a trace on a factor is unique only up to scale). For semiclassical states S = A/4G + S_out + S(g) − c*: the generalized entropy (G(vi)). [QG] Theorem G∞(ii)–(iii): III → II holds for every region of every spacetime; the trace is the A/4G density wherever the modular flow is geometric near the boundary.',
    refs: ['[BH] Theorem G(i)–(iii), (vi)', '[QG] Theorem G∞(ii)–(iii)'],
  },
  {
    how: 'Projecting onto an ADM-energy window Π = 1[E₁,E₂](H + X̂) cuts the II∞ factor down to ΠÂΠ, a type II₁ factor (Theorem G(iv)).',
    trace: 'The trace is finite: tr Π = ∫ exp(βX) dX over the window.',
    entropy: 'Every normal state has S ≤ ln tr Π, with equality only for the tracial state (G(v)). With the first law, ln tr Π = A(E₂)/4G − S_BH(M*) + ln[(1 − exp(−βΔE))/β]: the black hole seen from outside has capacity A/4G. With an observer in de Sitter space the same construction gives a II₁ algebra with maximal entropy 3π/GΛ ([QG] Theorem C).',
    refs: ['[BH] Theorem G(iv)–(v)', '[QG] Theorem C'],
  },
];

const TRANSITIONS = ['dress to a clock at infinity', 'fix an energy window'];

export async function mount(el: HTMLElement, ctx: PanelContext): Promise<void> {
  ctx.setHeaderExtras([badge('proved', '[BH] Theorem G'), badge('proved', '[QG] Theorems G∞, C, U'), badge('live')]);
  const [hor, gw, val] = await Promise.all([ctx.load('horizon.json'), ctx.load('gw.json'), ctx.load('validation.json')]);
  if (ctx.unmounted) return;
  const cap = hor.ok ? hor.data.capacity : null;

  // --- algebra ladder ---------------------------------------------------------------------------------
  const ladderP = panel({
    title: 'Algebra ladder',
    subtitle: 'What the gravitational constraint does to the exterior algebra, one step at a time',
    badges: [badge('proved', '[BH] Theorem G'), badge('proved', '[QG] Theorem G∞')],
    span: 7,
    stretch: true,
  });
  const ladder: Rung[] = cap?.ladder?.length ? cap.ladder : ALGEBRA_LADDER;
  let step = 0;
  const rungButtons: HTMLButtonElement[] = [];
  const arrow = (label: string): HTMLElement => {
    const ns = 'http://www.w3.org/2000/svg';
    const svg = document.createElementNS(ns, 'svg');
    svg.setAttribute('viewBox', '0 0 30 14');
    svg.setAttribute('aria-hidden', 'true');
    const path = document.createElementNS(ns, 'path');
    path.setAttribute('d', 'M1 7h26m-6-5l6 5-6 5');
    path.setAttribute('fill', 'none');
    path.setAttribute('stroke', 'currentColor');
    path.setAttribute('stroke-width', '1.6');
    svg.append(path);
    return h('div', { class: 'cap-arrow' }, svg, label);
  };
  const ladderRow = h('div', { class: 'cap-ladder', role: 'group', 'aria-label': 'Algebra ladder steps' });
  ladder.forEach((r, i) => {
    const btn = h(
      'button',
      { type: 'button', class: 'cap-rung', 'aria-pressed': String(i === 0), onclick: () => show(i) },
      h('span', { class: 'cap-rung-step' }, `Step ${i + 1}`),
      h('span', { class: 'cap-rung-type' }, prettyType(r.type)),
      h('span', { class: 'cap-rung-label' }, prettyText(r.label)),
      h('span', { class: 'cap-rung-row' }, h('span', { class: 'k' }, 'trace'), h('span', { class: 'v' }, prettyText(r.trace))),
      h('span', { class: 'cap-rung-row' }, h('span', { class: 'k' }, 'entropy'), h('span', { class: 'v' }, prettyText(r.entropy))),
    );
    rungButtons.push(btn);
    if (i > 0) ladderRow.append(arrow(TRANSITIONS[i - 1] ?? ''));
    ladderRow.append(btn);
  });
  const detail = h('div', { class: 'cap-detail', 'aria-live': 'polite' });
  const prev = h('button', { type: 'button', class: 'btn', onclick: () => show(step - 1) }, '← Previous');
  const next = h('button', { type: 'button', class: 'btn', onclick: () => show(step + 1) }, 'Next →');
  const stepLabel = h('span', { class: 'small muted' });
  function show(i: number): void {
    step = Math.max(0, Math.min(ladder.length - 1, i));
    rungButtons.forEach((b, k) => b.setAttribute('aria-pressed', String(k === step)));
    const r = ladder[step]!;
    const d = RUNG_DETAIL[step];
    clear(detail);
    append(
      detail,
      h('h3', null, `Type ${prettyType(r.type)} · ${prettyText(r.label)}`),
      d ? h('p', null, h('strong', null, step === 0 ? 'Start. ' : `What changes (${TRANSITIONS[step - 1]}). `), d.how) : null,
      d ? h('p', null, h('strong', null, 'Trace. '), d.trace) : null,
      d ? h('p', null, h('strong', null, 'Entropy. '), d.entropy) : null,
      h('p', { class: 'row', style: { gap: '6px' } }, ...(d?.refs ?? [r.source]).map((s) => badge('proved', s)), h('span', { class: 'small muted' }, `bundle source: ${r.source}`)),
    );
    prev.disabled = step === 0;
    next.disabled = step === ladder.length - 1;
    stepLabel.textContent = `Step ${step + 1} of ${ladder.length}`;
  }
  append(
    ladderP.body,
    ladderRow,
    h('div', { class: 'cap-stepper-nav' }, prev, next, stepLabel),
    detail,
    h(
      'div',
      { class: 'callout m2-callout' },
      h('div', { class: 'm2-section-label' }, 'The capacity in one line · [BH] Theorem G(v)'),
      h('div', { class: 'm2-formula' }, 'S ≤ ln ∫_window exp(S_BH(M)) dM − S_BH(M*)'),
      h('div', { class: 'm2-formula' }, '  = A(E₂)/4G − S_BH(M*) + ln[(1 − exp(−β(E₂ − E₁)))/β]'),
      h('div', { class: 'small muted', style: { marginTop: '6px' } }, 'Only entropy differences are canonical: the trace of a factor is unique up to scale, fixed by one state-independent constant c*.'),
    ),
    cap ? null : h('div', { class: 'm2-callout' }, missingBundle('horizon.json', 'The ladder above is the quoted structure of [BH] Theorem G; the computed capacity numbers need the bundle.')),
  );
  show(0);

  const provesP = panel({ title: 'What the paper proves', subtitle: 'Status as the papers’ own tables give it ([BH] Sec. 16, [QG] Sec. 20)', span: 5, stretch: true });
  provesP.body.append(
    proofList([
      { kind: 'proved', src: '[BH] Theorem G', text: 'Gravitational dressing turns the type III₁ exterior algebra into a type II∞ crossed product; an energy window makes it type II₁ with maximal entropy ln tr Π ≈ A/4G, and its entropy is the generalized entropy.' },
      { kind: 'proved', src: '[QG] Theorem G∞', text: 'Type III₁ → type II for every region of every globally hyperbolic spacetime, with no symmetry assumed. Symmetry enters only in identifying the trace with the A/4G density.' },
      { kind: 'proved', src: '[QG] Theorem C', text: 'A universe with Λ > 0 and an observer has a type II₁ algebra; its unique trace is the no-boundary state, with maximal entropy 3π/GΛ.' },
      { kind: 'proved', src: '[QG] Theorem U', text: 'Area has a discrete spectrum with smallest non-zero value a₀ = 4√3 πγ l_P² = 5.97 l_P², γ fixed by A/4G (SU(2) counting); horizon Hilbert spaces are finite.' },
      { kind: 'proved', src: '[BH] Theorem I', text: 'Capacity alone does not give the Page curve: the exterior obeys the bound at all times while the radiation can still follow Hawking’s curve. The remaining input is H* (Theorem J).' },
      { kind: 'assumed', src: '[BH] Sec. 16', text: 'H* for an evaporating hole is derived by the gravitational path integral (replica wormholes), not proved as a theorem in the paper.' },
      { kind: 'computed', src: 'Q-CHAS', text: 'Every number below is computed from these formulas, in the browser or by the Python pipeline. LIGO remnants come from aligned-spin fits to the tutorial’s component parameters (detector frame).' },
    ], true),
  );
  el.append(grid(ladderP.root, provesP.root));

  // --- horizon calculator ---------------------------------------------------------------------------------
  const events: GwEvent[] = gw.ok ? gw.data.events.filter((e) => Number.isFinite(e.mf) && Number.isFinite(e.af)) : [];
  const calcP = panel({
    title: 'Horizon capacity calculator',
    subtitle: 'Choose a black hole: area, capacity A/4G in nats and bits, Hawking temperature and the number of Planck-area quanta',
    badges: [badge('live', 'physics/capacity.ts'), badge('computed', 'Kerr formulas')],
  });
  const state = { mass: events[0]?.mf ?? 10, spin: events[0]?.af ?? 0, preset: events[0]?.name ?? 'custom', xCap: 10, mStarMode: 'pipeline' as 'pipeline' | 'hole' };
  const a0 = cap?.min_area_lp2 ?? minimalArea();
  const options = [{ value: 'custom', label: 'Custom mass and spin' }, ...events.map((e) => ({ value: e.name, label: `${e.name} remnant` }))];
  const presetSel = select({
    label: 'Black hole',
    options,
    value: state.preset,
    onChange: (v) => {
      state.preset = v;
      const e = events.find((x) => x.name === v);
      if (e) {
        state.mass = e.mf;
        state.spin = e.af;
        massSlider.set(Math.log10(e.mf));
        spinSlider.set(e.af);
        massInput.value = String(Number(e.mf.toPrecision(6)));
      }
      update();
    },
  });
  const toCustom = (): void => {
    state.preset = 'custom';
    presetSel.select.value = 'custom';
  };
  const massSlider = slider({
    label: 'Mass log₁₀(M / M☉)',
    min: -19,
    max: 10,
    step: 0.01,
    value: Math.log10(state.mass),
    format: (v) => `${fmtSci(10 ** v, 3)} M☉`,
    onInput: (v) => {
      state.mass = 10 ** v;
      massInput.value = String(Number(state.mass.toPrecision(6)));
      toCustom();
      update();
    },
  });
  const massInput = h('input', { type: 'number', min: 1e-19, max: 1e10, step: 'any', value: String(Number(state.mass.toPrecision(6))), 'aria-label': 'Mass in solar masses' });
  massInput.addEventListener('change', () => {
    const v = Number(massInput.value);
    if (!(v > 0) || !Number.isFinite(v)) return;
    state.mass = Math.min(1e10, Math.max(1e-19, v));
    massSlider.set(Math.log10(state.mass));
    toCustom();
    update();
  });
  const spinSlider = slider({
    label: 'Spin a* = Jc / GM²',
    min: 0,
    max: 0.998,
    step: 0.001,
    value: state.spin,
    format: (v) => v.toFixed(3),
    onInput: (v) => {
      state.spin = v;
      toCustom();
      update();
    },
  });
  calcP.body.append(
    h(
      'div',
      { class: 'm2-controls' },
      presetSel.root,
      massSlider.root,
      h('label', { class: 'control' }, h('span', { class: 'control-label' }, 'Mass (M☉), exact'), massInput),
      spinSlider.root,
    ),
  );
  const calcTiles = h('div');
  const calcFoot = h('div', { class: 'chart-caption' });
  calcP.body.append(calcTiles, calcFoot);
  if (!gw.ok) calcP.body.append(h('div', { class: 'm2-callout' }, loadFailure(gw)));

  // --- energy window chart ---------------------------------------------------------------------------------
  const winP = panel({
    title: 'Energy-window capacity ln tr Π',
    subtitle: 'Maximal entropy of the type II₁ algebra of an ADM window of width ΔE above M*, relative to S_BH(M*), against βΔE',
    badges: [badge('proved', '[BH] Theorem G(iv)–(v)'), badge('live')],
    span: 8,
  });
  const modeSeg = segmented({
    label: 'Reference mass M*',
    options: [
      { value: 'pipeline', label: `Pipeline M* = ${fmtNum(cap?.window_mass_planck ?? 10, 3)}` },
      { value: 'hole', label: 'Selected hole' },
    ],
    value: state.mStarMode,
    onChange: (v) => {
      state.mStarMode = v as 'pipeline' | 'hole';
      drawWindow();
    },
  });
  const capSlider = slider({
    label: 'Window width βΔE',
    min: -2,
    max: 2.5,
    step: 0.01,
    value: Math.log10(state.xCap),
    format: (v) => fmtNum(10 ** v, 3),
    onInput: (v) => {
      state.xCap = 10 ** v;
      drawWindowTiles();
    },
  });
  // the chart (and its reference line) follows on release; tiles follow every tick
  capSlider.input.addEventListener('change', () => drawWindowChart());
  const winControls = h('div', { class: 'm2-controls' }, h('div', { class: 'control' }, h('span', { class: 'control-label' }, 'Reference mass'), modeSeg.root), capSlider.root);
  const winTiles = h('div');
  const winChartHost = h('div');
  winP.body.append(winControls, winTiles, winChartHost);
  let winChart: ChartHandle | null = null;
  const X = logspace(1e-2, 10 ** 2.5, 151);
  const mStarNow = (): number => (state.mStarMode === 'pipeline' ? (cap?.window_mass_planck ?? 10) : blackHole(state.mass, 0).massPlanck);
  const drawWindow = (): void => {
    drawWindowTiles();
    drawWindowChart();
  };

  function drawWindowTiles(): void {
    const mStar = mStarNow();
    const one = windowCurve(mStar, [state.xCap]);
    const lnTr = one.lnTrace[0]!;
    const dS = one.deltaS[0]!;
    const sStar = 4 * Math.PI * mStar * mStar;
    clear(winTiles);
    winTiles.append(
      tilesN(
        4,
        statTile({ label: 'ln tr Π (relative to S_BH(M*))', value: fmtNum(lnTr, 5), unit: 'nats', badge: badge('live') }),
        statTile({ label: 'Area term A(E₂)/4G − S_BH(M*)', value: fmtNum(dS, 5), unit: 'nats', badge: badge('live') }),
        statTile({ label: 'Log correction ln tr Π − area term', value: fmtNum(lnTr - dS, 4), unit: 'nats', sub: `≈ −ln β = −ln(8πM*) = ${fmtNum(-Math.log(8 * Math.PI * mStar), 4)}` }),
        statTile({
          label: 'Type II₁ cap S ≤ S_BH(M*) + ln tr Π',
          value: sStar + lnTr > 1e6 ? fmtPow10(Math.log10(sStar + lnTr)) : fmtNum(sStar + lnTr, 6),
          unit: 'nats',
          sub: state.mStarMode === 'pipeline' ? `M* = ${fmtNum(mStar, 3)} m_P, S_BH(M*) = ${fmtNum(sStar, 6)}` : 'selected hole, Schwarzschild M*',
        }),
      ),
    );
  }

  function drawWindowChart(): void {
    const mStar = mStarNow();
    const sStar = 4 * Math.PI * mStar * mStar;
    const curve = windowCurve(mStar, X);
    winChart?.destroy();
    clear(winChartHost);
    const series: LineSeries[] = [
      { name: 'ln tr Π (exact, live)', x: curve.x, y: curve.lnTrace, color: 1 },
      { name: 'Theorem G(v) first-order form', x: curve.x, y: curve.theoremGv, color: 2, dash: 'dash' },
      { name: 'area term A(E₂)/4G − S_BH(M*)', x: curve.x, y: curve.deltaS, color: 'context' },
    ];
    if (state.mStarMode === 'pipeline' && cap?.window_examples?.length) {
      const ex = cap.window_examples.filter((w) => Number.isFinite(w.ln_trace));
      series.push({
        name: 'pipeline (qchas) ln tr Π',
        x: ex.map((w) => (w.beta ?? 8 * Math.PI * w.m_star) * w.width),
        y: ex.map((w) => w.ln_trace),
        color: 4,
        markers: true,
        width: 0,
      });
    }
    winChart = keep(
      ctx,
      lineChart(winChartHost, series, {
        xTitle: 'βΔE (window width in units of the Hawking temperature)',
        yTitle: 'entropy above S_BH(M*) (nats)',
        xType: 'log',
        height: 320,
        // label-less: on a log axis the shared refLine annotation would sit at 10^value
        refLines: [{ axis: 'x', value: state.xCap }],
        ariaLabel: 'Energy-window capacity against window width',
        layout: { xaxis: { dtick: 1 } } as Partial<Layout>,
        caption:
          state.mStarMode === 'pipeline'
            ? 'Orange markers: the Python pipeline’s window examples (horizon.json capacity.window_examples); lines: the TypeScript port, computed live. They coincide.'
            : `M* = ${fmtSci(mStar, 3)} m_P: the capacity is the area term minus ln β ≈ ${fmtNum(Math.log(8 * Math.PI * mStar), 3)} nats, negligible against S_BH(M*) = ${fmtPow10(Math.log10(sStar))}.`,
      }),
    );
  }

  // --- de Sitter and area gap ---------------------------------------------------------------------------------
  const dsP = panel({
    title: 'Capacity of our universe',
    subtitle: 'De Sitter entropy 3π/GΛ of the measured Λ, and the area gap',
    badges: [badge('proved', '[QG] Theorem C'), badge('computed')],
    span: 4,
    stretch: true,
  });
  dsP.root.classList.add('dense-table');
  if (hor.ok) {
    const ds = hor.data.capacity.de_sitter;
    const cap = hor.data.capacity;
    const live = deSitterEntropy();
    dsP.body.append(
      tilesN(
        2,
        statTile({
          label: 'De Sitter entropy S = 3π/GΛ',
          value: fmtPow10(ds.log10_entropy),
          unit: 'nats',
          badge: badge('computed', 'horizon.json'),
          sub: [`log₁₀ S = ${fmtNum(ds.log10_entropy, 6)}; live port ${fmtNum(live.log10_entropy, 6)}`],
        }),
        statTile({ label: 'Cosmological horizon √(3/Λ)', value: fmtNum(ds.horizon_radius_Gpc, 4), unit: 'Gpc', badge: badge('computed') }),
        statTile({
          label: 'Area gap a₀ = 4√3 πγ l_P²',
          value: fmtNum(cap.min_area_lp2, 5),
          unit: 'l_P²',
          badge: badge('computed', `γ = ${BARBERO_IMMIRZI_SU2.value}`),
          sub: [badge('quoted', MIN_AREA_LP2.source), ` paper ${MIN_AREA_LP2.value}`],
        }),
        statTile({ label: 'Cosmological constant Λ', value: fmtSci(ds.Lambda_m2 ?? live.Lambda_m2, 4), unit: 'm⁻²', badge: badge('computed') }),
      ),
      h(
        'div',
        { class: 'm2-callout' },
        kv([
          ['inputs', [`H₀ = ${PLANCK2018.H0} km/s/Mpc, Ω_Λ = ${fmtNum(ds.omega_lambda ?? live.omega_lambda, 4)} (Planck 2018) `, badge('assumed', 'Planck 2018 best fit')]],
          ['formula', [h('span', { class: 'm2-formula' }, 'S = 3π / (G Λ) = 3π / (Λ l_P²)'), ' ', badge('proved', '[QG] Theorem C')]],
        ]),
      ),
      h('p', { class: 'small muted', style: { margin: '10px 0 0' } }, 'The paper proves the formula; it prints no numerical value. The number is ours, from the Planck 2018 best fit.'),
    );
    const v = validationBlock(val, ['Crossed product ([QG] Theorem U)'], { compact: true });
    if (v) dsP.body.append(h('div', { class: 'm2-section-label', style: { marginTop: '14px' } }, 'Ours vs paper'), v);
  } else dsP.body.append(loadFailure(hor));

  el.append(grid(calcP.root), grid(winP.root, dsP.root));

  function update(): void {
    const e = events.find((x) => x.name === state.preset);
    const bh = blackHole(state.mass, state.spin, a0);
    const vsPipeline = (ours: number, theirs: number | null | undefined, fmt: (v: number) => string): Child =>
      e && theirs !== null && theirs !== undefined && Number.isFinite(theirs)
        ? h('span', { title: 'Value written by the Python pipeline (gw.json)' }, `gw.json ${fmt(theirs)}`, ' ', Math.abs(ours / theirs - 1) < 1e-6 ? statusPill('good', 'same') : statusPill('warning', 'differs'))
        : null;
    clear(calcTiles);
    calcTiles.append(
      tilesN(
        4,
        statTile({ label: 'Horizon radius r₊', value: fmtNum(bh.rPlusKm, 5), unit: 'km', badge: badge('live') }),
        statTile({ label: 'Horizon area A', value: fmtSci(bh.areaM2, 4), unit: 'm²', badge: badge('live') }),
        statTile({
          label: 'Capacity A/4G',
          value: fmtPow10(bh.log10EntropyNats, 4),
          unit: 'nats',
          badge: badge('live'),
          sub: vsPipeline(bh.entropyNats, e?.S_BH_nats ?? (e ? 10 ** e.S_BH_log10 : null), (v) => fmtPow10(Math.log10(v), 4)),
        }),
        statTile({ label: 'Capacity in bits', value: fmtPow10(bh.log10EntropyBits, 4), unit: 'bits', badge: badge('live'), sub: 'A/(4G ln 2)' }),
        statTile({
          label: 'Area quanta A/a₀',
          value: fmtPow10(bh.log10AreaQuanta, 4),
          badge: badge('live'),
          sub: [`a₀ = ${fmtNum(a0, 4)} l_P² `, badge('proved', '[QG] Theorem U')],
        }),
        statTile({
          label: 'Hawking temperature T_H',
          value: fmtSci(bh.temperatureK, 4),
          unit: 'K',
          badge: badge('live'),
          sub: vsPipeline(bh.temperatureK, e?.T_H_K, (v) => `${fmtSci(v, 4)} K`),
        }),
        statTile({
          label: 'Evaporation time t_evap',
          value: fmtPow10(bh.log10TEvapYr),
          unit: 'yr',
          badge: badge('assumed', 'massless scalar'),
          sub: vsPipeline(bh.log10TEvapYr, e?.t_evap_log10_yr, (v) => `${fmtPow10(v)} yr`),
        }),
        statTile({
          label: `Page time ${fmtNum(pageTimeFraction(), 3)} t_evap`,
          value: fmtPow10(bh.log10PageTimeYr),
          unit: 'yr',
          badge: badge('assumed', 'massless scalar'),
          sub: vsPipeline(bh.log10PageTimeYr, e?.page_time_log10_yr, (v) => `${fmtPow10(v)} yr`),
        }),
      ),
    );
    calcFoot.textContent =
      `Planck units: A = 8πM²(1 + √(1 − a*²)), S = A/4G (one quarter nat per Planck area), T_H = √(1 − a*²)/(4π r₊); ` +
      `lifetime t_evap = M³/(3α) with α = 7.44×10⁻⁵ ([BH] Sec. 17): massless-scalar Hawking emission of an isolated Schwarzschild hole (spin ignored). ` +
      (e
        ? `${e.name}: final mass and spin from aligned-spin remnant fits to the LIGO tutorial’s component parameters (${e.approximant ?? 'tutorial'} template), ${e.mass_frame ?? 'detector'} frame, so redshift is not removed.`
        : 'Custom hole: mass and spin as set above.');
    if (state.mStarMode === 'hole') {
      drawWindowTiles();
      // rebuilding the chart on every slider tick would block the main thread; debounce it
      clearTimeout(chartTimer);
      chartTimer = window.setTimeout(drawWindowChart, 180);
    }
  }
  let chartTimer = 0;
  ctx.onCleanup(() => clearTimeout(chartTimer));

  update();
  drawWindow();
}
