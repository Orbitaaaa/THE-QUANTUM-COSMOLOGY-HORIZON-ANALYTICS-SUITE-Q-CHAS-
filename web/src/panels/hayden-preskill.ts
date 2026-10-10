/**
 * 2.4 Hayden-Preskill Information Recovery Matrix ([QG] Theorem R, [BH] Theorem E(iv), Sec. 11).
 *
 * Decoding infidelity of a k = 2 qubit diary against the number s of further emitted qubits, for
 * Haar scrambling of an old hole, GUE dynamics e^{-iHt} at t = 1, 3, 10, 100 and a young hole,
 * with the 2-design bound 2^{1+k-s} and the paper's table; the same numbers as a recovery matrix;
 * a live decoding-budget calculator; and the diary mutual information of the 20-qubit run.
 * Bundle fields: horizon.decoding.{haar_old, eth, haar_young, paper}, horizon.evaporation
 * .{early_diary, late_diary, paper}; validation group "Hayden-Preskill decoding ([QG] Theorem R)".
 */

import './m2.css';
import './hayden-preskill.css';
import type { Layout } from 'plotly.js';
import { type ChartHandle, heatmapChart, scatterChart, type ScatterSeries } from '../lib/charts';
import type { DecodingRun, HorizonBundle } from '../lib/data';
import type { PanelContext } from '../lib/panel';
import { decodingInfidelityBound, haydenPreskillEmissions } from '../lib/physics/page';
import { badge, clear, fmtNum, fmtSci, grid, h, loadFailure, panel, segmented, slider, statTile } from '../lib/ui';
import { emissionsForBound, log10OrNull, positiveOrNull } from './horizon-calc';
import { diaryChart, keep, proofList, tilesN, validationBlock } from './m2-ui';

export const title = 'Hayden–Preskill Information Recovery Matrix';
export const subtitle =
  'Throw a diary into a black hole past its Page time: after a few more emitted qubits it can be decoded from the radiation. Exact decoding bounds for Haar and GUE scrambling against the theorem.';

interface Dyn {
  key: string;
  label: string;
  run: DecodingRun;
  paper: number[] | undefined;
  color: number;
}

function dynamics(d: HorizonBundle['decoding']): Dyn[] {
  const eth = Object.keys(d.eth).sort((a, b) => Number(a) - Number(b));
  return [
    { key: 'haar_old', label: 'Haar, old hole', run: d.haar_old, paper: d.paper?.haar_old, color: 1 },
    ...eth.map((t, i) => ({ key: `eth_${t}`, label: `GUE e^(−iHt), t = ${t}`, run: d.eth[t]!, paper: d.paper?.eth?.[t], color: 2 + i })),
    // slot 7 rather than 6: slot 6 (green) sits too close to slot 3 when every pair is on screen
    { key: 'haar_young', label: 'Haar, young hole', run: d.haar_young, paper: d.paper?.haar_young, color: 7 },
  ];
}

export async function mount(el: HTMLElement, ctx: PanelContext): Promise<void> {
  ctx.setHeaderExtras([badge('proved', '[QG] Theorem R'), badge('proved', '[BH] Theorem E(iv)'), badge('computed', 'qchas')]);
  const [hor, val] = await Promise.all([ctx.load('horizon.json'), ctx.load('validation.json')]);
  if (ctx.unmounted) return;

  const curveP = panel({
    title: 'Decoding infidelity vs emitted qubits',
    subtitle: 'k = 2 diary qubits thrown into a 7-qubit black hole; mean over 6 samples per entry',
    badges: [badge('computed', 'horizon.json decoding'), badge('quoted', '[QG] Sec. 15 table')],
    span: 8,
  });
  const provesP = panel({ title: 'What the paper proves', subtitle: '[QG] Secs. 15, 20 and [BH] Secs. 7, 11', span: 4, stretch: true });
  provesP.body.append(
    proofList(
      [
        { kind: 'proved', src: '[QG] Theorem R', text: 'If the dynamics is a unitary 2-design, the average decoding infidelity of a k-qubit diary after s further emissions obeys E[1 − F_dec] ≤ 2^(1+k−s): the interior comes out as soon as s exceeds k by a few qubits.' },
        { kind: 'proved', src: '[BH] Theorem E(iv)', text: 'A diary of k qudits thrown in after the Page time is decoded from the radiation after about k + log₂(1/ε) further emissions (decoupling for 2-designs).' },
        { kind: 'computed', src: '[QG] Sec. 15', text: 'Tested with GUE dynamics e^(−iHt): from t = 10 it is indistinguishable from Haar; at t ≤ 1 the diary has not yet spread (the scrambling time). A young hole (before its Page time) releases the diary only after about half of it has been emitted.' },
        { kind: 'computed', src: '[BH] Sec. 11', text: 'In the exactly unitary 20-qubit evaporation, a diary thrown in at step 14 is recovered at step 17, three emissions later; one thrown in at step 5 waits for the Page time.' },
      ],
      true,
    ),
  );
  el.append(grid(curveP.root, provesP.root));

  if (!hor.ok) {
    curveP.body.append(loadFailure(hor));
    return;
  }
  const d = hor.data.decoding;
  const dyn = dynamics(d);
  const k = d.haar_old.k;
  const nbh = d.haar_old.n_bh;
  const s = d.haar_old.s;
  const sPaper = d.paper?.s ?? s;

  // --- curves ------------------------------------------------------------------------------------------
  const series: ScatterSeries[] = [
    ...dyn.map((x) => ({
      name: x.label,
      x: x.run.s,
      y: x.run.infidelity.map((v) => positiveOrNull(v)),
      error: x.run.std?.map((e, i) => (positiveOrNull(x.run.infidelity[i]) === null ? 0 : e)),
      color: x.color,
      mode: 'lines+markers' as const,
      size: 7,
    })),
    { name: `bound 2^(1+k−s), k = ${k}`, x: s, y: s.map((v) => decodingInfidelityBound(k, v)), color: 'context', mode: 'lines+markers', size: 6, symbol: 'square-open' },
  ];
  const paperX: number[] = [];
  const paperY: (number | null)[] = [];
  const paperText: string[] = [];
  for (const x of dyn)
    x.paper?.forEach((v, i) => {
      const y = positiveOrNull(v);
      if (y === null) return;
      paperX.push(sPaper[i] ?? i);
      paperY.push(y);
      paperText.push(`${x.label}, paper`);
    });
  if (paperX.length) series.push({ name: '[QG] Sec. 15 table (all rows)', x: paperX, y: paperY, text: paperText, color: '--text-2', symbol: 'circle-open', size: 11 });
  keep(
    ctx,
    scatterChart(curveP.body, series, {
      xTitle: 'further emitted qubits s',
      yTitle: 'decoding infidelity 1 − F (upper bound on 1 − F_dec)',
      yType: 'log',
      yRange: [2e-5, 1.6],
      height: 400,
      layout: { yaxis: { dtick: 1 }, xaxis: { dtick: 1 } } as Partial<Layout>,
      format: (v) => fmtSci(v, 3),
      ariaLabel: 'Decoding infidelity against emitted qubits for each dynamics, log scale',
      caption: `At s = ${s[s.length - 1]} every qubit of the ${nbh}-qubit hole has been emitted and the infidelity is zero to machine precision (|value| < 10⁻¹⁵, ours and the paper’s), so that column is off this log scale. Error bars: sample standard deviation over ${d.haar_old.samples ?? 6} draws.`,
    }),
  );

  // --- recovery matrix ------------------------------------------------------------------------------------
  const matP = panel({
    title: 'Recovery matrix',
    subtitle: 'Rows: dynamics; columns: emitted qubits s; colour: recovery depth −log₁₀(1 − F), stronger colour = better recovered',
    badges: [badge('computed'), badge('quoted', '[QG] Sec. 15')],
    span: 7,
  });
  const sCols = s.filter((_, i) => dyn.some((x) => positiveOrNull(x.run.infidelity[i]) !== null));
  type Mode = 'ours' | 'paper' | 'ratio';
  let mode: Mode = 'ours';
  const seg = segmented({
    label: 'Matrix values',
    options: [
      { value: 'ours', label: 'Q-CHAS' },
      { value: 'paper', label: 'Paper' },
      { value: 'ratio', label: 'log₁₀ ours / paper' },
    ],
    value: mode,
    onChange: (v) => ((mode = v as Mode), drawMatrix()),
  });
  const matHost = h('div');
  matP.body.append(h('div', { class: 'm2-pill-row' }, seg.root), matHost);
  let mat: ChartHandle | null = null;
  function drawMatrix(): void {
    const z = dyn.map((x) =>
      sCols.map((sv) => {
        const i = s.indexOf(sv);
        const j = sPaper.indexOf(sv);
        const ours = log10OrNull(x.run.infidelity[i]);
        const paper = log10OrNull(x.paper?.[j]);
        // recovery depth -log10(1 - F): more is better recovered, so more colour
        if (mode === 'ours') return ours === null ? null : -ours;
        if (mode === 'paper') return paper === null ? null : -paper;
        return ours === null || paper === null ? null : ours - paper;
      }),
    );
    mat?.destroy();
    clear(matHost);
    mat = keep(
      ctx,
      heatmapChart(
        matHost,
        {
          z,
          x: sCols,
          y: dyn.map((x) => x.label),
          scale: mode === 'ratio' ? 'diverging' : 'sequential',
          zmid: mode === 'ratio' ? 0 : undefined,
          zmin: mode === 'ratio' ? undefined : 0,
          zmax: mode === 'ratio' ? undefined : 5,
          colorbarTitle: mode === 'ratio' ? 'log₁₀ ratio' : '−log₁₀(1 − F)',
          zLabel: mode === 'ratio' ? 'log₁₀ ours/paper' : '−log₁₀(1 − F)',
        },
        {
          xTitle: 'further emitted qubits s',
          yTitle: 'dynamics',
          height: 320,
          format: (v) => fmtNum(v, 3),
          // categorical rows: the shared axis defaults to 'linear', which cannot place strings
          layout: { yaxis: { type: 'category', autorange: 'reversed' }, xaxis: { dtick: 1 }, margin: { l: 150, r: 18, t: 12, b: 46, pad: 0 } } as Partial<Layout>,
          ariaLabel: 'Recovery matrix of log10 decoding infidelity by dynamics and emitted qubits',
          caption:
            mode === 'ratio'
              ? 'Ours against the printed table: within the sampling noise of 6 draws (~10–20 % at small infidelity), i.e. |log₁₀ ratio| ≲ 0.1.'
              : 'Recovery depth 5 means infidelity 10⁻⁵ (colour clamped there); hover or use Table for the values. Haar and GUE from t = 10 recover the diary within 2–3 qubits beyond k; t = 1 stalls near 0.55; a young hole needs about half of itself emitted.',
        },
      ),
    );
  }
  drawMatrix();

  // --- live budget ------------------------------------------------------------------------------------------
  const budP = panel({
    title: 'Decoding budget',
    subtitle: 'How many further emissions until a k-qubit diary is decodable to infidelity ε',
    badges: [badge('live', 'physics/page.ts')],
    span: 5,
    stretch: true,
  });
  const bud = { k: 2, eps: 0.01 };
  const budTiles = h('div');
  const kS = slider({ label: 'Diary size k (qubits)', min: 1, max: 40, step: 1, value: bud.k, format: (v) => String(v), onInput: (v) => ((bud.k = v), drawBudget()) });
  const eS = slider({ label: 'Target infidelity ε', min: -12, max: -1, step: 0.1, value: -2, format: (v) => fmtSci(10 ** v, 2), onInput: (v) => ((bud.eps = 10 ** v), drawBudget()) });
  budP.body.append(h('div', { class: 'm2-controls' }, kS.root, eS.root), budTiles);
  function drawBudget(): void {
    clear(budTiles);
    budTiles.append(
      tilesN(
        2,
        statTile({
          label: 'Hayden–Preskill estimate k + log₂(1/ε)',
          value: fmtNum(haydenPreskillEmissions(bud.k, bud.eps), 4),
          unit: 'emissions',
          badge: badge('proved', '[BH] Theorem E(iv)'),
        }),
        statTile({
          label: 'Theorem R: 2^(1+k−s) ≤ ε once s ≥',
          value: emissionsForBound(bud.k, bud.eps),
          unit: 'emissions',
          badge: badge('proved', '[QG] Theorem R'),
          sub: `bound there: ${fmtSci(decodingInfidelityBound(bud.k, emissionsForBound(bud.k, bud.eps)), 3)}`,
        }),
      ),
      h(
        'p',
        { class: 'small muted', style: { margin: '10px 0 0' } },
        'Both hold for a hole past its Page time with 2-design scrambling; the hole must hold at least that many more qubits. Before the Page time the diary waits for it (young hole, early diary).',
      ),
    );
  }
  drawBudget();
  const vb = validationBlock(val, ['Hayden-Preskill decoding ([QG] Theorem R)'], { summaryOnly: true });
  if (vb)
    budP.body.append(
      h('div', { class: 'm2-section-label', style: { marginTop: '16px' } }, 'Ours vs paper · Haar, old hole, s = 0–8'),
      vb,
      h('p', { class: 'small muted', style: { margin: '6px 0 0' } }, 'Each entry is a mean of 6 Haar samples, as in the paper, so the sampling noise is ~10–20 % at the smallest infidelities; the recovery matrix’s log₁₀ ours/paper view shows every entry.'),
    );
  el.append(grid(matP.root, budP.root));

  // --- diary in the black hole -------------------------------------------------------------------------------
  const ev = hor.data.evaporation;
  const diaryP = panel({
    title: 'A diary in the black hole',
    subtitle: 'Exactly unitary evaporation of 20 qubits: diary information in the radiation, I(Ref : R) in bits (2 = complete)',
    badges: [badge('computed', '[BH] Sec. 11'), badge('quoted', '[BH] Fig. 2')],
    span: 7,
  });
  keep(ctx, diaryChart(diaryP.body, ev, 300));
  diaryP.body.append(
    h(
      'div',
      { class: 'chart-caption' },
      `Ours: early diary (in at ${ev.early_diary.diary_step}) recovered at step ${ev.early_diary.recovered_step ?? '—'}, late diary (in at ${ev.late_diary.diary_step}) at ${ev.late_diary.recovered_step ?? '—'}; paper: ${ev.paper.early_recovered_step} and ${ev.paper.late_recovered_step}. Horizontal line: the recovery threshold of ${ev.paper.recovery_threshold_bits ?? 1.9} bits.`,
    ),
  );
  const expP = panel({ title: 'How the mirror works', subtitle: 'The protocol behind both computations', span: 5, stretch: true });
  const steps: [string, string][] = [
    ['Old hole', 'Past its Page time the black hole is maximally entangled with the radiation it has already emitted.'],
    ['Throw in a diary', `k qubits (here k = ${k}), each maximally entangled with an outside reference, fall in.`],
    ['Scramble', 'The hole’s dynamics acts as a unitary 2-design (Haar, or GUE e^(−iHt) after the scrambling time).'],
    ['Emit s more qubits', 'Decoupling: the reference becomes uncorrelated with the remaining hole, so the early radiation plus the s new qubits hold the diary, up to 2^(1+k−s).'],
    ['Young hole or early diary', 'Before the Page time the radiation is not yet entangled with the hole; the diary comes out only after about half of the hole has been emitted.'],
  ];
  expP.body.append(
    h(
      'ol',
      { class: 'hp-steps' },
      ...steps.map(([t, x]) => h('li', null, h('strong', null, `${t}. `), x)),
    ),
    h('p', { class: 'small muted', style: { margin: '10px 0 0' } }, 'This is the decoupling inequality for 2-designs. It assumes H* (no interior independent of the hole and the radiation), which is what makes the radiation able to hold the diary at all.'),
  );
  el.append(grid(diaryP.root, expP.root));
}
