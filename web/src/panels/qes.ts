/**
 * 2.3 Quantum Extremal Surface & Atmosphere Locator ([BH] Secs. 9, 10, 15).
 *
 * A live QES solve (physics/qes.ts solveQes) for any mass, central charge and collection
 * radius, the [BH] Sec. 9 fit laws, a schematic radial diagram, the 12-row table ours vs paper,
 * and the quantum atmosphere: the share of the entanglement beyond g cutoff lengths.
 * Bundle fields: horizon.qes.{rows, laws}, horizon.atmosphere.{gaps, share, paper_share,
 * area_coeff, continuum, n_sites, n_ball, ell_max}; validation groups "QES ([BH] Sec. 9)" and
 * "Quantum atmosphere ([BH] Sec. 10)".
 */

import './m2.css';
import './qes.css';
import { scatterChart } from '../lib/charts';
import type { QesRow } from '../lib/data';
import type { Layout } from 'plotly.js';
import type { PanelContext } from '../lib/panel';
import { QES_DEPTH_COEFF, QES_LAG_OFFSET_RS } from '../lib/physics/constants';
import { type QesResult, depthFit, lagFit, solveQes } from '../lib/physics/qes';
import { append, badge, clear, fmtNum, fmtSci, grid, h, loadFailure, panel, slider, statTile, statusPill, table } from '../lib/ui';
import { keep, proofList, tilesN, validationBlock } from './m2-ui';

export const title = 'Quantum Extremal Surface & Atmosphere Locator';
export const subtitle =
  'Where the island ends: the quantum extremal surface sits just inside the horizon and trails the radiation by a scrambling time; the entanglement it counts lives in a Planck-thin atmosphere.';

const NS = 'http://www.w3.org/2000/svg';
function svgEl<K extends keyof SVGElementTagNameMap>(tag: K, attrs: Record<string, string | number>, text?: string): SVGElementTagNameMap[K] {
  const e = document.createElementNS(NS, tag);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, String(v));
  if (text !== undefined) e.textContent = text;
  return e;
}

/** Schematic radial diagram: island [0, a], QES just inside r_s, gap, radiation [r_b, inf). */
function radialDiagram(r: QesResult): SVGSVGElement {
  const W = 680;
  const H = 190;
  const y = 82;
  const x0 = 34;
  const xs = 344; // horizon
  const xa = 326; // QES (depth exaggerated)
  const xb = 492; // r_b
  const xEnd = 650;
  const svg = svgEl('svg', { viewBox: `0 0 ${W} ${H}`, class: 'm2-svg qes-diagram', role: 'img', 'aria-label': 'Schematic radial diagram of the island, the QES, the horizon and the radiation region' });
  // regions
  svg.append(
    svgEl('rect', { x: x0, y: y - 26, width: xs - x0, height: 52, class: 'qes-interior', rx: 3 }),
    svgEl('rect', { x: x0, y: y - 26, width: xa - x0, height: 52, class: 'qes-island', rx: 3 }),
    svgEl('rect', { x: xb, y: y - 26, width: xEnd - xb, height: 52, class: 'qes-rad', rx: 3 }),
  );
  // axis with a break between the horizon and r_b
  svg.append(
    svgEl('line', { x1: x0, y1: y, x2: 410, y2: y, class: 'qes-axis' }),
    svgEl('line', { x1: 426, y1: y, x2: xEnd + 14, y2: y, class: 'qes-axis' }),
    svgEl('path', { d: `M404 ${y + 9} l8 -18 M420 ${y + 9} l8 -18`, class: 'qes-axis' }),
    svgEl('path', { d: `M${xEnd + 14} ${y} l-7 -4 v8 z`, class: 'qes-axis-head' }),
  );
  // markers
  svg.append(
    svgEl('line', { x1: xs, y1: y - 40, x2: xs, y2: y + 40, class: 'qes-horizon' }),
    svgEl('line', { x1: xa, y1: y - 40, x2: xa, y2: y + 40, class: 'qes-qes' }),
    svgEl('line', { x1: xb, y1: y - 40, x2: xb, y2: y + 40, class: 'qes-rb' }),
    svgEl('circle', { cx: x0, cy: y, r: 3.5, class: 'qes-origin' }),
  );
  // labels
  const t = (x: number, yy: number, s: string, cls = '', anchor = 'middle'): SVGTextElement =>
    svgEl('text', { x, y: yy, 'text-anchor': anchor, class: cls }, s);
  svg.append(
    t(x0, y + 44, 'r = 0', 'qes-small', 'middle'),
    t((x0 + xa) / 2, y + 5, 'island I = [0, a]', 'qes-region'),
    t(xa - 4, y - 48, 'QES  r = r_s(1 − δ)', 'qes-label qes-label-qes', 'end'),
    t(xs + 4, y - 48, 'horizon r_s', 'qes-label', 'start'),
    t(xb, y + 56, `r_b = ${fmtNum(r.rb_over_rs, 3)} r_s`, 'qes-label', 'middle'),
    t((xb + xEnd) / 2, y + 5, 'radiation R = [r_b, ∞)', 'qes-region'),
    t(xEnd + 10, y + 22, '∞', 'qes-small', 'end'),
  );
  // depth callout
  svg.append(
    svgEl('path', { d: `M${xa} ${y + 32} v8 h${xs - xa} v-8`, class: 'qes-bracket' }),
    t((xa + xs) / 2, y + 56, `δ = ${fmtSci(r.depth, 3)}`, 'qes-label mono', 'middle'),
    t(x0, H - 12, `Schematic, not to scale: δ is drawn ${fmtSci((xs - xa) / (xs - x0) / r.depth, 2)}× enlarged. The QES trails the radiation by u_b − v_QES = ${fmtNum(r.lag, 6)} (t_scr = ${fmtNum(r.t_scr, 6)}).`, 'qes-small', 'start'),
  );
  return svg;
}

export async function mount(el: HTMLElement, ctx: PanelContext): Promise<void> {
  ctx.setHeaderExtras([badge('computed', '[BH] Secs. 9, 10, 15'), badge('proved', '[BH] Theorem D'), badge('live')]);
  const [hor, val] = await Promise.all([ctx.load('horizon.json'), ctx.load('validation.json')]);
  if (ctx.unmounted) return;

  // --- live solver -----------------------------------------------------------------------------------
  const solverP = panel({
    title: 'Live QES solver',
    subtitle: 'Extremize S_gen(a) = πr_a²/G + (c/6) ln[(U_a − U_b)(v_b − v_a) e^{ρ_a+ρ_b}/ε²] over the island endpoint (Planck units, ε = l_P)',
    badges: [badge('live', 'physics/qes.ts'), badge('computed', '[BH] Sec. 9')],
    span: 7,
  });
  const st = { m: 10, c: 1, rb: 20 };
  const mSlider = slider({ label: 'Mass M / m_P (log scale)', min: 1, max: 3, step: 0.005, value: 1, format: (v) => fmtNum(10 ** v, 4), onInput: (v) => ((st.m = 10 ** v), schedule()) });
  const cSlider = slider({ label: 'Central charge c (s-wave matter)', min: 1, max: 30, step: 1, value: 1, format: (v) => String(v), onInput: (v) => ((st.c = v), schedule()) });
  const rbSlider = slider({ label: 'Collection radius r_b / r_s', min: 20, max: 60, step: 1, value: 20, format: (v) => String(v), onInput: (v) => ((st.rb = v), schedule()) });
  const solveTiles = h('div');
  const lawTiles = h('div', { style: { marginTop: '10px' } });
  const statusRow = h('div', { class: 'm2-pill-row', style: { marginTop: '10px' } });
  const diagramHost = h('div', { class: 'qes-diagram-host' });
  solverP.body.append(h('div', { class: 'm2-controls' }, mSlider.root, cSlider.root, rbSlider.root), solveTiles, lawTiles, statusRow, diagramHost);

  const provesP = panel({ title: 'What the paper proves', subtitle: '[BH] Secs. 6, 9, 10, 15 and the status table of Sec. 16', span: 5, stretch: true });
  provesP.body.append(
    proofList(
      [
        { kind: 'proved', src: '[BH] Theorem D', text: 'The island is forced once S_H > A/4G, with minimal size S_H − A/4G: the interior partners of the early radiation must be encoded in the radiation. Given capacity and H*.' },
        { kind: 'computed', src: '[BH] Sec. 9', text: 'Its location: the QES lies just inside the horizon, at a depth set by c/S_BH, and trails the radiation by the scrambling time plus ~4.24 r_s. S_gen = S_BH + O(c), so the island saddle reproduces the capacity branch of Theorem E.' },
        { kind: 'computed', src: '[BH] Secs. 10, 15', text: 'The atmosphere: for a free field on a radial lattice the entanglement across a sphere lives within a few cutoff lengths, with a universal profile f(g) in the continuum limit.' },
        { kind: 'assumed', src: '[BH] Sec. 9', text: 'Matter is modelled as an s-wave field of central charge c in the in-vacuum; the island rule itself is derived by the gravitational path integral (replica wormholes), not proved.' },
      ],
      true,
    ),
  );
  el.append(grid(solverP.root, provesP.root));

  let raf = 0;
  function schedule(): void {
    cancelAnimationFrame(raf);
    raf = requestAnimationFrame(render);
  }
  ctx.onCleanup(() => cancelAnimationFrame(raf));
  function render(): void {
    const r = solveQes(st.m, st.c, st.rb);
    const dFit = depthFit(st.c, r.s_bh);
    const lFit = lagFit(st.m, st.c);
    const rs = 2 * st.m;
    clear(solveTiles);
    solveTiles.append(
      tilesN(
        4,
        statTile({ label: 'Depth (r_s − r_QES)/r_s', value: fmtSci(r.depth, 4), badge: badge('live') }),
        statTile({ label: 'Lag u_b − v_QES', value: fmtNum(r.lag, 6), unit: 't_P', badge: badge('live'), sub: `= ${fmtNum(r.lag / rs, 4)} r_s` }),
        statTile({ label: 'Scrambling time t_scr', value: fmtNum(r.t_scr, 6), unit: 't_P', badge: badge('live'), sub: '(1/2πT) ln(S_BH/c)' }),
        statTile({ label: 'S_gen − S_BH', value: fmtNum(r.s_gen_minus_s_bh, 4), unit: 'nats', badge: badge('live'), sub: `S_BH = ${fmtNum(r.s_bh, 6)}` }),
      ),
    );
    clear(lawTiles);
    lawTiles.append(
      tilesN(
        2,
        statTile({
          label: `Depth law ${QES_DEPTH_COEFF.value} c / S_BH`,
          value: fmtSci(dFit, 4),
          badge: badge('quoted', QES_DEPTH_COEFF.source),
          sub: `solver / law = ${fmtNum(r.depth / dFit, 4)}; coefficient 0.0399 ± 0.0007`,
        }),
        statTile({
          label: `Lag law (1/2πT) ln(S_BH/c) + ${QES_LAG_OFFSET_RS.value} r_s`,
          value: fmtNum(lFit, 6),
          badge: badge('quoted', QES_LAG_OFFSET_RS.source),
          sub: `solver − law = ${fmtNum((r.lag - lFit) / rs, 3)} r_s (offset quoted at r_b = 20 r_s)`,
        }),
      ),
    );
    clear(statusRow);
    append(
      statusRow,
      r.converged ? statusPill('good', `converged · |∇S_gen| = ${fmtSci(r.residual, 2)}`) : statusPill('warning', `not converged · residual ${fmtSci(r.residual, 2)}`),
      h('span', { class: 'small muted' }, `${r.iterations} Levenberg–Marquardt iterations on the analytic gradient`),
      r.s_bh < 30 * st.c ? statusPill('warning', 'S_BH not ≫ c: the fit laws are asymptotic') : null,
    );
    clear(diagramHost);
    diagramHost.append(radialDiagram(r));
  }
  render();

  // --- paper table --------------------------------------------------------------------------------------
  const tableP = panel({
    title: 'The island located · [BH] Sec. 9 table',
    subtitle: 'Our solver (Python pipeline) against every printed digit; click a row to load it into the live solver',
    badges: [badge('computed', 'horizon.json qes.rows'), badge('quoted', '[BH] Sec. 9')],
    flush: true,
  });
  if (!hor.ok) {
    tableP.body.style.padding = '14px';
    tableP.body.append(loadFailure(hor));
  } else {
    const rows = hor.data.qes.rows;
    const pair = (ours: number, paper: number | undefined, f: (v: number) => string) => [
      h('span', { class: 'm2-inline-num' }, f(ours)),
      h('div', { class: 'small muted' }, paper === undefined ? '—' : `paper ${f(paper)}`),
    ];
    const f2 = (v: number): string => v.toFixed(2);
    const fe = (v: number): string => v.toExponential(2);
    const tbl = table<QesRow>(
      [
        { key: 'mass', label: 'M / m_P', align: 'num' },
        { key: 'c', label: 'c', align: 'num' },
        { key: 'rb_over_rs', label: 'r_b / r_s', align: 'num' },
        { key: 's_bh', label: 'S_BH', align: 'num', format: (r) => fmtNum(r.s_bh ?? 4 * Math.PI * r.mass ** 2, 6) },
        { key: 'depth', label: '(r_s − r_QES)/r_s', align: 'num', format: (r) => pair(r.depth, r.paper.depth, fe) },
        { key: 'lag', label: 'lag u_b − v_QES', align: 'num', format: (r) => pair(r.lag, r.paper.lag, f2) },
        { key: 't_scr', label: 't_scr', align: 'num', format: (r) => pair(r.t_scr, r.paper.t_scr, f2) },
        { key: 'ds', label: 'S_gen − S_BH', align: 'num', format: (r) => pair(r.s_gen_minus_s_bh, r.paper.s_gen_minus_s_bh, (v) => (v >= 0 ? '+' : '') + f2(v)) },
        {
          key: 'ok',
          label: 'Printed digits',
          format: (r) => {
            const ok =
              (r.paper.depth === undefined || fe(r.depth) === fe(r.paper.depth)) &&
              (r.paper.lag === undefined || f2(r.lag) === f2(r.paper.lag)) &&
              (r.paper.t_scr === undefined || f2(r.t_scr) === f2(r.paper.t_scr)) &&
              (r.paper.s_gen_minus_s_bh === undefined || f2(r.s_gen_minus_s_bh) === f2(r.paper.s_gen_minus_s_bh));
            return ok ? statusPill('good', 'all match') : statusPill('warning', 'last digit');
          },
        },
      ],
      rows,
    );
    tbl.querySelectorAll('tbody tr').forEach((tr, i) => {
      const r = rows[i];
      if (!r) return;
      tr.classList.add('qes-row');
      tr.setAttribute('tabindex', '0');
      tr.setAttribute('title', 'Load into the live solver');
      const load = (): void => {
        st.m = r.mass;
        st.c = r.c;
        st.rb = r.rb_over_rs;
        mSlider.set(Math.log10(r.mass));
        cSlider.set(Math.min(30, r.c));
        rbSlider.set(r.rb_over_rs);
        render();
        solverP.root.scrollIntoView({ behavior: 'smooth', block: 'start' });
      };
      tr.addEventListener('click', load);
      tr.addEventListener('keydown', (e) => {
        if ((e as KeyboardEvent).key === 'Enter') load();
      });
    });
    const v = validationBlock(val, ['QES ([BH] Sec. 9)'], { summaryOnly: true });
    append(tableP.body, v ? h('div', { style: { padding: '12px 14px 4px' } }, v) : null, tbl);
    const laws = hor.data.qes.laws;
    tableP.body.append(
      h(
        'div',
        { class: 'chart-caption', style: { padding: '0 14px 12px' } },
        `Fit laws across the table ([BH] Sec. 9): (r_s − r_QES)/r_s = ${laws.depth_coeff ?? QES_DEPTH_COEFF.value} (±${laws.depth_err ?? 0.0007}) c/S_BH and u_b − v_QES = (1/2πT) ln(S_BH/c) + ${laws.lag_offset_rs ?? QES_LAG_OFFSET_RS.value} r_s. "Printed digits" compares our value rounded to the paper’s printed precision.`,
      ),
    );
  }
  el.append(grid(tableP.root));

  // --- quantum atmosphere ---------------------------------------------------------------------------------
  const atmP = panel({
    title: 'Quantum atmosphere · [BH] Sec. 10',
    subtitle: 'Share of the entanglement 2S between a ball and the exterior beyond a gap of g cutoff lengths',
    badges: [badge('computed', 'Srednicki lattice'), badge('quoted', '[BH] Secs. 10, 15')],
    span: 7,
  });
  const atmSide = panel({ title: 'Ours vs paper vs continuum', subtitle: 'Radial lattice of 60 sites, ball of 20, partial waves ℓ ≤ 3000', span: 5, stretch: true });
  if (!hor.ok) {
    atmP.body.append(loadFailure(hor));
    atmSide.body.append(loadFailure(hor));
  } else {
    const a = hor.data.atmosphere;
    const cont = a.continuum;
    const contPts = [1, 2, 3].map((g) => ({ g, v: cont[`f${g}`] })).filter((p): p is { g: number; v: number } => typeof p.v === 'number');
    const paperPts = a.gaps.map((g, i) => ({ g, v: a.paper_share[i] })).filter((p): p is { g: number; v: number } => typeof p.v === 'number');
    keep(
      ctx,
      scatterChart(
        atmP.body,
        [
          { name: 'Q-CHAS lattice (computed)', x: a.gaps, y: a.share, color: 1, mode: 'lines+markers' },
          { name: '[BH] Sec. 10 table', x: paperPts.map((p) => p.g), y: paperPts.map((p) => p.v), color: '--text-2', symbol: 'circle-open', size: 11 },
          { name: 'continuum limit f(g), [BH] Sec. 15', x: contPts.map((p) => p.g), y: contPts.map((p) => p.v), color: 3, symbol: 'diamond', size: 10 },
        ],
        {
          xTitle: 'gap g (cutoff lengths)',
          yTitle: 'share of the entanglement beyond the gap',
          yType: 'log',
          height: 330,
          layout: { yaxis: { dtick: 1 } } as Partial<Layout>,
          format: (v) => fmtNum(v, 4),
          ariaLabel: 'Share of entanglement beyond a gap of g cutoff lengths, log scale',
          caption: 'The share falls by a factor ~6 within the first cutoff length and by ~1000 within twelve: the entanglement across the surface is concentrated in a layer a few cutoff lengths thick.',
        },
      ),
    );
    const rows = a.gaps.map((g, i) => ({ g, ours: a.share[i]!, paper: a.paper_share[i] ?? null, cont: g >= 1 && g <= 3 ? (cont[`f${g}`] as number | undefined) ?? null : null }));
    const pct = (v: number | null): string => (v === null ? '—' : `${(100 * v).toFixed(v >= 0.01 ? 2 : 3)} %`);
    atmSide.body.append(
      tilesN(
        2,
        statTile({
          label: 'Area-law coefficient S/(n + ½)²',
          value: fmtNum(a.area_coeff, 4),
          badge: badge('computed'),
          sub: [badge('quoted', '[BH] Sec. 10'), ' paper 0.295 (Srednicki 0.30)'],
        }),
        statTile({
          label: 'Continuum limit S/(R/a)²',
          value: typeof cont.area_coeff === 'number' ? fmtNum(cont.area_coeff, 4) : '—',
          badge: badge('quoted', String(cont.source ?? '[BH] Sec. 15')),
          sub: 'resolution-independent',
        }),
      ),
      h(
        'div',
        { style: { marginTop: '12px' } },
        table(
          [
            { key: 'g', label: 'Gap g', align: 'num' },
            { key: 'ours', label: 'Q-CHAS', align: 'num', format: (r) => pct(r.ours) },
            { key: 'paper', label: 'Paper Sec. 10', align: 'num', format: (r) => pct(r.paper) },
            { key: 'cont', label: 'f(g), n → ∞', align: 'num', format: (r) => pct(r.cont) },
          ],
          rows,
        ),
      ),
    );
    const v = validationBlock(val, ['Quantum atmosphere ([BH] Sec. 10)'], { summaryOnly: true });
    if (v) atmSide.body.append(h('div', { style: { marginTop: '12px' } }, v));
    atmSide.body.append(
      h(
        'div',
        { class: 'callout m2-callout m2-takeaway' },
        h('strong', null, 'A surface density: one quarter nat per Planck area. '),
        'Identifying the cutoff through S = A/4G (’t Hooft’s brick wall, renormalized by Susskind–Uglum) makes this layer Planck-thin. The horizon’s capacity is therefore set by its area, which is why the Page transfer law (Theorem E) is controlled by A/4G.',
        h('div', { class: 'muted', style: { marginTop: '6px' } }, 'At a fixed physical gap the mutual information converges while 2S grows as (R/a)²: the share beyond any fixed distance vanishes as a² ([BH] Sec. 15).'),
      ),
    );
  }
  el.append(grid(atmP.root, atmSide.root));
}
