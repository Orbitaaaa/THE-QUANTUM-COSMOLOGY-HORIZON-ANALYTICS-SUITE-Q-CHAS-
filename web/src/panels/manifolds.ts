/**
 * 1.1 Manifold Selector & Principle Filter ([TOP] Secs. 1-5, 9.2-9.4).
 *
 * The 18 complete flat 3-manifolds against principles P1-P4: toggling principles recomputes the
 * survivors live (m1-logic.survivors, a port of qchas.topology.manifolds.survivors), names the
 * principle that excludes each manifold, and shows the inclusion-minimal selecting sets. Theorem
 * B is made visible with the E1 / E2 vacuum-variance slices; the curved candidates, the Weeks
 * table and the Planck evidence tables are quoted with their sources.
 *
 * Bundle fields: topology.{manifolds, principles, curved, weeks, evidence_torus,
 * exclusion_torus, other_flat, vacuum.{E1, E2, E1_smoothed, E2_smoothed, E2_axes, note}};
 * cmb.likelihood (optional: our WMAP re-analysis next to Table 9a when it has been computed).
 */

import './m1.css';
import { heatmapChart } from '../lib/charts';
import type { CmbBundle, EvidenceTorusEntry, ExclusionTorusEntry, ManifoldEntry, TopologyBundle, VacuumSlice } from '../lib/data';
import type { PanelContext } from '../lib/panel';
import { cssVar } from '../lib/theme';
import {
  badge,
  clear,
  fmtFixed,
  fmtNum,
  grid,
  h,
  loadFailure,
  note,
  panel,
  segmented,
  statTile,
  statusPill,
} from '../lib/ui';
import {
  type PrincipleKey,
  PRINCIPLE_KEYS,
  excludedBy,
  failedPrinciples,
  failingExactly,
  fmtSet,
  minimalSelectingSets,
  survivors,
} from './m1-logic';
import { type RawColumn, keep, proofList, rawTable, sectionLabel, tilesN } from './m1-ui';

export const title = 'Manifold Selector & Principle Filter';
export const subtitle =
  'Which of the 18 flat 3-manifolds survive principles P1–P4, and by which theorem each one is excluded. Toggle the principles; the survivors are recomputed live.';

const P_HEAD: Record<PrincipleKey, string> = {
  P1: 'P1 orientable',
  P2: 'P2 closed',
  P3: 'P3 local homogeneity',
  P4: 'P4 parallel spinor',
};

/** Paper's status tag ([TOP] Sec. 3) rendered with the matching provenance colour. */
function statusTag(status: string): HTMLElement {
  const s = status.toLowerCase();
  const kind = s.startsWith('experimentally') ? 'measured' : s.startsWith('speculative') ? 'open' : 'assumed';
  return h('span', { class: `badge badge-${kind} m1-tag`, title: `[TOP] Sec. 3 status tag: ${status}` }, status);
}

function principleCell(status: string | undefined): HTMLElement {
  if (status === 'allowed') return statusPill('good', 'allowed');
  if (status === 'no spin structure') return statusPill('serious', 'no spin structure');
  return statusPill('critical', status ?? '—');
}

const yesNo = (b: boolean): HTMLElement => h('span', { class: b ? 'm1-yes' : 'm1-no' }, b ? 'yes' : 'no');

export async function mount(el: HTMLElement, ctx: PanelContext): Promise<void> {
  ctx.setHeaderExtras([badge('proved', '[TOP] Main Theorem'), badge('quoted', '[TOP] Sec. 9'), badge('live', 'm1-logic.ts')]);
  const [topo, cmb] = await Promise.all([ctx.load('topology.json'), ctx.load('cmb.json')]);
  if (ctx.unmounted) return;
  if (!topo.ok) {
    const p = panel({ title: 'Manifold Selector', span: 12 });
    p.body.append(loadFailure(topo));
    el.append(grid(p.root));
    return;
  }
  const T = topo.data;

  mountSelector(el, T);
  mountPrinciples(el, T);
  mountVacuum(el, ctx, T);
  mountCurved(el, T);
  mountEvidence(el, T, cmb.ok ? cmb.data : null);
}

// --- principle filter and the 18-manifold table ----------------------------------------------

function mountSelector(el: HTMLElement, T: TopologyBundle): void {
  const ms = T.manifolds;
  const active = new Set<PrincipleKey>(PRINCIPLE_KEYS);
  const byKey = new Map(T.principles.map((p) => [p.key as PrincipleKey, p]));

  const sel = panel({
    title: 'Principle filter · [TOP] Main Theorem',
    subtitle: 'Each principle excludes manifolds by a proved theorem (P1 Theorem A, P2 Theorem C′, P3 Theorem B-loc, P4 Theorem D). Survivors are highlighted.',
    badges: [badge('proved', '[TOP] Secs. 4–5'), badge('live')],
    span: 12,
  });

  const chips = PRINCIPLE_KEYS.map((k) => {
    const p = byKey.get(k);
    const btn = h(
      'button',
      { type: 'button', class: 'm1-pchip', 'aria-pressed': 'true', title: p ? `${p.statement} (${p.theorem})` : k },
      h('span', { class: 'm1-pchip-key' }, k),
      h('span', { class: 'm1-pchip-name' }, p?.name ?? k),
      p ? statusTag(p.status) : null,
    );
    btn.addEventListener('click', () => {
      if (active.has(k)) active.delete(k);
      else active.add(k);
      update();
    });
    return { k, btn };
  });

  const presets: { label: string; set: PrincipleKey[] }[] = [
    { label: '{P2, P3}', set: ['P2', 'P3'] },
    { label: '{P2, P4}', set: ['P2', 'P4'] },
    { label: 'All four', set: ['P1', 'P2', 'P3', 'P4'] },
    { label: 'None', set: [] },
  ];
  const presetRow = h(
    'div',
    { class: 'm1-presets', role: 'group', 'aria-label': 'Principle presets' },
    h('span', { class: 'small muted' }, 'Presets'),
    ...presets.map((pr) => {
      const b = h('button', { type: 'button', class: 'btn' }, pr.label);
      b.addEventListener('click', () => {
        active.clear();
        pr.set.forEach((k) => active.add(k));
        update();
      });
      return b;
    }),
  );

  const resultTiles = h('div');
  const cols: RawColumn<ManifoldEntry>[] = [
    { label: 'Label', cell: (m) => h('strong', { class: 'mono' }, m.label) },
    { label: 'Name', cell: (m) => [m.name, m.note ? h('span', { class: 'm1-note-mark', title: m.note }, ' *') : null], className: 'wrap' },
    { label: 'Closed', cell: (m) => yesNo(m.closed) },
    { label: 'Orientable', cell: (m) => yesNo(m.orientable) },
    { label: 'Holonomy', cell: (m) => m.holonomy, className: 'wrap' },
    { label: 'Translations only', cell: (m) => yesNo(m.translations_only) },
    ...PRINCIPLE_KEYS.map(
      (k): RawColumn<ManifoldEntry> => ({
        label: P_HEAD[k],
        cell: (m) => principleCell(m.principles[k]),
        className: `m1-pcol m1-pcol-${k}`,
      }),
    ),
    {
      label: 'Verdict under the active principles',
      cell: (m) => {
        const ex = excludedBy(m, active);
        if (!active.size) return h('span', { class: 'muted small' }, 'no principle active');
        return ex.length
          ? h('span', { class: 'm1-verdict' }, statusPill('critical', 'excluded'), h('span', { class: 'small text-2' }, ` by ${ex.join(', ')}`))
          : statusPill('good', 'survives');
      },
      className: 'wrap',
    },
  ];
  const tbl = rawTable(cols, ms, (m) => (active.size && !excludedBy(m, active).length ? 'm1-row-survivor' : active.size ? 'm1-row-out' : null));
  const e16 = ms.find((m) => m.note);
  sel.body.append(
    h('div', { class: 'm1-pchips', role: 'group', 'aria-label': 'Principles' }, ...chips.map((c) => c.btn)),
    presetRow,
    resultTiles,
    h('div', { class: 'm1-manifold-table' }, tbl.root),
    ...(e16 ? [h('div', { class: 'chart-caption' }, `* ${e16.label}: ${e16.note}.`)] : []),
  );
  el.append(grid(sel.root));

  const minimal = minimalSelectingSets(ms);
  const single = failingExactly(ms, 1);

  function update(): void {
    chips.forEach((c) => c.btn.setAttribute('aria-pressed', String(active.has(c.k))));
    sel.root.dataset.active = [...active].join(' ');
    for (const k of PRINCIPLE_KEYS) sel.root.classList.toggle(`m1-off-${k}`, !active.has(k));
    const surv = survivors(ms, [...active].sort());
    tbl.render(ms);
    clear(resultTiles);
    resultTiles.append(
      tilesN(
        3,
        statTile({
          label: `Survivors of ${active.size ? fmtSet([...active].sort()) : 'no principle'}`,
          value: surv.length <= 4 ? surv.join(' · ') : `${surv.length} manifolds`,
          sub: [badge('live'), `${surv.length} of ${ms.length}${surv.length > 4 ? `: ${surv.join(', ')}` : ''}`],
        }),
        statTile({
          label: 'Inclusion-minimal sets selecting E1',
          value: minimal.map(fmtSet).join(' · ') || '—',
          sub: [badge('computed', 'm1-logic'), 'paper: “either pair {P2, P3} or {P2, P4} already forces it”'],
        }),
        statTile({
          label: 'Manifolds failing exactly one principle',
          value: single.join(' · ') || 'none',
          sub: [badge('computed', '[TOP] Sec. 4 table'), single.length ? `each fails only ${[...new Set(single.flatMap((l) => failedPrinciples(ms.find((m) => m.label === l)!)))].join(', ')}` : ''],
        }),
      ),
    );
  }
  update();
}

// --- the four principles, what is proved, independence, the false sentence ---------------------

function mountPrinciples(el: HTMLElement, T: TopologyBundle): void {
  const ms = T.manifolds;
  const pr = panel({
    title: 'The four principles · [TOP] Sec. 3',
    subtitle: 'Each tag gives the physical status the paper assigns; each principle enters through one proved theorem.',
    badges: [badge('quoted', '[TOP] Sec. 3')],
    span: 7,
    stretch: true,
  });
  pr.body.append(
    h(
      'div',
      { class: 'm1-principles' },
      ...T.principles.map((p) =>
        h(
          'article',
          { class: 'm1-principle' },
          h('div', { class: 'm1-principle-head' }, h('span', { class: 'm1-pchip-key' }, p.key), h('strong', null, p.name), statusTag(p.status)),
          h('p', null, p.statement),
          h('div', { class: 'row small' }, badge('proved', `[TOP] ${p.theorem}`), h('span', { class: 'muted' }, p.source ?? '[TOP] Sec. 3')),
        ),
      ),
    ),
  );

  const by = (l: string): ManifoldEntry | undefined => ms.find((m) => m.label === l);
  const sat = (m: ManifoldEntry): string => {
    const f = failedPrinciples(m);
    const ok = PRINCIPLE_KEYS.filter((k) => !f.includes(k));
    return `satisfies ${ok.length ? ok.join(', ') : 'none'}${f.length ? `, fails ${f.join(', ')}` : ''}`;
  };
  const single = failingExactly(ms, 1);
  const s2 = survivors(ms, ['P2']);
  const s3 = survivors(ms, ['P3']);

  const pv = panel({
    title: 'What the paper proves',
    subtitle: '[TOP] Lemma 1, Theorems A, B, B-loc, C′, D and the Main Theorem; status table of Sec. 13',
    span: 5,
    stretch: true,
  });
  pv.body.append(
    proofList(
      [
        { kind: 'proved', src: '[TOP] Lemma 1', text: 'All 18 flat 3-manifolds solve Einstein’s equations with the same scale factor: local physics cannot decide the topology, so the selection must be global.' },
        { kind: 'proved', src: '[TOP] Theorem A', text: 'Under P1 (the measured CP violation, J ≠ 0) space is orientable: excludes E7–E10, E13–E15, E17.' },
        { kind: 'proved', src: '[TOP] Theorem C′', text: 'A finite partition function at every temperature ⇔ finite volume ⇔ closed ⇔ no-boundary creation. P2 excludes exactly E11–E18.' },
        { kind: 'proved', src: '[TOP] Theorems B, B-loc', text: 'A homogeneous vacuum on a closed flat space exists only on the torus; locally homogeneous observables on any flat space ⇔ every deck transformation is a translation (E1, E11, E16 without screw, E18).' },
        { kind: 'proved', src: '[TOP] Theorem D', text: 'A parallel spinor exists ⇔ translations only and the periodic spin structure (P4).' },
        { kind: 'proved', src: '[TOP] Main Theorem', text: 'P2 with P3, or P2 with P4, forces the flat 3-torus T³ = ℝ³/Λ; P1 then holds automatically.' },
        { kind: 'assumed', src: '[TOP] Sec. 3', text: 'P2 and P3 are hypotheses, P4 is speculative and optional; only P1 is experimentally grounded. The selection is a theorem from these premises, not an observation.' },
        { kind: 'quoted', src: '[TOP] Sec. 9', text: 'Observation does not select a topology: Planck shows no evidence for a torus of any allowed size, nor for E2–E10 or the spherical candidates.' },
      ],
      true,
    ),
    sectionLabel('Independence, from the table (computed)'),
    h(
      'ul',
      { class: 'm1-list' },
      ...['E18', 'E2', 'E7'].map((l) => {
        const m = by(l);
        return m ? h('li', null, h('strong', { class: 'mono' }, l), ` ${m.name}: ${sat(m)}`) : null;
      }),
      h('li', null, `P2 alone leaves ${s2.length} (${s2[0]}–${s2[s2.length - 1]}); P3 alone leaves ${s3.length} (${s3.join(', ')}); together only ${survivors(ms, ['P2', 'P3']).join(', ')}.`),
    ),
    h(
      'div',
      { class: 'callout m1-callout' },
      statusPill('warning', 'Paper sentence contradicted by its own table'),
      h(
        'p',
        null,
        '[TOP] Secs. 4 and 5.6 state that “every manifold except E1 fails at least two principles”. Reading the Sec. 4 status table, ',
        h('strong', null, single.join(', ')),
        ' fail only P2 (they are open but orientable, translation-only and admit the periodic spin structure). The Main Theorem is unaffected: both selecting pairs contain P2.',
      ),
      h('div', { class: 'row small' }, badge('computed', 'm1-logic.failingExactly'), badge('quoted', '[TOP] Secs. 4, 5.6')),
    ),
  );
  el.append(grid(pr.root, pv.root));
}

// --- Theorem B made visible: E1 vs E2 vacuum variance --------------------------------------------

function mountVacuum(el: HTMLElement, ctx: PanelContext, T: TopologyBundle): void {
  const vac = T.vacuum;
  const p = panel({
    title: 'Theorem B made visible · E1 vs E2 vacuum variance',
    subtitle: 'Renormalized local variance ⟨φ²⟩ − ⟨φ²⟩_Minkowski of a free scalar in its ground state, minus its minimum over the slice z = 0 (perpendicular to the E2 screw axes)',
    badges: [badge('proved', '[TOP] Theorems B, B-loc'), badge('computed', 'qchas.topology.vacuum')],
    span: 12,
  });
  el.append(grid(p.root));
  if (!vac?.E1 || !vac.E2) {
    p.body.append(h('div', { class: 'empty-state' }, 'topology.json carries no vacuum slices (built before qchas.topology.vacuum existed): run qchas build.'));
    return;
  }
  const state = { smoothed: !!(vac.E1_smoothed && vac.E2_smoothed) };
  const ctl = segmented({
    label: 'Smoothing',
    options: [
      { value: 'raw', label: 'Image sum' },
      { value: 'smooth', label: 'Smoothed σ = 0.35 (as Fig. 1)' },
    ],
    value: state.smoothed ? 'smooth' : 'raw',
    onChange: (v) => {
      state.smoothed = v === 'smooth';
      render();
    },
  });
  if (!vac.E1_smoothed || !vac.E2_smoothed) ctl.root.querySelectorAll('button')[1]?.setAttribute('disabled', '');
  const tilesHost = h('div');
  const left = h('div', { class: 'm1-heat' });
  const right = h('div', { class: 'm1-heat' });
  p.body.append(
    h('div', { class: 'm1-toolbar' }, h('span', { class: 'control-label' }, 'Field'), ctl.root),
    tilesHost,
    h('div', { class: 'm1-two' }, left, right),
    h(
      'div',
      { class: 'chart-caption' },
      `Shared colour scale. Grid: one lattice cell, ${vac.E1.x.length} × ${vac.E1.y.length} points, image sum with cutoff ${fmtNum(vac.E1.cutoff_mL, 3)}/m. ${vac.note ?? ''}. `,
      'Rings mark where the half-turn screw axes cross the slice. [TOP] Fig. 1 shows the same comparison over 2 × 2 cells.',
    ),
  );
  let charts: { destroy(): void }[] = [];

  function slice(s: VacuumSlice, zmax: number, host: HTMLElement, label: string, axes: [number, number][]): void {
    const ink = cssVar('--text-1');
    const r = 0.025;
    charts.push(
      keep(
        ctx,
        heatmapChart(
          host,
          { z: s.excess_minus_min, x: s.x, y: s.y, scale: 'sequential', zmin: 0, zmax, colorbarTitle: 'excess − min', zLabel: 'excess − min' },
          {
            xTitle: 'x (lattice units)',
            yTitle: 'y (lattice units)',
            height: 360,
            ariaLabel: `${label}: vacuum variance excess over the slice`,
            format: (v) => fmtNum(v, 4),
            layout: {
              title: { text: label, font: { size: 12.5 }, x: 0.02, xanchor: 'left' } as never,
              margin: { t: 34 } as never,
              yaxis: { scaleanchor: 'x', scaleratio: 1 } as never,
              shapes: axes.map(([x, y]) => ({ type: 'circle', xref: 'x', yref: 'y', x0: x - r, x1: x + r, y0: y - r, y1: y + r, line: { color: ink, width: 1.5 } })) as never,
            },
          },
        ),
      ),
    );
  }

  function render(): void {
    charts.forEach((c) => c.destroy());
    charts = [];
    clear(left);
    clear(right);
    const e1 = (state.smoothed ? vac!.E1_smoothed : vac!.E1)!;
    const e2 = (state.smoothed ? vac!.E2_smoothed : vac!.E2)!;
    const zmax = Math.max(e1.max - e1.min, e2.max - e2.min) || 1;
    const [x0, x1] = [Math.min(...e2.x), Math.max(...e2.x)];
    const [y0, y1] = [Math.min(...e2.y), Math.max(...e2.y)];
    const axes = (vac!.E2_axes ?? []).filter(([x, y]) => x >= x0 && x <= x1 && y >= y0 && y <= y1);
    slice(e1, zmax, left, 'E1 3-torus: constant', []);
    slice(e2, zmax, right, 'E2 half-turn space: louder on the screw axes', axes);
    clear(tilesHost);
    tilesHost.append(
      tilesN(
        4,
        statTile({ label: 'E1: spread over the slice (max − min)', value: fmtNum(e1.max - e1.min, 3), sub: [badge('computed'), 'Theorem B-loc: exactly constant'] }),
        statTile({ label: 'E2: spread over the slice (max − min)', value: fmtNum(e2.max - e2.min, 3), sub: [badge('computed'), `range ${fmtNum(e2.min, 4)} – ${fmtNum(e2.max, 4)}`] }),
        statTile({ label: 'E1: excess everywhere', value: fmtNum(e1.min, 5), sub: [badge('computed'), 'units of m², image sum'] }),
        statTile({ label: 'Mass m·L', value: fmtNum(e1.m_L, 3), sub: [badge('assumed', 'not stated in [TOP]'), state.smoothed ? `smoothing σ = ${fmtNum(e1.smoothing_sigma, 2)} (Fig. 1)` : 'unsmoothed field'] }),
      ),
    );
  }
  render();
  ctx.onThemeChange(() => render());
}

// --- curved candidates and the Weeks manifold ------------------------------------------------------

function mountCurved(el: HTMLElement, T: TopologyBundle): void {
  const sp = panel({
    title: 'Curved candidates · spherical spaces',
    subtitle: 'Homogeneous quotients S³/G tested against Planck at Ω_k on a pre-registered grid; ln B against simply connected S³ at the same curvature',
    badges: [badge('quoted', '[TOP] Sec. 9.3')],
    span: 6,
    stretch: true,
  });
  sp.body.append(
    rawTable(
      [
        { label: 'Space', cell: (r) => r.name, className: 'wrap' },
        { label: '|G|', num: true, cell: (r) => String(r.order) },
        { label: 'Shortest loop / R_c', cell: (r) => r.shortest_loop },
        { label: 'Allowed Ω_k', num: true, cell: (r) => `[${r.omega_k_range.map((v) => v.toFixed(4)).join(', ')}]` },
        { label: 'ln B vs S³', num: true, cell: (r) => (r.lnB_vs_S3 > 0 ? '+' : '') + r.lnB_vs_S3.toFixed(2) },
        { label: 'Notes', cell: (r) => r.note, className: 'wrap' },
      ],
      T.curved,
    ).root,
    h(
      'p',
      { class: 'm1-takeaway' },
      'No spherical space shows evidence. Near the circle limit the data disfavour them: the dodecahedral configuration once proposed to explain the low quadrupole is rejected with ln B ≈ −4. ',
      badge('quoted', '[TOP] Sec. 9.3'),
    ),
  );

  const wk = panel({
    title: 'Weeks manifold · Prop. 9.4',
    subtitle: 'Smallest closed hyperbolic space: matched circles need 2χ*/R_c above the systole 0.58460',
    badges: [badge('proved', '[TOP] Prop. 9.4'), badge('quoted', '[TOP] Sec. 9.4')],
    span: 6,
    stretch: true,
  });
  wk.body.append(
    rawTable(
      [
        { label: 'Ω_k', num: true, cell: (r) => r.omega_k.toFixed(4) },
        { label: 'R_c (Gpc)', num: true, cell: (r) => fmtFixed(r.Rc_Gpc, 1) },
        { label: 'χ*/R_c', num: true, cell: (r) => fmtFixed(r.chi_over_Rc, 3) },
        { label: 'Shortest loop / D', num: true, cell: (r) => fmtFixed(r.shortest_loop_over_D, 2) },
        { label: 'V_Weeks / V_LSS', num: true, cell: (r) => fmtFixed(r.volume_ratio, 1) },
        { label: 'Circles', cell: (r) => (r.circles.startsWith('no') ? h('span', { class: 'muted' }, r.circles) : h('span', null, r.circles)), className: 'wrap' },
      ],
      T.weeks,
      (r) => (r.circles.startsWith('no') ? null : 'm1-row-flag'),
    ).root,
    h(
      'p',
      { class: 'm1-takeaway' },
      'Circles need Ω_k ≥ 0.0080, 3.8σ above the Planck + BAO value Ω_k = 0.0007 ± 0.0019. Below that a Weeks universe could show only statistically, which the paper estimates (not proves) to be undetectable. ',
      badge('quoted', '[TOP] Prop. 9.4'),
    ),
  );
  el.append(grid(sp.root, wk.root));
}

// --- Planck evidence tables ------------------------------------------------------------------------

function mountEvidence(el: HTMLElement, T: TopologyBundle, cmb: CmbBundle | null): void {
  const state = { shape: 'all' };
  const filter = segmented({
    label: 'Lattice shape',
    options: [
      { value: 'all', label: 'All shapes' },
      { value: 'FCC', label: 'FCC' },
      { value: 'cubic', label: 'Cubic' },
      { value: 'BCC', label: 'BCC' },
    ],
    value: 'all',
    onChange: (v) => {
      state.shape = v;
      render();
    },
  });
  const head = panel({
    title: 'Confronting the Planck 2018 sky · [TOP] Sec. 9',
    subtitle: 'Quoted from the paper. Planck shows no evidence for a torus of any allowed size; small tori are excluded at 95 %.',
    badges: [badge('quoted', '[TOP] Sec. 9'), badge('measured', 'Planck 2018 Commander')],
    span: 12,
  });
  const lik = cmb?.likelihood ?? null;
  head.body.append(
    tilesN(
      4,
      statTile({ label: 'Pre-registered FCC family', value: 'ln B = −0.23', sub: [badge('quoted', '[TOP] Sec. 9'), 'p = 0.25: no evidence'] }),
      statTile({ label: 'Averaged over all 27 torus models', value: 'ln B = −0.06', sub: [badge('quoted', '[TOP] Sec. 9'), 'p = 0.15 (look-elsewhere built in)'] }),
      statTile({ label: '95 % lower bounds on the shortest loop', value: '1.04 · 1.00 · 0.985', unit: 'D', sub: [badge('quoted', '[TOP] Table 9b'), 'FCC · BCC · cubic; D = 27.74 Gpc'] }),
      statTile({ label: 'E2–E10, prior-averaged', value: '−0.73 … +0.23', unit: 'ln B', sub: [badge('quoted', '[TOP] Sec. 9.2'), 'none shows evidence'] }),
    ),
    h(
      'div',
      { class: 'm1-toolbar' },
      h('span', { class: 'control-label' }, 'Shape (scopes Tables 9a and 9b)'),
      filter.root,
      h(
        'span',
        { class: 'small muted m1-grow' },
        lik
          ? 'Our orientation-scanned likelihood on the build map is shown next to Table 9a: an independent re-analysis of a different map, not a reproduction of the Planck numbers.'
          : 'Q-CHAS’s own orientation-scanned likelihood is not in this build (cmb.json likelihood = null, see Data & Provenance), so only the paper’s Planck values are shown.',
      ),
    ),
  );
  el.append(grid(head.root));

  const pa = panel({
    title: 'Table 9a · evidence for a torus',
    subtitle: 'Optimal orientation-averaged statistic, calibrated on 141 simulated skies with no topology',
    badges: [badge('quoted', '[TOP] Table 9a')],
    span: 6,
    stretch: true,
  });
  const pb = panel({
    title: 'Table 9b · exclusion of tori',
    subtitle: 'Calibrated on simulated skies in which that torus is true (orientation uniformly random)',
    badges: [badge('quoted', '[TOP] Table 9b')],
    span: 6,
    stretch: true,
  });
  const hostA = h('div');
  const hostB = h('div');
  pa.body.append(hostA);
  pb.body.append(hostB);
  el.append(grid(pa.root, pb.root));

  const src = cmb?.source ?? 'map';
  const ours = (r: EvidenceTorusEntry) =>
    lik?.find((x) => x.shape.toLowerCase() === r.shape.toLowerCase() && Math.abs(x.d_over_D - r.d_over_D) < 1e-6) ?? null;
  const sign = (v: number, d = 2): string => (v > 0 ? '+' : v < 0 ? '−' : '') + Math.abs(v).toFixed(d);

  function render(): void {
    const keepRow = (s: string): boolean => state.shape === 'all' || s.toLowerCase() === state.shape.toLowerCase();
    const colsA: RawColumn<EvidenceTorusEntry>[] = [
      { label: 'Shape', cell: (r) => r.shape },
      { label: 'd / D', num: true, cell: (r) => r.d_over_D.toFixed(2) },
      { label: 'max Δln L', num: true, cell: (r) => sign(r.max_dlnL) },
      { label: 'ln B (optimal)', num: true, cell: (r) => sign(r.lnB) },
      { label: '95 % of ln B, no topology', num: true, cell: (r) => sign(r.null95) },
      { label: 'p (no topology)', num: true, cell: (r) => r.p.toFixed(2) },
    ];
    if (lik)
      colsA.push(
        { label: `ln B · ${src} (ours)`, num: true, cell: (r) => (ours(r) ? sign(ours(r)!.lnB) : '—') },
        { label: `p · ${src} (ours)`, num: true, cell: (r) => (ours(r) ? ours(r)!.p.toFixed(2) : '—') },
      );
    clear(hostA);
    hostA.append(
      rawTable(colsA, T.evidence_torus.filter((r) => keepRow(r.shape))).root,
      h('p', { class: 'm1-takeaway' }, 'Largest single-model value ln B = +1.12 (BCC, 1.05 D) has look-elsewhere-corrected p = 0.11. None reaches “positive” on the Kass–Raftery scale. ', badge('quoted', '[TOP] Sec. 9')),
    );
    clear(hostB);
    hostB.append(
      rawTable<ExclusionTorusEntry>(
        [
          { label: 'Shape', cell: (r) => r.shape },
          { label: 'd / D', num: true, cell: (r) => r.d_over_D.toFixed(3).replace(/0$/, '') },
          { label: 'ln B, Planck', num: true, cell: (r) => sign(r.lnB_planck) },
          { label: 'median if true', num: true, cell: (r) => sign(r.median_if_true) },
          { label: '5 % if true', num: true, cell: (r) => sign(r.q05_if_true) },
          { label: 'P(ln B ≤ Planck | true)', num: true, cell: (r) => r.p_le_planck.replace('<=', '≤ ') },
          { label: 'Sims', num: true, cell: (r) => String(r.sims) },
          { label: 'Verdict', cell: (r) => (r.verdict.startsWith('excluded') ? statusPill('critical', r.verdict) : statusPill('good', r.verdict)) },
        ],
        T.exclusion_torus.filter((r) => keepRow(r.shape)),
        (r) => (r.verdict.startsWith('excluded') ? 'm1-row-excl' : null),
      ).root,
    );
  }
  render();

  const o = panel({
    title: 'All other compact flat topologies · E2–E10',
    subtitle: 'Pre-registered priors: size uniform in [1, 2] D, observer uniform in the cell, Haar orientation; 2 ≤ ℓ ≤ 20',
    badges: [badge('quoted', '[TOP] Sec. 9.2')],
    span: 12,
  });
  o.body.append(
    rawTable(
      [
        { label: 'Space', cell: (r) => r.space },
        { label: 'Orientable', cell: (r) => yesNo(r.orientable) },
        { label: 'Prior samples', num: true, cell: (r) => String(r.samples) },
        { label: 'ln B (prior average)', num: true, cell: (r) => r.lnB.replace('-', '−') },
        { label: 'Largest single configuration', num: true, cell: (r) => sign(r.largest) },
      ],
      T.other_flat,
    ).root,
    h(
      'p',
      { class: 'm1-takeaway' },
      'None of the nine spaces shows evidence. E5 is almost completely excluded before any likelihood (1.9 × 10⁻⁴ of its prior volume survives the circle-search limit). The non-orientable E7–E10, forbidden by Theorem A, are not favoured either. Every observational test is consistent with the Main Theorem, but none selects a topology. ',
      badge('quoted', '[TOP] Sec. 9.2'),
    ),
    note(
      [
        'The likelihood is exact for Gaussian multipoles: Δln L(R) = −½ aᵀ(C(R)⁻¹ − C_ΛCDM⁻¹)a − ½ ln det C₀ / det C_ΛCDM, maximized over all of SO(3). ([TOP] Sec. 9)',
        'At d ≥ 1.05 D even a real torus would typically give ln B ≲ 1: the unknown orientation consumes almost all of the available information, and no rotation-invariant law can supply it (Proposition G).',
        'The paper’s numbers are for the Planck 2018 Commander map. Analyses run in this dashboard use the build map named in the status bar (WMAP 9-year ILC when Planck is unreachable): a different map, so they are an independent re-analysis, not a reproduction.',
      ],
      'How the evidence is computed',
    ),
  );
  el.append(grid(o.root));
}
