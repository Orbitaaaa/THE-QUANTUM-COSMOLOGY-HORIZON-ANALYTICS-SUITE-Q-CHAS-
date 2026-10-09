/**
 * Mission overview: pipeline status, the two analysis modules with key numbers from the bundles,
 * a sky monitor, the onboard physics self-test (live TypeScript physics against printed paper
 * values) and the three source papers.
 */

import { type BundleName, BUNDLE_NAMES, type LoadResult, loadSkyMap } from '../lib/data';
import { lineChart } from '../lib/charts';
import type { GlobeView, SkyView } from '../lib/sky';
import type { PanelContext, RouteInfo } from '../lib/panel';
import { FCC_PREREGISTERED_AXES_GAL, FIG4_PAPER, MIN_AREA_LP2, PAGE_MASS_FRACTION_PAPER, PAGE_TIME_PAPER, Z2_CERTIFIED } from '../lib/physics/constants';
import { circlePairs } from '../lib/physics/circles';
import { minimalArea } from '../lib/physics/capacity';
import { cubicBasis, epsteinZ2, rhombohedralCurve } from '../lib/physics/lattice';
import { curves, pageMassFraction, pageTimeFraction } from '../lib/physics/page';
import { PAPER_TABLE, solveQes } from '../lib/physics/qes';
import {
  type ProvenanceKind,
  badge,
  fmtBytes,
  fmtNum,
  fmtSci,
  fmtUtc,
  grid,
  h,
  missingBundle,
  note,
  panel,
  segmented,
  statTile,
  statusPill,
  table,
  tiles,
} from '../lib/ui';
import { ROUTES } from '../routes';

export const title = 'Mission Overview';
export const subtitle =
  'Q-CHAS confronts three new results in mathematical physics with real public survey data: WMAP and Planck CMB maps, SPHEREx spectral images and LIGO-Virgo strain.';

const BUILT = BUNDLE_NAMES.filter((n) => n !== 'manifest.json');

export async function mount(el: HTMLElement, ctx: PanelContext): Promise<void> {
  const m = ctx.manifest;

  // --- mission status ------------------------------------------------------------------------
  const status = panel({
    title: 'Pipeline status',
    subtitle: m ? `Bundles written by qchas ${m.qchas_version}` : 'No pipeline output found under ./data/',
    badges: [badge('measured', m?.cmb_map_source ?? undefined)].filter(() => !!m),
  });
  const built = new Set(m?.bundles ?? []);
  const totalBytes = m ? Object.values(m.datasets ?? {}).reduce((s, d) => s + (d.bytes || 0), 0) : null;
  const nWarn = m?.warnings?.length ?? 0;
  const nErr = m?.errors?.length ?? 0;
  status.body.append(
    tiles(
      statTile({ label: 'CMB map in use', value: m?.cmb_map_source ?? null, mono: true, sub: 'Planck PR3 preferred, WMAP 9-yr ILC fallback' }),
      statTile({ label: 'Pipeline build', value: m ? fmtUtc(m.generated_utc) : null, mono: true, sub: m ? `qchas ${m.qchas_version}` : 'not run' }),
      statTile({ label: 'Datasets on record', value: m ? Object.keys(m.datasets ?? {}).length : null, sub: m ? `${fmtBytes(totalBytes)} with SHA-256` : 'see Data & Provenance' }),
      statTile({ label: 'Bundles built', value: m ? `${BUILT.filter((b) => built.has(b)).length} / ${BUILT.length}` : null, sub: 'web/public/data/' }),
      statTile({
        label: 'Build warnings',
        value: m ? nWarn : null,
        sub: m
          ? nErr
            ? statusPill('serious', `${nErr} builder${nErr === 1 ? '' : 's'} failed`)
            : statusPill(nWarn ? 'warning' : 'good', nWarn ? 'review in Data & Provenance' : 'none')
          : 'no build yet',
      }),
    ),
  );
  if (!m)
    status.body.append(
      h('div', { style: { marginTop: '12px' } }, missingBundle('manifest.json', 'Every analysis view stays empty until the Python pipeline has fetched the real datasets and exported its bundles. The live physics below runs in your browser regardless.')),
    );
  status.body.append(
    h(
      'div',
      { style: { marginTop: '12px' } },
      note(
        [
          'The Planck 2018 sky shows no evidence for a torus of any allowed size and excludes small tori at 95 %: it requires d > 1.04 D for FCC, d > 1.00 D for BCC and d > 0.985 D for the cubic torus. [TOP] Sec. 9',
          'A pre-registered search finds no topological resonance (p = 0.78). Resonance amplifies the topological signal up to 64-fold only at known orientation. [TOP] Sec. 11',
          'A torus as large as 3 D could be decided only by a multi-field primordial line mapped in three dimensions, not by the CMB alone. [TOP] Sec. 11',
          "Field theory on a fixed background forces Hawking's curve; finite capacity A/4G, derived from the crossed product, forces the island and the Page transfer law. [BH] Theorems B, D, E, G",
          'Q-CHAS reproduces these statements and the numbers behind them. It claims no detection.',
        ],
        'What the papers prove, and what the data say',
      ),
    ),
  );
  el.append(grid(status.root));

  // --- modules ------------------------------------------------------------------------------------
  const names: BundleName[] = ['lattice.json', 'resonance.json', 'cmb.json', 'horizon.json', 'bounce.json', 'gw.json'];
  const loaded = Object.fromEntries(await Promise.all(names.map(async (n) => [n, await ctx.load(n)] as const))) as {
    [K in BundleName]?: LoadResult<unknown>;
  };
  if (ctx.unmounted) return;
  const get = <T>(n: BundleName): T | null => {
    const r = loaded[n];
    return r && r.ok ? (r.data as T) : null;
  };
  type Num = number | null | undefined;
  const lattice = get<{
    critical?: { FCC?: { Z: number } };
    flows?: { n?: number; counts?: Record<string, number>; paper_counts?: Record<string, number>; source?: string };
  }>('lattice.json');
  const resonance = get<{ search?: { p_max?: Num } | null }>('resonance.json');
  const cmb = get<{ stats?: { S12_obs?: Num }; source?: string }>('cmb.json');
  const horizon = get<{ page?: { page_time?: Num }; capacity?: { de_sitter?: { log10_entropy?: Num } } }>('horizon.json');
  const bounce = get<{ runs?: { a_max?: number; a_min?: Num; paper?: { a_min?: number | string | null } }[] }>('bounce.json');
  // compare with a run the paper tabulates (a_max = 2 in [QG] Theorem Lambda), else the first run
  const run0 = bounce?.runs?.find((r) => r.paper?.a_min !== undefined && r.paper.a_min !== null) ?? bounce?.runs?.[0];
  const gw = get<{ events?: unknown[] }>('gw.json');

  interface KeyNumber {
    label: string;
    value: string | null;
    bundle: BundleName;
    kind: ProvenanceKind;
    paper?: string;
  }
  const counts = lattice?.flows?.counts;
  const paperCounts = lattice?.flows?.paper_counts;
  // runs of different size are compared by their FCC share, not by raw counts
  const share = (c?: Record<string, number>): string => {
    const fcc = c?.FCC ?? 0;
    const total = Object.values(c ?? {}).reduce((a, b) => a + b, 0);
    return total ? `${fcc} : ${c?.BCC ?? 0} (FCC ${Math.round((100 * fcc) / total)} % of ${total})` : '—';
  };
  const m1: KeyNumber[] = [
    { label: 'Casimir Z(2) at the FCC minimum', value: lattice?.critical?.FCC ? fmtNum(lattice.critical.FCC.Z, 8) : null, bundle: 'lattice.json', kind: 'computed', paper: `${Z2_CERTIFIED.value.FCC.toFixed(8)} · ${Z2_CERTIFIED.source}` },
    {
      label: 'Shape-flow basins FCC : BCC',
      value: counts ? share(counts) : null,
      bundle: 'lattice.json',
      kind: 'computed',
      paper: paperCounts ? `${share(paperCounts)} · [TOP] Remark after Prop. E.2` : undefined,
    },
    { label: 'Resonance search p (max statistic)', value: resonance?.search ? fmtNum(resonance.search.p_max, 2) : null, bundle: 'resonance.json', kind: 'measured', paper: '0.78 · [TOP] Sec. 11.6' },
    { label: 'S1/2 of the map in use', value: cmb?.stats ? fmtNum(cmb.stats.S12_obs, 4) : null, bundle: 'cmb.json', kind: 'measured' },
  ];
  const m2: KeyNumber[] = [
    { label: 'Page time t_P / t_evap', value: horizon?.page ? fmtNum(horizon.page.page_time, 3) : null, bundle: 'horizon.json', kind: 'computed', paper: `${PAGE_TIME_PAPER.value} · ${PAGE_TIME_PAPER.source}` },
    { label: 'De Sitter entropy log₁₀ S', value: horizon?.capacity?.de_sitter ? fmtNum(horizon.capacity.de_sitter.log10_entropy, 5) : null, bundle: 'horizon.json', kind: 'computed' },
    {
      label: `Bounce ⟨a⟩_min (ℓ units${run0?.a_max ? `, a_max = ${run0.a_max}` : ''})`,
      value: run0 ? fmtNum(run0.a_min, 4) : null,
      bundle: 'bounce.json',
      kind: 'computed',
      paper: run0?.paper?.a_min !== undefined && run0.paper.a_min !== null ? `${run0.paper.a_min} · [QG] Theorem Q` : undefined,
    },
    { label: 'LIGO events analysed', value: gw?.events ? String(gw.events.length) : null, bundle: 'gw.json', kind: 'measured' },
  ];

  const moduleCard = (code: string, name: string, text: string, group: RouteInfo['group'], keys: KeyNumber[], span: 6) => {
    const p = panel({ title: `Module ${code}`, subtitle: name, span, stretch: true });
    const routes = ROUTES.filter((r) => r.group === group);
    p.body.append(
      h('p', { class: 'text-2', style: { margin: '0' } }, text),
      h(
        'ul',
        { class: 'tool-list' },
        ...routes.map((r) =>
          h(
            'li',
            null,
            h(
              'a',
              { href: `#/${r.id}` },
              h('span', { class: 'nav-code' }, r.code),
              h('span', null, r.label),
              h('span', { class: 'tool-desc' }, r.live ? badge('live', undefined, 'Computes closed-form physics live in the browser') : null),
            ),
          ),
        ),
      ),
      h('div', { class: 'muted small', style: { margin: '14px 0 6px' } }, 'Key numbers'),
      table(
        [
          { key: 'label', label: 'Quantity', className: 'wrap' },
          {
            key: 'value',
            label: 'Q-CHAS',
            align: 'num',
            format: (r) => (r.value ?? h('span', { class: 'muted', title: `${r.bundle} not built yet` }, 'not built')),
          },
          { key: 'paper', label: 'Paper', className: 'wrap small muted', format: (r) => r.paper ?? '—' },
          { key: 'kind', label: '', format: (r) => badge(r.kind) },
        ],
        keys,
      ),
    );
    return p.root;
  };
  el.append(
    grid(
      moduleCard(
        '1',
        'Global Topology & All-Sky Survey Analyst',
        'Selects the shape of flat space from four principles, relaxes the torus under Standard-Model Casimir stress, draws circles-in-the-sky and ghost images, bounds what any survey can learn, and runs the topological-resonance test on the real CMB map.',
        'm1',
        m1,
        6,
      ),
      moduleCard(
        '2',
        'Dressed Horizon & Evaporation Simulator',
        "Turns the type III₁ → II algebra ladder into a horizon capacity A/4G, plots Hawking's curve against the Page transfer law, locates the quantum extremal surface, simulates exactly unitary evaporation and Hayden-Preskill recovery, and reads real LIGO ringdowns.",
        'm2',
        m2,
        6,
      ),
    ),
  );

  // --- sky monitor + self-test ----------------------------------------------------------------------------
  const skyPanel = panel({
    title: 'Sky monitor',
    subtitle: 'Galactic frame · RING HEALPix · drag to rotate, click to read a pixel',
    badges: [badge('quoted', '[TOP] Sec. 9')],
    span: 6,
    stretch: true,
    flush: true,
  });
  const skyHost = h('div', { style: { padding: '0' } });
  const skyNote = h('div', { class: 'panel-foot' }, 'Points A1-A6: best-fit FCC lattice axes at d = 1.05 D, not significant ([TOP] Sec. 9).');
  skyPanel.body.append(skyHost);
  skyPanel.root.append(skyNote);

  const selfTest = panel({
    title: 'Onboard physics self-test',
    subtitle: 'Live TypeScript physics against the numbers printed in the papers',
    badges: [badge('live')],
    span: 6,
    stretch: true,
    flush: true,
  });
  el.append(grid(skyPanel.root, selfTest.root));
  selfTest.body.append(selfTestTable());

  void mountSky(skyHost, skyPanel.head, skyNote, ctx);

  // --- live previews of two paper figures ---------------------------------------------------------
  const pagePrev = panel({
    title: 'Page transfer law, live',
    subtitle: "Hawking's curve, the capacity A/4G and S(R) ≤ min{S_H, A/4G}, in units of S_BH(M₀)",
    badges: [badge('computed', '[BH] Fig. 1')],
    actions: [h('a', { class: 'btn', href: '#/page-curve' }, 'Open 2.2')],
    span: 6,
    stretch: true,
  });
  const zPrev = panel({
    title: 'Casimir energy of the torus shape, live',
    subtitle: 'Z(2) at unit covolume along the rhombohedral family (Ewald sum in this browser)',
    badges: [badge('computed', '[TOP] Fig. 2')],
    actions: [h('a', { class: 'btn', href: '#/lattice' }, 'Open 1.2')],
    span: 6,
    stretch: true,
  });
  el.append(grid(pagePrev.root, zPrev.root));
  const pc = curves(201);
  const h1 = lineChart(
    pagePrev.body,
    [
      { name: 'Page curve S(R)', x: pc.t, y: pc.page, color: 1 },
      { name: "Hawking's curve S_H", x: pc.t, y: pc.hawking, color: 2 },
      { name: 'Capacity A/4G', x: pc.t, y: pc.bh_entropy, color: 'context', dash: 'dash' },
    ],
    {
      xTitle: 'time / lifetime',
      yTitle: 'entropy / S_BH(M₀)',
      height: 280,
      refLines: [{ axis: 'x', value: pc.page_time, label: `Page time ${pc.page_time.toFixed(3)}` }],
      ariaLabel: 'Page curve, Hawking curve and black-hole entropy against time',
    },
  );
  const alpha = Array.from({ length: 141 }, (_, i) => 45 + 0.5 * i);
  const rc = rhombohedralCurve(alpha);
  const h2 = lineChart(
    zPrev.body,
    [{ name: 'Z(2), rhombohedral α', x: rc.alpha_deg, y: rc.Z, color: 1, markers: false }],
    {
      xTitle: 'rhombohedral angle α (deg)',
      yTitle: 'Z(2) at unit covolume',
      height: 280,
      refLines: [
        { axis: 'x', value: 60, label: 'FCC' },
        { axis: 'x', value: 90, label: 'SC' },
        { axis: 'x', value: (Math.acos(-1 / 3) * 180) / Math.PI, label: 'BCC' },
      ],
      ariaLabel: 'Casimir functional Z(2) along the rhombohedral family',
    },
  );
  ctx.onCleanup(() => {
    h1.destroy();
    h2.destroy();
  });

  // --- papers --------------------------------------------------------------------------------------------
  const papers = panel({ title: 'Source papers', subtitle: 'Keys used in every source badge and in the code' });
  papers.body.append(
    h(
      'div',
      { class: 'card-list' },
      paperCard('[TOP]', 'Cosmic topology · research note · 27 Sep 2026', 'Selecting the shape of flat space', 'Four principles about quantum fields at the Big Bang single out the 3-torus among the 18 flat 3-manifolds; the note derives its exact CMB signatures and a causal law limiting what any observation can learn about topology.', 'Status table: [TOP] Sec. 13'),
      paperCard('[QG]', 'Quantum gravity · research note v5 · 29 Sep 2026', 'Quantum gravity from the crossed product', 'Imposing a gravitational constraint on a type III₁ field algebra is a crossed product; time, the state of the universe, a unitary quantum Big Bang and black-hole information follow from that one operation.', 'Status table: [QG] Sec. 20'),
      paperCard('[BH]', 'Black-hole information · research note · 28 Sep 2026', 'What field theory forbids, and what the horizon must do', "Field theory on a fixed background forces Hawking's curve; the finite capacity A/4G of the dressed horizon forces the island and the Page transfer law S(R) ≤ min{S_H, A/4G}.", 'Status table: [BH] Sec. 16'),
    ),
  );
  el.append(grid(papers.root));
}

function paperCard(key: string, eyebrow: string, titleText: string, summary: string, status: string): HTMLElement {
  return h(
    'article',
    { class: 'card' },
    h('div', { class: 'row', style: { justifyContent: 'space-between' } }, h('span', { class: 'card-eyebrow' }, eyebrow), h('code', null, key)),
    h('h3', { class: 'card-title' }, titleText),
    h('p', { class: 'card-text' }, summary),
    h('div', { class: 'muted small mono' }, status),
  );
}

interface Check {
  name: string;
  ours: string;
  paper: string;
  ok: boolean;
  source: string;
}

/** Half a unit in the last printed digit of a decimal or scientific literal. */
export function halfUlp(printed: string): number {
  const [mant, exp] = printed.toLowerCase().split('e');
  const decimals = (mant!.split('.')[1] ?? '').length;
  return 0.5 * 10 ** (-decimals + Number(exp ?? 0));
}

/** A numeric check that passes when ours agrees with the printed value to its last digit. */
function printedCheck(name: string, ours: number, printed: string, source: string, digits = 4, tol = halfUlp(printed)): Check {
  const sci = /e/i.test(printed);
  const show = (x: number, d: number): string => (sci || Math.abs(x) < 1e-3 ? fmtSci(x, d) : String(Number(x.toPrecision(d))));
  const paper = sci ? fmtSci(Number(printed), printed.split(/e/i)[0]!.replace(/[-.]/g, '').length) : printed;
  return { name, ours: show(ours, digits), paper, ok: Math.abs(ours - Number(printed)) <= tol * (1 + 1e-12), source };
}

function selfTestTable(): HTMLElement {
  const t0 = performance.now();
  const qes = solveQes(10, 1, 20);
  const qesBig = solveQes(100, 10, 60);
  const pairs = (s: string) => circlePairs(cubicBasis(s), FIG4_PAPER.value.d_over_D);
  const fcc = pairs('FCC');
  const bcc = pairs('BCC');
  const checks: Check[] = [
    ...(['SC', 'BCC', 'FCC'] as const).map((s) =>
      printedCheck(`Casimir Z(2), ${s} torus`, epsteinZ2(cubicBasis(s)), String(Z2_CERTIFIED.value[s]), Z2_CERTIFIED.source, 10, 1e-9),
    ),
    printedCheck('Page time t_P / t_evap', pageTimeFraction(), String(PAGE_TIME_PAPER.value), PAGE_TIME_PAPER.source),
    printedCheck('Mass fraction at the Page time', pageMassFraction(), String(PAGE_MASS_FRACTION_PAPER.value), PAGE_MASS_FRACTION_PAPER.source),
    printedCheck('QES depth, M = 10, c = 1, r_b = 20 r_s', qes.depth, '3.11e-5', PAPER_TABLE.source),
    printedCheck('QES lag, M = 100, c = 10, r_b = 60 r_s', qesBig.lag, '4636.30', PAPER_TABLE.source, 7),
    printedCheck('Circle radius at d = 0.8 D (deg)', fcc[0]?.radius_deg ?? NaN, String(FIG4_PAPER.value.radius_deg), FIG4_PAPER.source),
    {
      name: 'Circle pairs at d = 0.8 D, FCC / BCC',
      ours: `${fcc.length} / ${bcc.length}`,
      paper: `${FIG4_PAPER.value.pairs.FCC} / ${FIG4_PAPER.value.pairs.BCC}`,
      ok: fcc.length === FIG4_PAPER.value.pairs.FCC && bcc.length === FIG4_PAPER.value.pairs.BCC,
      source: FIG4_PAPER.source,
    },
    printedCheck('Area gap a₀ (l_P²)', minimalArea(), String(MIN_AREA_LP2.value), MIN_AREA_LP2.source),
  ];
  const ms = performance.now() - t0;
  const passed = checks.filter((c) => c.ok).length;
  return h(
    'div',
    null,
    table(
      [
        {
          key: 'name',
          label: 'Check',
          className: 'wrap',
          format: (c) => h('div', null, c.name, h('div', { class: 'small muted mono' }, c.source)),
        },
        { key: 'ours', label: 'Browser', align: 'num' },
        { key: 'paper', label: 'Paper (printed)', align: 'num' },
        { key: 'ok', label: 'Result', format: (c) => (c.ok ? statusPill('good', 'match') : statusPill('warning', 'differs')) },
      ],
      checks,
    ),
    h(
      'div',
      { class: 'panel-foot' },
      `${passed} / ${checks.length} agree to the last printed digit (Z(2) to 1e-9) · computed in ${ms.toFixed(0)} ms in this browser`,
    ),
  );
}

async function mountSky(host: HTMLElement, head: HTMLElement, foot: HTMLElement, ctx: PanelContext): Promise<void> {
  const sky = await import('../lib/sky');
  if (ctx.unmounted) return;
  const mapRes = await loadSkyMap();
  if (ctx.unmounted) return;
  let view: SkyView | null = null;
  let builds = 0;
  let camera: 'outer' | 'inner' = 'outer';
  ctx.onCleanup(() => view?.dispose());
  const mapInput = mapRes.ok
    ? {
        map: mapRes.data.map,
        mask: mapRes.data.mask,
        // the exported 99.5 % display range; min/max are Galactic-plane extremes that wash out the CMB
        limit: mapRes.data.entry.range ? Math.abs(mapRes.data.entry.range[1]) || undefined : undefined,
        unit: mapRes.data.meta.unit === 'uK' ? 'μK' : mapRes.data.meta.unit,
      }
    : null;
  const cameraSwitch = segmented({
    label: 'Camera',
    value: camera,
    options: [
      { value: 'outer', label: 'Outside' },
      { value: 'inner', label: 'Observer' },
    ],
    onChange: (v) => {
      camera = v as 'outer' | 'inner';
      if (view?.kind === 'globe') (view as GlobeView).setView(camera);
    },
  }).root;
  const build = async (mode: 'globe' | 'mollweide'): Promise<void> => {
    const my = ++builds;
    view?.dispose();
    view?.el.remove();
    host.replaceChildren();
    const next = await sky.createSkyView(host, { mode, height: 460, view: camera });
    if (ctx.unmounted || my !== builds) {
      next.dispose();
      next.el.remove();
      return;
    }
    view = next;
    cameraSwitch.hidden = next.kind !== 'globe';
    next.setMap(mapInput);
    FCC_PREREGISTERED_AXES_GAL.value.forEach((lb, i) => next.addPoint(lb, '--series-2', `A${i + 1}`));
  };
  const actions = head.querySelector('.panel-actions');
  if (sky.webglAvailable()) {
    const projection = segmented({
      label: 'Projection',
      value: 'globe',
      options: [
        { value: 'globe', label: 'Globe' },
        { value: 'mollweide', label: 'Mollweide' },
      ],
      onChange: (v) => void build(v as 'globe' | 'mollweide'),
    }).root;
    actions?.prepend(projection, cameraSwitch);
  }
  if (mapRes.ok) {
    foot.prepend(`${mapRes.data.entry.description} (${mapRes.data.meta.source}). `);
  } else {
    foot.prepend(mapRes.reason === 'missing' ? 'sky_map.json not built yet (run qchas build): graticule only. ' : `${mapRes.message}. `);
  }
  await build(sky.webglAvailable() ? 'globe' : 'mollweide');
}
