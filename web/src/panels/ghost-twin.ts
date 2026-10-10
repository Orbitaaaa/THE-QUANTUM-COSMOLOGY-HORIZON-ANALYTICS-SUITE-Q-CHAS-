/**
 * 1.3 Antipodal Ghost-Twin Tracker ([TOP] Sec. 7, Theorem F, Corollary F.1; Secs. 8-9 for the data).
 *
 * The real CMB sky (WMAP 9-yr ILC: the 2 <= l <= 30 reconstruction by default, or the full map) on
 * the Three.js globe or in Mollweide. For any cubic/BCC/FCC torus of size d/D and orientation (the
 * paper's Sec. 9 FCC orientation, Euler angles, or Haar-random) every matched-circle pair is drawn
 * live (physics/circles.ts). Clicking the sky gives the matched twin n' = n + g/chi* when the point
 * lies on a circle, and the ghost images of a source at a chosen distance in that direction, drawn
 * as chords through the observer, with an info card. A live key-fact box makes explicit why, for
 * every size the data allow, topology can only show statistically. SPHEREx QR3 footprints and the
 * six Level-2 thumbnails can be overlaid; clicking a thumbnail turns the globe to its pointing.
 *
 * Bundle fields: sky_map.json (+ binaries), topology.circles.fcc_preregistered,
 * topology.{evidence_torus, exclusion_torus}, cmb.likelihood (our WMAP re-analysis, optional),
 * spherex.{bands, frames, thumbnails} (optional).
 */

import './m1.css';
import './ghost-twin.css';
import { type Theme, sequentialLut } from '../lib/colormap';
import { type LikelihoodRow, type SkyMapData, type SpherexBundle, type SpherexThumbnail, dataUrl, loadSkyMap } from '../lib/data';
import { lb2pix } from '../lib/healpix';
import type { PanelContext } from '../lib/panel';
import { CIRCLE_SEARCH_BOUND, CHI_STAR_MPC, D_LSS_GPC, LIKELIHOOD_BOUNDS } from '../lib/physics/constants';
import { type CirclePair, type GhostImage, CHI_STAR, angularSeparation, circlePairs, ghostImages, lbToVec, vecToLb } from '../lib/physics/circles';
import { cubicBasis } from '../lib/physics/lattice';
import { identity3 } from '../lib/physics/vec';
import type { GlobeView, LB, OverlayId, SkyMapInput, SkyPick, SkyView } from '../lib/sky';
import { currentTheme } from '../lib/theme';
import {
  append,
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
  slider,
  statusPill,
  toggle,
} from '../lib/ui';
import { asMat3, bandKey, comovingDistanceMpc, eulerFromMatrix, eulerZYZ, haarRotation, imageThresholdChi, redshiftForDistance, snapToCircle } from './m1-logic';
import { SHAPE_SLOT, kv, proofList, rawTable } from './m1-ui';

export const title = 'Ghost-Twin Tracker';
export const subtitle =
  'Matched circles, twin points and ghost images of a cubic, BCC or FCC torus drawn live on the real CMB sky, with the sizes the data already exclude flagged.';

type Shape = 'SC' | 'BCC' | 'FCC';
const SHAPE_LABEL: Record<Shape, string> = { SC: 'cubic', BCC: 'BCC', FCC: 'FCC' };
const TABLE_SHAPE: Record<Shape, string> = { SC: 'cubic', BCC: 'BCC', FCC: 'FCC' };
const D_MIN = 0.5;
const D_MAX = 1.5;
const SNAP_TOL_DEG = 2;
const MAX_CHORDS = 6;
const CHI_GPC = CHI_STAR_MPC.value / 1e3;

const likelihoodBound = (s: Shape): number => LIKELIHOOD_BOUNDS.value[TABLE_SHAPE[s] as keyof typeof LIKELIHOOD_BOUNDS.value];

interface Selection {
  l: number;
  b: number;
  /** the pick snapped onto a circle, with its twin */
  hit: ReturnType<typeof snapToCircle>;
}

export async function mount(el: HTMLElement, ctx: PanelContext): Promise<void> {
  ctx.setHeaderExtras([badge('proved', '[TOP] Theorem F'), badge('measured', ctx.manifest?.cmb_map_source ?? 'CMB map'), badge('live', 'physics/circles.ts')]);
  const [topoRes, mapRes, sxRes, cmbRes] = await Promise.all([
    ctx.load('topology.json'),
    loadSkyMap('lowl'),
    ctx.load('spherex.json'),
    ctx.load('cmb.json'),
  ]);
  if (ctx.unmounted) return;
  if (!mapRes.ok) {
    const p = panel({ title: 'Ghost-Twin Tracker', span: 12 });
    p.body.append(loadFailure(mapRes));
    el.append(grid(p.root));
    return;
  }
  const topo = topoRes.ok ? topoRes.data : null;
  const spherex = sxRes.ok ? sxRes.data : null;
  const likelihood = cmbRes.ok ? (cmbRes.data.likelihood ?? []) : [];
  // the map our own likelihood ran on (WMAP9-ILC in the sandbox build, Planck PR3 when reachable)
  const ourMap = (cmbRes.ok ? (cmbRes.data.likelihood_meta?.map ?? cmbRes.data.source) : null) ?? ctx.manifest?.cmb_map_source ?? 'the build map';
  const lnBText = (x: number | null | undefined): string => (x === null || x === undefined || !Number.isFinite(x) ? '—' : `${x < 0 ? '−' : '+'}${Math.abs(x).toFixed(2)}`);
  const preR = asMat3(topo?.circles?.fcc_preregistered?.R);

  // --- state ---------------------------------------------------------------------------------------
  const st = {
    shape: 'FCC' as Shape,
    d: 1.05,
    orient: (preR ? 'paper' : 'euler') as 'paper' | 'euler' | 'random',
    R: preR ?? identity3(),
    euler: preR ? eulerFromMatrix(preR) : ([0, 0, 0] as [number, number, number]),
    frac: 1,
    mapName: 'lowl',
    mode: 'outer' as 'outer' | 'inner' | 'mollweide',
    showCircles: true,
    showFoot: false,
    sel: null as Selection | null,
    thumb: null as SpherexThumbnail | null,
  };
  const maps = new Map<string, SkyMapData>([['lowl', mapRes.data]]);
  let pairs: CirclePair[] = [];
  let pairsError: string | null = null;

  // --- layout ---------------------------------------------------------------------------------------
  const skyPanel = panel({
    title: 'The real CMB sky with a torus overlaid',
    subtitle: 'Drag to rotate, click to place a point. Circles, twins and images follow every control below in real time.',
    badges: [badge('measured', mapRes.data.meta.source), badge('live')],
    span: 8,
  });
  const side = panel({ title: 'Torus', subtitle: 'Lattice, size and orientation of the fundamental domain.', badges: [badge('live')], span: 4 });
  side.root.classList.add('m1-ghost-side');
  el.append(grid(skyPanel.root, side.root));

  const infoPanel = panel({ title: 'Clicked point', subtitle: 'Twin on a matched circle and the ghost images of a source in that direction.', badges: [badge('live', 'circles.ts')], span: 7 });
  const factPanel = panel({ title: 'Why topology can only show statistically', subtitle: 'Computed live for the current torus.', badges: [badge('proved', '[TOP] Cor. F.1'), badge('live')], span: 5 });
  factPanel.root.classList.add('m1-keyfact');
  el.append(grid(infoPanel.root, factPanel.root));

  // sky toolbar
  const viewSel = segmented({
    label: 'View',
    options: [
      { value: 'outer', label: 'Globe' },
      { value: 'inner', label: 'Observer inside' },
      { value: 'mollweide', label: 'Mollweide' },
    ],
    value: st.mode,
    onChange: (v) => void setMode(v as typeof st.mode),
  });
  const mapSel = segmented({
    label: 'Map',
    options: mapRes.data.meta.maps.map((m) => ({ value: m.name, label: m.name === 'lowl' ? '2 ≤ ℓ ≤ 30' : m.name === 'full' ? 'Full map' : m.name })),
    value: st.mapName,
    onChange: (v) => void setMapName(v),
  });
  const circTog = toggle({ label: 'Circles', checked: st.showCircles, onChange: (v) => ((st.showCircles = v), drawCircles()) });
  const footTog = toggle({
    label: spherex ? `SPHEREx footprints (${spherex.frames.length})` : 'SPHEREx footprints (bundle missing)',
    checked: false,
    onChange: (v) => ((st.showFoot = v), drawFootprints()),
  });
  footTog.input.disabled = !spherex;
  const skyHost = h('div', { class: 'm1-sky-host' });
  const mapLabel = h('div', { class: 'chart-caption' });
  const legend = h('div', { class: 'm1-legend' });
  skyPanel.body.append(
    h('div', { class: 'm1-toolbar' }, viewSel.root, mapSel.root, h('div', { class: 'm1-layer-row' }, circTog.root, footTog.root)),
    skyHost,
    legend,
    mapLabel,
  );

  // side controls
  const shapeSel = segmented({
    label: 'Lattice',
    options: (['SC', 'BCC', 'FCC'] as Shape[]).map((s) => ({ value: s, label: SHAPE_LABEL[s] })),
    value: st.shape,
    onChange: (v) => {
      st.shape = v as Shape;
      update();
    },
  });
  const dSlider = slider({
    label: 'Size d/D (shortest loop / last-scattering diameter)',
    min: D_MIN,
    max: D_MAX,
    step: 0.005,
    value: st.d,
    format: (v) => v.toFixed(3),
    onInput: (v) => {
      st.d = v;
      update();
    },
  });
  const ruler = h('div', { class: 'm1-ruler', 'aria-hidden': 'true' });
  const sizePresets = h(
    'div',
    { class: 'm1-jump', role: 'group', 'aria-label': 'Size presets' },
    ...[
      { label: '0.8 D (Fig. 4)', v: 0.8 },
      { label: '0.985 D (circle bound)', v: 0.985 },
      { label: '1.05 D (Sec. 9 fit)', v: 1.05 },
    ].map((p) => {
      const b = h('button', { type: 'button', class: 'btn', title: p.v < 1 ? 'An excluded size, shown so the circles are visible' : 'An allowed size: no circles' }, p.label);
      b.addEventListener('click', () => dSlider.set(p.v, true));
      return b;
    }),
  );
  const sizeStatus = h('div', { class: 'row small gt-status' });
  const orientSel = segmented({
    label: 'Orientation',
    options: [
      ...(preR ? [{ value: 'paper', label: 'Paper (Sec. 9)' }] : []),
      { value: 'euler', label: 'Euler angles' },
      { value: 'random', label: 'Random' },
    ],
    value: st.orient,
    onChange: (v) => setOrient(v as typeof st.orient),
  });
  const eulerSliders = (['α', 'β', 'γ'] as const).map((name, i) =>
    slider({
      label: `${name} (ZYZ, degrees)`,
      min: i === 1 ? 0 : 0,
      max: i === 1 ? 180 : 360,
      step: 0.5,
      value: st.euler[i]!,
      format: (v) => `${v.toFixed(1)}°`,
      onInput: (v) => {
        st.euler[i] = v;
        st.R = eulerZYZ(...st.euler);
        update();
      },
    }),
  );
  const eulerBox = h('div', { class: 'm1-euler', hidden: st.orient !== 'euler' }, ...eulerSliders.map((s) => s.root));
  const reroll = h('button', { type: 'button', class: 'btn', hidden: st.orient !== 'random' }, 'Draw another random orientation');
  reroll.addEventListener('click', () => setOrient('random'));
  const orientNote = h('div', { class: 'small muted' });
  const distSlider = slider({
    label: 'Source distance (fraction of χ*)',
    min: 0.05,
    max: 1,
    step: 0.005,
    value: st.frac,
    format: (v) => `${v.toFixed(3)} χ*`,
    onInput: (v) => {
      st.frac = v;
      updateSelection();
    },
  });
  const zGalaxy = comovingDistanceMpc(1) / CHI_STAR_MPC.value;
  const distPresets = h(
    'div',
    { class: 'm1-jump', role: 'group', 'aria-label': 'Source distance presets' },
    ...[
      { label: 'Last scattering', v: 1 },
      { label: `z ≈ 1 galaxy (${zGalaxy.toFixed(2)})`, v: Math.round(zGalaxy / 0.005) * 0.005 },
      { label: `z ≈ 0.1 (${(comovingDistanceMpc(0.1) / CHI_STAR_MPC.value).toFixed(3)})`, v: Math.max(0.05, Math.round(comovingDistanceMpc(0.1) / CHI_STAR_MPC.value / 0.005) * 0.005) },
    ].map((p) => {
      const b = h('button', { type: 'button', class: 'btn' }, p.label);
      b.addEventListener('click', () => distSlider.set(p.v, true));
      return b;
    }),
  );
  const distInfo = h('div', { class: 'small muted' });
  side.body.append(
    h('div', { class: 'm1-field' }, h('span', { class: 'control-label' }, 'Lattice'), shapeSel.root),
    h('div', { class: 'm1-field' }, dSlider.root, ruler, sizePresets, sizeStatus),
    h('div', { class: 'm1-field' }, h('span', { class: 'control-label' }, 'Orientation'), orientSel.root, eulerBox, reroll, orientNote),
    h('div', { class: 'm1-field' }, distSlider.root, distPresets, distInfo),
  );

  // --- sky view ---------------------------------------------------------------------------------------
  let sky: SkyView | null = null;
  let offPick: (() => void) | null = null;
  const layers: Record<'circles' | 'foot' | 'sel' | 'thumb', OverlayId[]> = { circles: [], foot: [], sel: [], thumb: [] };
  let building = 0;

  async function buildView(): Promise<void> {
    const token = ++building;
    offPick?.();
    sky?.dispose();
    sky = null;
    clear(skyHost);
    (Object.keys(layers) as (keyof typeof layers)[]).forEach((k) => (layers[k] = []));
    const { createSkyView } = await import('../lib/sky');
    if (ctx.unmounted || token !== building) return;
    const view = await createSkyView(skyHost, {
      mode: st.mode === 'mollweide' ? 'mollweide' : 'globe',
      view: st.mode === 'inner' ? 'inner' : 'outer',
      height: st.mode === 'mollweide' ? undefined : 540,
    });
    if (ctx.unmounted || token !== building) return view.dispose();
    sky = view;
    sky.setMap(mapInput());
    offPick = sky.onPick(onPick);
    drawCircles();
    drawFootprints();
    drawSelection();
    drawThumb();
  }
  ctx.onCleanup(() => {
    building++;
    offPick?.();
    sky?.dispose();
    sky = null;
  });

  function mapInput(): SkyMapInput {
    const m = maps.get(st.mapName)!;
    return {
      map: m.map,
      // the low-l map is a full-sky reconstruction from the a_lm; the mask matters for the raw map
      mask: st.mapName === 'lowl' ? null : m.mask,
      limit: m.entry.range?.[1],
      unit: 'μK',
    };
  }

  function setLayer(name: keyof typeof layers, adds: (() => OverlayId)[]): void {
    if (!sky) return;
    for (const id of layers[name]) sky.removeOverlay(id);
    layers[name] = adds.map((f) => f());
  }

  async function setMode(mode: typeof st.mode): Promise<void> {
    const wasGlobe = st.mode !== 'mollweide';
    st.mode = mode;
    if (sky && wasGlobe && mode !== 'mollweide' && sky.kind === 'globe') {
      (sky as GlobeView).setView(mode === 'inner' ? 'inner' : 'outer');
      return;
    }
    await buildView();
  }

  async function setMapName(name: string): Promise<void> {
    if (!maps.has(name)) {
      skyHost.classList.add('m1-busy');
      const r = await loadSkyMap(name);
      skyHost.classList.remove('m1-busy');
      if (ctx.unmounted) return;
      if (!r.ok) {
        mapLabel.textContent = `Could not load the ${name} map: ${r.message ?? r.reason}`;
        return;
      }
      maps.set(name, r.data);
    }
    st.mapName = name;
    sky?.setMap(mapInput());
    renderMapLabel();
    renderInfo();
  }

  function renderMapLabel(): void {
    const m = maps.get(st.mapName)!;
    const mf = ctx.manifest?.cmb_map_source;
    mapLabel.replaceChildren(
      h('strong', null, `${m.meta.source}: `),
      `${m.entry.description}. HEALPix Nside ${m.meta.nside}, Galactic, colour limit ±${fmtNum(m.entry.range?.[1] ?? null, 3)} μK. `,
      st.mapName !== 'lowl' && m.mask ? `Hatched: ${m.meta.mask?.name ?? 'mask'} (f_sky ${fmtNum(m.meta.mask?.fsky ?? null, 3)}). ` : '',
      !mf
        ? ''
        : /^planck/i.test(mf)
          ? `Our CMB numbers are an independent recomputation on ${mf}; the paper's Planck results are quoted next to them.`
          : `Our CMB numbers are an independent re-analysis on ${mf}; Planck PR3 maps were not reachable from the build host, so the paper's Planck results are quoted, not reproduced.`,
    );
  }

  // --- overlays -------------------------------------------------------------------------------------
  const shapeColor = (): string => `--series-${SHAPE_SLOT[st.shape] ?? 3}`;

  function drawCircles(): void {
    if (!sky) return;
    const s = sky;
    const labelled = pairs.length <= 7;
    setLayer(
      'circles',
      st.showCircles
        ? pairs.flatMap((p, k) => [
            () => s.addSmallCircle(p.centre_minus_lb, p.radius_deg, shapeColor(), { label: labelled ? `${k + 1}` : undefined }),
            () => s.addSmallCircle(p.centre_plus_lb, p.radius_deg, shapeColor(), { label: labelled ? `${k + 1}′` : undefined }),
          ])
        : [],
    );
  }

  function bandColors(theme: Theme): Map<string, string> {
    const lut = sequentialLut(theme);
    const keys = spherex ? Object.keys(spherex.bands).sort() : [];
    const out = new Map<string, string>();
    keys.forEach((k, i) => {
      // skip the extremes, which sit too close to the surface colours
      const idx = Math.round(70 + (i / Math.max(1, keys.length - 1)) * 175);
      out.set(k, `rgb(${lut[3 * idx]}, ${lut[3 * idx + 1]}, ${lut[3 * idx + 2]})`);
    });
    return out;
  }

  function drawFootprints(): void {
    if (!sky || !spherex) return;
    const s = sky;
    const colors = bandColors(currentTheme());
    const t0 = performance.now();
    setLayer(
      'foot',
      st.showFoot
        ? spherex.frames.map((f) => () => s.addPolygon(f.corners_gal, colors.get(bandKey(f.band_um, spherex.bands) ?? '') ?? '--series-other'))
        : [],
    );
    // build time of the footprint layer, for profiling
    if (st.showFoot) skyHost.dataset.footMs = (performance.now() - t0).toFixed(0);
    renderLegend();
  }

  function drawSelection(): void {
    if (!sky) return;
    const s = sky;
    const sel = st.sel;
    if (!sel) return setLayer('sel', []);
    const adds: (() => OverlayId)[] = [];
    const src: LB = sel.hit ? vecToLb(sel.hit.point) : [sel.l, sel.b];
    adds.push(() => s.addPoint(src, '--text-1', 'n'));
    if (sel.hit) {
      const tw = vecToLb(sel.hit.twin);
      // chords carry their label at the first end, so start at the twin
      adds.push(() => s.addChordThroughOrigin(tw, src, '--series-1', { label: 'twin n′' }));
    }
    const imgs = images();
    imgs.slice(0, MAX_CHORDS).forEach((im, k) => {
      adds.push(() => s.addChordThroughOrigin([im.l, im.b], src, im.inside_horizon ? '--series-4' : '--series-other', { label: `${k + 1}` }));
    });
    setLayer('sel', adds);
  }

  function drawThumb(): void {
    if (!sky) return;
    const s = sky;
    const t = st.thumb;
    if (!t) return setLayer('thumb', []);
    const frame = spherex?.frames.find((f) => f.obsid === t.obsid && f.detector === t.detector);
    setLayer('thumb', [
      ...(frame ? [() => s.addPolygon(frame.corners_gal, '--text-1', { label: `D${t.detector}` })] : []),
      () => s.addPoint([t.l, t.b], '--text-1'),
    ]);
  }

  function renderLegend(): void {
    const sw = (color: string, box = false): HTMLElement => h('span', { class: `m1-swatch${box ? ' box' : ''}`, style: `background: ${color.startsWith('--') ? `var(${color})` : color}` });
    const items: HTMLElement[] = [
      h('span', null, sw(`var(${shapeColor()})`), `matched circles, ${SHAPE_LABEL[st.shape]}`),
      h('span', null, sw('var(--series-1)'), 'point and twin'),
      h('span', null, sw('var(--series-4)'), 'ghost image inside the sphere'),
      h('span', null, sw('var(--series-other)'), 'ghost image beyond last scattering'),
    ];
    if (st.showFoot && spherex) {
      const colors = bandColors(currentTheme());
      for (const [k, v] of Object.entries(spherex.bands).sort()) items.push(h('span', null, sw(colors.get(k) ?? '', true), `band ${k}: ${v[0]}–${v[1]} μm`));
    }
    legend.replaceChildren(...items);
  }

  // --- interaction ---------------------------------------------------------------------------------------
  function onPick(p: SkyPick): void {
    const n = lbToVec(p.l, p.b);
    st.sel = { l: p.l, b: p.b, hit: snapToCircle(n, pairs, SNAP_TOL_DEG) };
    drawSelection();
    renderInfo();
  }

  function setOrient(o: typeof st.orient): void {
    st.orient = o;
    orientSel.set(o);
    if (o === 'paper' && preR) st.R = preR;
    if (o === 'random') st.R = haarRotation();
    if (o === 'euler') {
      st.euler = eulerFromMatrix(st.R);
      eulerSliders.forEach((s, i) => s.set(Number(st.euler[i]!.toFixed(1))));
      st.R = eulerZYZ(...st.euler);
    }
    eulerBox.hidden = o !== 'euler';
    reroll.hidden = o !== 'random';
    update();
  }

  // --- recompute ------------------------------------------------------------------------------------------
  let raf = 0;
  let timer = 0;
  function update(): void {
    renderRuler();
    // Every overlay change rebuilds the globe's whole overlay layer, so with the 240 SPHEREx
    // footprints on, circle updates wait until the control has been still for a moment.
    if (st.showFoot && sky?.kind === 'globe') {
      clearTimeout(timer);
      timer = window.setTimeout(recompute, 220);
      return;
    }
    if (raf) return;
    raf = requestAnimationFrame(() => {
      raf = 0;
      recompute();
    });
  }
  ctx.onCleanup(() => {
    cancelAnimationFrame(raf);
    clearTimeout(timer);
  });

  function recompute(): void {
    try {
      pairs = circlePairs(cubicBasis(st.shape), st.d, st.R, CHI_STAR);
      pairsError = null;
    } catch (e) {
      pairs = [];
      pairsError = String(e instanceof Error ? e.message : e);
    }
    if (st.sel) st.sel.hit = snapToCircle(lbToVec(st.sel.l, st.sel.b), pairs, SNAP_TOL_DEG);
    renderRuler();
    renderOrientNote();
    drawCircles();
    drawSelection();
    renderLegend();
    renderFacts();
    renderInfo();
    renderStatusTable();
  }

  function updateSelection(): void {
    renderDistInfo();
    drawSelection();
    renderInfo();
    renderFacts();
  }

  // --- side readouts ------------------------------------------------------------------------------------
  function renderRuler(): void {
    const pos = (v: number): string => `${((v - D_MIN) / (D_MAX - D_MIN)) * 100}%`;
    const bound = Math.max(CIRCLE_SEARCH_BOUND.value, likelihoodBound(st.shape));
    ruler.replaceChildren(
      h('div', { class: 'm1-ruler-track' }, h('div', { class: 'm1-ruler-excl', style: `width: ${pos(bound)}` })),
      ...[0.5, 1, 1.5].map((v) => h('div', { class: 'm1-ruler-mark', style: `left: ${pos(v)}` }, h('span', null, v === 1 ? 'D' : v.toFixed(1)))),
      h('div', { class: 'm1-ruler-mark gt-tick', style: `left: ${pos(CIRCLE_SEARCH_BOUND.value)}`, title: 'Circle searches: d ≥ 0.985 D' }),
      h('div', { class: 'm1-ruler-mark gt-tick', style: `left: ${pos(likelihoodBound(st.shape))}`, title: `[TOP] likelihood bound for ${SHAPE_LABEL[st.shape]}` }),
      h('div', { class: 'm1-ruler-now', style: `left: ${pos(st.d)}` }),
    );
    const lb = likelihoodBound(st.shape);
    let pill: HTMLElement;
    if (st.d < CIRCLE_SEARCH_BOUND.value) pill = statusPill('critical', 'excluded: circle searches and likelihood');
    else if (st.d < lb) pill = statusPill('critical', `excluded: ${SHAPE_LABEL[st.shape]} likelihood bound ${lb} D`);
    else pill = statusPill('good', 'allowed by current data');
    sizeStatus.replaceChildren(
      pill,
      h(
        'span',
        { class: 'muted' },
        `Hatched: excluded at 95 %: circle searches d < ${CIRCLE_SEARCH_BOUND.value} D; [TOP] likelihood d < ${LIKELIHOOD_BOUNDS.value.FCC} D (FCC), ${LIKELIHOOD_BOUNDS.value.BCC} D (BCC), ${LIKELIHOOD_BOUNDS.value.cubic} D (cubic). `,
      ),
      badge('quoted', '[TOP] Sec. 9'),
    );
  }

  function renderOrientNote(): void {
    const e = eulerFromMatrix(st.R);
    const txt =
      st.orient === 'paper'
        ? `The paper's best-fit FCC orientation at 1.05 D on Planck ([TOP] Sec. 9). It is not significant; the paper proposes it as a fixed target for future tests${st.shape !== 'FCC' ? `. Here it is applied to the ${SHAPE_LABEL[st.shape]} lattice` : ''}.`
        : st.orient === 'random'
          ? 'A Haar-random rotation: no rotation-invariant theory can predict the orientation ([TOP] Prop. G).'
          : 'R = Rz(α) Ry(β) Rz(γ) takes the lattice frame to Galactic coordinates.';
    orientNote.replaceChildren(txt, h('br'), `ZYZ: ${e.map((x) => `${x.toFixed(1)}°`).join(', ')}`);
  }

  function renderDistInfo(): void {
    const r = st.frac * CHI_STAR_MPC.value;
    const z = redshiftForDistance(r);
    distInfo.textContent = `${fmtFixed(r / 1e3, 2)} Gpc comoving${st.frac >= 0.999 ? ' (the last-scattering surface, z ≈ 1090)' : `, z ≈ ${fmtNum(z, 3)}`} (Planck 2018 ΛCDM).`;
  }

  // --- info card --------------------------------------------------------------------------------------------
  function images(): GhostImage[] {
    const sel = st.sel;
    if (!sel) return [];
    const dir = sel.hit ? sel.hit.point : lbToVec(sel.l, sel.b);
    try {
      return ghostImages(dir, st.frac * CHI_STAR, cubicBasis(st.shape), st.d, st.R, CHI_STAR, 2 * CHI_STAR).filter((g) => !g.is_source);
    } catch {
      return [];
    }
  }

  const valueAt = (lb: LB): string => {
    const m = maps.get(st.mapName)!;
    const v = m.map[lb2pix(m.meta.nside, lb[0], lb[1])];
    return v === undefined || !Number.isFinite(v) ? '—' : `${fmtNum(v, 3)} μK`;
  };

  function renderInfo(): void {
    const sel = st.sel;
    clear(infoPanel.body);
    if (!sel) {
      infoPanel.body.append(
        h(
          'p',
          { class: 'm1-info-empty' },
          'Click anywhere on the sky. On (or within 2° of) a circle you get the matched twin n′ = n + g/χ* (Corollary F.1); anywhere, the ghost images of a source at the chosen distance in that direction.',
        ),
      );
      return;
    }
    const src: LB = sel.hit ? vecToLb(sel.hit.point) : [sel.l, sel.b];
    const fmtLb = (lb: LB): string => `l ${lb[0].toFixed(1)}°, b ${lb[1] >= 0 ? '+' : '−'}${Math.abs(lb[1]).toFixed(1)}°`;
    const rows: [string, Node | string | (Node | string)[]][] = [['Direction n', [fmtLb(src), ' · ', valueAt(src), ' ', badge('measured', maps.get(st.mapName)!.meta.source)]]];
    if (sel.hit) {
      const p = pairs[sel.hit.index]!;
      const tw = vecToLb(sel.hit.twin);
      rows.push(
        ['On circle', [`pair ${sel.hit.index + 1} of ${pairs.length}, radius ${p.radius_deg.toFixed(2)}°, |g| = ${p.length_over_D.toFixed(3)} D (snapped ${sel.hit.offsetDeg.toFixed(2)}°)`]],
        ['Matched twin n′', [fmtLb(tw), ' · ', valueAt(tw), ' ', badge('live', 'Cor. F.1')]],
        ['Angle n to n′', `${angularSeparation(src, tw).toFixed(2)}°: the same point of the last-scattering surface seen in two directions`],
      );
    } else rows.push(['On a circle?', pairs.length ? `no: the nearest circle is more than ${SNAP_TOL_DEG}° away` : 'no circles at this size']);
    const imgs = images();
    const inside = imgs.filter((g) => g.inside_horizon);
    rows.push([
      'Source',
      `${(st.frac).toFixed(3)} χ* = ${fmtFixed(st.frac * CHI_GPC, 2)} Gpc in direction n; ${imgs.length} ghost images within 2χ*, ${inside.length} inside the observable sphere`,
    ]);
    infoPanel.body.append(kv(rows));
    if (sel.hit)
      infoPanel.body.append(
        h(
          'p',
          { class: 'm1-takeaway' },
          `In a torus of this size and orientation the two temperatures would be equal up to Doppler and integrated Sachs-Wolfe terms. One pair on one map proves nothing either way; [TOP] Sec. 9 tests all pairs and orientations at once and finds no evidence for a torus.`,
        ),
      );
    if (imgs.length) {
      const tbl = rawTable(
        [
          { label: '#', num: true, cell: (g: GhostImage) => String(imgs.indexOf(g) + 1) },
          { label: 'Direction (l, b)', cell: (g) => `${g.l.toFixed(1)}°, ${g.b >= 0 ? '+' : '−'}${Math.abs(g.b).toFixed(1)}°` },
          { label: 'Distance', num: true, cell: (g) => `${(g.distance_over_D / CHI_STAR).toFixed(3)} χ*` },
          { label: 'Gpc', num: true, cell: (g) => fmtFixed(g.distance_over_D * D_LSS_GPC.value, 2) },
          { label: 'Lattice vector n', cell: (g) => h('span', { class: 'm1-mono' }, `(${g.n.join(', ')})`) },
          {
            label: 'Observable?',
            cell: (g) =>
              g.inside_horizon ? statusPill('good', 'inside the sphere') : h('span', { class: 'muted' }, 'beyond last scattering'),
          },
        ],
        imgs.slice(0, 10),
        (g) => (g.inside_horizon ? null : 'm1-row-out'),
      );
      infoPanel.body.append(
        tbl.root,
        h(
          'div',
          { class: 'chart-caption' },
          `Nearest ${Math.min(10, imgs.length)} of ${imgs.length}; the ${Math.min(MAX_CHORDS, imgs.length)} nearest are drawn as chords from n through the observer. An image beyond last scattering (> χ*) cannot be seen in any electromagnetic survey.`,
        ),
      );
    }
  }

  // --- key facts ----------------------------------------------------------------------------------------------
  function renderFacts(): void {
    const thr = imageThresholdChi(st.d);
    const nearest = Math.max(0, 2 * st.d - st.frac);
    const lb = likelihoodBound(st.shape);
    const allowed = st.d >= Math.max(CIRCLE_SEARCH_BOUND.value, lb);
    const radii = pairs.map((p) => p.radius_deg);
    clear(factPanel.body);
    append(
      factPanel.body,
      h(
        'ul',
        null,
        h(
          'li',
          null,
          h('strong', null, 'Circles: '),
          pairsError
            ? `cannot enumerate (${pairsError}).`
            : pairs.length
              ? `${pairs.length} matched pair${pairs.length > 1 ? 's' : ''}, radii ${fmtFixed(Math.min(...radii), 1)}°–${fmtFixed(Math.max(...radii), 1)}° (one pair per ±g with |g| < D; radius arccos(|g|/D)).`
              : `none: the shortest loop d = ${st.d.toFixed(3)} D is at least the last-scattering diameter, so the sphere does not meet any of its images.`,
          ' ',
          badge('live', 'circlePairs'),
        ),
        h(
          'li',
          null,
          h('strong', null, 'Ghost images inside the observable sphere '),
          `need a source at r ≥ (2d/D − 1) χ* = ${thr.toFixed(3)} χ*`,
          thr >= 1 ? ': no source inside the sphere has any image inside it.' : thr > 0.9 ? ': only sources essentially at last scattering.' : '.',
          ' ',
          badge('live', 'imageThresholdChi'),
        ),
        h(
          'li',
          null,
          h('strong', null, `A source at ${st.frac.toFixed(2)} χ*: `),
          `its nearest possible image is at least (2d/D − r/χ*) χ* = ${nearest.toFixed(3)} χ* away`,
          nearest > 1 ? ', beyond the last-scattering surface.' : ', which can lie inside the sphere for the right direction.',
          ` A z ≈ 1 galaxy sits at ${zGalaxy.toFixed(2)} χ*.`,
        ),
        h(
          'li',
          null,
          h('strong', null, 'This size: '),
          allowed ? 'allowed by current data.' : `excluded at 95 % (${st.d < CIRCLE_SEARCH_BOUND.value ? 'circle searches' : `[TOP] likelihood, ${SHAPE_LABEL[st.shape]} ≥ ${lb} D`}).`,
          ' ',
          badge('quoted', '[TOP] Sec. 9'),
        ),
      ),
      h(
        'p',
        { class: 'm1-takeaway' },
        h('strong', null, 'So: '),
        'for every allowed size (d ≥ 0.985 D, and ≥ 1.04 D for FCC) there are no circles once d ≥ D, and the images of nearby sources lie beyond last scattering. Topology can only show statistically, in the correlations of the low multipoles ([TOP] Secs. 8–9), and Planck shows no evidence for a torus of any allowed size. ',
        'Topological resonance would give up to 64× more information, but only at known orientation, and the pre-registered search of Planck finds no sign of it (p = 0.78). Outside these exits strong evidence is possible only for d ≲ 1.1 D; deciding a torus as large as 3 D needs a multi-field primordial line mapped in three dimensions, not the CMB alone.',
      ),
      h('div', { class: 'row small' }, badge('proved', '[TOP] Theorem F, Cor. F.1'), badge('quoted', '[TOP] Secs. 9, 11.6–11.8')),
    );
  }

  // --- data status at this lattice ------------------------------------------------------------------------------
  const statusPanel = panel({
    title: 'What the data say at each size',
    subtitle: `Quoted Planck results ([TOP] Tables 9a, 9b) next to our independent re-analysis on ${ourMap}. The row nearest the slider is highlighted.`,
    badges: [badge('quoted', '[TOP] Table 9'), badge('computed', 'cmb.likelihood')],
    span: 12,
  });
  const statusHost = h('div');
  statusPanel.body.append(statusHost);
  el.append(grid(statusPanel.root));

  interface SizeRow {
    d: number;
    paperLnB: number | null;
    paperVerdict: string | null;
    paperP: string | null;
    ours: LikelihoodRow | null;
  }
  function renderStatusTable(): void {
    const shape = TABLE_SHAPE[st.shape];
    const ex = (topo?.exclusion_torus ?? []).filter((r) => r.shape === shape);
    const ev = (topo?.evidence_torus ?? []).filter((r) => r.shape === shape);
    const ours = likelihood.filter((r) => r.shape === shape);
    const sizes = [...new Set([...ex.map((r) => r.d_over_D), ...ev.map((r) => r.d_over_D), ...ours.map((r) => r.d_over_D)])].sort((a, b) => a - b);
    const rows: SizeRow[] = sizes.map((d) => {
      const e = ex.find((r) => r.d_over_D === d);
      const v = ev.find((r) => r.d_over_D === d);
      return {
        d,
        paperLnB: e?.lnB_planck ?? v?.lnB ?? null,
        paperVerdict: e?.verdict ?? null,
        paperP: v ? String(v.p) : null,
        ours: ours.find((r) => r.d_over_D === d) ?? null,
      };
    });
    const nearest = rows.reduce<SizeRow | null>((b, r) => (!b || Math.abs(r.d - st.d) < Math.abs(b.d - st.d) ? r : b), null);
    const verdictCell = (v: string | null | undefined): HTMLElement | string =>
      !v ? '—' : v.startsWith('excluded') ? statusPill('critical', v) : statusPill('good', v);
    const tbl = rawTable<SizeRow>(
      [
        { label: 'd/D', num: true, cell: (r) => r.d.toFixed(3) },
        { label: 'Planck ln B (paper)', num: true, cell: (r) => lnBText(r.paperLnB) },
        { label: 'p (no torus, paper)', num: true, title: 'Table 9a: fraction of simulated skies without topology with ln B at least as large', cell: (r) => r.paperP ?? '—' },
        { label: 'Verdict (paper)', cell: (r) => verdictCell(r.paperVerdict) },
        { label: `${ourMap} ln B (ours)`, num: true, cell: (r) => lnBText(r.ours?.lnB) },
        { label: 'p (no torus, ours)', num: true, cell: (r) => (r.ours && Number.isFinite(r.ours.p) ? r.ours.p.toFixed(2) : '—') },
        { label: 'Verdict (ours)', cell: (r) => verdictCell((r.ours as unknown as { exclusion?: { verdict?: string } } | null)?.exclusion?.verdict) },
      ],
      rows,
      (r) => (r === nearest ? 'm1-row-flag' : null),
    );
    statusHost.replaceChildren(
      rows.length ? tbl.root : h('div', { class: 'empty-state' }, 'No tabulated sizes for this lattice in the bundles.'),
      h(
        'div',
        { class: 'chart-caption' },
        `Lattice ${SHAPE_LABEL[st.shape]}. ln B is the orientation-averaged Bayes factor against infinite flat space; negative values disfavour the torus. Neither analysis finds evidence for a torus at any allowed size; both exclude the small ones. Our numbers come from a different map (${ourMap}, not the paper's Planck 2018 Commander map) and are not expected to agree in detail.`,
      ),
    );
  }

  // --- SPHEREx gallery --------------------------------------------------------------------------------------------
  if (spherex?.thumbnails.length) mountGallery(el, spherex, (t) => {
    st.thumb = st.thumb === t ? null : t;
    drawThumb();
    if (st.thumb && sky?.kind === 'globe') (sky as GlobeView).lookAt([t.l, t.b]);
    return st.thumb === t;
  });

  // --- what is proved ---------------------------------------------------------------------------------------------
  const proofs = panel({ title: 'What the paper proves, measures and assumes', span: 12 });
  proofs.body.append(
    proofList([
      { kind: 'proved', src: '[TOP] Theorem F, Cor. F.1', text: 'For a 3-torus every lattice vector g with |g| < D gives one pair of antipodal circles of angular radius arccos(|g|/D), matched with zero twist by n′ = n + g/χ*; the pattern is the same for every observer.' },
      { kind: 'proved', src: '[TOP] Prop. G', text: 'No rotation-invariant theory can predict the orientation; size and orientation must be measured.' },
      { kind: 'measured', src: '[TOP] Sec. 9', text: 'Planck 2018: no evidence for a torus of any allowed size; exclusion of small tori at 95 % (FCC d ≥ 1.04 D, BCC ≥ 1.00 D, cubic ≥ 0.985 D). The FCC orientation shown is a best fit that is not significant.' },
      { kind: 'computed', src: 'cmb.likelihood', text: `Our likelihood on ${ourMap} (same method, different map) also finds no evidence.` },
      { kind: 'assumed', src: '[TOP] Sec. 7', text: `Distances use χ* = ${CHI_STAR_MPC.value} Mpc (Planck 2018 ΛCDM), D = 2χ* = ${D_LSS_GPC.value} Gpc.` },
    ]),
    note(
      [
        'Circles exist only for d < D; the circle-search bound d ≥ 0.985 D therefore leaves a sliver 0.985 ≤ d < 1 with circles of radius below 9.9°, and the likelihood bounds close it for FCC and BCC.',
        'Images: |x + g| ≥ |g| − |x| ≥ d − r, so an image of a source at r lies inside χ* only if r ≥ d − χ*; at d ≥ D that needs r ≥ χ*.',
        'Resonance (Sec. 11.6) needs a narrow primordial line on a dual-lattice shell; single-field inflation cannot make it.',
      ],
      'Derivations used in the live box',
    ),
  );
  el.append(grid(proofs.root));

  // --- go -------------------------------------------------------------------------------------------------------------
  renderMapLabel();
  renderDistInfo();
  recompute();
  await buildView();
}

// --- SPHEREx Level-2 gallery ------------------------------------------------------------------------------------------

function mountGallery(el: HTMLElement, S: SpherexBundle, onSelect: (t: SpherexThumbnail) => boolean): void {
  const p = panel({
    title: `SPHEREx Quick-Release Level-2 frames (${S.release.toUpperCase()})`,
    subtitle: `${S.thumbnails.length} real frames read from ${S.bucket} (${S.source ?? 'NASA/IPAC IRSA'}). Click one to turn the globe to its pointing; turn on the footprint layer to see all ${S.n_frames} sampled frames.`,
    badges: [badge('measured', 'SPHEREx QR3')],
    span: 12,
  });
  const buttons: HTMLButtonElement[] = [];
  const cards = S.thumbnails.map((t) => {
    const band = bandKey(t.band_um, S.bands);
    const b = h(
      'button',
      { type: 'button', class: 'm1-thumb', 'aria-pressed': 'false', title: t.key ?? t.obsid },
      h('img', { src: dataUrl(t.file), alt: `SPHEREx frame ${t.obsid}, detector ${t.detector}, ${t.band_um[0]}–${t.band_um[1]} μm`, loading: 'lazy' }),
      h('span', { class: 'm1-thumb-title' }, `${t.obsid} · D${t.detector}`),
      h(
        'span',
        { class: 'm1-thumb-meta' },
        `band ${band ?? '?'}: ${t.band_um[0]}–${t.band_um[1]} μm`,
        h('br'),
        t.date ? `${t.date.slice(0, 16).replace('T', ' ')} UTC` : '',
        h('br'),
        `l ${t.l.toFixed(1)}°, b ${t.b.toFixed(1)}°`,
        h('br'),
        `stretch ${fmtNum(t.p1 as number, 2)}–${fmtNum(t.p995 as number, 2)} ${t.unit}`,
      ),
    );
    b.addEventListener('click', () => {
      const on = onSelect(t);
      buttons.forEach((x) => x.setAttribute('aria-pressed', 'false'));
      b.setAttribute('aria-pressed', String(on));
    });
    buttons.push(b);
    return b;
  });
  p.body.append(
    h('div', { class: 'm1-gallery' }, ...cards),
    h('div', { class: 'chart-caption' }, `${S.note ?? ''} Images are shown with a linear stretch from the 1st to the 99.5th percentile of each frame, in the unit its FITS header gives (BUNIT: ${[...new Set((S.thumbnails ?? []).map((t) => t.unit).filter(Boolean))].join(', ') || 'not stated'}).`),
  );
  el.append(grid(p.root));
}
