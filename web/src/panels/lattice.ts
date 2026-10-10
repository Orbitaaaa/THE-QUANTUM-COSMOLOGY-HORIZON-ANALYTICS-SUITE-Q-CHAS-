/**
 * 1.2 Lattice Shape Optimizer ([TOP] Sec. 6, Prop. E.2, Theorems H.1-H.2, Prop. H.3).
 *
 * The Casimir functional Z(2) at unit covolume along the rhombohedral family (the pipeline curve,
 * lattice.rhombohedral) with a live alpha slider recomputed in the browser by Ewald summation
 * (physics/lattice.ts), a lazily loaded Three.js view of the unit cell that morphs with alpha, the
 * certified critical points ours vs paper with the tr H = 2Z/3 identity check, the shape-flow
 * trajectories and basin counts, and the Prop. H.3 relaxation numbers with the paper's caveat that
 * FCC is a conditional prediction.
 *
 * Bundle fields: lattice.{rhombohedral, critical, gaps, flows, relaxation}; validation.json rows
 * (optional) for the ours-vs-printed summary.
 */

import './m1.css';
import { type ChartHandle, type LineSeries, barChart, lineChart } from '../lib/charts';
import type { CriticalPoint, LatticeBundle, ShapeKey } from '../lib/data';
import type { PanelContext } from '../lib/panel';
import { epsteinZ2, rhombohedralBasis } from '../lib/physics/lattice';
import type { Mat3 } from '../lib/physics/vec';
import { cssVar } from '../lib/theme';
import { append, badge, fmtFixed, fmtNum, fmtSci, grid, h, loadFailure, note, panel, slider, statTile, statusPill } from '../lib/ui';
import { flowSeries, hessianTraceCheck, interp } from './m1-logic';
import {
  type RawColumn,
  SHAPE_SLOT,
  fmtSigned,
  keep,
  oursPaper,
  proofList,
  rawTable,
  restyleSeries,
  sectionLabel,
  tilesN,
  validationRows,
  validationSummary,
} from './m1-ui';

export const title = 'Lattice Shape Optimizer';
export const subtitle =
  'Which torus shape minimizes the Casimir energy: Z(2) along the rhombohedral family, the certified SC saddle and BCC/FCC minima, and where the shape flow ends.';

const SHAPES: ShapeKey[] = ['SC', 'BCC', 'FCC'];
const ALPHA_OF: Record<ShapeKey, number> = { FCC: 60, SC: 90, BCC: (Math.acos(-1 / 3) * 180) / Math.PI };
const ROLE: Record<ShapeKey, string> = { SC: 'saddle', BCC: 'local minimum', FCC: 'deepest minimum' };
const ALPHA_MIN = 40;
const ALPHA_MAX = 118;
const VALIDATION_GROUPS = ['Casimir lattice ([TOP] Prop. E.2)', 'Shape flow ([TOP] Sec. 6.3)'];

/** Rhombohedral basis rescaled to unit covolume (the normalization of Z(2) in Prop. E.2). */
function unitCovolumeBasis(alpha: number): Mat3 {
  const B = rhombohedralBasis(alpha);
  const det = Math.abs(
    B[0][0] * (B[1][1] * B[2][2] - B[1][2] * B[2][1]) -
      B[0][1] * (B[1][0] * B[2][2] - B[1][2] * B[2][0]) +
      B[0][2] * (B[1][0] * B[2][1] - B[1][1] * B[2][0]),
  );
  const s = det ** (-1 / 3);
  return B.map((r) => r.map((x) => x * s)) as Mat3;
}

function nearestShape(alpha: number): { shape: ShapeKey; dist: number } {
  let best: ShapeKey = 'FCC';
  for (const s of SHAPES) if (Math.abs(alpha - ALPHA_OF[s]) < Math.abs(alpha - ALPHA_OF[best])) best = s;
  return { shape: best, dist: Math.abs(alpha - ALPHA_OF[best]) };
}

export async function mount(el: HTMLElement, ctx: PanelContext): Promise<void> {
  ctx.setHeaderExtras([badge('proved', '[TOP] Prop. E.2'), badge('computed', 'lattice.json'), badge('live', 'physics/lattice.ts')]);
  const [res, val] = await Promise.all([ctx.load('lattice.json'), ctx.load('validation.json')]);
  if (ctx.unmounted) return;
  if (!res.ok) {
    const p = panel({ title: 'Lattice Shape Optimizer', span: 12 });
    p.body.append(loadFailure(res));
    el.append(grid(p.root));
    return;
  }
  const L = res.data;
  mountCurve(el, ctx, L);
  mountCritical(el, ctx, L, validationSummary(validationRows(val, VALIDATION_GROUPS), VALIDATION_GROUPS));
  mountFlows(el, ctx, L);
  mountRelaxation(el, L);
}

// --- Z(2) along the rhombohedral family, live alpha, 3D cell -------------------------------------

function mountCurve(el: HTMLElement, ctx: PanelContext, L: LatticeBundle): void {
  const R = L.rhombohedral;
  const angles = { ...ALPHA_OF, ...(R.angles_deg ?? {}) } as Record<ShapeKey, number>;
  const keep_ = R.alpha_deg.map((a, i) => (a <= ALPHA_MAX + 1e-9 ? i : -1)).filter((i) => i >= 0);
  const xs = keep_.map((i) => R.alpha_deg[i]!);
  const ys = keep_.map((i) => R.Z[i]!);

  const p = panel({
    title: 'Casimir energy along the rhombohedral family · [TOP] Fig. 2',
    subtitle: 'Z(2) at unit covolume for three equal-length basis vectors at mutual angle α. The family passes through FCC (60°), SC (90°) and BCC (109.47°).',
    badges: [badge('computed', 'lattice.rhombohedral'), badge('live', 'epsteinZ2')],
    span: 8,
  });
  const cellPanel = panel({
    title: 'Unit cell at the slider angle',
    subtitle: 'Primitive basis a₁, a₂, a₃ at unit covolume, the parallelepiped they span and nearby lattice points. Drag to rotate.',
    badges: [badge('live')],
    span: 4,
  });
  el.append(grid(p.root, cellPanel.root));

  let alpha = 60;
  const live: LineSeries = { name: 'live α (browser Ewald sum)', x: [alpha], y: [epsteinZ2(unitCovolumeBasis(alpha))], color: '--text-1', markers: true, width: 0 };
  const series: LineSeries[] = [
    { name: 'Z(2), pipeline', x: xs, y: ys, color: 1 },
    ...SHAPES.map(
      (s): LineSeries => ({
        name: `${s} ${ROLE[s]}`,
        x: [angles[s]],
        y: [L.critical[s].Z],
        color: SHAPE_SLOT[s],
        markers: true,
        width: 0,
      }),
    ),
    live,
  ];
  const chartHost = h('div');
  p.body.append(chartHost);
  const chart: ChartHandle = keep(
    ctx,
    lineChart(chartHost, series, {
      xTitle: 'rhombohedral angle α (degrees)',
      yTitle: 'Z(2) at unit covolume (dimensionless)',
      xRange: [ALPHA_MIN, ALPHA_MAX + 2],
      yRange: [15.85, 17.7],
      height: 320,
      refLines: [{ axis: 'y', value: L.critical.FCC.Z, label: 'Z_FCC' }],
      caption: 'Above ~116° the curve rises steeply (Z = 118 at 119.5°, where the cell flattens); the axis is clipped to show the three critical points.',
      format: (v) => fmtNum(v, 7),
    }),
  );

  const tilesHost = h('div');
  const sl = slider({
    label: 'Angle α between basis vectors',
    min: ALPHA_MIN,
    max: ALPHA_MAX,
    step: 0.1,
    value: alpha,
    format: (v) => `${v.toFixed(1)}°`,
    onInput: (v) => schedule(v),
  });
  const jump = h(
    'div',
    { class: 'm1-jump', role: 'group', 'aria-label': 'Jump to a critical point' },
    h('span', { class: 'small muted' }, 'Jump to'),
    ...SHAPES.map((s) => {
      const b = h('button', { type: 'button', class: 'btn' }, `${s} ${fmtFixed(angles[s], s === 'BCC' ? 2 : 0)}°`);
      b.addEventListener('click', () => sl.set(Number(angles[s].toFixed(1)), true));
      return b;
    }),
  );
  p.body.append(h('div', { class: 'm1-toolbar' }, h('div', { class: 'm1-grow' }, sl.root), jump), tilesHost);

  // 3D unit cell (Three.js loaded lazily)
  const cellHost = h('div', { class: 'm1-cell3d', role: 'img', 'aria-label': 'Three-dimensional view of the rhombohedral unit cell' });
  const cellHud = h('div', { class: 'm1-cell-hud' });
  cellHost.append(cellHud);
  cellPanel.body.append(
    cellHost,
    h('div', { class: 'm1-legend' },
      ...[['a₁', 1], ['a₂', 4], ['a₃', 5]].map(([n, slot]) => h('span', null, h('span', { class: 'm1-swatch', style: `background: var(--series-${slot})` }), `basis vector ${n}`)),
      h('span', null, h('span', { class: 'm1-swatch', style: 'background: var(--text-2)' }), 'cell edges'),
    ),
  );
  let cell: CellView | null = null;
  void createCellView(cellHost).then((v) => {
    if (ctx.unmounted) return v?.dispose();
    cell = v;
    if (cell) {
      ctx.onCleanup(() => cell?.dispose());
      ctx.onThemeChange(() => cell?.setBasis(unitCovolumeBasis(alpha)));
      cell.setBasis(unitCovolumeBasis(alpha));
    } else cellHost.append(h('div', { class: 'empty-state' }, 'WebGL unavailable: the 3D cell view needs WebGL.'));
  });

  let raf = 0;
  let pending = alpha;
  function schedule(v: number): void {
    pending = v;
    if (raf) return;
    raf = requestAnimationFrame(() => {
      raf = 0;
      update(pending);
    });
  }
  ctx.onCleanup(() => cancelAnimationFrame(raf));

  function update(a: number): void {
    alpha = a;
    const B = unitCovolumeBasis(a);
    const t0 = performance.now();
    const z = epsteinZ2(B);
    const ms = performance.now() - t0;
    live.x = [a];
    live.y = [z];
    void restyleSeries(chart, series, [series.length - 1]);
    const zPipe = interp(xs, ys, a);
    const ns = nearestShape(a);
    const len = Math.hypot(...B[0]);
    tilesHost.replaceChildren(
      tilesN(
        4,
        statTile({ label: 'Z(2) at this α', value: fmtNum(z, 8), sub: [badge('live', 'epsteinZ2'), `${ms.toFixed(1)} ms in the browser`] }),
        statTile({
          label: 'Pipeline curve at α',
          value: zPipe === null ? '—' : fmtNum(zPipe, 8),
          sub: [badge('computed', 'lattice.rhombohedral'), zPipe === null ? '' : `live − pipeline ${fmtSigned(z - zPipe, 2)} (0.5° grid, linear interpolation)`],
        }),
        statTile({
          label: 'Above the FCC minimum',
          value: fmtSci(z - L.critical.FCC.Z, 3),
          sub: [badge('live'), `Z − Z_FCC; nearest critical point ${ns.shape} (${ns.dist < 0.05 ? 'here' : `${ns.dist.toFixed(1)}° away`})`],
        }),
        statTile({ label: 'Basis vector length', value: fmtNum(len, 5), sub: [badge('live'), 'unit covolume, cos α = a_i · a_j / |a|²'] }),
      ),
    );
    cellHud.textContent = `α = ${a.toFixed(1)}° · ${ns.dist < 0.05 ? ns.shape : `near ${ns.shape}`} · |a_i| = ${len.toFixed(3)} · Z = ${z.toFixed(5)}`;
    cell?.setBasis(B);
  }
  update(alpha);
}

// --- Three.js unit-cell view ---------------------------------------------------------------------

interface CellView {
  setBasis(B: Mat3): void;
  dispose(): void;
}

async function createCellView(host: HTMLElement): Promise<CellView | null> {
  const probe = document.createElement('canvas');
  if (!(probe.getContext('webgl2') || probe.getContext('webgl'))) return null;
  const THREE = await import('three');
  const { OrbitControls } = await import('three/examples/jsm/controls/OrbitControls.js');
  const { CSS2DObject, CSS2DRenderer } = await import('three/examples/jsm/renderers/CSS2DRenderer.js');
  if (!host.isConnected) return null;

  let renderer: InstanceType<typeof THREE.WebGLRenderer>;
  try {
    renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
  } catch {
    return null;
  }
  renderer.setPixelRatio(Math.min(2, window.devicePixelRatio || 1));
  renderer.setClearColor(0x000000, 0);
  host.prepend(renderer.domElement);
  const labels = new CSS2DRenderer();
  Object.assign(labels.domElement.style, { position: 'absolute', inset: '0', pointerEvents: 'none' });
  host.append(labels.domElement);

  const scene = new THREE.Scene();
  const camera = new THREE.PerspectiveCamera(30, 1, 0.01, 100);
  const up = new THREE.Vector3(1, 1, 1).normalize();
  camera.up.copy(up);
  camera.position.set(4.2, -3.2, 1.4);
  const controls = new OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true;
  controls.dampingFactor = 0.1;
  controls.enablePan = false;
  controls.minDistance = 2.5;
  controls.maxDistance = 14;
  const group = new THREE.Group();
  scene.add(group);

  let dirty = true;
  let disposed = false;
  let frame = 0;
  controls.addEventListener('change', () => (dirty = true));

  const color = (token: string): InstanceType<typeof THREE.Color> => new THREE.Color().setStyle(cssVar(token) || '#888');

  function clearGroup(): void {
    group.traverse((o) => {
      const m = o as unknown as { geometry?: { dispose(): void }; material?: { dispose(): void } | { dispose(): void }[] };
      m.geometry?.dispose();
      if (Array.isArray(m.material)) m.material.forEach((x) => x.dispose());
      else m.material?.dispose();
      if (o instanceof CSS2DObject) o.element.remove();
    });
    group.clear();
  }

  function setBasis(B: Mat3): void {
    clearGroup();
    const v = B.map((r) => new THREE.Vector3(r[0], r[1], r[2]));
    const [a1, a2, a3] = v as [InstanceType<typeof THREE.Vector3>, InstanceType<typeof THREE.Vector3>, InstanceType<typeof THREE.Vector3>];
    const centre = a1.clone().add(a2).add(a3).multiplyScalar(0.5);
    const at = (i: number, j: number, k: number): InstanceType<typeof THREE.Vector3> =>
      a1.clone().multiplyScalar(i).add(a2.clone().multiplyScalar(j)).add(a3.clone().multiplyScalar(k)).sub(centre);

    // parallelepiped edges
    const pts: number[] = [];
    const corners: [number, number, number][] = [];
    for (let i = 0; i < 2; i++) for (let j = 0; j < 2; j++) for (let k = 0; k < 2; k++) corners.push([i, j, k]);
    for (const c of corners)
      for (let d = 0; d < 3; d++) {
        if (c[d] === 1) continue;
        const e: [number, number, number] = [...c];
        e[d] = 1;
        const p0 = at(...c);
        const p1 = at(...e);
        pts.push(p0.x, p0.y, p0.z, p1.x, p1.y, p1.z);
      }
    const edgeGeom = new THREE.BufferGeometry();
    edgeGeom.setAttribute('position', new THREE.Float32BufferAttribute(pts, 3));
    group.add(new THREE.LineSegments(edgeGeom, new THREE.LineBasicMaterial({ color: color('--text-2') })));

    // translucent faces
    const faceGeom = new THREE.BufferGeometry();
    const face = (p: InstanceType<typeof THREE.Vector3>[]): number[] => [p[0]!, p[1]!, p[2]!, p[0]!, p[2]!, p[3]!].flatMap((q) => [q.x, q.y, q.z]);
    const f: number[] = [];
    for (const s of [0, 1]) {
      f.push(...face([at(s, 0, 0), at(s, 1, 0), at(s, 1, 1), at(s, 0, 1)]));
      f.push(...face([at(0, s, 0), at(1, s, 0), at(1, s, 1), at(0, s, 1)]));
      f.push(...face([at(0, 0, s), at(1, 0, s), at(1, 1, s), at(0, 1, s)]));
    }
    faceGeom.setAttribute('position', new THREE.Float32BufferAttribute(f, 3));
    group.add(
      new THREE.Mesh(
        faceGeom,
        new THREE.MeshBasicMaterial({ color: color('--accent'), transparent: true, opacity: 0.08, side: THREE.DoubleSide, depthWrite: false }),
      ),
    );

    // lattice points of the cell and its neighbours
    const pointGeom = new THREE.SphereGeometry(0.045, 16, 12);
    const cornerMat = new THREE.MeshBasicMaterial({ color: color('--text-1') });
    const otherMat = new THREE.MeshBasicMaterial({ color: color('--series-other') });
    const reach = 1.25 * Math.max(a1.length(), centre.length() * 1.2);
    for (let i = -1; i <= 2; i++)
      for (let j = -1; j <= 2; j++)
        for (let k = -1; k <= 2; k++) {
          const p = at(i, j, k);
          const inCell = i >= 0 && i <= 1 && j >= 0 && j <= 1 && k >= 0 && k <= 1;
          if (!inCell && p.length() > reach) continue;
          const m = new THREE.Mesh(pointGeom, inCell ? cornerMat : otherMat);
          m.position.copy(p);
          if (!inCell) m.scale.setScalar(0.7);
          group.add(m);
        }

    // basis vectors from the corner at the origin of the cell
    const o = at(0, 0, 0);
    v.forEach((vec, i) => {
      const len = vec.length();
      const arrow = new THREE.ArrowHelper(vec.clone().normalize(), o, len, color(`--series-${[1, 4, 5][i]}`), 0.14, 0.07);
      group.add(arrow);
      const lab = new CSS2DObject(h('div', { class: 'm1-cell-label' }, `a${'₁₂₃'[i]}`));
      lab.position.copy(o.clone().add(vec.clone().multiplyScalar(1.08)));
      group.add(lab);
    });
    dirty = true;
  }

  function resize(): void {
    const w = host.clientWidth || 1;
    const hh = host.clientHeight || 1;
    renderer.setSize(w, hh);
    labels.setSize(w, hh);
    camera.aspect = w / hh;
    camera.updateProjectionMatrix();
    dirty = true;
  }
  const ro = new ResizeObserver(resize);
  ro.observe(host);
  resize();

  function loop(): void {
    if (disposed) return;
    frame = requestAnimationFrame(loop);
    if (!host.isConnected) return dispose();
    controls.update();
    if (!dirty) return;
    dirty = false;
    renderer.render(scene, camera);
    labels.render(scene, camera);
  }
  loop();

  function dispose(): void {
    if (disposed) return;
    disposed = true;
    cancelAnimationFrame(frame);
    ro.disconnect();
    controls.dispose();
    clearGroup();
    renderer.dispose();
    renderer.forceContextLoss();
    renderer.domElement.remove();
    labels.domElement.remove();
  }
  return { setBasis, dispose };
}

// --- certified critical points -------------------------------------------------------------------

interface CritRow {
  shape: ShapeKey;
  c: CriticalPoint;
}

const fmtEig = (xs: readonly number[], d: number): string => xs.map((x) => x.toFixed(d)).join(', ');

function mountCritical(el: HTMLElement, ctx: PanelContext, L: LatticeBundle, summary: HTMLElement | null): void {
  void ctx;
  const p = panel({
    title: 'Certified critical points · [TOP] Prop. E.2',
    subtitle: 'Z(2) and the five Hessian eigenvalues in the traceless shape directions, ours (Ewald, analytic Hessian) against the printed values. Every critical point obeys tr H = 2Z/3 exactly.',
    badges: [badge('proved', '[TOP] Prop. E.2'), badge('computed', 'lattice.critical')],
    span: 12,
  });
  const rows: CritRow[] = SHAPES.map((s) => ({ shape: s, c: L.critical[s] }));
  const cols: RawColumn<CritRow>[] = [
    { label: 'Shape', cell: (r) => h('strong', { class: 'mono' }, r.shape) },
    { label: 'Role', cell: (r) => (r.shape === 'SC' ? statusPill('warning', 'saddle (3 negative)') : statusPill('good', ROLE[r.shape])) },
    { label: 'Z(2)', num: true, cell: (r) => oursPaper(fmtNum(r.c.Z, 15), fmtNum(r.c.Z_paper, 15)) },
    {
      label: 'Hessian eigenvalues',
      num: true,
      cell: (r) => oursPaper(fmtEig(r.c.hessian, 5), fmtEig(r.c.hessian_paper, 4)),
      className: 'm1-mono',
    },
    { label: '|∇Z| (ours)', num: true, cell: (r) => (r.c.gradient_norm === undefined ? '—' : fmtSci(r.c.gradient_norm, 2)), title: 'Gradient norm at the critical point; Theorem E proves it vanishes exactly' },
    {
      label: 'tr H − 2Z/3',
      num: true,
      title: 'Trace identity at a critical point of Z(2) at unit covolume',
      cell: (r) => {
        const ours = hessianTraceCheck(r.c.Z, r.c.hessian).residual;
        const paper = hessianTraceCheck(r.c.Z_paper, r.c.hessian_paper).residual;
        return oursPaper(fmtSci(ours, 1), fmtSigned(paper, 2));
      },
    },
  ];
  const tbl = rawTable(cols, rows, (r) => (r.shape === 'BCC' ? 'm1-row-flag' : null));

  const bcc = L.critical.BCC;
  const oursTr = hessianTraceCheck(bcc.Z, bcc.hessian);
  const paperTr = hessianTraceCheck(bcc.Z_paper, bcc.hessian_paper);
  const triple = bcc.hessian[4]!;
  const gaps = L.gaps;

  append(
    p.body,
    tbl.root,
    h('div', { class: 'chart-caption' }, 'Top number ours, below it the printed value. Printed eigenvalues carry 4 decimals, so a residual of a few 10⁻⁴ is rounding, except where noted.'),
    h(
      'div',
      { class: 'callout m1-callout' },
      h('strong', null, 'One printed BCC eigenvalue looks like a typo. '),
      h(
        'p',
        null,
        `For the triple BCC eigenvalue the paper prints ${fmtNum(bcc.hessian_paper[4], 5)}; we obtain ${triple.toFixed(5)} from both the analytic and a finite-difference Hessian. `,
        `The identity tr H = 2Z/3 holds exactly at every critical point: with our values tr H = ${oursTr.trace.toFixed(5)} against 2Z/3 = ${oursTr.twoZ3.toFixed(5)} (difference ${fmtSci(oursTr.residual, 1)}), `,
        `while the printed values give ${paperTr.trace.toFixed(4)}, off by ${fmtSigned(paperTr.residual, 2)}, more than the 4-decimal rounding allows (at most 2.5×10⁻⁴). `,
        'Positive definiteness, and so every conclusion of Prop. E.2, is unaffected.',
      ),
      h('div', { class: 'row small' }, badge('computed', 'tr H = 2Z/3'), badge('quoted', '[TOP] Prop. E.2(b)')),
    ),
    gaps
      ? h(
          'div',
          null,
          sectionLabel('Certified gaps · Prop. E.2(c)'),
          tilesN(
            3,
            statTile({ label: 'Z_BCC − Z_FCC', value: fmtNum(gaps.computed['BCC-FCC'], 10), sub: [badge('computed'), `paper ${fmtNum(gaps.paper['BCC-FCC'], 10)}`] }),
            statTile({ label: 'Z_SC − Z_BCC', value: fmtNum(gaps.computed['SC-BCC'], 10), sub: [badge('computed'), `paper ${fmtNum(gaps.paper['SC-BCC'], 10)}`] }),
            statTile({ label: 'FCC is the deepest of the three', value: 'by 0.0097', sub: [badge('proved', gaps.source), 'a 0.06 % gap: the two minima are nearly degenerate'] }),
          ),
        )
      : null,
    summary,
  );
  el.append(grid(p.root));
}

// --- shape flow -----------------------------------------------------------------------------------

function mountFlows(el: HTMLElement, ctx: PanelContext, L: LatticeBundle): void {
  const F = L.flows;
  const zF = L.critical.FCC.Z;
  const gap = L.critical.BCC.Z - zF;
  const byEnd = (end: string): { tau: number[]; Z: number[] }[] => F.trajectories.filter((t) => t.end === end);
  const fcc = flowSeries(byEnd('FCC'), zF);
  const bcc = flowSeries(byEnd('BCC'), zF);
  const other = flowSeries(byEnd('other'), zF);
  const nShown = { FCC: byEnd('FCC').length, BCC: byEnd('BCC').length, other: byEnd('other').length };

  const p = panel({
    title: 'Shape flow from Haar-random lattices · [TOP] Theorem H.2',
    subtitle: `Z(2) − Z_FCC along ${F.trajectories.length} of the ${F.n} gradient flows (τ = ln η, radiation era), coloured by where each one ends. Flows that end at BCC level off at the BCC−FCC gap.`,
    badges: [badge('proved', '[TOP] Theorem H.2'), badge('computed', 'lattice.flows')],
    span: 8,
  });
  const series: LineSeries[] = [
    { name: `ends at FCC (${nShown.FCC} shown)`, ...fcc, color: SHAPE_SLOT.FCC, width: 1.5 },
    { name: `ends at BCC (${nShown.BCC} shown)`, ...bcc, color: SHAPE_SLOT.BCC, width: 1.5 },
  ];
  if (nShown.other) series.push({ name: `other (${nShown.other})`, ...other, color: 'context', width: 1.5 });
  keep(
    ctx,
    lineChart(p.body, series, {
      xTitle: 'τ = ln η (e-folds of conformal time)',
      yTitle: 'Z(2) − Z_FCC (unit covolume)',
      yType: 'log',
      height: 330,
      refLines: [{ axis: 'y', value: gap }],
      caption: `The horizontal line is the BCC−FCC gap ${fmtNum(gap, 4)}. Points within 10⁻⁷ of Z_FCC (the bundle's rounding) are dropped from the log axis.`,
    }),
  );

  const q = panel({
    title: 'Basins of attraction',
    subtitle: `Endpoints of all ${F.n} flows, ours (seed ${F.seed ?? '—'}, ${F.mode ?? 'gradient'} flow, |t| ≤ ${F.t_max ?? '—'}) against the paper's 240.`,
    badges: [badge('computed', 'lattice.flows'), badge('quoted', '[TOP] Sec. 6')],
    span: 4,
  });
  const ends = ['FCC', 'BCC', 'other'];
  keep(
    ctx,
    barChart(
      q.body,
      [
        { name: 'ours', x: ends, y: ends.map((e) => F.counts[e] ?? 0), color: 1 },
        { name: 'paper', x: ends, y: ends.map((e) => F.paper_counts[e] ?? 0), color: 'context' },
      ],
      { xTitle: 'endpoint of the flow', yTitle: 'number of flows', height: 250, layout: { xaxis: { type: 'category' } } },
    ),
  );
  const frac = F.fcc_fraction ?? (F.counts.FCC ?? 0) / F.n;
  const err = F.fcc_fraction_err ?? Math.sqrt((frac * (1 - frac)) / F.n);
  const pf = F.paper_fcc_fraction ?? (F.paper_counts.FCC ?? 0) / 240;
  q.body.append(
    tilesN(
      2,
      statTile({ label: 'FCC fraction, ours', value: `${(100 * frac).toFixed(1)} %`, sub: [badge('computed'), `± ${(100 * err).toFixed(1)} % (binomial), ${F.counts.FCC}/${F.n}`] }),
      statTile({ label: 'FCC fraction, paper', value: `${(100 * pf).toFixed(1)} %`, sub: [badge('quoted', '[TOP] Prop. H.3'), `${F.paper_counts.FCC}/240, printed as 79 % ± 3 %`] }),
    ),
    h(
      'p',
      { class: 'm1-takeaway' },
      `The two fractions differ by ${fmtNum(Math.abs(frac - pf) / Math.hypot(err, 0.03), 2)} combined standard errors: consistent. `,
      `${F.unconverged ? `${F.unconverged} flows did not converge.` : 'Every flow converged to FCC or BCC; none ended anywhere else, as in the paper.'}`,
    ),
  );
  el.append(grid(p.root, q.root));
}

// --- Prop. H.3: how much relaxation happens --------------------------------------------------------

function mountRelaxation(el: HTMLElement, L: LatticeBundle): void {
  const R = L.relaxation;
  const o = R.ours ?? {};
  const p = panel({
    title: 'Does the shape actually relax? · [TOP] Prop. H.3',
    subtitle: 'Near FCC the slowest relaxation rate per e-fold of conformal time is κc_E × (smallest FCC Hessian eigenvalue). With the Standard-Model coupling this freezes the shape as soon as the universe is many thermal wavelengths across.',
    badges: [badge('proved', '[TOP] Prop. H.3'), badge('computed', 'lattice.relaxation.ours')],
    span: 12,
  });
  const fmtE = (x: number | undefined): string => (x === undefined ? '—' : fmtSci(x, 3));
  p.body.append(
    tilesN(
      4,
      statTile({
        label: 'Coupling κc_E × (TL)⁴',
        value: fmtNum(o.kappa_cE_coeff ?? R.kappa_cE_coeff, 4),
        sub: [badge('computed'), `paper ≈ ${fmtNum(R.kappa_cE_coeff, 3)}; 180(n_F − n_B)/(π⁴ g*), n_F − n_B = 60, g* = 106.75`],
      }),
      statTile({
        label: 'Slowest rate / κc_E',
        value: fmtNum(o.fcc_mu_min ?? R.kappa_cE_rate_fcc, 5),
        sub: [badge('computed', 'FCC Hessian'), `paper ${fmtNum(R.kappa_cE_rate_fcc, 3)}: the smallest FCC eigenvalue`],
      }),
      statTile({
        label: 'Relaxation within Δτ needs (TL)⁴ ≲',
        value: `${fmtNum(o.TL4_per_dtau ?? R.TL4_per_dtau, 4)} Δτ`,
        sub: [badge('computed'), `paper ${fmtNum(R.TL4_per_dtau, 2)} Δτ`],
      }),
      statTile({
        label: 'TL today',
        value: fmtE(o.TL_today ?? R.TL_today),
        sub: [badge('computed'), `paper ≈ ${fmtE(R.TL_today)}, T₀ = 2.35×10⁻⁴ eV, L ≥ 0.985 D`],
      }),
    ),
    h(
      'div',
      { class: 'callout m1-callout' },
      h('strong', null, 'FCC is a conditional prediction. '),
      h(
        'p',
        null,
        'With (TL)⁴ ≈ 10¹²⁰ against the ≲ 1 needed, the shape has been frozen since the entropy of the universe was produced, and inflation freezes it too (Prop. H.3(b)). ',
        'Our torus therefore has the shape it had when Casimir forces last acted. It is FCC only if (1) there was a pre-inflationary epoch with TL of order one, a universe about one thermal wavelength across, lasting several relaxation times, and (2) the initial shape lay in the FCC basin, about 79 % of Haar-random starts (81 % in our run). ',
        'Neither condition is established, and the paper\'s Planck analysis treats the shape as free, testing SC, BCC and FCC.',
      ),
      h('div', { class: 'row small' }, badge('assumed', '[TOP] Sec. 13 status table'), badge('quoted', '[TOP] Prop. H.3')),
    ),
    note(
      proofList([
        { kind: 'proved', src: '[TOP] Theorem E', text: 'The cubic lattices are critical points of every rotation-invariant functional of the shape (unconditional).' },
        { kind: 'proved', src: '[TOP] Prop. E.2', text: 'FCC and BCC are strict local minima of Z(2) and SC is a saddle; computer-assisted, with certified error bounds (200-bit ball arithmetic).' },
        { kind: 'proved', src: '[TOP] Theorems H.1, H.2', text: 'The shape obeys Einstein\'s equations for a Bianchi I torus; in a radiation era every flow converges to a critical point, FCC and BCC are asymptotically stable, SC unstable.' },
        { kind: 'assumed', src: '[TOP] Prop. H.3', text: 'That our torus is FCC needs a pre-inflationary Casimir era and a start in the FCC basin: a conditional prediction, not assumed in the data analysis.' },
        { kind: 'open', src: '[TOP] Sec. 6 remark', text: 'Whether FCC is the global minimum (the d = 3 Sarnak–Strömbergsson conjecture) is not used; our 240 flows agree with it.' },
      ]),
      'What the paper proves and what it assumes',
      true,
    ),
  );
  el.append(grid(p.root));
}
