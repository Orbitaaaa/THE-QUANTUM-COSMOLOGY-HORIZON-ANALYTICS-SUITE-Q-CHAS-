/**
 * Act II · Experiment — "Hunt for the seams" (docs/SHOW.md).
 *
 * Left: the real WMAP sky with the exact matched circles of the chosen torus. Right: the
 * wave-space crystal, the dual-lattice wavevectors kappa with |kappa| chi_* <= 14 that a torus
 * allows. Judges change the size d/D, the shape, and switch between the standard search
 * (smooth spectrum, unknown orientation) and a resonance (a narrow primordial line on one
 * dual-lattice shell, known orientation). The meter shows the information K in nats from
 * cmb.json (smooth spectrum) and resonance.json (cubic torus, |m|^2 = 16 shell, f = 0.5).
 */

import * as THREE from 'three';
import { LineMaterial } from 'three/examples/jsm/lines/LineMaterial.js';
import { LineSegments2 } from 'three/examples/jsm/lines/LineSegments2.js';
import { LineSegmentsGeometry } from 'three/examples/jsm/lines/LineSegmentsGeometry.js';
import type { CmbBundle, ResonanceBundle } from '../../lib/data';
import { circlePairs } from '../../lib/physics/circles';
import { CIRCLE_SEARCH_BOUND } from '../../lib/physics/constants';
import { type CubicShape, cubicBasis, dualBasis, latticeVectors } from '../../lib/physics/lattice';
import { type Mat3, identity3, matVec } from '../../lib/physics/vec';
import { scaledBasis } from '../../lib/physics/circles';
import { h } from '../../lib/ui';
import { Act, type ActContext, type ActModule, DragOrbit } from '../act';
import { EXPERIMENT_TAKEAWAY } from '../script';
import { POINTS_FRAG, cmbValueTexture, easeInOut, galToThree, lerp, pointScale, ramp, rng, skyDome } from '../gfx';
import { CMB_COLOR_GLSL, GALACTIC_GLSL } from '../shaders/common';
import { fmt, syncFill } from '../ui';

const R_GLOBE = 4.2;
const R_CRYSTAL = 4.4; // world radius of |kappa| chi_* = 14
const KCHI_MAX = 14;
const MAX_PTS = 2600;
const MAX_RINGS = 96;
const MAX_SHELL = 12;
const SHAPE_LABEL: Record<CubicShape, string> = { SC: 'cubic', BCC: 'BCC', FCC: 'FCC' };

interface State {
  shape: CubicShape;
  d: number;
  res: number; // 0 standard, 1 resonance (blended)
}

const GLOBE_FRAG = /* glsl */ `
${GALACTIC_GLSL}
${CMB_COLOR_GLSL}
uniform sampler2D uMap;
uniform float uRes;
uniform vec3 uK[${MAX_SHELL}];
uniform float uPhi[${MAX_SHELL}];
uniform int uN;
varying vec3 vLocal;
varying vec3 vN;
varying vec3 vView;
void main() {
  vec3 n = normalize(vLocal);
  float v = (texture2D(uMap, galUv(threeToGal(n))).r - 0.5) * 2.0;
  vec3 col = cmbColor(v * 1.5) * 0.85;
  float p = 0.0;
  for (int i = 0; i < ${MAX_SHELL}; i++) { if (i >= uN) break; p += cos(dot(uK[i], n) + uPhi[i]); }
  p /= sqrt(max(1.0, float(uN)) * 0.5);
  vec3 pc = cmbColor(p * 0.5) * 0.85;
  float crest = smoothstep(0.9, 1.0, abs(sin(p * 1.9)));
  col = mix(col, pc, 0.75 * uRes) + vec3(1.0, 0.7, 0.3) * crest * 0.16 * uRes;
  float f = pow(1.0 - abs(dot(normalize(vN), normalize(vView))), 3.0);
  col += mix(vec3(0.25, 0.75, 1.0), vec3(1.0, 0.65, 0.25), uRes) * f * 0.45;
  gl_FragColor = vec4(col, 1.0);
}`;

const GLOBE_VERT = /* glsl */ `
varying vec3 vLocal;
varying vec3 vN;
varying vec3 vView;
void main() {
  vLocal = position;
  vec4 wp = modelMatrix * vec4(position, 1.0);
  vN = normalize(mat3(modelMatrix) * normal);
  vView = cameraPosition - wp.xyz;
  gl_Position = projectionMatrix * viewMatrix * wp;
}`;

const PTS_VERT = /* glsl */ `
attribute float aSize;
attribute vec3 aColor;
attribute float aAlpha;
attribute float aShell;
uniform float uScale;
uniform float uRes;
uniform float uTime;
varying vec3 vColor;
varying float vAlpha;
void main() {
  float flare = aShell * uRes;
  vColor = mix(aColor, vec3(1.9, 1.25, 0.45), flare);
  vAlpha = aAlpha * mix(1.0, 0.35, uRes * (1.0 - aShell));
  vec4 mv = modelViewMatrix * vec4(position, 1.0);
  float pulse = 1.0 + flare * (1.4 + 0.35 * sin(uTime * 5.0));
  gl_PointSize = min(aSize * pulse * uScale / max(0.05, -mv.z), 48.0);
  gl_Position = projectionMatrix * mv;
}`;

function lookupR(topoR: number[][] | undefined): Mat3 {
  if (topoR && topoR.length === 3 && topoR.every((r) => r.length === 3 && r.every(Number.isFinite))) return topoR.map((r) => [r[0]!, r[1]!, r[2]!]) as Mat3;
  return identity3();
}

/** log-linear interpolation of y(x) on a sorted grid; null outside */
function interpLog(xs: number[], ys: number[], x: number): number | null {
  if (!xs.length || x < xs[0]! - 1e-9 || x > xs[xs.length - 1]! + 1e-9) return null;
  for (let i = 0; i < xs.length - 1; i++) {
    const a = xs[i]!;
    const b = xs[i + 1]!;
    if (x >= a - 1e-9 && x <= b + 1e-9) {
      const f = (x - a) / (b - a);
      return Math.exp(Math.log(ys[i]!) * (1 - f) + Math.log(ys[i + 1]!) * f);
    }
  }
  return ys[ys.length - 1]!;
}

class Experiment extends Act {
  readonly id = 'experiment' as const;
  readonly duration = 13;
  private readonly R: Mat3;
  private readonly globe = new THREE.Group();
  private readonly crystal = new THREE.Group();
  private readonly globeMat: THREE.ShaderMaterial;
  private readonly rings: THREE.InstancedMesh;
  private readonly ringMat: THREE.ShaderMaterial;
  private readonly pts: THREE.Points;
  private readonly ptsMat: THREE.ShaderMaterial;
  private readonly beamGeo = new LineSegmentsGeometry();
  private readonly beamMat: LineMaterial;
  private readonly shellMat: THREE.ShaderMaterial;
  private readonly shellSphere: THREE.Mesh;
  private readonly orbit: DragOrbit;
  private user: State | null = null;
  private resBlend = 0;
  private key = '';
  private nShell = 0;
  private shellKChi = 0;
  private readonly kl: { smooth: (s: CubicShape, d: number) => number | null; line: (d: number) => number | null };
  // HUD
  private readonly labelSky: HTMLElement;
  private readonly labelK: HTMLElement;
  private readonly panel: HTMLElement;
  private readonly modeBtns: HTMLButtonElement[];
  private readonly shapeBtns: Map<CubicShape, HTMLButtonElement>;
  private readonly dSlider: HTMLInputElement;
  private readonly dVal: HTMLElement;
  private readonly dNote: HTMLElement;
  private readonly barSmooth: HTMLElement;
  private readonly valSmooth: HTMLElement;
  private readonly barLine: HTMLElement;
  private readonly valLine: HTMLElement;
  private readonly lossBand: HTMLElement;
  private readonly rowLoss: HTMLElement;
  private readonly gainEl: HTMLElement;
  private readonly gainNote: HTMLElement;
  private readonly shellNote: HTMLElement;
  private readonly takeaway: HTMLElement;

  constructor(ctx: ActContext, tex: THREE.Texture, cmb: CmbBundle | null, res: ResonanceBundle | null, topoR: number[][] | undefined) {
    super(ctx);
    this.R = lookupR(topoR);
    this.camera.fov = 36;
    this.camera.position.set(0, 1.2, 34);
    this.camera.lookAt(0, 0.6, 0);
    this.bloom = { strength: 0.85, radius: 0.5, threshold: 0.5 };
    this.scene.add(skyDome(1500, 0.55));
    this.scene.add(this.globe, this.crystal);

    // --- information tables ---------------------------------------------------------------
    const klx = cmb?.kl?.d_over_D ?? [];
    const lineRows = (res?.kl_table ?? []).filter((r) => r.shape === 'cubic' && r.shell === 16 && Math.abs(r.f - 0.5) < 1e-9).sort((a, b) => a.d_over_D - b.d_over_D);
    this.kl = {
      smooth: (s, d) => {
        const key = s === 'SC' ? 'cubic' : s;
        const ys = (cmb?.kl?.[key] as number[] | undefined) ?? [];
        return ys.length === klx.length ? interpLog(klx, ys, d) : null;
      },
      line: (d) =>
        interpLog(
          lineRows.map((r) => r.d_over_D),
          lineRows.map((r) => r.kl),
          d,
        ),
    };

    // --- globe ----------------------------------------------------------------------------
    const kU: THREE.Vector3[] = Array.from({ length: MAX_SHELL }, () => new THREE.Vector3());
    const r0 = rng(11);
    this.globeMat = new THREE.ShaderMaterial({
      uniforms: { uMap: { value: tex }, uRes: { value: 0 }, uK: { value: kU }, uPhi: { value: Array.from({ length: MAX_SHELL }, () => 2 * Math.PI * r0()) }, uN: { value: 0 } },
      vertexShader: GLOBE_VERT,
      fragmentShader: GLOBE_FRAG,
    });
    this.globe.add(new THREE.Mesh(new THREE.SphereGeometry(R_GLOBE, 128, 64), this.globeMat));
    const halo = new THREE.Mesh(
      new THREE.SphereGeometry(R_GLOBE * 1.12, 64, 32),
      new THREE.ShaderMaterial({
        uniforms: { uColor: { value: new THREE.Color(0.2, 0.7, 1.0) } },
        vertexShader: GLOBE_VERT,
        fragmentShader: /* glsl */ `
          uniform vec3 uColor; varying vec3 vN; varying vec3 vView;
          void main() { float f = abs(dot(normalize(vN), normalize(vView))); float g = pow(smoothstep(0.0, 0.62, f), 2.0) * (1.0 - smoothstep(0.62, 1.0, f));
            gl_FragColor = vec4(uColor * g * 0.55, 1.0); }`,
        side: THREE.BackSide,
        transparent: true,
        blending: THREE.AdditiveBlending,
        depthWrite: false,
      }),
    );
    this.globe.add(halo);
    this.ringMat = new THREE.ShaderMaterial({
      uniforms: { uTime: { value: 0 }, uOpacity: { value: 1 } },
      vertexShader: /* glsl */ `
        varying float vU;
        void main() { vU = uv.x; gl_Position = projectionMatrix * viewMatrix * modelMatrix * instanceMatrix * vec4(position, 1.0); }`,
      fragmentShader: /* glsl */ `
        uniform float uTime; uniform float uOpacity; varying float vU;
        void main() { float flow = 0.75 + 0.25 * sin(vU * 6.2831853 * 5.0 - uTime * 2.5);
          gl_FragColor = vec4(vec3(0.35, 0.95, 1.0) * 1.35 * flow * uOpacity, 1.0); }`,
      transparent: true,
      blending: THREE.AdditiveBlending,
      depthWrite: false,
    });
    this.rings = new THREE.InstancedMesh(new THREE.TorusGeometry(1, 0.012, 6, 160), this.ringMat, MAX_RINGS);
    this.rings.count = 0;
    this.rings.frustumCulled = false;
    this.globe.add(this.rings);

    // --- crystal ----------------------------------------------------------------------------
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.BufferAttribute(new Float32Array(MAX_PTS * 3), 3));
    g.setAttribute('aSize', new THREE.BufferAttribute(new Float32Array(MAX_PTS), 1));
    g.setAttribute('aColor', new THREE.BufferAttribute(new Float32Array(MAX_PTS * 3), 3));
    g.setAttribute('aAlpha', new THREE.BufferAttribute(new Float32Array(MAX_PTS), 1));
    g.setAttribute('aShell', new THREE.BufferAttribute(new Float32Array(MAX_PTS), 1));
    this.ptsMat = new THREE.ShaderMaterial({
      uniforms: { uScale: { value: 800 }, uRes: { value: 0 }, uTime: { value: 0 } },
      vertexShader: PTS_VERT,
      fragmentShader: POINTS_FRAG,
      transparent: true,
      blending: THREE.AdditiveBlending,
      depthWrite: false,
    });
    this.pts = new THREE.Points(g, this.ptsMat);
    this.pts.frustumCulled = false;
    this.crystal.add(this.pts);
    this.beamMat = new LineMaterial({ color: 0xffc35a, linewidth: 2.2, transparent: true, opacity: 0, blending: THREE.AdditiveBlending, depthWrite: false });
    const beams = new LineSegments2(this.beamGeo, this.beamMat);
    beams.frustumCulled = false;
    this.crystal.add(beams);
    this.shellMat = new THREE.ShaderMaterial({
      uniforms: { uOpacity: { value: 0 } },
      vertexShader: GLOBE_VERT,
      fragmentShader: /* glsl */ `
        uniform float uOpacity; varying vec3 vN; varying vec3 vView;
        void main() { float f = pow(1.0 - abs(dot(normalize(vN), normalize(vView))), 3.0);
          gl_FragColor = vec4(vec3(1.0, 0.62, 0.22) * (0.015 + f * f) * uOpacity, 1.0); }`,
      transparent: true,
      blending: THREE.AdditiveBlending,
      depthWrite: false,
      side: THREE.DoubleSide,
    });
    this.shellSphere = new THREE.Mesh(new THREE.SphereGeometry(1, 64, 32), this.shellMat);
    this.crystal.add(this.shellSphere);
    const boundary = new THREE.Mesh(
      new THREE.SphereGeometry(R_CRYSTAL, 64, 32),
      new THREE.ShaderMaterial({
        vertexShader: GLOBE_VERT,
        fragmentShader: /* glsl */ `
          varying vec3 vN; varying vec3 vView;
          void main() { float f = pow(1.0 - abs(dot(normalize(vN), normalize(vView))), 4.0);
            gl_FragColor = vec4(vec3(0.3, 0.55, 1.0) * f * 0.22, 1.0); }`,
        transparent: true,
        blending: THREE.AdditiveBlending,
        depthWrite: false,
      }),
    );
    this.crystal.add(boundary);

    // --- HUD --------------------------------------------------------------------------------
    this.labelSky = h('div', { class: 'obj-label' }, h('b', null, 'Real sky'), ' · WMAP 9-year, with the exact circles of this torus');
    this.labelK = h('div', { class: 'obj-label' }, h('b', null, 'Wave-space crystal'), ' · every wave a torus allows, |κ|χ* ≤ 14');
    this.modeBtns = [
      h('button', { type: 'button', class: 'standard', 'aria-pressed': 'true' }, 'Standard search') as HTMLButtonElement,
      h('button', { type: 'button', class: 'resonance', 'aria-pressed': 'false' }, 'Resonance shell') as HTMLButtonElement,
    ];
    this.shapeBtns = new Map((['SC', 'BCC', 'FCC'] as CubicShape[]).map((s) => [s, h('button', { type: 'button', 'aria-pressed': 'false' }, SHAPE_LABEL[s]) as HTMLButtonElement]));
    this.dSlider = h('input', { type: 'range', min: '0.6', max: '1.5', step: '0.005', value: '1.3', 'aria-label': 'Universe size d over D' }) as HTMLInputElement;
    this.dVal = h('span', { class: 'mono dval' }, '1.30');
    this.dNote = h('span', { class: 'dnote2' });
    const bar = (cls: string): [HTMLElement, HTMLElement] => {
      const b = h('div', { class: `kbar ${cls}` });
      return [b, h('div', { class: 'ktrack' }, b)];
    };
    const [bs, ts] = bar('smooth');
    const [bl, tl] = bar('line');
    this.barSmooth = bs;
    this.barLine = bl;
    this.valSmooth = h('span', { class: 'kval mono' });
    this.valLine = h('span', { class: 'kval mono' });
    this.lossBand = h('div', { class: 'kloss' });
    this.rowLoss = h('div', { class: 'krow loss' }, h('div', { class: 'klab' }, 'Resonance · ', h('b', null, 'unknown'), ' orientation'), h('div', { class: 'ktrack' }, this.lossBand), h('span', { class: 'kval mono' }, '—'));
    this.gainEl = h('div', { class: 'gain' }, '×1');
    this.gainNote = h('div', { class: 'gain-note' }, 'more information');
    this.shellNote = h('div', { class: 'shell-note' });
    const tick = (k: number, label: string): HTMLElement => h('span', { class: 'ktick', style: { left: `${kPos(k)}%` } }, label);
    this.panel = h(
      'div',
      { class: 'show-panel exp-panel show-ui' },
      h('div', { class: 'seg res', role: 'group', 'aria-label': 'Search method' }, ...this.modeBtns),
      h(
        'div',
        { class: 'exp-row' },
        h('div', { class: 'seg', role: 'group', 'aria-label': 'Torus shape' }, ...this.shapeBtns.values()),
        h('label', { class: 'dock-field grow' }, h('span', { class: 'dock-k' }, 'Size d/D ', this.dVal, ' ', this.dNote), h('span', { class: 'dock-scale' }, this.dSlider, h('span', { class: 'dock-excl', style: { width: `${((0.985 - 0.6) / 0.9) * 100}%` } }))),
      ),
      h('div', { class: 'panel-k' }, h('span', null, 'Information K (nats, log scale)'), h('span', { class: 'kone-key' }, 'K = 1 detectable')),
      h(
        'div',
        { class: 'kmeter' },
        h('div', { class: 'krow' }, h('div', { class: 'klab' }, 'Standard search · smooth spectrum'), ts, this.valSmooth),
        h('div', { class: 'krow' }, h('div', { class: 'klab' }, 'Resonance · ', h('b', null, 'known'), ' orientation'), tl, this.valLine),
        this.rowLoss,
        h('div', { class: 'kaxis' }, tick(0.01, '0.01'), tick(0.1, '0.1'), h('span', { class: 'ktick one', style: { left: `${kPos(1)}%` } }, '1'), tick(10, '10'), tick(100, '100'), h('span', { class: 'kone', style: { left: `${kPos(1)}%` } })),
      ),
      h('div', { class: 'gain-wrap' }, this.gainEl, h('div', null, this.gainNote, this.shellNote)),
      h('div', { class: 'exp-foot' }, 'Real sky: no line found — p = 0.78 (Planck) · 0.79 (WMAP, ours).'),
    );
    this.takeaway = h('div', { class: 'exp-takeaway' }, EXPERIMENT_TAKEAWAY);
    ctx.hud.append(this.labelSky, this.labelK, this.panel, this.takeaway);

    const touch = (): State => {
      if (!this.user) this.user = { ...this.scripted(this.time) };
      this.pause();
      ctx.userInput();
      return this.user;
    };
    this.modeBtns.forEach((b, i) => this.listen(b, 'click', () => (touch().res = i)));
    for (const [s, b] of this.shapeBtns) this.listen(b, 'click', () => (touch().shape = s));
    this.listen(this.dSlider, 'input', () => (touch().d = Number(this.dSlider.value)));

    this.orbit = new DragOrbit(ctx.stage.canvas.parentElement ?? ctx.stage.canvas, () => ctx.userInput(), 1.0);
    this.own(() => this.orbit.dispose());
    this.update(0);
  }

  override play(): void {
    this.user = null;
    super.play();
  }

  override seek(t: number): void {
    this.user = null;
    super.seek(t);
  }

  private scripted(t: number): State {
    return { shape: 'SC', d: lerp(1.3, 1.05, easeInOut((t - 7.2) / 1.8)), res: ramp(t, 3.4, 4.6) };
  }

  protected update(dt: number): void {
    const t = this.time;
    const target = this.user ?? this.scripted(t);
    if (this.user) this.resBlend += (target.res - this.resBlend) * Math.min(1, dt * 5);
    else this.resBlend = target.res;
    const st: State = { ...target, res: this.resBlend };
    this.orbit.step(dt);

    // layout: globe left, crystal right (stacked on narrow screens)
    const aspect = this.width / Math.max(1, this.height);
    const halfH = 34 * Math.tan((this.camera.fov * Math.PI) / 360);
    const halfW = halfH * aspect;
    if (aspect > 1.15) {
      this.globe.position.set(-0.4 * halfW, 1.6, 0);
      this.crystal.position.set(0.12 * halfW, 1.6, 0);
    } else {
      this.globe.position.set(0, 0.42 * halfH, 0);
      this.crystal.position.set(0, -0.05 * halfH, 0);
    }
    const spin = 0.35 + 0.08 * (this.wall + 0) + this.orbit.yaw;
    const q = new THREE.Quaternion().setFromEuler(new THREE.Euler(0.32 + this.orbit.pitch, spin, 0, 'YXZ'));
    this.globe.quaternion.copy(q);
    this.crystal.quaternion.copy(q);

    this.rebuild(st.shape, st.d);
    this.globeMat.uniforms.uRes!.value = st.res;
    this.ptsMat.uniforms.uRes!.value = st.res;
    this.ptsMat.uniforms.uTime!.value = this.wall;
    this.ptsMat.uniforms.uScale!.value = pointScale(this.camera, this.height, this.pixelRatio);
    this.beamMat.opacity = 0.85 * st.res;
    this.beamMat.resolution.set(this.width, this.height);
    this.shellMat.uniforms.uOpacity!.value = 0.6 * st.res;
    this.ringMat.uniforms.uTime!.value = this.wall;

    this.paintHud(st, t);
  }

  /** rebuild circles, crystal points and the resonance shell when shape or size change */
  private rebuild(shape: CubicShape, d: number): void {
    const key = `${shape}:${d.toFixed(4)}`;
    if (key === this.key) return;
    this.key = key;
    const basis = cubicBasis(shape);
    // circles on the globe
    const pairs = d < 1 ? circlePairs(basis, d, this.R) : [];
    const m = new THREE.Matrix4();
    const qq = new THREE.Quaternion();
    const z = new THREE.Vector3(0, 0, 1);
    let k = 0;
    for (const p of pairs) {
      const ca = p.length_over_D;
      const sa = Math.sqrt(Math.max(0, 1 - ca * ca));
      for (const c of [p.centre_minus, p.centre_plus]) {
        if (k >= MAX_RINGS) break;
        const n = galToThree(c);
        qq.setFromUnitVectors(z, n);
        m.compose(n.clone().multiplyScalar(R_GLOBE * ca * 1.002), qq, new THREE.Vector3().setScalar(R_GLOBE * sa * 1.002));
        this.rings.setMatrixAt(k++, m);
      }
    }
    this.rings.count = k;
    this.rings.instanceMatrix.needsUpdate = true;

    // dual lattice: kappa = 2 pi inv(B)^T n, shown out to |kappa| chi_* = 14
    const B = scaledBasis(basis, d);
    const K = dualBasis(B);
    const vecs = latticeVectors(K, (2 * KCHI_MAX) * (1 + 1e-9)).slice(0, MAX_PTS);
    const kmin = vecs[0]?.length ?? 1;
    // resonance shell: |m|^2 = 16 (6 vectors along the cube axes) for cubic, the first shell otherwise
    const shellLen = shape === 'SC' ? 4 * kmin : kmin;
    const g = this.pts.geometry;
    const pos = g.getAttribute('position') as THREE.BufferAttribute;
    const size = g.getAttribute('aSize') as THREE.BufferAttribute;
    const col = g.getAttribute('aColor') as THREE.BufferAttribute;
    const alpha = g.getAttribute('aAlpha') as THREE.BufferAttribute;
    const shell = g.getAttribute('aShell') as THREE.BufferAttribute;
    const kVecs: THREE.Vector3[] = [];
    const beam: number[] = [];
    const scale = R_CRYSTAL / KCHI_MAX;
    vecs.forEach((v, i) => {
      const kchi = v.length * 0.5;
      const w = galToThree(matVec(this.R, v.g), 0.5 * scale);
      pos.setXYZ(i, w.x, w.y, w.z);
      const x = kchi / KCHI_MAX;
      size.setX(i, 0.05 + 0.32 * Math.pow(Math.max(0.12, kchi / 2.5), -0.9));
      const c = new THREE.Color().setHSL(0.53 + 0.17 * x, 0.85, 0.62 - 0.12 * x);
      col.setXYZ(i, c.r * 1.25, c.g * 1.25, c.b * 1.25);
      alpha.setX(i, 0.95 - 0.55 * x);
      const on = Math.abs(v.length - shellLen) < 1e-6 * shellLen ? 1 : 0;
      shell.setX(i, on);
      if (on) {
        const kv = galToThree(matVec(this.R, v.g), 0.5); // kappa chi_* in sky (three) coordinates
        if (kVecs.length < MAX_SHELL) kVecs.push(kv);
        beam.push(-w.x * 1.35, -w.y * 1.35, -w.z * 1.35, w.x * 1.35, w.y * 1.35, w.z * 1.35);
      }
    });
    g.setDrawRange(0, vecs.length);
    for (const a of [pos, size, col, alpha, shell]) a.needsUpdate = true;
    this.nShell = kVecs.length;
    this.shellKChi = shellLen * 0.5;
    const uK = this.globeMat.uniforms.uK!.value as THREE.Vector3[];
    kVecs.forEach((kv, i) => uK[i]!.copy(kv));
    this.globeMat.uniforms.uN!.value = kVecs.length;
    if (beam.length) this.beamGeo.setPositions(beam);
    this.shellSphere.scale.setScalar(Math.max(0.01, this.shellKChi * scale));
  }

  private paintHud(st: State, t: number): void {
    // labels under the objects
    const place = (el: HTMLElement, obj: THREE.Object3D, radius: number): void => {
      const p = obj.position.clone().add(new THREE.Vector3(0, -radius * 1.18, 0)).project(this.camera);
      el.style.transform = `translate(-50%, 0) translate(${((p.x + 1) / 2) * this.width}px, ${((1 - p.y) / 2) * this.height}px)`;
    };
    place(this.labelSky, this.globe, R_GLOBE);
    place(this.labelK, this.crystal, R_CRYSTAL);

    const resOn = st.res > 0.5;
    this.modeBtns[0]!.setAttribute('aria-pressed', String(!resOn));
    this.modeBtns[1]!.setAttribute('aria-pressed', String(resOn));
    for (const [s, b] of this.shapeBtns) b.setAttribute('aria-pressed', String(s === st.shape));
    if (document.activeElement !== this.dSlider) this.dSlider.value = st.d.toFixed(3);
    syncFill(this.dSlider);
    this.dVal.textContent = st.d.toFixed(2);
    const excluded = st.d < CIRCLE_SEARCH_BOUND.value;
    this.dNote.textContent = excluded ? 'excluded by data' : '';
    this.dNote.classList.toggle('bad', excluded);

    const ks = this.kl.smooth(st.shape, st.d);
    const kl = st.shape === 'SC' ? this.kl.line(st.d) : null;
    const big = st.d < 1;
    setBar(this.barSmooth, big ? 100 : ks);
    this.valSmooth.textContent = big ? '≫ 19' : ks === null ? '< 0.03' : fmt(ks);
    const lineShown = kl !== null ? lerp(ks ?? kl, kl, st.res) : null;
    setBar(this.barLine, lineShown);
    this.valLine.textContent = kl === null ? '—' : fmt(lineShown ?? kl);
    const gain = kl !== null && ks ? kl / ks : null;
    const g = gain !== null ? lerp(1, gain, st.res) : null;
    this.gainEl.textContent = g !== null ? `×${g < 10 ? g.toFixed(1) : g.toFixed(0)}` : '—';
    this.gainNote.textContent = g !== null ? 'more information with a resonance' : big ? 'circles already visible (size excluded)' : 'resonance K computed for the cubic torus, 1.05–1.3 D';
    this.shellNote.textContent = this.nShell ? `shell: ${this.nShell} wavevectors, |κ|χ* = ${this.shellKChi.toFixed(2)}${st.shape === 'SC' ? ' (|m|² = 16, f = 0.5)' : ''}` : '';
    // the loss band spans the gain that an orientation search consumes
    if (kl !== null && ks) {
      const a = kPos(ks);
      const b = kPos(lerp(ks, kl, st.res));
      this.lossBand.style.left = `${a}%`;
      this.lossBand.style.width = `${Math.max(0, b - a)}%`;
    } else this.lossBand.style.width = '0';
    this.rowLoss.classList.toggle('hot', !this.user && t >= 10);
    const tk = !this.user && this.isPlaying === false && t >= this.duration - 0.01 ? 1 : !this.user ? ramp(t, 11.2, 11.8) : 0;
    this.takeaway.style.opacity = String(tk);
    this.panel.classList.toggle('dim-line', st.res < 0.02);
  }

  override dispose(): void {
    this.beamGeo.dispose();
    super.dispose();
  }
}

/** position on the meter's log axis, 0.01 .. 100 nats -> 0 .. 100 % */
function kPos(k: number): number {
  return Math.max(0, Math.min(100, ((Math.log10(Math.max(k, 1e-3)) + 2) / 4) * 100));
}

function setBar(el: HTMLElement, k: number | null): void {
  el.style.width = k === null ? '0' : `${kPos(k)}%`;
}

export const experiment: ActModule = {
  async mount(ctx) {
    const [sky, cmb, res, topo] = await Promise.all([ctx.data.sky(), ctx.data.cmb(), ctx.data.resonance(), ctx.data.topology()]);
    const tex = sky ? cmbValueTexture(sky) : new THREE.DataTexture(new Uint8Array(16).fill(128), 2, 2);
    tex.needsUpdate = true;
    return new Experiment(ctx, tex, cmb, res, topo?.circles?.fcc_preregistered?.R as number[][] | undefined);
  },
};
