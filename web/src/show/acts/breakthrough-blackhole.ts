/**
 * Act III·A · Breakthrough — "A black hole gives its information back" (docs/SHOW.md).
 *
 * An evaporating Schwarzschild hole (M0 = 1, r_s = 2m, shadow b_c = 3 sqrt(3) m with
 * m(t) = (1 - t)^(1/3) of [BH] Theorem C) in front of a procedural deep sky that is lensed per
 * pixel (strong-deflection log term + weak-field 4m/b), with a photon ring. The horizon image is
 * tiled with Planck cells (ILLUSTRATIVE size, fixed area each, so their number tracks A/4G).
 * Hawking pairs are emitted at the rate of Hawking's curve; each escaping quantum is tied to an
 * infalling partner by a gold thread. After the Page time the number of threads crossing the
 * horizon follows min{S_H, A/4G} exactly: the earliest partners migrate into the radiation and
 * their threads turn cyan, until at t = 1 every thread lives in the radiation. The diary uses the
 * exact 20-qubit run (horizon.json -> evaporation): thrown in after 5 emissions it is recovered at
 * step 13; after 14, at step 17. The HUD draws Hawking's curve and the Page curve live from
 * lib/physics/page.ts.
 */

import * as THREE from 'three';
import type { EvaporationRun, HorizonBundle } from '../../lib/data';
import { bhEntropy, hawkingEntropy, massFraction, pageCurve, pageTimeFraction } from '../../lib/physics/page';
import { BETA_RADIATION } from '../../lib/physics/constants';
import { h } from '../../lib/ui';
import { Act, type ActContext, type ActModule, DragOrbit } from '../act';
import { POINTS_FRAG, clamp, lerp, pointScale, ramp, rng } from '../gfx';
import { NOISE_GLSL, SKY_GLSL } from '../shaders/common';
import { button } from '../ui';

const BC0 = 3 * Math.sqrt(3); // shadow radius for m = 1
const LIFE = 17.5; // act seconds per lifetime
const N_PAIRS = 280;
const N_DUST = 9000;
const N_STEPS = 21; // emissions in the 20-qubit run
const T_PAGE = pageTimeFraction();
const BETA = BETA_RADIATION.value;

const PAIR_GLSL = /* glsl */ `
float mfrac(float t) { return pow(max(1.0 - t, 0.0), 1.0 / 3.0); }
float bcAt(float t) { return ${BC0.toFixed(6)} * mfrac(t); }
vec3 outPos(vec3 dir, float te, float sp, float t) {
  float re = bcAt(te) * 1.04;
  return dir * (re + sp * max(t - te, 0.0));
}
vec3 inPos(vec3 dir, float te, float t) {
  float k = 1.0 - exp(-max(t - te, 0.0) * 7.0);
  return dir * bcAt(min(t, 0.999)) * (0.9 - 0.5 * k);
}`;

// --- lensed sky -------------------------------------------------------------------------------
const LENS_VERT = /* glsl */ `
varying vec2 vNdc;
void main() { vNdc = position.xy; gl_Position = vec4(position.xy, 0.0, 1.0); }`;

const LENS_FRAG = /* glsl */ `
${NOISE_GLSL}
${SKY_GLSL}
uniform mat4 uInvProj;
uniform mat4 uCamWorld;
uniform vec3 uCamPos;
uniform float uM;
uniform float uFlash;
varying vec2 vNdc;
void main() {
  vec4 v = uInvProj * vec4(vNdc, 1.0, 1.0);
  vec3 D = normalize((uCamWorld * vec4(v.xyz / v.w, 0.0)).xyz);
  vec3 O = uCamPos;
  float tc = -dot(O, D);
  vec3 P = O + D * max(tc, 0.0);
  float b = length(P);
  float bc = ${BC0.toFixed(6)} * uM;
  vec3 col;
  vec3 band = normalize(vec3(0.3, 1.0, 0.25));
  if (uM < 0.002) {
    col = deepSky(D, band);
  } else if (tc > 0.0 && b < bc) {
    col = vec3(0.0);
  } else {
    float x = max(b / bc - 1.0, 1e-4);
    float weak = 4.0 * uM / b * (1.0 + 0.75 * 3.14159 * uM / b);
    float strong = -log(x) - 0.4002;
    float alpha = tc > 0.0 ? max(weak, strong) : 4.0 * uM / max(length(O), b);
    vec3 H = normalize(-O);
    vec3 E = H - dot(H, D) * D;
    E = length(E) > 1e-5 ? normalize(E) : vec3(0.0, 1.0, 0.0);
    vec3 Dp = cos(alpha) * D + sin(alpha) * E;
    col = deepSky(Dp, band);
    // photon ring: the light that circles the hole piles up just outside b_c
    float w = 0.016 * bc;
    float ring = exp(-pow((b - bc * 1.004) / w, 2.0)) * 1.7 + exp(-pow((b - bc * 1.03) / (3.0 * w), 2.0)) * 0.28;
    col += vec3(1.0, 0.82, 0.58) * ring * (tc > 0.0 ? 1.0 : 0.0);
  }
  col += vec3(1.0, 0.9, 0.8) * uFlash;
  gl_FragColor = vec4(col, 1.0);
}`;

// --- horizon with Planck cells and the island --------------------------------------------------
const HOR_VERT = /* glsl */ `
varying vec3 vObj;
varying vec3 vNv;
void main() {
  vObj = position;
  vNv = normalize(normalMatrix * normal);
  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
}`;

const HOR_FRAG = /* glsl */ `
${NOISE_GLSL}
uniform float uR;
uniform float uCell;
uniform float uTime;
uniform float uIsland;
uniform float uIslandFlash;
uniform float uHit;
uniform vec3 uHitDir;
varying vec3 vObj;
varying vec3 vNv;
vec3 worley(vec3 p, out vec3 cid) {
  vec3 c = floor(p);
  float f1 = 9.0, f2 = 9.0;
  cid = c;
  for (int i = -1; i <= 1; i++) for (int j = -1; j <= 1; j++) for (int k = -1; k <= 1; k++) {
    vec3 cc = c + vec3(float(i), float(j), float(k));
    vec3 q = cc + 0.15 + 0.7 * hash33(cc);
    float d = length(p - q);
    if (d < f1) { f2 = f1; f1 = d; cid = cc; } else if (d < f2) { f2 = d; }
  }
  return vec3(f1, f2, 0.0);
}
void main() {
  vec3 n = normalize(vObj);
  vec3 cid;
  vec3 w = worley(n * uR / uCell, cid);
  float edge = 1.0 - smoothstep(0.0, 0.06, w.y - w.x);
  float ph = hash13(cid * 1.7 + 3.1);
  float flick = pow(0.5 + 0.5 * sin(uTime * (0.8 + 2.2 * ph) + ph * 40.0), 40.0);
  float rho = length(vNv.xy);              // 0 at the centre of the disk, 1 at the limb
  float limb = 1.0 - smoothstep(0.78, 0.98, rho);
  vec3 gold = vec3(1.0, 0.66, 0.22);
  vec3 col = gold * (edge * 0.15 + flick * 0.2) * limb * (0.55 + 0.45 * (1.0 - rho));
  // island: a thin magenta shell just inside, seen through the dark horizon
  float s = rho / 0.84;
  float shell = exp(-pow((s - 1.0) / 0.018, 2.0)) * 0.9 + (s < 1.0 ? 0.012 : 0.0);
  col += vec3(1.0, 0.31, 0.85) * shell * uIsland * (0.8 + 0.2 * sin(uTime * 2.0));
  col += vec3(1.0, 0.45, 0.9) * uIslandFlash * exp(-pow((s - 1.0) / 0.12, 2.0)) * 1.2;
  // the diary's impact ripples across the cells
  float hd = acos(clamp(dot(n, normalize(uHitDir)), -1.0, 1.0));
  col += vec3(0.6, 0.95, 1.0) * uHit * exp(-pow((hd - (1.0 - uHit) * 2.4) * 6.0, 2.0)) * 1.4;
  gl_FragColor = vec4(col, 1.0);
}`;

// --- particles and threads ----------------------------------------------------------------------
const OUT_VERT = /* glsl */ `
${PAIR_GLSL}
attribute vec3 aDir;
attribute float aTe;
attribute float aSp;
attribute float aSize;
attribute float aBright;
uniform float uT;
uniform float uScale;
varying vec3 vColor;
varying float vAlpha;
void main() {
  vec3 p = outPos(aDir, aTe, aSp, uT);
  float age = uT - aTe;
  float r = length(p);
  vAlpha = age < 0.0 ? 0.0 : smoothstep(0.0, 0.01, age) * (1.0 - smoothstep(16.0, 30.0, r));
  vColor = vec3(0.75, 0.88, 1.25) * aBright;
  vec4 mv = modelViewMatrix * vec4(p, 1.0);
  gl_PointSize = min(aSize * uScale / max(0.05, -mv.z), 14.0);
  gl_Position = projectionMatrix * mv;
}`;

const PARTNER_VERT = /* glsl */ `
${PAIR_GLSL}
attribute vec3 aDir;
attribute float aTe;
attribute float aTau;
attribute vec3 aDirJ;
attribute float aTeJ;
attribute float aSpJ;
uniform float uT;
uniform float uScale;
uniform float uXray;
varying vec3 vColor;
varying float vAlpha;
void main() {
  float mig = smoothstep(aTau, aTau + 0.035, uT);
  vec3 p = mix(inPos(aDir, aTe, uT), outPos(aDirJ, aTeJ, aSpJ, uT), mig);
  vAlpha = (uT < aTe ? 0.0 : 1.0) * (1.0 - mig) * uXray;
  vColor = vec3(1.0, 0.22, 0.12) * 1.2;
  vec4 mv = modelViewMatrix * vec4(p, 1.0);
  gl_PointSize = min(0.22 * uScale / max(0.05, -mv.z), 9.0);
  gl_Position = projectionMatrix * mv;
}`;

const THREAD_VERT = /* glsl */ `
${PAIR_GLSL}
attribute float aEnd;
attribute vec3 aDir;
attribute float aTe;
attribute float aSp;
attribute float aTau;
attribute vec3 aDirJ;
attribute float aTeJ;
attribute float aSpJ;
uniform float uT;
varying vec3 vColor;
varying float vAlpha;
void main() {
  float mig = smoothstep(aTau, aTau + 0.035, uT);
  vec3 outer = outPos(aDir, aTe, aSp, uT);
  vec3 inner = mix(inPos(aDir, aTe, uT), outPos(aDirJ, aTeJ, aSpJ, uT), mig);
  vec3 p = aEnd < 0.5 ? outer : inner;
  float age = uT - aTe;
  float span = length(outer - inner);
  vAlpha = age < 0.0 ? 0.0 : smoothstep(0.0, 0.015, age) * (1.0 - smoothstep(12.0, 24.0, length(outer))) * (1.0 - 0.7 * smoothstep(4.0, 14.0, span));
  vColor = mix(vec3(1.0, 0.70, 0.28), vec3(0.35, 0.95, 1.0), mig) * mix(mix(0.3, 0.8, mig), 1.0, aEnd);
  gl_Position = projectionMatrix * modelViewMatrix * vec4(p, 1.0);
}`;

const THREAD_FRAG = /* glsl */ `
uniform float uOpacity;
varying vec3 vColor;
varying float vAlpha;
void main() { gl_FragColor = vec4(vColor * vAlpha * uOpacity, 1.0); }`;

type Scenario = 'early' | 'late';

/** Hawking emission time of the k-th of n quanta: (1 - m^2) = (k + 1/2) / n. */
const emissionTime = (u: number): number => 1 - Math.pow(1 - clamp(u), 1.5);

/** First time at which the number of re-routed threads N (S_H - S(R)) / beta reaches k. */
function rerouteTime(k: number, n: number): number {
  const rerouted = (t: number): number => (n * (hawkingEntropy(t) - pageCurve(t))) / BETA;
  if (rerouted(1) < k) return 2;
  let lo = T_PAGE;
  let hi = 1;
  for (let i = 0; i < 50; i++) {
    const mid = 0.5 * (lo + hi);
    if (rerouted(mid) >= k) hi = mid;
    else lo = mid;
  }
  return hi;
}

function unitDir(r: () => number): THREE.Vector3 {
  const z = 2 * r() - 1;
  const ph = 2 * Math.PI * r();
  const s = Math.sqrt(1 - z * z);
  return new THREE.Vector3(s * Math.cos(ph), z, s * Math.sin(ph));
}

class BlackHole extends Act {
  readonly id = 'breakthrough' as const;
  readonly duration = 19;
  private readonly lensMat: THREE.ShaderMaterial;
  private readonly horizon: THREE.Mesh<THREE.SphereGeometry, THREE.ShaderMaterial>;
  private readonly outMat: THREE.ShaderMaterial;
  private readonly dustMat: THREE.ShaderMaterial;
  private readonly partnerMat: THREE.ShaderMaterial;
  private readonly threadMat: THREE.ShaderMaterial;
  private readonly threadXray: THREE.ShaderMaterial;
  private readonly flash: THREE.Sprite;
  private readonly diary: THREE.Group;
  private readonly diaryBack: THREE.Group;
  private readonly orbit: DragOrbit;
  private scenario: Scenario = 'early';
  private readonly runs: Record<Scenario, EvaporationRun | null>;
  // HUD
  private readonly hudCanvas: HTMLCanvasElement;
  private readonly hudRead: HTMLElement;
  private readonly diaryLabel: HTMLElement;
  private readonly dock: HTMLElement;
  private readonly scenBtns: Record<Scenario, HTMLButtonElement>;
  private hudKey = '';

  constructor(ctx: ActContext, horizon: HorizonBundle | null) {
    super(ctx);
    this.camera.fov = 42;
    this.bloom = { strength: 0.8, radius: 0.28, threshold: 0.72 };
    this.runs = { early: horizon?.evaporation?.early_diary ?? null, late: horizon?.evaporation?.late_diary ?? null };

    // lensed sky (full-screen, drawn first)
    this.lensMat = new THREE.ShaderMaterial({
      uniforms: {
        uInvProj: { value: new THREE.Matrix4() },
        uCamWorld: { value: new THREE.Matrix4() },
        uCamPos: { value: new THREE.Vector3() },
        uM: { value: 1 },
        uFlash: { value: 0 },
      },
      vertexShader: LENS_VERT,
      fragmentShader: LENS_FRAG,
      depthTest: false,
      depthWrite: false,
    });
    const quad = new THREE.Mesh(new THREE.PlaneGeometry(2, 2), this.lensMat);
    quad.frustumCulled = false;
    quad.renderOrder = -100;
    this.scene.add(quad);
    this.scene.background = null;

    // horizon image
    this.horizon = new THREE.Mesh(
      new THREE.SphereGeometry(1, 128, 64),
      new THREE.ShaderMaterial({
        uniforms: {
          uR: { value: BC0 },
          uCell: { value: 0.95 },
          uTime: { value: 0 },
          uIsland: { value: 0 },
          uIslandFlash: { value: 0 },
          uHit: { value: 0 },
          uHitDir: { value: new THREE.Vector3(0.4, 0.3, 1).normalize() },
        },
        vertexShader: HOR_VERT,
        fragmentShader: HOR_FRAG,
      }),
    );
    this.scene.add(this.horizon);

    // Hawking pairs (deterministic)
    const r = rng(2024);
    const dirs: THREE.Vector3[] = [];
    const te: number[] = [];
    const sp: number[] = [];
    for (let i = 0; i < N_PAIRS; i++) {
      dirs.push(unitDir(r));
      te.push(emissionTime((i + 0.5) / N_PAIRS));
      sp.push(20 + 8 * r());
    }
    const tau = Array.from({ length: N_PAIRS }, (_, i) => rerouteTime(i + 1, N_PAIRS));
    // the partner of an early quantum migrates to a quantum emitted after its re-route time
    const used = new Set<number>();
    const target = tau.map((tt, i) => {
      let best = -1;
      let bestDot = -2;
      for (let k = 0; k < N_PAIRS; k++) {
        if (te[k]! < tt + 0.004 || te[k]! > tt + 0.25 || used.has(k)) continue;
        const d = dirs[k]!.dot(dirs[i]!);
        if (d > bestDot) [best, bestDot] = [k, d];
      }
      if (best < 0) best = N_PAIRS - 1 - (used.size % 12);
      used.add(best);
      return best;
    });

    const outGeo = new THREE.BufferGeometry();
    const f3 = (fn: (i: number) => THREE.Vector3, n: number): Float32Array => {
      const a = new Float32Array(3 * n);
      for (let i = 0; i < n; i++) fn(i).toArray(a, 3 * i);
      return a;
    };
    const f1 = (fn: (i: number) => number, n: number): Float32Array => Float32Array.from({ length: n }, (_, i) => fn(i));
    outGeo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(3 * N_PAIRS), 3));
    outGeo.setAttribute('aDir', new THREE.BufferAttribute(f3((i) => dirs[i]!, N_PAIRS), 3));
    outGeo.setAttribute('aTe', new THREE.BufferAttribute(f1((i) => te[i]!, N_PAIRS), 1));
    outGeo.setAttribute('aSp', new THREE.BufferAttribute(f1((i) => sp[i]!, N_PAIRS), 1));
    outGeo.setAttribute('aSize', new THREE.BufferAttribute(f1(() => 0.2, N_PAIRS), 1));
    outGeo.setAttribute('aBright', new THREE.BufferAttribute(f1(() => 1.35, N_PAIRS), 1));
    this.outMat = pointsMat(OUT_VERT, { uT: 0, uScale: 800 });
    this.scene.add(noCull(new THREE.Points(outGeo, this.outMat)));

    // radiation dust (no threads drawn)
    const dGeo = new THREE.BufferGeometry();
    const rd = rng(77);
    const dDirs = Array.from({ length: N_DUST }, () => unitDir(rd));
    dGeo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(3 * N_DUST), 3));
    dGeo.setAttribute('aDir', new THREE.BufferAttribute(f3((i) => dDirs[i]!, N_DUST), 3));
    dGeo.setAttribute('aTe', new THREE.BufferAttribute(f1((i) => emissionTime((i + rd()) / N_DUST), N_DUST), 1));
    dGeo.setAttribute('aSp', new THREE.BufferAttribute(f1(() => 16 + 14 * rd(), N_DUST), 1));
    dGeo.setAttribute('aSize', new THREE.BufferAttribute(f1(() => 0.035 + 0.08 * rd() ** 3, N_DUST), 1));
    dGeo.setAttribute('aBright', new THREE.BufferAttribute(f1(() => 0.45 + 0.5 * rd(), N_DUST), 1));
    this.dustMat = pointsMat(OUT_VERT, { uT: 0, uScale: 800 });
    this.scene.add(noCull(new THREE.Points(dGeo, this.dustMat)));

    // partners (drawn through the horizon, dim) and threads
    const pairAttrs = (geo: THREE.BufferGeometry, idx: (v: number) => number, n: number): void => {
      geo.setAttribute('aDir', new THREE.BufferAttribute(f3((v) => dirs[idx(v)]!, n), 3));
      geo.setAttribute('aTe', new THREE.BufferAttribute(f1((v) => te[idx(v)]!, n), 1));
      geo.setAttribute('aSp', new THREE.BufferAttribute(f1((v) => sp[idx(v)]!, n), 1));
      geo.setAttribute('aTau', new THREE.BufferAttribute(f1((v) => tau[idx(v)]!, n), 1));
      geo.setAttribute('aDirJ', new THREE.BufferAttribute(f3((v) => dirs[target[idx(v)]!]!, n), 3));
      geo.setAttribute('aTeJ', new THREE.BufferAttribute(f1((v) => te[target[idx(v)]!]!, n), 1));
      geo.setAttribute('aSpJ', new THREE.BufferAttribute(f1((v) => sp[target[idx(v)]!]!, n), 1));
    };
    const pGeo = new THREE.BufferGeometry();
    pGeo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(3 * N_PAIRS), 3));
    pairAttrs(pGeo, (v) => v, N_PAIRS);
    this.partnerMat = pointsMat(PARTNER_VERT, { uT: 0, uScale: 800, uXray: 0.85 });
    this.partnerMat.depthTest = false;
    this.scene.add(noCull(new THREE.Points(pGeo, this.partnerMat)));

    const tGeo = new THREE.BufferGeometry();
    tGeo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(6 * N_PAIRS), 3));
    tGeo.setAttribute('aEnd', new THREE.BufferAttribute(f1((v) => v % 2, 2 * N_PAIRS), 1));
    pairAttrs(tGeo, (v) => v >> 1, 2 * N_PAIRS);
    const threadMat = (opacity: number, depthTest: boolean): THREE.ShaderMaterial =>
      new THREE.ShaderMaterial({
        uniforms: { uT: { value: 0 }, uOpacity: { value: opacity } },
        vertexShader: THREAD_VERT,
        fragmentShader: THREAD_FRAG,
        transparent: true,
        blending: THREE.AdditiveBlending,
        depthWrite: false,
        depthTest,
      });
    this.threadMat = threadMat(0.55, true);
    this.threadXray = threadMat(0.14, false);
    this.scene.add(noCull(new THREE.LineSegments(tGeo, this.threadMat)), noCull(new THREE.LineSegments(tGeo, this.threadXray)));

    // final flash
    this.flash = new THREE.Sprite(new THREE.SpriteMaterial({ map: radialTexture(), color: 0xffffff, blending: THREE.AdditiveBlending, depthWrite: false, depthTest: false, transparent: true }));
    this.scene.add(this.flash);

    // diary: a glowing cube that falls in, and its recovered copy in the radiation
    this.diary = diaryCube(0xbff8ff, 1.0);
    this.diaryBack = diaryCube(0x59f3ff, 1.2);
    this.scene.add(this.diary, this.diaryBack);

    // --- HUD ---------------------------------------------------------------------------------
    this.hudCanvas = h('canvas', { class: 'page-canvas', width: '380', height: '210', 'aria-label': 'Page curve' }) as HTMLCanvasElement;
    this.hudRead = h('div', { class: 'page-read' });
    const panel = h(
      'div',
      { class: 'show-panel bh-panel show-ui' },
      h('div', { class: 'panel-k' }, h('span', null, 'Radiation entropy S(R)'), h('b', null, 'live')),
      this.hudCanvas,
      h(
        'div',
        { class: 'page-legend' },
        h('span', { class: 'lg lg-hawking' }, 'Hawking · field theory'),
        h('span', { class: 'lg lg-page' }, 'Page curve · min{S_H, A/4G}'),
        h('span', { class: 'lg lg-cap' }, 'capacity A/4G'),
      ),
      this.hudRead,
    );
    this.diaryLabel = h('div', { class: 'diary-label' });
    this.scenBtns = {
      early: button('Early · after 5 emissions', () => this.throwDiary('early'), 'show-btn'),
      late: button('Late · after 14 emissions', () => this.throwDiary('late'), 'show-btn'),
    };
    this.dock = h('div', { class: 'show-dock bh-dock show-ui' }, h('span', { class: 'dock-k' }, 'Throw in a diary'), this.scenBtns.early, this.scenBtns.late);
    ctx.hud.append(panel, this.diaryLabel, this.dock);

    this.orbit = new DragOrbit(ctx.stage.canvas.parentElement ?? ctx.stage.canvas, () => ctx.userInput(), 0.9);
    this.own(() => this.orbit.dispose());
    this.update(0);
  }

  private throwDiary(s: Scenario): void {
    this.ctx.userInput();
    this.scenario = s;
    const run = this.runs[s];
    const step = run?.diary_step ?? (s === 'early' ? 5 : 14);
    this.time = Math.max(0, (step / N_STEPS - 0.05) * LIFE);
    this.isPlaying = true;
  }

  /** lifetime fraction at act time t */
  private life(t: number): number {
    return clamp(t / LIFE);
  }

  protected update(dt: number): void {
    const t = this.time;
    const L = this.life(t);
    const m = massFraction(L);
    this.orbit.step(dt);

    // camera: slow orbit, closing in as the hole shrinks
    const yaw = 0.25 + 0.022 * t + (this.isPlaying ? 0 : 0.02 * this.wall) + this.orbit.yaw;
    const pitch = 0.16 + this.orbit.pitch;
    const dist = 17 + 11 * Math.pow(m, 0.8) + 4 * ramp(t, LIFE, this.duration);
    this.camera.position.set(dist * Math.cos(pitch) * Math.sin(yaw), dist * Math.sin(pitch), dist * Math.cos(pitch) * Math.cos(yaw));
    this.camera.lookAt(0, 0, 0);
    this.camera.updateMatrixWorld();
    this.camera.updateProjectionMatrix();
    const lu = this.lensMat.uniforms;
    (lu.uInvProj!.value as THREE.Matrix4).copy(this.camera.projectionMatrixInverse);
    (lu.uCamWorld!.value as THREE.Matrix4).copy(this.camera.matrixWorld);
    (lu.uCamPos!.value as THREE.Vector3).copy(this.camera.position);
    lu.uM!.value = L >= 0.9995 ? 0 : m;

    // horizon image
    const R = BC0 * m;
    this.horizon.visible = R > 0.02;
    this.horizon.scale.setScalar(Math.max(R, 1e-3));
    const hu = this.horizon.material.uniforms;
    hu.uR!.value = R;
    hu.uTime!.value = this.wall;
    const island = ramp(L, T_PAGE, T_PAGE + 0.03);
    hu.uIsland!.value = island;
    hu.uIslandFlash!.value = Math.exp(-Math.max(0, L - T_PAGE) * 60) * (L >= T_PAGE ? 1 : 0);

    // particles
    const sc = pointScale(this.camera, this.height, this.pixelRatio);
    // the radiation keeps flying (and the last partners finish migrating) after the hole is gone
    const Lr = t / LIFE;
    for (const mat of [this.outMat, this.dustMat, this.partnerMat]) {
      mat.uniforms.uT!.value = Lr;
      mat.uniforms.uScale!.value = sc;
    }
    this.threadMat.uniforms.uT!.value = Lr;
    this.threadXray.uniforms.uT!.value = Lr;

    // final flash and fade of the shadow
    const fl = Math.exp(-Math.pow((t - LIFE) / 0.3, 2)) * (t < LIFE ? 0.6 : 1);
    this.flash.scale.setScalar(2 + 22 * fl);
    (this.flash.material as THREE.SpriteMaterial).opacity = Math.min(1, fl);
    lu.uFlash!.value = 0.12 * fl;

    this.updateDiary(L);
    this.paintHud(L, m, island);
  }

  private updateDiary(L: number): void {
    const run = this.runs[this.scenario];
    const step = run?.diary_step ?? (this.scenario === 'early' ? 5 : 14);
    const rec = run?.recovered_step ?? (this.scenario === 'early' ? 13 : 17);
    const t0 = step / N_STEPS;
    const tIn = t0 + 0.035;
    const tRec = rec / N_STEPS;
    // fall: from the camera side into the hole along a curving path
    const fall = clamp((L - (t0 - 0.03)) / (tIn - (t0 - 0.03)));
    const on = L >= t0 - 0.03 && L < tIn;
    this.diary.visible = on;
    const camDir = this.camera.position.clone().normalize();
    const side = new THREE.Vector3(0, 1, 0).cross(camDir).normalize();
    const start = camDir.clone().multiplyScalar(12).add(side.clone().multiplyScalar(-7)).add(new THREE.Vector3(0, 3.5, 0));
    const R = BC0 * massFraction(L);
    if (on) {
      const e = fall * fall;
      const p = start.clone().multiplyScalar(1 - e).add(camDir.clone().multiplyScalar(R * 0.95).multiplyScalar(e));
      this.diary.position.copy(p);
      this.diary.scale.setScalar(lerp(1.0, 0.2, e));
      this.diary.rotation.set(L * 40, L * 55, 0);
    }
    // the impact ripple on the cells
    const hu = this.horizon.material.uniforms;
    const hit = L >= tIn ? Math.exp(-(L - tIn) * 45) : 0;
    hu.uHit!.value = hit;
    (hu.uHitDir!.value as THREE.Vector3).copy(camDir);
    // recovery: the diary re-forms in the radiation cloud
    const back = ramp(L, tRec - 0.012, tRec + 0.012) * (1 - ramp(L, tRec + 0.2, tRec + 0.26));
    this.diaryBack.visible = back > 0.001;
    const outDir = camDir.clone().multiplyScalar(0.35).add(side.clone().multiplyScalar(-0.85)).add(new THREE.Vector3(0, -0.2, 0)).normalize();
    const rr = 8.5 + 6 * Math.max(0, L - tRec);
    this.diaryBack.position.copy(outDir.multiplyScalar(rr));
    this.diaryBack.scale.setScalar(0.9 * back + 0.6 * Math.exp(-Math.max(0, L - tRec) * 40) * back);
    this.diaryBack.rotation.set(L * 20, L * 30, 0);
    for (const c of this.diaryBack.children) ((c as THREE.LineSegments).material as THREE.LineBasicMaterial).opacity = back;

    // label next to the recovered diary
    const shown = back > 0.5 ? 1 : on ? 1 : 0;
    const anchor = back > 0.5 ? this.diaryBack.position : this.diary.position;
    const p = anchor.clone().project(this.camera);
    this.diaryLabel.style.opacity = String(shown);
    this.diaryLabel.style.transform = `translate(${((p.x + 1) / 2) * this.width + 18}px, ${((1 - p.y) / 2) * this.height - 10}px)`;
    const mi = run ? interpStep(run.diary_mi, L * N_STEPS) : 0;
    this.diaryLabel.textContent =
      back > 0.5 ? `diary recovered · step ${rec} · ${mi.toFixed(2)} / 2 bits in the radiation` : `diary · thrown in after ${step} emissions`;
    this.diaryLabel.classList.toggle('back', back > 0.5);
  }

  private paintHud(L: number, m: number, island: number): void {
    const key = `${L.toFixed(4)}:${this.scenario}:${this.width}:${this.isPlaying}`;
    if (key === this.hudKey) return;
    this.hudKey = key;
    drawPageCurve(this.hudCanvas, L);
    const run = this.runs[this.scenario];
    const bits = run ? interpStep(run.s_rad, L * N_STEPS) : null;
    const sR = pageCurve(L);
    this.hudRead.replaceChildren(
      readout('mass left', `${(m * 100).toFixed(1)}%`),
      readout('S(R) / S_BH(0)', sR.toFixed(3)),
      readout('island', island > 0.5 ? 'on' : 'off', island > 0.5 ? 'island' : ''),
      readout('20-qubit run', bits === null ? '—' : `${bits.toFixed(2)} bits`),
    );
    for (const s of ['early', 'late'] as const) this.scenBtns[s].setAttribute('aria-pressed', String(this.scenario === s));
    this.dock.classList.toggle('on', !this.isPlaying || this.time > 13);
  }
}

function readout(k: string, v: string, cls = ''): HTMLElement {
  return h('div', { class: `rd ${cls}`.trim() }, h('span', { class: 'rk' }, k), h('span', { class: 'rv mono' }, v));
}

/** linear interpolation of a per-step array (steps 1..n) at a fractional step */
function interpStep(arr: number[], s: number): number {
  if (!arr.length) return 0;
  const x = clamp(s - 1, 0, arr.length - 1);
  const i = Math.floor(x);
  const f = x - i;
  return arr[i]! * (1 - f) + (arr[Math.min(i + 1, arr.length - 1)] ?? arr[i]!) * f;
}

function drawPageCurve(c: HTMLCanvasElement, L: number): void {
  const dpr = Math.min(2, window.devicePixelRatio || 1);
  const W = 380;
  const H = 210;
  if (c.width !== W * dpr) {
    c.width = W * dpr;
    c.height = H * dpr;
  }
  const g = c.getContext('2d');
  if (!g) return;
  g.setTransform(dpr, 0, 0, dpr, 0, 0);
  g.clearRect(0, 0, W, H);
  const pad = { l: 30, r: 10, t: 12, b: 24 };
  const X = (t: number): number => pad.l + t * (W - pad.l - pad.r);
  const Y = (s: number): number => H - pad.b - (s / 1.95) * (H - pad.t - pad.b);
  g.strokeStyle = 'rgba(160,190,255,0.14)';
  g.lineWidth = 1;
  for (const s of [0.5, 1, 1.5]) {
    g.beginPath();
    g.moveTo(X(0), Y(s));
    g.lineTo(X(1), Y(s));
    g.stroke();
  }
  g.fillStyle = 'rgba(139,151,173,0.9)';
  g.font = '10px "IBM Plex Mono", monospace';
  g.fillText('0', X(0) - 4, H - 8);
  g.fillText('lifetime →', X(1) - 62, H - 8);
  for (const s of [1, 1.855]) g.fillText(s === 1 ? '1' : '1.86', 2, Y(s) + 3);
  // Page time
  const tp = pageTimeFraction();
  g.setLineDash([3, 4]);
  g.strokeStyle = 'rgba(255,79,216,0.55)';
  g.beginPath();
  g.moveTo(X(tp), Y(0));
  g.lineTo(X(tp), Y(1.9));
  g.stroke();
  g.setLineDash([]);
  g.fillStyle = 'rgba(255,79,216,0.85)';
  g.fillText(`Page time ${tp.toFixed(3)}`, X(tp) + 4, Y(1.86));
  const curve = (f: (t: number) => number, color: string, width: number, dash: number[], upTo: number, glow = 0): void => {
    g.save();
    g.setLineDash(dash);
    g.strokeStyle = color;
    g.lineWidth = width;
    if (glow) {
      g.shadowColor = color;
      g.shadowBlur = glow;
    }
    g.beginPath();
    const n = 160;
    for (let i = 0; i <= n; i++) {
      const t = (i / n) * upTo;
      if (i === 0) g.moveTo(X(t), Y(f(t)));
      else g.lineTo(X(t), Y(f(t)));
    }
    g.stroke();
    g.restore();
  };
  // full curves faint, live part bright
  curve(hawkingEntropy, 'rgba(255,122,58,0.22)', 1.5, [4, 4], 1);
  curve(bhEntropy, 'rgba(255,195,90,0.22)', 1.2, [], 1);
  curve(pageCurve, 'rgba(89,243,255,0.2)', 2, [], 1);
  if (L > 0) {
    curve(hawkingEntropy, 'rgba(255,122,58,0.95)', 1.6, [4, 4], L);
    curve(bhEntropy, 'rgba(255,195,90,0.75)', 1.3, [], L);
    curve(pageCurve, '#59f3ff', 2.6, [], L, 10);
  }
  // markers
  const dot = (x: number, y: number, color: string, r: number): void => {
    g.beginPath();
    g.fillStyle = color;
    g.shadowColor = color;
    g.shadowBlur = 12;
    g.arc(x, y, r, 0, Math.PI * 2);
    g.fill();
    g.shadowBlur = 0;
  };
  g.strokeStyle = 'rgba(233,238,248,0.35)';
  g.beginPath();
  g.moveTo(X(L), Y(0));
  g.lineTo(X(L), Y(1.9));
  g.stroke();
  dot(X(L), Y(hawkingEntropy(L)), '#ff7a3a', 3);
  dot(X(L), Y(pageCurve(L)), '#bff8ff', 4);
}

function pointsMat(vert: string, u: Record<string, number>): THREE.ShaderMaterial {
  return new THREE.ShaderMaterial({
    uniforms: Object.fromEntries(Object.entries(u).map(([k, v]) => [k, { value: v }])),
    vertexShader: vert,
    fragmentShader: POINTS_FRAG,
    transparent: true,
    blending: THREE.AdditiveBlending,
    depthWrite: false,
  });
}

function noCull<T extends THREE.Object3D>(o: T): T {
  o.frustumCulled = false;
  return o;
}

function radialTexture(): THREE.CanvasTexture {
  const c = document.createElement('canvas');
  c.width = c.height = 256;
  const g = c.getContext('2d')!;
  const grd = g.createRadialGradient(128, 128, 0, 128, 128, 128);
  grd.addColorStop(0, 'rgba(255,255,255,1)');
  grd.addColorStop(0.1, 'rgba(220,245,255,0.8)');
  grd.addColorStop(0.35, 'rgba(120,200,255,0.25)');
  grd.addColorStop(1, 'rgba(0,0,0,0)');
  g.fillStyle = grd;
  g.fillRect(0, 0, 256, 256);
  return new THREE.CanvasTexture(c);
}

function diaryCube(color: number, size: number): THREE.Group {
  const grp = new THREE.Group();
  const edges = new THREE.EdgesGeometry(new THREE.BoxGeometry(size, size * 1.3, size * 0.35));
  const mat = new THREE.LineBasicMaterial({ color: new THREE.Color(color).multiplyScalar(2.2), transparent: true, blending: THREE.AdditiveBlending, depthWrite: false });
  grp.add(new THREE.LineSegments(edges, mat));
  const inner = new THREE.LineSegments(new THREE.EdgesGeometry(new THREE.BoxGeometry(size * 0.6, size * 0.06, 0.01)), mat);
  inner.position.y = size * 0.25;
  const inner2 = inner.clone();
  inner2.position.y = 0;
  const inner3 = inner.clone();
  inner3.position.y = -size * 0.25;
  grp.add(inner, inner2, inner3);
  grp.visible = false;
  return grp;
}

export const blackhole: ActModule = {
  async mount(ctx) {
    const horizon = await ctx.data.horizon();
    return new BlackHole(ctx, horizon);
  },
};
