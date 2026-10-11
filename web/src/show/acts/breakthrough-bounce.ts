/**
 * Act III·B · Breakthrough — "The Big Bang bounces" (docs/SHOW.md).
 *
 * A luminous landscape of the exact |Psi(a, T)|^2 for a_max = 2 ([QG] Theorems Q, Lambda; the odd
 * combination of two coherent states of the half-line oscillator, qchas/horizon/bounce.py):
 *
 *   |Psi|^2 = N^2 sqrt(12/pi) [e^{-12 (a - |x|)^2} (1 - e^{-24 a |x|})^2 + 4 e^{-12 (a^2 + x^2)} sin^2(p a)],
 *   x = a_max cos T, p = -12 a_max sin T,
 *
 * revealed along the radiation time T, with the interference fringes at the bounce. Beside it a
 * universe bubble of radius <a>(T) (bounce.json, exact) bounces at <a>_min = 0.165 ell (a_max = 2;
 * 0.163 ell for large universes, 0.26 l_P), while a red ghost of the classical a_max |cos T|
 * crashes to zero. Graviton modes ride across the bubble as integer-frequency oscillators.
 */

import * as THREE from 'three';
import type { BounceBundle, BounceRun } from '../../lib/data';
import { h } from '../../lib/ui';
import { Act, type ActContext, type ActModule, DragOrbit } from '../act';
import { clamp, ramp, skyDome } from '../gfx';

const A_MAX = 2;
const A_SPAN = 2.8;
const NA = 240;
const NT = 320;
const SX = 12 / Math.PI; // world x per unit T
const SZ = 2.2; // world z per unit a
const HY = 3.0; // world height of the tallest peak
const Z0 = 3.0; // world z of a = 0
const SWEEP = 8.2; // seconds for T: 0 -> pi
const C2 = Math.sqrt(12 / Math.PI);

/** Exact |Psi(a, T)|^2 on a > 0 (module docstring). */
export function bounceDensity(a: number, T: number, amax = A_MAX): number {
  const x = amax * Math.cos(T);
  const p = -12 * amax * Math.sin(T);
  const ax = Math.abs(x);
  const norm2 = 1 / -Math.expm1(-12 * amax * amax);
  return norm2 * C2 * (Math.exp(-12 * (a - ax) ** 2) * Math.expm1(-24 * a * ax) ** 2 + 4 * Math.exp(-12 * (a * a + x * x)) * Math.sin(p * a) ** 2);
}

const SURF_VERT = /* glsl */ `
attribute float aRho;
attribute float aT;
varying float vRho;
varying float vT;
varying vec2 vGrid;
varying vec3 vPos;
void main() {
  vRho = aRho;
  vT = aT;
  vGrid = uv;
  vPos = position;
  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
}`;

const SURF_FRAG = /* glsl */ `
uniform float uReveal;
uniform float uTime;
varying float vRho;
varying float vT;
varying vec2 vGrid;
varying vec3 vPos;
vec3 ramp4(float x) {
  vec3 c0 = vec3(0.01, 0.02, 0.07), c1 = vec3(0.05, 0.22, 0.75), c2 = vec3(0.25, 0.85, 1.1), c3 = vec3(1.25, 1.1, 0.75);
  return x < 0.33 ? mix(c0, c1, x / 0.33) : x < 0.66 ? mix(c1, c2, (x - 0.33) / 0.33) : mix(c2, c3, (x - 0.66) / 0.34);
}
void main() {
  if (vT > uReveal + 0.002) discard;
  float x = pow(clamp(vRho, 0.0, 1.0), 0.5);
  vec3 col = ramp4(x) * (0.18 + 1.05 * x);
  // faint survey grid on the landscape
  vec2 g = abs(fract(vGrid * vec2(24.0, 14.0)) - 0.5);
  float line = 1.0 - smoothstep(0.0, 0.035, min(g.x, g.y));
  col += vec3(0.2, 0.45, 0.9) * line * 0.12;
  // the reveal front
  float front = exp(-pow((vT - uReveal) / 0.012, 2.0));
  col += vec3(0.6, 0.95, 1.2) * front * (0.4 + 2.0 * x);
  gl_FragColor = vec4(col, 1.0);
}`;

const BUBBLE_VERT = /* glsl */ `
uniform float uR;
uniform float uPhase;
uniform float uAmp;
varying vec3 vN;
varying vec3 vView;
varying float vWave;
void main() {
  vec3 n = normalize(position);
  float th = acos(clamp(n.y, -1.0, 1.0));
  float ph = atan(n.z, n.x);
  // graviton modes n = 3, 4, 5: integer-frequency oscillators that run straight through the bounce
  float w = 0.0;
  w += cos(3.0 * uPhase) * pow(sin(th), 3.0) * cos(3.0 * ph);
  w += 0.8 * cos(4.0 * uPhase + 1.0) * pow(sin(th), 4.0) * cos(4.0 * ph + 0.7);
  w += 0.6 * cos(5.0 * uPhase + 2.0) * pow(sin(th), 5.0) * sin(5.0 * ph);
  vWave = w;
  vec3 p = n * uR * (1.0 + uAmp * w);
  vec4 wp = modelMatrix * vec4(p, 1.0);
  vN = normalize(mat3(modelMatrix) * n);
  vView = cameraPosition - wp.xyz;
  gl_Position = projectionMatrix * viewMatrix * wp;
}`;

const BUBBLE_FRAG = /* glsl */ `
uniform vec3 uColor;
uniform float uGlow;
uniform float uWire;
varying vec3 vN;
varying vec3 vView;
varying float vWave;
void main() {
  float f = pow(1.0 - abs(dot(normalize(vN), normalize(vView))), 2.2);
  vec3 col = uColor * (0.06 + 1.2 * f) * uGlow;
  col += vec3(1.0, 0.75, 0.35) * smoothstep(0.35, 0.9, abs(vWave)) * 0.35 * uGlow * uWire;
  gl_FragColor = vec4(col, 1.0);
}`;

function lineMesh(points: THREE.Vector3[], color: THREE.Color, opacity = 1): THREE.Line {
  const g = new THREE.BufferGeometry().setFromPoints(points);
  const m = new THREE.LineBasicMaterial({ color, transparent: true, opacity, blending: THREE.AdditiveBlending, depthWrite: false });
  const l = new THREE.Line(g, m);
  l.frustumCulled = false;
  return l;
}

class Bounce extends Act {
  readonly id = 'bounce' as const;
  readonly duration = 9;
  private readonly surfMat: THREE.ShaderMaterial;
  private readonly bubble: THREE.Mesh<THREE.SphereGeometry, THREE.ShaderMaterial>;
  private readonly ghost: THREE.Mesh<THREE.SphereGeometry, THREE.ShaderMaterial>;
  private readonly marker: THREE.Mesh;
  private readonly ghostMarker: THREE.Mesh;
  private readonly flash: THREE.Mesh<THREE.RingGeometry, THREE.MeshBasicMaterial>;
  private readonly world = new THREE.Group();
  private readonly orbit: DragOrbit;
  private readonly run: BounceRun | null;
  private readonly aMin: number;
  private readonly readT: HTMLElement;
  private readonly readA: HTMLElement;
  private readonly readC: HTMLElement;
  private readonly bounceTag: HTMLElement;
  private readonly bubbleLabel: HTMLElement;
  private readonly ghostLabel: HTMLElement;

  constructor(ctx: ActContext, data: BounceBundle | null) {
    super(ctx);
    this.camera.fov = 38;
    this.bloom = { strength: 0.6, radius: 0.35, threshold: 0.85 };
    this.scene.add(skyDome(1500, 0.5));
    this.scene.add(this.world);
    this.run = data?.runs.find((r) => Math.abs(r.a_max - A_MAX) < 1e-9) ?? null;
    this.aMin = this.run?.a_min ?? 0.16462;

    // --- the |Psi|^2 landscape ------------------------------------------------------------
    const pos = new Float32Array(NA * NT * 3);
    const rho = new Float32Array(NA * NT);
    const tt = new Float32Array(NA * NT);
    const uv = new Float32Array(NA * NT * 2);
    let rmax = 0;
    const raw = new Float32Array(NA * NT);
    for (let j = 0; j < NT; j++) {
      const T = (j / (NT - 1)) * Math.PI;
      for (let i = 0; i < NA; i++) {
        const a = 1e-4 + (i / (NA - 1)) * A_SPAN;
        const d = bounceDensity(a, T);
        raw[j * NA + i] = d;
        rmax = Math.max(rmax, d);
      }
    }
    for (let j = 0; j < NT; j++) {
      const T = (j / (NT - 1)) * Math.PI;
      for (let i = 0; i < NA; i++) {
        const k = j * NA + i;
        const a = (i / (NA - 1)) * A_SPAN;
        const r = raw[k]! / rmax;
        pos.set([(T - Math.PI / 2) * SX, HY * Math.pow(r, 0.75), Z0 - a * SZ], 3 * k);
        rho[k] = r;
        tt[k] = T / Math.PI;
        uv.set([a / A_SPAN, T / Math.PI], 2 * k);
      }
    }
    const idx: number[] = [];
    for (let j = 0; j < NT - 1; j++)
      for (let i = 0; i < NA - 1; i++) {
        const k = j * NA + i;
        idx.push(k, k + 1, k + NA, k + 1, k + NA + 1, k + NA);
      }
    const geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.BufferAttribute(pos, 3));
    geo.setAttribute('aRho', new THREE.BufferAttribute(rho, 1));
    geo.setAttribute('aT', new THREE.BufferAttribute(tt, 1));
    geo.setAttribute('uv', new THREE.BufferAttribute(uv, 2));
    geo.setIndex(idx);
    this.surfMat = new THREE.ShaderMaterial({ uniforms: { uReveal: { value: 0 }, uTime: { value: 0 } }, vertexShader: SURF_VERT, fragmentShader: SURF_FRAG, side: THREE.DoubleSide });
    const surf = new THREE.Mesh(geo, this.surfMat);
    surf.frustumCulled = false;
    this.world.add(surf);

    // floor traces: classical a_max |cos T| (red) and the exact <a>(T) (cyan)
    const floorY = 0.02;
    const cls: THREE.Vector3[] = [];
    const qm: THREE.Vector3[] = [];
    for (let j = 0; j <= 400; j++) {
      const T = (j / 400) * Math.PI;
      cls.push(new THREE.Vector3((T - Math.PI / 2) * SX, floorY, Z0 - A_MAX * Math.abs(Math.cos(T)) * SZ));
      qm.push(new THREE.Vector3((T - Math.PI / 2) * SX, floorY + 0.01, Z0 - this.aMean(T) * SZ));
    }
    const overlay = (l: THREE.Line): THREE.Line => {
      (l.material as THREE.LineBasicMaterial).depthTest = false;
      l.renderOrder = 5;
      return l;
    };
    this.world.add(overlay(lineMesh(cls, new THREE.Color(1.6, 0.25, 0.18), 0.9)), overlay(lineMesh(qm, new THREE.Color(0.4, 1.4, 1.8), 1)));
    // axes: T along x, a into the screen
    this.world.add(
      lineMesh([new THREE.Vector3(-6, 0, Z0), new THREE.Vector3(6.4, 0, Z0)], new THREE.Color(0.35, 0.45, 0.7), 0.6),
      lineMesh([new THREE.Vector3(-6, 0, Z0), new THREE.Vector3(-6, 0, Z0 - A_SPAN * SZ)], new THREE.Color(0.35, 0.45, 0.7), 0.6),
    );
    // current-time markers on the floor
    this.marker = new THREE.Mesh(new THREE.SphereGeometry(0.09, 16, 8), new THREE.MeshBasicMaterial({ color: new THREE.Color(0.8, 2.2, 2.6) }));
    this.ghostMarker = new THREE.Mesh(new THREE.SphereGeometry(0.075, 16, 8), new THREE.MeshBasicMaterial({ color: new THREE.Color(2.4, 0.4, 0.3) }));
    this.world.add(this.marker, this.ghostMarker);

    // --- the universe bubble and its classical ghost -----------------------------------------
    const bubbleMat = (color: THREE.Color, glow: number, wire: number): THREE.ShaderMaterial =>
      new THREE.ShaderMaterial({
        uniforms: { uR: { value: 1 }, uPhase: { value: 0 }, uAmp: { value: 0.06 }, uColor: { value: color }, uGlow: { value: glow }, uWire: { value: wire } },
        vertexShader: BUBBLE_VERT,
        fragmentShader: BUBBLE_FRAG,
        transparent: true,
        blending: THREE.AdditiveBlending,
        depthWrite: false,
        side: THREE.DoubleSide,
      });
    this.bubble = new THREE.Mesh(new THREE.SphereGeometry(1, 96, 48), bubbleMat(new THREE.Color(0.35, 0.9, 1.2), 0.42, 1));
    this.ghost = new THREE.Mesh(new THREE.SphereGeometry(1, 32, 16), bubbleMat(new THREE.Color(1.3, 0.18, 0.12), 0.3, 0));
    (this.ghost.material as THREE.ShaderMaterial).wireframe = true;
    this.ghost.material.uniforms.uAmp!.value = 0;
    this.scene.add(this.bubble, this.ghost);
    this.flash = new THREE.Mesh(new THREE.RingGeometry(0.9, 1, 96), new THREE.MeshBasicMaterial({ color: new THREE.Color(1.5, 1.3, 2.2), transparent: true, opacity: 0, blending: THREE.AdditiveBlending, depthWrite: false, side: THREE.DoubleSide }));
    this.scene.add(this.flash);

    // --- HUD ----------------------------------------------------------------------------------
    this.readT = h('span', { class: 'rv mono' });
    this.readA = h('span', { class: 'rv mono' });
    this.readC = h('span', { class: 'rv mono' });
    this.bounceTag = h('div', { class: 'bounce-tag' }, h('span', { class: 'bt-k' }, 'bounce'), h('span', { class: 'bt-v' }, `⟨a⟩min = ${this.aMin.toFixed(3)} ℓ ≈ ${(this.aMin * 1.5958).toFixed(2)} ℓP`));
    const panel = h(
      'div',
      { class: 'show-panel bounce-panel show-ui' },
      h('div', { class: 'panel-k' }, h('span', null, 'Closed radiation universe · a_max = 2'), h('b', null, 'exact')),
      h(
        'div',
        { class: 'page-read three' },
        h('div', { class: 'rd' }, h('span', { class: 'rk' }, 'radiation time T'), this.readT),
        h('div', { class: 'rd' }, h('span', { class: 'rk' }, 'quantum ⟨a⟩'), this.readA),
        h('div', { class: 'rd ghost' }, h('span', { class: 'rk' }, 'classical a'), this.readC),
      ),
      h(
        'div',
        { class: 'bounce-legend' },
        h('span', { class: 'lg lg-q' }, '|Ψ(a, T)|² landscape'),
        h('span', { class: 'lg lg-qm' }, '⟨a⟩(T), exact'),
        h('span', { class: 'lg lg-cl' }, 'classical a_max |cos T|'),
      ),
    );
    this.bubbleLabel = h('div', { class: 'obj-label' }, h('b', null, 'Universe'), ' · radius ⟨a⟩(T)');
    this.ghostLabel = h('div', { class: 'obj-label ghost-label' }, 'classical: crushed to a = 0');
    ctx.hud.append(panel, this.bounceTag, this.bubbleLabel, this.ghostLabel);

    this.orbit = new DragOrbit(ctx.stage.canvas.parentElement ?? ctx.stage.canvas, () => ctx.userInput(), 0.8);
    this.own(() => this.orbit.dispose());
    this.update(0);
  }

  /** exact <a>(T) for a_max = 2 from bounce.json (periodic, period pi), classical fallback */
  private aMean(T: number): number {
    const run = this.run;
    if (!run || run.T.length < 2) return Math.max(A_MAX * Math.abs(Math.cos(T)), 0.165);
    const per = Math.PI;
    let x = ((T % per) + per) % per;
    const Ts = run.T;
    if (x > Ts[Ts.length - 1]!) x = Ts[Ts.length - 1]!;
    let lo = 0;
    let hi = Ts.length - 1;
    while (hi - lo > 1) {
      const mid = (lo + hi) >> 1;
      if (Ts[mid]! <= x) lo = mid;
      else hi = mid;
    }
    const f = (x - Ts[lo]!) / Math.max(1e-12, Ts[hi]! - Ts[lo]!);
    return run.a_mean[lo]! * (1 - f) + run.a_mean[hi]! * f;
  }

  protected update(dt: number): void {
    const t = this.time;
    this.orbit.step(dt);
    const scripted = this.isPlaying || t < this.duration - 1e-3;
    // radiation time: one sweep 0 -> pi in the script, then the universe keeps breathing
    const T = scripted ? Math.PI * clamp(t / SWEEP) : Math.PI * (1 + ((this.wall * 0.38) % 1));
    const reveal = scripted ? clamp(t / SWEEP) : 1;
    this.surfMat.uniforms.uReveal!.value = reveal;
    const Tm = T % Math.PI;

    // camera
    const yaw = 0.1 + 0.012 * t + (this.isPlaying ? 0 : 0.2 * Math.sin(0.12 * this.wall)) + this.orbit.yaw;
    const pitch = 0.6 + this.orbit.pitch * 0.6;
    const dist = 18.5;
    const target = new THREE.Vector3(1.4, 0.4, -0.4);
    this.camera.position.set(target.x + dist * Math.cos(pitch) * Math.sin(yaw), target.y + dist * Math.sin(pitch), target.z + dist * Math.cos(pitch) * Math.cos(yaw));
    this.camera.lookAt(target);

    // floor markers
    const x = (Tm - Math.PI / 2) * SX;
    const aq = this.aMean(Tm);
    const ac = A_MAX * Math.abs(Math.cos(Tm));
    this.marker.position.set(x, 0.05, Z0 - aq * SZ);
    this.ghostMarker.position.set(x, 0.05, Z0 - ac * SZ);

    // bubble beside the landscape
    const bpos = new THREE.Vector3(8.4, 4.6, -3.4);
    const R = 0.95 * aq;
    this.bubble.position.copy(bpos);
    this.ghost.position.copy(bpos);
    const bu = this.bubble.material.uniforms;
    bu.uR!.value = R;
    bu.uPhase!.value = T * 2;
    bu.uAmp!.value = 0.05;
    const gu = this.ghost.material.uniforms;
    gu.uR!.value = Math.max(1e-4, 0.95 * ac);
    gu.uGlow!.value = 0.3 * (ac > 0.02 ? 1 : 0);
    // the bounce flash
    const nearB = Math.exp(-Math.pow((Tm - Math.PI / 2) / 0.09, 2));
    this.flash.position.copy(bpos);
    this.flash.lookAt(this.camera.position);
    const after = Tm - Math.PI / 2;
    const wave = after >= 0 && after < 0.6 ? 1 - after / 0.6 : 0;
    this.flash.scale.setScalar(0.35 + 7 * Math.max(0, after));
    this.flash.material.opacity = wave * 0.9;

    // HUD
    this.readT.textContent = `${(Tm / Math.PI).toFixed(2)} π`;
    this.readA.textContent = `${aq.toFixed(3)} ℓ`;
    this.readC.textContent = `${ac.toFixed(3)} ℓ`;
    const tag = scripted ? ramp(t, SWEEP / 2 - 0.4, SWEEP / 2) : 1;
    this.bounceTag.style.opacity = String(tag);
    this.bounceTag.classList.toggle('hot', nearB > 0.5);
    const place = (el: HTMLElement, p: THREE.Vector3, dy: number): void => {
      const q = p.clone().add(new THREE.Vector3(0, dy, 0)).project(this.camera);
      el.style.transform = `translate(-50%, 0) translate(${((q.x + 1) / 2) * this.width}px, ${((1 - q.y) / 2) * this.height}px)`;
    };
    place(this.bubbleLabel, bpos, -2.3);
    place(this.ghostLabel, bpos, -2.75);
    this.ghostLabel.style.opacity = String(clamp(1 - ac / 0.4));
    place(this.bounceTag, new THREE.Vector3(0, 0, Z0 - this.aMin * SZ), -0.1);
  }
}

export const bounce: ActModule = {
  async mount(ctx: ActContext) {
    const data = await ctx.data.bounce();
    return new Bounce(ctx, data);
  },
};
