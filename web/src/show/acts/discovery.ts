/**
 * Act I · Discovery — "Space may be a crystal" (docs/SHOW.md).
 *
 * Inside the real WMAP sky; pull back through the last-scattering sphere; copies of the same sky
 * ignite at FCC lattice translations and their exact intersection circles light up; the
 * fundamental domain (rhombic dodecahedron) and its honeycomb appear; fly back to one globe whose
 * matched circles shrink and vanish as d grows from 0.8 D to 1.05 D.
 *
 * Units: 10 world units per D, so the last-scattering sphere has radius chi_* = 5.
 * Geometry is exact: lattice vectors from lib/physics/lattice, circles from lib/physics/circles,
 * orientation = the paper's best-fit (not significant) FCC axes (topology.json).
 */

import * as THREE from 'three';
import { LineMaterial } from 'three/examples/jsm/lines/LineMaterial.js';
import { LineSegments2 } from 'three/examples/jsm/lines/LineSegments2.js';
import { LineSegmentsGeometry } from 'three/examples/jsm/lines/LineSegmentsGeometry.js';
import { circlePairs } from '../../lib/physics/circles';
import { cubicBasis, latticeVectors } from '../../lib/physics/lattice';
import { type Mat3, type Vec3, cross, identity3, matVec, scale, unit, vecMat } from '../../lib/physics/vec';
import { CIRCLE_SEARCH_BOUND } from '../../lib/physics/constants';
import { h } from '../../lib/ui';
import { Act, type ActContext, type ActModule, DragOrbit } from '../act';
import { syncFill } from '../ui';
import { POINTS_FRAG, cmbValueTexture, easeInOut, galToThree, lerp, pointScale, ramp, rng, skyDome, window4 } from '../gfx';
import { CMB_COLOR_GLSL, GALACTIC_GLSL } from '../shaders/common';

const S = 10; // world units per D
const R_LSS = 0.5 * S; // chi_* = D/2
const D_CRYSTAL = 0.8; // illustrative size (excluded by data)
const MAX_GLOBE_CIRCLES = 96;

const SPHERE_VERT = /* glsl */ `
#ifdef USE_INSTANCING
attribute float aIgnite;
attribute float aDim;
#endif
uniform float uTime;
uniform float uIgnite;
varying vec3 vLocal;
varying vec3 vN;
varying vec3 vView;
varying float vI;
varying float vFlash;
void main() {
  vLocal = position;
#ifdef USE_INSTANCING
  mat4 m = modelMatrix * instanceMatrix;
  float ign = aIgnite;
  float dim = aDim;
#else
  mat4 m = modelMatrix;
  float ign = uIgnite;
  float dim = 1.0;
#endif
  vec4 wp = m * vec4(position, 1.0);
  vN = normalize(mat3(m) * normal);
  vView = cameraPosition - wp.xyz;
  float x = clamp((uTime - ign) / 0.9, 0.0, 1.0);
  vI = x * x * (3.0 - 2.0 * x) * dim;
  vFlash = uTime > ign ? exp(-(uTime - ign) * 3.0) : 0.0;
  gl_Position = projectionMatrix * viewMatrix * wp;
}`;

const SPHERE_FRAG = /* glsl */ `
${GALACTIC_GLSL}
${CMB_COLOR_GLSL}
uniform sampler2D uMap;
uniform float uGain;
uniform float uFace;
uniform float uRimC;
uniform float uRim;
uniform float uBack;
uniform float uTime;
varying vec3 vLocal;
varying vec3 vN;
varying vec3 vView;
varying float vI;
varying float vFlash;
void main() {
  vec4 tx = texture2D(uMap, galUv(threeToGal(vLocal)));
  float v = (tx.r - 0.5) * 2.0;
  vec3 c = cmbColor(v * 1.55);
  float f = pow(1.0 - abs(dot(normalize(vN), normalize(vView))), 2.6);
  vec3 col = c * (uFace + uRimC * f) * uGain;
  col += vec3(0.30, 0.85, 1.0) * f * (uRim + 1.5 * vFlash);
  col += c * vFlash * 0.6;
  if (!gl_FrontFacing) col *= uBack;
  gl_FragColor = vec4(col * vI, 1.0);
}`;

const RING_VERT = /* glsl */ `
attribute float aIgnite;
attribute float aGain;
uniform float uTime;
varying float vI;
varying float vU;
varying float vFlash;
void main() {
  vU = uv.x;
  float x = clamp((uTime - aIgnite) / 0.7, 0.0, 1.0);
  vI = x * x * (3.0 - 2.0 * x) * aGain;
  vFlash = uTime > aIgnite ? exp(-(uTime - aIgnite) * 2.5) : 0.0;
  gl_Position = projectionMatrix * viewMatrix * modelMatrix * instanceMatrix * vec4(position, 1.0);
}`;

const RING_FRAG = /* glsl */ `
uniform float uTime;
uniform vec3 uColor;
uniform float uOpacity;
varying float vI;
varying float vU;
varying float vFlash;
void main() {
  float flow = 0.7 + 0.3 * sin(vU * 6.2831853 * 6.0 - uTime * 3.0);
  vec3 col = uColor * (1.0 * flow + 4.0 * vFlash);
  gl_FragColor = vec4(col * vI * uOpacity, 1.0);
}`;

const STREAM_VERT = /* glsl */ `
${GALACTIC_GLSL}
${CMB_COLOR_GLSL}
attribute vec3 aDir;
attribute float aPhase;
attribute float aSize;
uniform float uTime;
uniform float uScale;
uniform float uR;
uniform sampler2D uMap;
uniform float uOpacity;
varying vec3 vColor;
varying float vAlpha;
void main() {
  float f = fract(aPhase + uTime * 0.11);
  float r = uR * (1.0 - f) + 0.35;
  vec3 p = aDir * r;
  float v = (texture2D(uMap, galUv(threeToGal(aDir))).r - 0.5) * 2.0;
  vColor = cmbColor(v * 1.6) * 1.3 + vec3(0.05);
  vAlpha = smoothstep(0.0, 0.12, f) * smoothstep(0.7, 2.0, r) * uOpacity;
  vec4 mv = modelViewMatrix * vec4(p, 1.0);
  gl_PointSize = min(aSize * uScale / max(0.05, -mv.z), 7.0);
  gl_Position = projectionMatrix * mv;
}`;

interface Pair {
  i: number;
  j: number;
}

function lookupR(topoR: number[][] | undefined): Mat3 {
  if (topoR && topoR.length === 3 && topoR.every((r) => r.length === 3 && r.every(Number.isFinite))) return topoR.map((r) => [r[0]!, r[1]!, r[2]!]) as Mat3;
  return identity3();
}

/** Rotation taking +Z to the unit vector n. */
function quatZTo(n: THREE.Vector3): THREE.Quaternion {
  return new THREE.Quaternion().setFromUnitVectors(new THREE.Vector3(0, 0, 1), n.clone().normalize());
}

/** Edges of the rhombic dodecahedron (FCC Voronoi cell) for NN distance d, lattice frame. */
function rhombicDodecahedronEdges(d: number): [Vec3, Vec3][] {
  const c = d * Math.SQRT2;
  const edges: [Vec3, Vec3][] = [];
  for (const sx of [-1, 1])
    for (const sy of [-1, 1])
      for (const sz of [-1, 1]) {
        const v: Vec3 = [(sx * c) / 4, (sy * c) / 4, (sz * c) / 4];
        edges.push([v, [(sx * c) / 2, 0, 0]], [v, [0, (sy * c) / 2, 0]], [v, [0, 0, (sz * c) / 2]]);
      }
  return edges;
}

class Discovery extends Act {
  readonly id = 'discovery' as const;
  readonly duration = 19;
  private readonly R: Mat3;
  private readonly basis: Mat3;
  private readonly sphereMat: THREE.ShaderMaterial;
  private readonly copyMat: THREE.ShaderMaterial;
  private readonly ringMat: THREE.ShaderMaterial;
  private readonly globeRingMat: THREE.ShaderMaterial;
  private readonly streamMat: THREE.ShaderMaterial;
  private readonly sky: THREE.Mesh<THREE.SphereGeometry, THREE.ShaderMaterial>;
  private readonly central: THREE.Mesh;
  private readonly copies: THREE.InstancedMesh;
  private readonly globeRings: THREE.InstancedMesh;
  private readonly cellMat: LineMaterial;
  private readonly combMat: LineMaterial;
  private readonly chordMat: LineMaterial;
  private readonly chordGeo: LineSegmentsGeometry;
  private readonly twinPoints: THREE.Points;
  private readonly twinPos: Float32Array;
  private readonly orbit: DragOrbit;
  private readonly hudD: HTMLElement;
  private readonly hudDVal: HTMLElement;
  private readonly hudDNote: HTMLElement;
  private readonly counter: HTMLElement;
  private readonly counterNum: HTMLElement;
  private readonly slider: HTMLInputElement;
  private readonly freeCtl: HTMLElement;
  private userD: number | null = null;

  constructor(ctx: ActContext, tex: THREE.Texture, topoR: number[][] | undefined) {
    super(ctx);
    this.R = lookupR(topoR);
    this.basis = cubicBasis('FCC');
    this.camera.fov = 75;
    this.bloom = { strength: 0.95, radius: 0.55, threshold: 0.55 };

    this.sky = skyDome(1500, 0);
    this.scene.add(this.sky);

    const common = {
      uMap: { value: tex },
      uTime: { value: 0 },
      uGain: { value: 1 },
      uFace: { value: 0.35 },
      uRimC: { value: 1 },
      uRim: { value: 0.4 },
      uBack: { value: 1 },
      uIgnite: { value: -10 },
    };
    this.sphereMat = new THREE.ShaderMaterial({
      uniforms: THREE.UniformsUtils.clone(common),
      vertexShader: SPHERE_VERT,
      fragmentShader: SPHERE_FRAG,
      side: THREE.DoubleSide,
      transparent: true,
      blending: THREE.AdditiveBlending,
      depthWrite: false,
    });
    this.sphereMat.uniforms.uMap!.value = tex;
    this.central = new THREE.Mesh(new THREE.SphereGeometry(R_LSS, 160, 80), this.sphereMat);
    this.scene.add(this.central);

    // --- copies at FCC translations (d = 0.8 D, illustrative) --------------------------------
    const B = scaledFcc(this.basis, D_CRYSTAL);
    const vecs = latticeVectors(B, D_CRYSTAL * Math.SQRT2 * 1.001);
    const centres: THREE.Vector3[] = [new THREE.Vector3()];
    const shells: number[] = [0];
    for (const v of vecs) {
      centres.push(galToThree(matVec(this.R, v.g), S));
      shells.push(Math.round((v.length / D_CRYSTAL) ** 2));
    }
    const rand = rng(7);
    const ignite: number[] = shells.map((sh) => (sh === 0 ? -10 : 5.6 + 0.55 * (sh - 1) + 0.5 * rand()));
    const nCopies = centres.length - 1;
    this.copyMat = new THREE.ShaderMaterial({
      uniforms: THREE.UniformsUtils.clone(common),
      vertexShader: SPHERE_VERT,
      fragmentShader: SPHERE_FRAG,
      side: THREE.DoubleSide,
      transparent: true,
      blending: THREE.AdditiveBlending,
      depthWrite: false,
    });
    this.copyMat.uniforms.uMap!.value = tex;
    const copyGeo = new THREE.SphereGeometry(R_LSS, 72, 36);
    const aIgn = new Float32Array(nCopies);
    const aDim = new Float32Array(nCopies);
    this.copies = new THREE.InstancedMesh(copyGeo, this.copyMat, nCopies);
    const m = new THREE.Matrix4();
    for (let k = 0; k < nCopies; k++) {
      m.makeTranslation(centres[k + 1]!);
      this.copies.setMatrixAt(k, m);
      aIgn[k] = ignite[k + 1]!;
      aDim[k] = shells[k + 1] === 1 ? 1 : 0.55;
    }
    copyGeo.setAttribute('aIgnite', new THREE.InstancedBufferAttribute(aIgn, 1));
    copyGeo.setAttribute('aDim', new THREE.InstancedBufferAttribute(aDim, 1));
    this.copies.frustumCulled = false;
    this.scene.add(this.copies);

    // --- exact intersection circles between neighbouring copies (not the central sphere) ------
    const pairs: Pair[] = [];
    for (let i = 1; i < centres.length; i++)
      for (let j = i + 1; j < centres.length; j++) if (centres[i]!.distanceTo(centres[j]!) < 2 * R_LSS * 0.999) pairs.push({ i, j });
    const ringGeo = new THREE.TorusGeometry(1, 0.011, 6, 160);
    const ringIgn = new Float32Array(pairs.length);
    const ringGain = new Float32Array(pairs.length);
    this.ringMat = ringMaterial(0.35);
    const rings = new THREE.InstancedMesh(ringGeo, this.ringMat, pairs.length);
    const q = new THREE.Quaternion();
    const sc = new THREE.Vector3();
    pairs.forEach((p, k) => {
      const a = centres[p.i]!;
      const b = centres[p.j]!;
      const gvec = b.clone().sub(a);
      const half = gvec.length() / 2;
      const r = Math.sqrt(R_LSS * R_LSS - half * half);
      q.copy(quatZTo(gvec));
      sc.set(r, r, r);
      m.compose(a.clone().add(b).multiplyScalar(0.5), q, sc);
      rings.setMatrixAt(k, m);
      const mid = a.clone().add(b).multiplyScalar(0.5).length() / S;
      ringIgn[k] = Math.max(ignite[p.i]!, ignite[p.j]!) + 1.4 + 0.9 * mid;
      ringGain[k] = Math.max(0.2, 0.85 - 0.5 * mid);
    });
    ringGeo.setAttribute('aIgnite', new THREE.InstancedBufferAttribute(ringIgn, 1));
    ringGeo.setAttribute('aGain', new THREE.InstancedBufferAttribute(ringGain, 1));
    rings.frustumCulled = false;
    this.scene.add(rings);

    // --- circles on the central globe: exact matched pairs for the current d ------------------
    const gGeo = new THREE.TorusGeometry(1, 0.012, 6, 200);
    this.globeRingMat = ringMaterial(1);
    this.globeRings = new THREE.InstancedMesh(gGeo, this.globeRingMat, MAX_GLOBE_CIRCLES);
    gGeo.setAttribute('aIgnite', new THREE.InstancedBufferAttribute(new Float32Array(MAX_GLOBE_CIRCLES).fill(7.9), 1));
    gGeo.setAttribute('aGain', new THREE.InstancedBufferAttribute(new Float32Array(MAX_GLOBE_CIRCLES).fill(1.25), 1));
    this.globeRings.frustumCulled = false;
    this.scene.add(this.globeRings);

    // --- fundamental domain and honeycomb -------------------------------------------------------
    const edges = rhombicDodecahedronEdges(D_CRYSTAL);
    const toWorld = (v: Vec3, off: THREE.Vector3): number[] => {
      const p = galToThree(matVec(this.R, v), S).add(off);
      return [p.x, p.y, p.z];
    };
    const cellPos: number[] = [];
    for (const [a, b] of edges) cellPos.push(...toWorld(a, new THREE.Vector3()), ...toWorld(b, new THREE.Vector3()));
    const combPos: number[] = [];
    for (let k = 1; k < centres.length; k++) for (const [a, b] of edges) combPos.push(...toWorld(a, centres[k]!), ...toWorld(b, centres[k]!));
    this.cellMat = new LineMaterial({ color: 0xfff1d6, linewidth: 3.6, transparent: true, opacity: 0, blending: THREE.AdditiveBlending, depthWrite: false });
    this.combMat = new LineMaterial({ color: 0xffc35a, linewidth: 1.4, transparent: true, opacity: 0, blending: THREE.AdditiveBlending, depthWrite: false });
    const cellGeo = new LineSegmentsGeometry().setPositions(cellPos);
    const combGeo = new LineSegmentsGeometry().setPositions(combPos);
    const cell = new LineSegments2(cellGeo, this.cellMat);
    const comb = new LineSegments2(combGeo, this.combMat);
    cell.frustumCulled = comb.frustumCulled = false;
    this.scene.add(cell, comb);

    // --- light streaming from the last-scattering surface to the observer ----------------------
    const nStream = 7000;
    const sDir = new Float32Array(nStream * 3);
    const sPhase = new Float32Array(nStream);
    const sSize = new Float32Array(nStream);
    const r2 = rng(3);
    for (let k = 0; k < nStream; k++) {
      const z = 2 * r2() - 1;
      const ph = 2 * Math.PI * r2();
      const s = Math.sqrt(1 - z * z);
      sDir.set([s * Math.cos(ph), z, s * Math.sin(ph)], 3 * k);
      sPhase[k] = r2();
      sSize[k] = 0.018 + 0.05 * r2() ** 3;
    }
    const sGeo = new THREE.BufferGeometry();
    sGeo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(nStream * 3), 3));
    sGeo.setAttribute('aDir', new THREE.BufferAttribute(sDir, 3));
    sGeo.setAttribute('aPhase', new THREE.BufferAttribute(sPhase, 1));
    sGeo.setAttribute('aSize', new THREE.BufferAttribute(sSize, 1));
    this.streamMat = new THREE.ShaderMaterial({
      uniforms: { uTime: { value: 0 }, uScale: { value: 800 }, uR: { value: R_LSS * 0.985 }, uMap: { value: tex }, uOpacity: { value: 1 } },
      vertexShader: STREAM_VERT,
      fragmentShader: POINTS_FRAG,
      transparent: true,
      blending: THREE.AdditiveBlending,
      depthWrite: false,
    });
    const stream = new THREE.Points(sGeo, this.streamMat);
    stream.frustumCulled = false;
    this.scene.add(stream);

    // --- twin points and chords through the observer --------------------------------------------
    this.twinPos = new Float32Array(MAX_GLOBE_CIRCLES * 3);
    const tGeo = new THREE.BufferGeometry();
    tGeo.setAttribute('position', new THREE.BufferAttribute(this.twinPos, 3));
    tGeo.setAttribute('aSize', new THREE.BufferAttribute(new Float32Array(MAX_GLOBE_CIRCLES).fill(0.55), 1));
    tGeo.setAttribute('aColor', new THREE.BufferAttribute(new Float32Array(MAX_GLOBE_CIRCLES * 3).fill(1.6), 3));
    tGeo.setAttribute('aAlpha', new THREE.BufferAttribute(new Float32Array(MAX_GLOBE_CIRCLES).fill(1), 1));
    this.twinPoints = new THREE.Points(
      tGeo,
      new THREE.ShaderMaterial({
        uniforms: { uScale: { value: 800 } },
        vertexShader: /* glsl */ `
          attribute float aSize; attribute vec3 aColor; attribute float aAlpha;
          uniform float uScale; varying vec3 vColor; varying float vAlpha;
          void main() { vColor = aColor; vAlpha = aAlpha; vec4 mv = modelViewMatrix * vec4(position, 1.0);
            gl_PointSize = aSize * uScale / max(0.05, -mv.z); gl_Position = projectionMatrix * mv; }`,
        fragmentShader: POINTS_FRAG,
        transparent: true,
        blending: THREE.AdditiveBlending,
        depthWrite: false,
      }),
    );
    this.twinPoints.frustumCulled = false;
    this.scene.add(this.twinPoints);
    this.chordGeo = new LineSegmentsGeometry();
    this.chordMat = new LineMaterial({ color: 0xbff6ff, linewidth: 1.6, transparent: true, opacity: 0, blending: THREE.AdditiveBlending, depthWrite: false });
    const chords = new LineSegments2(this.chordGeo, this.chordMat);
    chords.frustumCulled = false;
    this.scene.add(chords);

    // observer
    const you = new THREE.Points(
      new THREE.BufferGeometry()
        .setAttribute('position', new THREE.BufferAttribute(new Float32Array(3), 3))
        .setAttribute('aSize', new THREE.BufferAttribute(new Float32Array([0.5]), 1))
        .setAttribute('aColor', new THREE.BufferAttribute(new Float32Array([2.2, 2.0, 1.6]), 3))
        .setAttribute('aAlpha', new THREE.BufferAttribute(new Float32Array([1]), 1)),
      this.twinPoints.material,
    );
    you.frustumCulled = false;
    this.scene.add(you);

    // --- HUD -------------------------------------------------------------------------------------
    this.counterNum = h('div', { class: 'counter-num' }, '18');
    this.counter = h(
      'div',
      { class: 'show-counter' },
      this.counterNum,
      h('div', { class: 'counter-label' }, 'flat shapes of space', h('br'), 'consistent with P1–P3'),
    );
    this.hudDVal = h('span', { class: 'mono' }, '0.80');
    this.hudDNote = h('span', { class: 'dnote' }, '');
    this.hudD = h('div', { class: 'show-dread' }, h('span', { class: 'dk' }, 'universe size'), h('span', null, 'd = ', this.hudDVal, ' D'), this.hudDNote);
    this.slider = h('input', { type: 'range', min: '0.6', max: '1.15', step: '0.005', value: '1.05', 'aria-label': 'Universe size d over D' }) as HTMLInputElement;
    this.freeCtl = h(
      'div',
      { class: 'show-dock show-ui' },
      h('label', { class: 'dock-field' }, h('span', { class: 'dock-k' }, 'Universe size d/D'), h('span', { class: 'dock-scale' }, this.slider, h('span', { class: 'dock-excl', style: { width: `${((0.985 - 0.6) / 0.55) * 100}%` } }))),
      h('span', { class: 'dock-hint' }, 'drag to orbit'),
    );
    ctx.hud.append(this.counter, this.hudD, this.freeCtl);
    this.listen(this.slider, 'input', () => {
      this.userD = Number(this.slider.value);
      ctx.userInput();
    });

    this.orbit = new DragOrbit(ctx.stage.canvas.parentElement ?? ctx.stage.canvas, () => ctx.userInput());
    this.own(() => this.orbit.dispose());
    this.update(0);
  }

  /** universe size at act time t (scripted) or set by the judge */
  private dAt(t: number): number {
    if (this.userD !== null && !this.isPlaying) return this.userD;
    return lerp(D_CRYSTAL, 1.05, easeInOut((t - 15.4) / 3.0));
  }

  protected update(dt: number): void {
    const t = this.time;
    if (this.isPlaying) this.userD = null;
    this.orbit.step(dt);
    const idle = this.isPlaying ? 0 : this.wall;

    // --- camera: yaw/pitch/dist, looking along dir (inside) or at the origin (outside) ---------
    const free = !this.isPlaying && t >= this.duration - 1e-3;
    const yaw = 0.6 + 0.05 * t + 0.9 * easeInOut((t - 6) / 10) + this.orbit.yaw + (free ? 0.03 * idle : 0);
    const pitch = Math.max(-1.2, Math.min(1.2, 0.08 + 0.32 * ramp(t, 4.6, 9) - 0.18 * ramp(t, 14.4, 17) + this.orbit.pitch));
    let dist = 46 * easeInOut((t - 4.6) / 3.8);
    if (t > 8.4) dist = lerp(46, 38, ramp(t, 8.4, 14.4));
    if (t > 14.4) dist = lerp(38, 17.5, easeInOut((t - 14.4) / 2.4));
    const dir = new THREE.Vector3(Math.cos(pitch) * Math.sin(yaw), Math.sin(pitch), Math.cos(pitch) * Math.cos(yaw));
    this.camera.position.copy(dir).multiplyScalar(-dist);
    this.camera.lookAt(this.camera.position.clone().add(dir));
    const fov = lerp(75, 48, ramp(t, 4.6, 8.4));
    if (Math.abs(this.camera.fov - fov) > 1e-3) {
      this.camera.fov = fov;
      this.camera.updateProjectionMatrix();
    }
    const scaleNow = pointScale(this.camera, this.height, this.pixelRatio);
    this.bloom.strength = lerp(0.55, 0.95, outside(dist));

    // --- sky, sphere and copies -------------------------------------------------------------------
    const outsideF = outside(dist);
    this.sky.material.uniforms.uIntensity!.value = outsideF * 0.85;
    const su = this.sphereMat.uniforms;
    su.uTime!.value = t;
    const globe = ramp(t, 15, 16.5);
    su.uGain!.value = (1 - 0.35 * ramp(t, 11, 13) * (1 - globe));
    su.uFace!.value = lerp(0.95, lerp(0.06, 0.62, globe), outsideF);
    su.uRimC!.value = lerp(0, lerp(0.9, 0.45, globe), outsideF);
    su.uRim!.value = lerp(0, lerp(0.18, 0.08, globe), outsideF);
    su.uBack!.value = lerp(1, lerp(0.5, 0.18, globe), outsideF);
    const cu = this.copyMat.uniforms;
    cu.uTime!.value = t;
    cu.uGain!.value = (1 - 0.6 * ramp(t, 11, 12.5)) * (1 - ramp(t, 14.4, 15.6));
    cu.uFace!.value = 0.025;
    cu.uRimC!.value = 0.55;
    cu.uRim!.value = 0.1;
    cu.uBack!.value = 0.45;
    this.copies.visible = t > 5.5 && t < 15.7;

    this.ringMat.uniforms.uTime!.value = t;
    this.ringMat.uniforms.uOpacity!.value = 1 - ramp(t, 14.2, 15.4);
    this.globeRingMat.uniforms.uTime!.value = t;
    this.globeRingMat.uniforms.uOpacity!.value = lerp(0.55, 0.62, ramp(t, 14.4, 15.6));

    this.streamMat.uniforms.uTime!.value = t + idle;
    this.streamMat.uniforms.uOpacity!.value = 1 - ramp(t, 5.0, 6.4);
    this.streamMat.uniforms.uScale!.value = scaleNow;
    (this.twinPoints.material as THREE.ShaderMaterial).uniforms.uScale!.value = scaleNow;

    this.cellMat.opacity = ramp(t, 11, 12) * (1 - ramp(t, 15.2, 16.4));
    this.combMat.opacity = 0.55 * ramp(t, 11.3, 12.6) * (1 - ramp(t, 14.2, 15.2));
    for (const mat of [this.cellMat, this.combMat, this.chordMat]) mat.resolution.set(this.width, this.height);

    // --- globe circles for the current d --------------------------------------------------------
    const d = this.dAt(t);
    this.updateGlobe(d, t + idle, window4(t, 15.2, 16.0, 99, 100));

    // --- HUD ----------------------------------------------------------------------------------------
    const cOn = window4(t, 11.2, 11.7, 14.2, 14.7);
    this.counter.style.opacity = String(cOn);
    this.counter.style.visibility = cOn > 0.01 ? 'visible' : 'hidden';
    const n = Math.round(lerp(18, 1, easeInOut((t - 11.6) / 1.8)));
    if (this.counterNum.textContent !== String(n)) this.counterNum.textContent = String(n);
    this.counter.classList.toggle('done', n === 1);
    const dOn = ramp(t, 15.0, 15.6);
    this.hudD.style.opacity = String(dOn);
    this.hudD.style.visibility = dOn > 0.01 ? 'visible' : 'hidden';
    this.hudDVal.textContent = d.toFixed(2);
    const excluded = d < CIRCLE_SEARCH_BOUND.value;
    this.hudDNote.textContent = excluded ? 'illustrative size · excluded by data' : d < 1 ? 'allowed by data' : 'circles vanish at d ≥ D';
    this.hudD.classList.toggle('excluded', excluded);
    this.freeCtl.classList.toggle('on', free);
    if (!free) this.slider.value = d.toFixed(3);
    syncFill(this.slider);
  }

  private updateGlobe(d: number, t: number, twinsOn: number): void {
    const pairs = d < 1 ? circlePairs(this.basis, d, this.R) : [];
    const m = new THREE.Matrix4();
    const q = new THREE.Quaternion();
    const sc = new THREE.Vector3();
    let k = 0;
    const chord: number[] = [];
    let tw = 0;
    pairs.forEach((p, pi) => {
      const ca = p.length_over_D; // cos(alpha) = |g| / (2 chi_*), chi_* = 1/2
      const sa = Math.sqrt(Math.max(0, 1 - ca * ca));
      for (const c of [p.centre_minus, p.centre_plus]) {
        if (k >= MAX_GLOBE_CIRCLES) break;
        const n = galToThree(c);
        q.copy(quatZTo(n));
        sc.setScalar(R_LSS * sa);
        m.compose(n.clone().multiplyScalar(R_LSS * ca), q, sc);
        this.globeRings.setMatrixAt(k++, m);
      }
      // a twin pair travelling round the matched circles (same azimuth on both)
      if (twinsOn > 0 && tw + 2 <= MAX_GLOBE_CIRCLES) {
        const u = unit(p.g_sky);
        const e1 = unit(cross(u, Math.abs(u[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0]));
        const e2 = cross(u, e1);
        const phi = 0.9 * t + pi * 1.7;
        const ring = scale(e1, sa * Math.cos(phi)).map((x, i) => x + sa * Math.sin(phi) * e2[i]!) as Vec3;
        const a = galToThree(ring.map((x, i) => x - ca * u[i]!), R_LSS);
        const b = galToThree(ring.map((x, i) => x + ca * u[i]!), R_LSS);
        this.twinPos.set([a.x, a.y, a.z], 3 * tw++);
        this.twinPos.set([b.x, b.y, b.z], 3 * tw++);
        chord.push(a.x, a.y, a.z, 0, 0, 0, 0, 0, 0, b.x, b.y, b.z);
      }
    });
    this.globeRings.count = k;
    this.globeRings.instanceMatrix.needsUpdate = true;
    const tg = this.twinPoints.geometry;
    tg.setDrawRange(0, twinsOn > 0 ? tw : 0);
    (tg.getAttribute('position') as THREE.BufferAttribute).needsUpdate = true;
    const alpha = tg.getAttribute('aAlpha') as THREE.BufferAttribute;
    for (let i = 0; i < tw; i++) alpha.setX(i, twinsOn * (0.75 + 0.25 * Math.sin(6 * t + i)));
    alpha.needsUpdate = true;
    if (chord.length) this.chordGeo.setPositions(chord);
    this.chordMat.opacity = chord.length ? 0.32 * twinsOn : 0;
  }

  override dispose(): void {
    this.chordGeo.dispose();
    super.dispose();
  }
}

/** 0 inside the last-scattering sphere, 1 once the camera is well outside it */
const outside = (dist: number): number => ramp(dist, R_LSS * 0.85, R_LSS * 1.6);

function ringMaterial(opacity: number): THREE.ShaderMaterial {
  return new THREE.ShaderMaterial({
    uniforms: { uTime: { value: 0 }, uColor: { value: new THREE.Color(0.35, 0.95, 1.0) }, uOpacity: { value: opacity } },
    vertexShader: RING_VERT,
    fragmentShader: RING_FRAG,
    transparent: true,
    blending: THREE.AdditiveBlending,
    depthWrite: false,
  });
}

/** FCC row basis scaled so its shortest vector has length d (D units). */
function scaledFcc(basis: Mat3, d: number): Mat3 {
  const s = d / Math.SQRT1_2;
  return basis.map((r) => vecMat(r, [[s, 0, 0], [0, s, 0], [0, 0, s]])) as Mat3;
}

export const discovery: ActModule = {
  async mount(ctx) {
    const [sky, topo] = await Promise.all([ctx.data.sky(), ctx.data.topology()]);
    const tex = sky ? cmbValueTexture(sky) : fallbackTexture();
    return new Discovery(ctx, tex, topo?.circles?.fcc_preregistered?.R as number[][] | undefined);
  },
};

/** If the map bundle is missing: a flat mid-grey value texture (the act says so in the HUD). */
function fallbackTexture(): THREE.DataTexture {
  const data = new Uint8Array(4 * 4 * 4).fill(128);
  const tex = new THREE.DataTexture(data, 4, 4);
  tex.needsUpdate = true;
  return tex;
}
