/**
 * Three.js globe for HEALPix maps (loaded lazily by `createSkyView` in sky.ts).
 *
 * The sphere is textured by an equirectangular image of the map; the fragment shader derives
 * (l, b) from the fragment position, so there is no UV seam. Galactic axes map to three.js as
 * (x, y, z)_gal -> (x, z, -y)_three, so Galactic north is "up". The outer view orbits the sphere;
 * the inner view puts the camera at the observer (the centre) and shows the sky as seen from
 * inside. Overlays are 2-px screen-space lines (Line2) and HTML labels (CSS2DRenderer).
 */

import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';
import { Line2 } from 'three/examples/jsm/lines/Line2.js';
import { LineGeometry } from 'three/examples/jsm/lines/LineGeometry.js';
import { LineMaterial } from 'three/examples/jsm/lines/LineMaterial.js';
import { CSS2DObject, CSS2DRenderer } from 'three/examples/jsm/renderers/CSS2DRenderer.js';
import { lbToVec, vecToLb } from './healpix';
import {
  type GlobeView,
  type LB,
  OverlayStore,
  type SkyMapInput,
  type SkyPick,
  SkyChrome,
  type SkyViewOptions,
  equirectTexture,
  formatPick,
  graticule,
  pickValue,
  polygonPoints,
  resolveCssColor,
  smallCirclePoints,
} from './sky';
import { cssVar, currentTheme, onThemeChange } from './theme';
import { h } from './ui';

type V3 = [number, number, number];

const toThree = (v: readonly number[], r = 1): THREE.Vector3 => new THREE.Vector3(v[0]! * r, v[2]! * r, -v[1]! * r);
const fromThree = (p: THREE.Vector3): V3 => [p.x, -p.z, p.y];

const VERT = /* glsl */ `
varying vec3 vPos;
void main() {
  vPos = position;
  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
}`;

const FRAG = /* glsl */ `
uniform sampler2D map;
uniform float hasMap;
uniform vec3 blank;
uniform float opacity;
varying vec3 vPos;
const float PI = 3.141592653589793;
void main() {
  vec3 g = normalize(vec3(vPos.x, -vPos.z, vPos.y));
  float l = atan(g.y, g.x);
  if (l < 0.0) l += 2.0 * PI;
  float b = asin(clamp(g.z, -1.0, 1.0));
  vec3 c = hasMap > 0.5 ? texture2D(map, vec2(l / (2.0 * PI), b / PI + 0.5)).rgb : blank;
  gl_FragColor = vec4(c, opacity);
}`;

const TEX_W = 1024;
const TEX_H = 512;

export function createGlobe(el: HTMLElement, opts: SkyViewOptions = {}): GlobeView {
  // create the WebGL context first: if it fails, nothing is left in `el` for the Mollweide fallback
  const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true, preserveDrawingBuffer: false });
  const root = h('div', { class: 'sky-view' });
  if (opts.height) root.style.height = `${opts.height}px`;
  else root.style.height = '420px';
  el.append(root);

  renderer.setPixelRatio(Math.min(2, window.devicePixelRatio || 1));
  renderer.setClearColor(0x000000, 0);
  root.append(renderer.domElement);
  renderer.domElement.setAttribute('role', 'img');
  renderer.domElement.setAttribute('aria-label', 'Sky globe in Galactic coordinates; drag to rotate');

  const labels = new CSS2DRenderer();
  Object.assign(labels.domElement.style, { position: 'absolute', inset: '0', pointerEvents: 'none' });
  root.append(labels.domElement);
  const chrome = new SkyChrome(root, opts);

  const scene = new THREE.Scene();
  const camera = new THREE.PerspectiveCamera(32, 1, 0.001, 100);
  const controls = new OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true;
  controls.dampingFactor = 0.08;
  controls.enablePan = false;

  let theme = currentTheme();
  let mode: 'outer' | 'inner' = opts.view ?? 'outer';
  let input: SkyMapInput | null = null;
  let opacity = 1;
  let userOpacity: number | null = null;
  let disposed = false;
  let dirty = true;
  let visible = true;

  const texture = new THREE.DataTexture(new Uint8Array(TEX_W * TEX_H * 4), TEX_W, TEX_H, THREE.RGBAFormat);
  texture.wrapS = THREE.RepeatWrapping;
  texture.wrapT = THREE.ClampToEdgeWrapping;
  texture.minFilter = THREE.LinearFilter;
  texture.magFilter = THREE.LinearFilter;
  texture.generateMipmaps = false;

  const blankColor = (): THREE.Color => new THREE.Color().setStyle(theme === 'dark' ? '#141d2e' : '#eef2f7', THREE.NoColorSpace);
  const sphereMat = new THREE.ShaderMaterial({
    vertexShader: VERT,
    fragmentShader: FRAG,
    uniforms: {
      map: { value: texture },
      hasMap: { value: 0 },
      blank: { value: blankColor() },
      opacity: { value: 1 },
    },
    side: THREE.DoubleSide,
  });
  const sphere = new THREE.Mesh(new THREE.SphereGeometry(1, 128, 64), sphereMat);
  scene.add(sphere);

  const gridGroup = new THREE.Group();
  const overlayGroup = new THREE.Group();
  scene.add(gridGroup, overlayGroup);
  const lineMaterials = new Set<LineMaterial>();

  // Overlay edits only mark the overlays stale; the next animation frame rebuilds them once, so
  // adding hundreds of footprints (or moving a slider) costs one rebuild per frame, not one per call.
  let overlaysStale = false;
  const store = new OverlayStore(() => {
    overlaysStale = true;
    dirty = true;
  });
  const pickers = new Set<(p: SkyPick) => void>();

  // --- helpers --------------------------------------------------------------------------------

  const overlayRadius = (): number => (mode === 'outer' ? 1.004 : 0.996);

  function line(points: V3[] | THREE.Vector3[], color: string, width = 2, alpha = 1, r = overlayRadius()): Line2 {
    const pos: number[] = [];
    for (const p of points) {
      const v = p instanceof THREE.Vector3 ? p : toThree(p, r);
      pos.push(v.x, v.y, v.z);
    }
    const geom = new LineGeometry();
    geom.setPositions(pos);
    const mat = new LineMaterial({
      color: new THREE.Color().setStyle(resolveCssColor(color)),
      linewidth: width,
      transparent: alpha < 1,
      opacity: alpha,
      worldUnits: false,
    });
    mat.resolution.set(root.clientWidth || 1, root.clientHeight || 1);
    lineMaterials.add(mat);
    const ln = new Line2(geom, mat);
    ln.computeLineDistances();
    return ln;
  }

  function label(text: string, pos: THREE.Vector3, cls = 'sky-label'): CSS2DObject {
    const obj = new CSS2DObject(h('div', { class: cls }, text));
    obj.position.copy(pos);
    return obj;
  }

  function disposeGroup(g: THREE.Group): void {
    g.traverse((o) => {
      if (o instanceof Line2) {
        o.geometry.dispose();
        lineMaterials.delete(o.material);
        o.material.dispose();
      } else if (o instanceof THREE.Mesh) {
        o.geometry.dispose();
        (o.material as THREE.Material).dispose();
      }
      if (o instanceof CSS2DObject) o.element.remove();
    });
    g.clear();
  }

  function rebuildGrid(): void {
    disposeGroup(gridGroup);
    if (opts.grid === false) return;
    const g = graticule(opts.gridStepDeg ?? 30);
    const c = cssVar('--text-3');
    for (const ln of g.lines) gridGroup.add(line(ln, c, 1, 0.35));
    gridGroup.add(line(g.equator, c, 1.2, 0.7));
    const r = overlayRadius();
    for (const l of [0, 90, 180, 270]) gridGroup.add(label(`${l}°`, toThree(lbToVec(l, 0), r * 1.0), 'sky-label sky-grid-label'));
  }

  function rebuildOverlays(): void {
    overlaysStale = false;
    disposeGroup(overlayGroup);
    const r = overlayRadius();
    let chords = 0;
    for (const o of store.items) {
      if (o.kind === 'circle') {
        overlayGroup.add(line(smallCirclePoints(o.center, o.radiusDeg), o.color, o.width ?? 2));
        if (o.label) overlayGroup.add(label(o.label, toThree(lbToVec(o.center[0], o.center[1]), r)));
      } else if (o.kind === 'polygon') {
        overlayGroup.add(line(polygonPoints(o.vertices), o.color, 2));
        if (o.label) overlayGroup.add(label(o.label, toThree(lbToVec(o.vertices[0]![0], o.vertices[0]![1]), r)));
      } else if (o.kind === 'point') {
        const pos = toThree(lbToVec(o.lb[0], o.lb[1]), r);
        const dot = new THREE.Mesh(
          new THREE.SphereGeometry(0.014, 16, 12),
          new THREE.MeshBasicMaterial({ color: new THREE.Color().setStyle(resolveCssColor(o.color)) }),
        );
        dot.position.copy(pos);
        overlayGroup.add(dot);
        if (o.label) overlayGroup.add(label(o.label, pos));
      } else if (o.kind === 'chord') {
        chords++;
        const a = toThree(lbToVec(o.a[0], o.a[1]), 1);
        const b = toThree(lbToVec(o.b[0], o.b[1]), 1);
        overlayGroup.add(line([a, new THREE.Vector3(0, 0, 0), b], o.color, 2));
        for (const p of [a, b]) {
          const dot = new THREE.Mesh(
            new THREE.SphereGeometry(0.014, 16, 12),
            new THREE.MeshBasicMaterial({ color: new THREE.Color().setStyle(resolveCssColor(o.color)) }),
          );
          dot.position.copy(p.clone().multiplyScalar(r));
          overlayGroup.add(dot);
        }
        if (o.label) overlayGroup.add(label(o.label, a.clone().multiplyScalar(r)));
      }
    }
    if (chords && mode === 'outer') {
      const obs = new THREE.Mesh(
        new THREE.SphereGeometry(0.02, 16, 12),
        new THREE.MeshBasicMaterial({ color: new THREE.Color().setStyle(cssVar('--text-1')) }),
      );
      overlayGroup.add(obs, label('observer', new THREE.Vector3(0, 0, 0)));
    }
    applyOpacity(chords > 0);
    dirty = true;
  }

  function applyOpacity(hasChords: boolean): void {
    opacity = userOpacity ?? (hasChords && mode === 'outer' ? 0.5 : 1);
    sphereMat.uniforms.opacity!.value = opacity;
    sphereMat.transparent = opacity < 1;
    sphereMat.depthWrite = opacity >= 1;
    sphereMat.needsUpdate = true;
  }

  function uploadTexture(): void {
    if (input) {
      const data = equirectTexture(input, theme, TEX_W, TEX_H);
      (texture.image.data as Uint8Array).set(data);
      texture.needsUpdate = true;
    }
    sphereMat.uniforms.hasMap!.value = input ? 1 : 0;
    sphereMat.uniforms.blank!.value = blankColor();
    chrome.setColorbar(input, theme);
    dirty = true;
  }

  function placeCamera(): void {
    if (mode === 'outer') {
      camera.fov = 32;
      camera.position.copy(toThree([1, 0, 0], 4.0));
      controls.target.set(0, 0, 0);
      controls.enableZoom = true;
      controls.minDistance = 1.6;
      controls.maxDistance = 9;
      controls.rotateSpeed = 0.6;
    } else {
      camera.fov = 75;
      // the camera sits (almost) at the observer and looks towards l = 0, b = 0
      camera.position.copy(toThree([-1, 0, 0], 1e-3));
      controls.target.set(0, 0, 0);
      controls.enableZoom = false;
      controls.minDistance = 0;
      controls.maxDistance = Infinity;
      controls.rotateSpeed = -0.3;
    }
    camera.updateProjectionMatrix();
    controls.update();
    dirty = true;
  }

  function resize(): void {
    const w = root.clientWidth || 1;
    const hh = root.clientHeight || 1;
    renderer.setSize(w, hh, false);
    renderer.domElement.style.width = `${w}px`;
    renderer.domElement.style.height = `${hh}px`;
    labels.setSize(w, hh);
    camera.aspect = w / hh;
    camera.updateProjectionMatrix();
    for (const m of lineMaterials) m.resolution.set(w, hh);
    dirty = true;
  }

  // label visibility: hide labels on the far side of the globe
  const tmp = new THREE.Vector3();
  function updateLabels(): void {
    const camDir = camera.getWorldDirection(new THREE.Vector3());
    for (const g of [gridGroup, overlayGroup])
      for (const o of g.children) {
        if (!(o instanceof CSS2DObject)) continue;
        if (o.position.lengthSq() < 1e-9) {
          o.visible = mode === 'outer';
          continue;
        }
        tmp.copy(o.position);
        o.visible = mode === 'outer' ? tmp.sub(camera.position).dot(o.position) < 0 : o.position.dot(camDir) > 0;
      }
  }

  // --- interaction -------------------------------------------------------------------------------

  const raycaster = new THREE.Raycaster();
  const ndc = new THREE.Vector2();
  function pickAt(ev: PointerEvent | MouseEvent): SkyPick | null {
    const r = renderer.domElement.getBoundingClientRect();
    ndc.set(((ev.clientX - r.left) / r.width) * 2 - 1, -((ev.clientY - r.top) / r.height) * 2 + 1);
    raycaster.setFromCamera(ndc, camera);
    const hit = raycaster.intersectObject(sphere, false)[0];
    if (!hit) return null;
    const [l, b] = vecToLb(fromThree(hit.point));
    return pickValue(input, l, b);
  }
  let downAt: [number, number] | null = null;
  renderer.domElement.addEventListener('pointerdown', (ev) => (downAt = [ev.clientX, ev.clientY]));
  renderer.domElement.addEventListener('pointerup', (ev) => {
    if (downAt && Math.hypot(ev.clientX - downAt[0], ev.clientY - downAt[1]) < 4) {
      const p = pickAt(ev);
      if (p) pickers.forEach((cb) => cb(p));
    }
    downAt = null;
  });
  renderer.domElement.addEventListener('pointermove', (ev) => {
    if (ev.buttons) return;
    const p = pickAt(ev);
    chrome.setHud(p ? formatPick(p, input?.unit) : '');
  });
  renderer.domElement.addEventListener('pointerleave', () => chrome.setHud(''));
  renderer.domElement.addEventListener(
    'wheel',
    (ev) => {
      if (mode !== 'inner') return;
      ev.preventDefault();
      camera.fov = Math.min(110, Math.max(15, camera.fov * (1 + Math.sign(ev.deltaY) * 0.08)));
      camera.updateProjectionMatrix();
      dirty = true;
    },
    { passive: false },
  );
  controls.addEventListener('change', () => (dirty = true));

  // --- loop and lifecycle ---------------------------------------------------------------------------

  let frame = 0;
  function loop(): void {
    if (disposed) return;
    frame = requestAnimationFrame(loop);
    if (!root.isConnected) return view.dispose();
    if (!visible) return;
    controls.update();
    if (overlaysStale) rebuildOverlays();
    if (!dirty) return;
    dirty = false;
    updateLabels();
    renderer.render(scene, camera);
    labels.render(scene, camera);
  }

  const ro = new ResizeObserver(resize);
  ro.observe(root);
  const io = new IntersectionObserver((entries) => {
    visible = entries.some((e) => e.isIntersecting);
    dirty = true;
  });
  io.observe(root);
  const offTheme = onThemeChange((t) => {
    theme = t;
    uploadTexture();
    rebuildGrid();
    rebuildOverlays();
  });

  resize();
  placeCamera();
  rebuildGrid();
  uploadTexture();
  loop();

  const view: GlobeView = {
    kind: 'globe',
    el: root,
    ...store.api(),
    setMap(next) {
      input = next;
      uploadTexture();
    },
    onPick(cb) {
      pickers.add(cb);
      return () => pickers.delete(cb);
    },
    resize,
    setView(next) {
      if (next === mode) return;
      mode = next;
      placeCamera();
      rebuildGrid();
      rebuildOverlays();
    },
    setOpacity(alpha) {
      userOpacity = Math.max(0, Math.min(1, alpha));
      applyOpacity(store.items.some((o) => o.kind === 'chord'));
      dirty = true;
    },
    lookAt(lb: LB) {
      const dir = toThree(lbToVec(lb[0], lb[1]), 1);
      if (mode === 'outer') camera.position.copy(dir.multiplyScalar(camera.position.length()));
      else camera.position.copy(dir.multiplyScalar(-1e-3));
      controls.update();
      dirty = true;
    },
    dispose() {
      if (disposed) return;
      disposed = true;
      cancelAnimationFrame(frame);
      ro.disconnect();
      io.disconnect();
      offTheme();
      pickers.clear();
      controls.dispose();
      disposeGroup(gridGroup);
      disposeGroup(overlayGroup);
      sphere.geometry.dispose();
      sphereMat.dispose();
      texture.dispose();
      renderer.dispose();
      renderer.forceContextLoss();
      root.remove();
    },
  };
  return view;
}
