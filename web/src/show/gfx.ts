/** Small graphics helpers shared by the acts: easing, frames, sky dome, CMB texture, glow sprites. */

import * as THREE from 'three';
import type { SkyMapData } from '../lib/data';
import { lb2pix } from '../lib/healpix';
import { NOISE_GLSL, SKY_GLSL } from './shaders/common';

export const clamp = (x: number, a = 0, b = 1): number => Math.min(b, Math.max(a, x));
export const lerp = (a: number, b: number, f: number): number => a + (b - a) * f;
/** 0 before a, 1 after b, smooth in between. */
export const ramp = (t: number, a: number, b: number): number => {
  const x = clamp((t - a) / (b - a));
  return x * x * (3 - 2 * x);
};
export const easeInOut = (x: number): number => {
  const u = clamp(x);
  return u < 0.5 ? 4 * u * u * u : 1 - (-2 * u + 2) ** 3 / 2;
};
/** a bump that rises over [a, b] and falls over [c, d] */
export const window4 = (t: number, a: number, b: number, c: number, d: number): number => ramp(t, a, b) * (1 - ramp(t, c, d));

/** Galactic (x, y, z) -> three.js (x, z, -y). */
export const galToThree = (v: readonly number[], r = 1): THREE.Vector3 => new THREE.Vector3(v[0]! * r, v[2]! * r, -v[1]! * r);

/** A deterministic PRNG (mulberry32) so every frame of the show is reproducible. */
export function rng(seed: number): () => number {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/** Background sphere with the procedural deep sky (stars + dusty band). */
export function skyDome(radius = 900, intensity = 1): THREE.Mesh<THREE.SphereGeometry, THREE.ShaderMaterial> {
  const mat = new THREE.ShaderMaterial({
    side: THREE.BackSide,
    depthWrite: false,
    uniforms: { uIntensity: { value: intensity }, uBand: { value: new THREE.Vector3(0.25, 1, 0.35).normalize() } },
    vertexShader: /* glsl */ `
      varying vec3 vDir;
      void main() { vDir = position; gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0); }`,
    fragmentShader: /* glsl */ `
      ${NOISE_GLSL}
      ${SKY_GLSL}
      uniform float uIntensity; uniform vec3 uBand;
      varying vec3 vDir;
      void main() { gl_FragColor = vec4(deepSky(vDir, uBand) * uIntensity, 1.0); }`,
  });
  const mesh = new THREE.Mesh(new THREE.SphereGeometry(radius, 48, 24), mat);
  mesh.renderOrder = -10;
  mesh.frustumCulled = false;
  return mesh;
}

/**
 * Equirectangular "value" texture of a HEALPix map for the show's own colormap: R = value /
 * range mapped to [0, 1] (0.5 = zero), G = 1 outside the mask. Column i <-> l, row 0 <-> b = -90.
 */
export function cmbValueTexture(sky: SkyMapData, width = 1024, height = 512, smoothDeg = 1.1): THREE.DataTexture {
  const { map, mask, meta, entry } = sky;
  const range = entry.range?.[1] ?? Math.max(Math.abs(entry.min), Math.abs(entry.max));
  const nside = meta.nside;
  const val = new Float32Array(width * height);
  const msk = new Uint8Array(width * height);
  for (let j = 0; j < height; j++) {
    const b = -90 + ((j + 0.5) * 180) / height;
    for (let i = 0; i < width; i++) {
      const l = ((i + 0.5) * 360) / width;
      const p = lb2pix(nside, l, b);
      val[j * width + i] = map[p] ?? 0;
      msk[j * width + i] = mask ? (mask[p] ? 255 : 0) : 255;
    }
  }
  // the low-l map is band-limited: a ~1 deg Gaussian removes the HEALPix pixel pattern
  smoothSphere(val, width, height, smoothDeg);
  const data = new Uint8Array(width * height * 4);
  for (let k = 0; k < width * height; k++) {
    data[4 * k] = Math.round(255 * clamp(0.5 + (0.5 * val[k]!) / range));
    data[4 * k + 1] = msk[k]!;
    data[4 * k + 3] = 255;
  }
  const tex = new THREE.DataTexture(data, width, height, THREE.RGBAFormat, THREE.UnsignedByteType);
  tex.wrapS = THREE.RepeatWrapping;
  tex.wrapT = THREE.ClampToEdgeWrapping;
  tex.magFilter = THREE.LinearFilter;
  tex.minFilter = THREE.LinearFilter;
  tex.needsUpdate = true;
  return tex;
}

/** Approximate Gaussian smoothing of an equirectangular map (3 box passes in longitude, widened by 1/cos b; one in latitude). */
function smoothSphere(v: Float32Array, w: number, hgt: number, sigmaDeg: number): void {
  const texDeg = 180 / hgt;
  const row = new Float32Array(w);
  for (let j = 0; j < hgt; j++) {
    const b = (-90 + ((j + 0.5) * 180) / hgt) * (Math.PI / 180);
    const sig = Math.min(w / 4, sigmaDeg / texDeg / Math.max(0.02, Math.cos(b)));
    const r = Math.max(0, Math.round(Math.sqrt((12 * sig * sig) / 3 + 1) / 2 - 0.5));
    if (r === 0) continue;
    for (let pass = 0; pass < 3; pass++) {
      let acc = 0;
      for (let k = -r; k <= r; k++) acc += v[j * w + ((k % w) + w) % w]!;
      for (let i = 0; i < w; i++) {
        row[i] = acc / (2 * r + 1);
        acc += v[j * w + ((i + r + 1) % w)]! - v[j * w + ((i - r) % w + w) % w]!;
      }
      v.set(row, j * w);
    }
  }
  const sigY = sigmaDeg / texDeg;
  const ry = Math.ceil(3 * sigY);
  const ker = Array.from({ length: 2 * ry + 1 }, (_, k) => Math.exp(-((k - ry) ** 2) / (2 * sigY * sigY)));
  const col = new Float32Array(hgt);
  for (let i = 0; i < w; i++) {
    for (let j = 0; j < hgt; j++) {
      let acc = 0;
      let wsum = 0;
      for (let k = -ry; k <= ry; k++) {
        const jj = j + k;
        if (jj < 0 || jj >= hgt) continue;
        const kw = ker[k + ry]!;
        acc += kw * v[jj * w + i]!;
        wsum += kw;
      }
      col[j] = acc / wsum;
    }
    for (let j = 0; j < hgt; j++) v[j * w + i] = col[j]!;
  }
}

/** A soft radial glow (white, alpha falloff) for sprites and point sprites. */
export function glowTexture(size = 128): THREE.CanvasTexture {
  const c = document.createElement('canvas');
  c.width = c.height = size;
  const g = c.getContext('2d')!;
  const grd = g.createRadialGradient(size / 2, size / 2, 0, size / 2, size / 2, size / 2);
  grd.addColorStop(0, 'rgba(255,255,255,1)');
  grd.addColorStop(0.18, 'rgba(255,255,255,0.55)');
  grd.addColorStop(0.45, 'rgba(255,255,255,0.12)');
  grd.addColorStop(1, 'rgba(255,255,255,0)');
  g.fillStyle = grd;
  g.fillRect(0, 0, size, size);
  const tex = new THREE.CanvasTexture(c);
  tex.colorSpace = THREE.SRGBColorSpace;
  return tex;
}

/** Standard vertex/fragment pair for additive soft points with per-point size and colour. */
export const POINTS_VERT = /* glsl */ `
attribute float aSize;
attribute vec3 aColor;
attribute float aAlpha;
uniform float uScale;
varying vec3 vColor;
varying float vAlpha;
void main() {
  vColor = aColor;
  vAlpha = aAlpha;
  vec4 mv = modelViewMatrix * vec4(position, 1.0);
  gl_PointSize = aSize * uScale / max(0.001, -mv.z);
  gl_Position = projectionMatrix * mv;
}`;

export const POINTS_FRAG = /* glsl */ `
varying vec3 vColor;
varying float vAlpha;
void main() {
  vec2 c = gl_PointCoord - 0.5;
  float r2 = dot(c, c) * 4.0;
  if (r2 > 1.0) discard;
  float core = exp(-r2 * 14.0);
  float halo = exp(-r2 * 3.2) * 0.35;
  gl_FragColor = vec4(vColor * (core * 1.6 + halo) * vAlpha, 1.0);
}`;

/** Pixel scale for POINTS_VERT so that aSize is a world-space diameter. */
export function pointScale(camera: THREE.PerspectiveCamera, heightPx: number, pixelRatio: number): number {
  return (heightPx * pixelRatio) / (2 * Math.tan((camera.fov * Math.PI) / 360));
}
