/**
 * Sky views for HEALPix RING maps in Galactic coordinates: a Three.js globe (outer or inner
 * view, orbit controls, click picking) and a 2D Mollweide canvas, sharing one overlay API:
 *
 *     const sky = await createSkyView(el, { mode: 'globe' });
 *     sky.setMap({ map, mask, limit: 300, unit: 'μK' });
 *     sky.addSmallCircle([l, b], radiusDeg, 'var(--series-1)');
 *     sky.addPoint([l, b], '--series-2', 'image 1');
 *     sky.addChordThroughOrigin(lbA, lbB, '--series-3');   // a point and its twin through the observer
 *     sky.addPolygon(cornersGal, '--series-4');               // survey footprint
 *     sky.onPick(({ l, b, value }) => ...);
 *
 * Colours may be CSS colours or palette tokens ('--series-1' / 'var(--series-1)'); tokens are
 * re-resolved when the theme changes. When WebGL is unavailable `createSkyView` falls back to
 * the Mollweide view and says so on screen. Three.js is loaded only for the globe.
 */

import { type Theme, divergingLut, sequentialLut } from './colormap';
import { lb2pix, lbToVec, mollweideForward, mollweideInverse, nsideFromNpix, vecToLb } from './healpix';
import { cssVar, currentTheme, onThemeChange } from './theme';
import { fmtNum, h } from './ui';

/** Galactic (l, b) in degrees. */
export type LB = readonly [number, number];
type V3 = [number, number, number];

export interface SkyMapInput {
  /** RING-ordered values, length 12 nside^2 */
  map: Float32Array;
  /** 1 = kept, 0 = masked (drawn as neutral) */
  mask?: Uint8Array | null;
  /** symmetric colour limit (default: 99.5th percentile of |value| over unmasked pixels) */
  limit?: number;
  unit?: string;
  /** 'diverging' (default, signed maps such as CMB temperature) or 'sequential' (magnitude) */
  scale?: 'diverging' | 'sequential';
}

export interface SkyPick {
  l: number;
  b: number;
  pix: number | null;
  value: number | null;
  masked: boolean;
}

export type OverlayId = number;

export interface SkyOverlayApi {
  /** Small circle of angular radius `radiusDeg` about `center` (e.g. a matched-circle). */
  addSmallCircle(center: LB, radiusDeg: number, color: string, opts?: { label?: string; width?: number }): OverlayId;
  /** A marked direction with an optional label. */
  addPoint(lb: LB, color: string, label?: string): OverlayId;
  /** Straight segment from the sky point A through the observer at the centre to B. */
  addChordThroughOrigin(a: LB, b: LB, color: string, opts?: { label?: string }): OverlayId;
  /** Closed polygon with great-circle edges (survey footprints). */
  addPolygon(vertices: LB[], color: string, opts?: { fill?: boolean; label?: string }): OverlayId;
  removeOverlay(id: OverlayId): void;
  clearOverlays(): void;
}

export interface SkyView extends SkyOverlayApi {
  readonly kind: 'globe' | 'mollweide';
  /** the container element the view renders into */
  readonly el: HTMLElement;
  setMap(input: SkyMapInput | null): void;
  /** Click picking; returns an unsubscribe function. */
  onPick(cb: (p: SkyPick) => void): () => void;
  resize(): void;
  dispose(): void;
}

export interface GlobeView extends SkyView {
  readonly kind: 'globe';
  /** 'outer': the sphere seen from outside; 'inner': the sky as seen by the observer at the centre. */
  setView(mode: 'outer' | 'inner'): void;
  /** Sphere opacity in [0, 1]; < 1 reveals chords through the observer in the outer view. */
  setOpacity(alpha: number): void;
  /** Turn the camera towards a direction. */
  lookAt(lb: LB): void;
}

export interface SkyViewOptions {
  mode?: 'globe' | 'mollweide';
  view?: 'outer' | 'inner';
  /** CSS height of the view in px (Mollweide default: half the width) */
  height?: number;
  /** Galactic graticule (default true) and its spacing in degrees (default 30) */
  grid?: boolean;
  gridStepDeg?: number;
  /** pointer readout of (l, b, value) (default true) */
  hud?: boolean;
  /** colour bar when a map is shown (default true) */
  colorbar?: boolean;
}

// --- overlays (shared model) ------------------------------------------------------------------

export type Overlay =
  | { id: number; kind: 'circle'; center: LB; radiusDeg: number; color: string; label?: string; width?: number }
  | { id: number; kind: 'point'; lb: LB; color: string; label?: string }
  | { id: number; kind: 'chord'; a: LB; b: LB; color: string; label?: string }
  | { id: number; kind: 'polygon'; vertices: LB[]; color: string; fill?: boolean; label?: string };

type DistributiveOmit<T, K extends PropertyKey> = T extends unknown ? Omit<T, K> : never;

export class OverlayStore {
  private next = 1;
  items: Overlay[] = [];
  constructor(private readonly changed: () => void) {}

  add(o: DistributiveOmit<Overlay, 'id'>): OverlayId {
    const id = this.next++;
    this.items.push({ ...o, id } as Overlay);
    this.changed();
    return id;
  }

  remove(id: OverlayId): void {
    this.items = this.items.filter((o) => o.id !== id);
    this.changed();
  }

  clear(): void {
    this.items = [];
    this.changed();
  }

  api(): SkyOverlayApi {
    return {
      addSmallCircle: (center, radiusDeg, color, opts) =>
        this.add({ kind: 'circle', center, radiusDeg, color, label: opts?.label, width: opts?.width }),
      addPoint: (lb, color, label) => this.add({ kind: 'point', lb, color, label }),
      addChordThroughOrigin: (a, b, color, opts) => this.add({ kind: 'chord', a, b, color, label: opts?.label }),
      addPolygon: (vertices, color, opts) => this.add({ kind: 'polygon', vertices, color, fill: opts?.fill, label: opts?.label }),
      removeOverlay: (id) => this.remove(id),
      clearOverlays: () => this.clear(),
    };
  }
}

// --- geometry ----------------------------------------------------------------------------------

const DEG = Math.PI / 180;
const norm3 = (v: V3): V3 => {
  const n = Math.hypot(v[0], v[1], v[2]);
  return [v[0] / n, v[1] / n, v[2] / n];
};
const cross3 = (a: V3, b: V3): V3 => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];

/** Unit vectors on the small circle of radius `radiusDeg` about `center` (closed: last = first). */
export function smallCirclePoints(center: LB, radiusDeg: number, n = 180): V3[] {
  const c = lbToVec(center[0], center[1]);
  const e1 = norm3(cross3(c, Math.abs(c[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0]));
  const e2 = cross3(c, e1);
  const cr = Math.cos(radiusDeg * DEG);
  const sr = Math.sin(radiusDeg * DEG);
  const out: V3[] = [];
  for (let i = 0; i <= n; i++) {
    const t = (2 * Math.PI * i) / n;
    const ct = Math.cos(t);
    const st = Math.sin(t);
    out.push([
      cr * c[0] + sr * (ct * e1[0] + st * e2[0]),
      cr * c[1] + sr * (ct * e1[1] + st * e2[1]),
      cr * c[2] + sr * (ct * e1[2] + st * e2[2]),
    ]);
  }
  return out;
}

/** Points along the shorter great-circle arc from a to b (unit vectors), endpoints included. */
export function greatArcPoints(a: V3, b: V3, n?: number): V3[] {
  const d = Math.max(-1, Math.min(1, a[0] * b[0] + a[1] * b[1] + a[2] * b[2]));
  const omega = Math.acos(d);
  const steps = n ?? Math.max(2, Math.ceil(omega / (2 * DEG)));
  if (omega < 1e-9) return [a, b];
  const s = Math.sin(omega);
  const out: V3[] = [];
  for (let i = 0; i <= steps; i++) {
    const t = i / steps;
    const fa = Math.sin((1 - t) * omega) / s;
    const fb = Math.sin(t * omega) / s;
    out.push([fa * a[0] + fb * b[0], fa * a[1] + fb * b[1], fa * a[2] + fb * b[2]]);
  }
  return out;
}

/** Closed outline of a polygon with great-circle edges. */
export function polygonPoints(vertices: LB[]): V3[] {
  const vs = vertices.map((v) => lbToVec(v[0], v[1]));
  const out: V3[] = [];
  vs.forEach((v, i) => {
    const arc = greatArcPoints(v, vs[(i + 1) % vs.length]!);
    out.push(...(i ? arc.slice(1) : arc));
  });
  return out;
}

/** Graticule lines (unit vectors) every `step` degrees: meridians then parallels. */
export function graticule(step = 30): { lines: V3[][]; equator: V3[] } {
  const lines: V3[][] = [];
  for (let l = 0; l < 360; l += step) {
    const m: V3[] = [];
    for (let b = -90; b <= 90; b += 2) m.push(lbToVec(l, b));
    lines.push(m);
  }
  for (let b = -90 + step; b < 90; b += step) {
    if (b === 0) continue;
    const p: V3[] = [];
    for (let l = 0; l <= 360; l += 2) p.push(lbToVec(l, b));
    lines.push(p);
  }
  const equator: V3[] = [];
  for (let l = 0; l <= 360; l += 2) equator.push(lbToVec(l, 0));
  return { lines, equator };
}

// --- colouring ------------------------------------------------------------------------------------

/** Percentile of |value| over unmasked finite pixels (healpy UNSEEN counts as masked). */
export function robustLimit(map: Float32Array, mask?: Uint8Array | null, q = 0.995): number {
  const vals: number[] = [];
  const stride = map.length > 200_000 ? Math.ceil(map.length / 200_000) : 1;
  for (let i = 0; i < map.length; i += stride) {
    const v = map[i]!;
    if ((mask && !mask[i]) || !Number.isFinite(v) || v < -1e30) continue;
    vals.push(Math.abs(v));
  }
  if (!vals.length) return 1;
  vals.sort((a, b) => a - b);
  return vals[Math.min(vals.length - 1, Math.floor(q * vals.length))] || 1;
}

export interface Colorizer {
  limit: number;
  /** writes RGB of pixel `pix` into out[offset..offset+2]; returns false for masked pixels */
  color(pix: number, out: Uint8ClampedArray, offset: number): boolean;
  lut: Uint8ClampedArray;
}

/** Masked pixels: a neutral base with a diagonal hatch, so "excluded" never reads as "zero". */
const MASK_RGB: Record<Theme, [number, number, number]> = { dark: [27, 34, 47], light: [226, 230, 236] };
const HATCH_RGB: Record<Theme, [number, number, number]> = { dark: [58, 68, 90], light: [190, 198, 211] };

/** Paint the hatch over a masked pixel at raster position (i, j). */
export function hatchMasked(out: Uint8ClampedArray, o: number, i: number, j: number, theme: Theme, period = 6): void {
  if ((i + j) % period === 0) out.set(HATCH_RGB[theme], o);
}

export function makeColorizer(input: SkyMapInput, theme: Theme): Colorizer {
  const sequential = input.scale === 'sequential';
  const lut = sequential ? sequentialLut(theme) : divergingLut(theme);
  const limit = input.limit ?? robustLimit(input.map, input.mask);
  const map = input.map;
  const mask = input.mask ?? null;
  const maskRgb = MASK_RGB[theme];
  let min = 0;
  if (sequential) {
    min = Infinity;
    for (let i = 0; i < map.length; i++) if (Number.isFinite(map[i]!) && map[i]! > -1e30) min = Math.min(min, map[i]!);
  }
  return {
    limit,
    lut,
    color(pix, out, o) {
      const v = map[pix]!;
      if ((mask && !mask[pix]) || !Number.isFinite(v) || v < -1e30) {
        out[o] = maskRgb[0];
        out[o + 1] = maskRgb[1];
        out[o + 2] = maskRgb[2];
        return false;
      }
      const t = sequential ? (v - min) / (limit - min || 1) : (v / limit + 1) / 2;
      const k = 3 * Math.max(0, Math.min(255, Math.round(t * 255)));
      out[o] = lut[k]!;
      out[o + 1] = lut[k + 1]!;
      out[o + 2] = lut[k + 2]!;
      return true;
    },
  };
}

/**
 * Equirectangular RGBA texture of a HEALPix map: column i <-> l = (i + 1/2) 360/W,
 * row j <-> b = -90 + (j + 1/2) 180/H (row 0 at the south pole, as sampled by the globe shader).
 */
export function equirectTexture(input: SkyMapInput, theme: Theme, width = 1024, height = 512): Uint8ClampedArray {
  const nside = nsideFromNpix(input.map.length);
  const col = makeColorizer(input, theme);
  const data = new Uint8ClampedArray(width * height * 4);
  for (let j = 0; j < height; j++) {
    const b = -90 + ((j + 0.5) * 180) / height;
    for (let i = 0; i < width; i++) {
      const l = ((i + 0.5) * 360) / width;
      const o = 4 * (j * width + i);
      if (!col.color(lb2pix(nside, l, b), data, o)) hatchMasked(data, o, i, j, theme);
      data[o + 3] = 255;
    }
  }
  return data;
}

/** Look up (l, b) in a map. */
export function pickValue(input: SkyMapInput | null, l: number, b: number): SkyPick {
  if (!input) return { l, b, pix: null, value: null, masked: false };
  const nside = nsideFromNpix(input.map.length);
  const pix = lb2pix(nside, l, b);
  const v = input.map[pix]!;
  const masked = !!(input.mask && !input.mask[pix]) || !Number.isFinite(v) || v < -1e30;
  return { l, b, pix, value: Number.isFinite(v) && v > -1e30 ? v : null, masked };
}

/** Resolve '--token' / 'var(--token)' / any CSS colour to a concrete colour string. */
export function resolveCssColor(color: string): string {
  const m = /^var\((--[\w-]+)\)$/.exec(color.trim());
  const token = m ? m[1]! : color.trim().startsWith('--') ? color.trim() : null;
  return token ? cssVar(token) || '#888888' : color;
}

export function formatPick(p: SkyPick, unit?: string): string {
  const lb = `l ${p.l.toFixed(1)}°  b ${p.b >= 0 ? '+' : '−'}${Math.abs(p.b).toFixed(1)}°`;
  if (p.value === null) return lb;
  return `${lb}  ·  ${fmtNum(p.value, 4)}${unit ? ` ${unit}` : ''}${p.masked ? ' (masked)' : ''}`;
}

function colorbar(lut: Uint8ClampedArray, lo: string, hi: string): HTMLElement {
  const c = document.createElement('canvas');
  c.width = 256;
  c.height = 1;
  const ctx = c.getContext('2d');
  if (ctx) {
    const img = ctx.createImageData(256, 1);
    for (let i = 0; i < 256; i++) {
      img.data.set([lut[3 * i]!, lut[3 * i + 1]!, lut[3 * i + 2]!, 255], 4 * i);
    }
    ctx.putImageData(img, 0, 0);
  }
  return h('div', { class: 'sky-colorbar', 'aria-hidden': 'true' }, h('span', null, lo), c, h('span', null, hi));
}

/** Shared chrome: HUD readout and colour bar for both renderers. */
export class SkyChrome {
  readonly hud: HTMLElement;
  private bar: HTMLElement | null = null;
  constructor(
    private readonly root: HTMLElement,
    private readonly opts: SkyViewOptions,
  ) {
    this.hud = h('div', { class: 'sky-hud', 'aria-live': 'off' });
    if (opts.hud !== false) root.append(this.hud);
  }
  setHud(text: string): void {
    this.hud.textContent = text;
  }
  setColorbar(input: SkyMapInput | null, theme: Theme): void {
    this.bar?.remove();
    this.bar = null;
    if (!input || this.opts.colorbar === false) return;
    const col = makeColorizer(input, theme);
    const unit = input.unit ? ` ${input.unit}` : '';
    const lo = input.scale === 'sequential' ? 'min' : `−${fmtNum(col.limit, 3)}`;
    this.bar = colorbar(col.lut, lo, `${input.scale === 'sequential' ? '' : '+'}${fmtNum(col.limit, 3)}${unit}`);
    this.root.append(this.bar);
  }
}

export function webglAvailable(): boolean {
  try {
    const c = document.createElement('canvas');
    return !!(window.WebGLRenderingContext && (c.getContext('webgl2') || c.getContext('webgl')));
  } catch {
    return false;
  }
}

/**
 * Create a sky view in `el`. A globe needs WebGL; without it (or if the 3D renderer fails)
 * the Mollweide view is returned instead, with an on-screen note.
 */
export async function createSkyView(el: HTMLElement, opts: SkyViewOptions = {}): Promise<SkyView> {
  if (opts.mode !== 'mollweide') {
    if (webglAvailable()) {
      try {
        const { createGlobe } = await import('./sky-globe');
        return createGlobe(el, opts);
      } catch (err) {
        console.warn('Q-CHAS: globe unavailable, using Mollweide', err);
      }
    }
    const view = createMollweide(el, opts);
    view.el.append(h('div', { class: 'sky-fallback-note' }, 'WebGL unavailable: Mollweide projection'));
    return view;
  }
  return createMollweide(el, opts);
}

// --- Mollweide canvas renderer --------------------------------------------------------------------

const SQRT2 = Math.SQRT2;

/** 2D Mollweide view (Galactic centre in the middle, l increasing to the left). */
export function createMollweide(el: HTMLElement, opts: SkyViewOptions = {}): SkyView {
  const canvas = h('canvas', { role: 'img', 'aria-label': 'All-sky map, Mollweide projection, Galactic coordinates' });
  const root = h('div', { class: 'sky-view' }, canvas);
  el.append(root);
  const chrome = new SkyChrome(root, opts);
  let input: SkyMapInput | null = null;
  let theme = currentTheme();
  let mapLayer: HTMLCanvasElement | null = null;
  // canvas size, device-pixel ratio and the 2:1 box (x0, y0, ew, eh) holding the ellipse
  let dims = { w: 0, h: 0, dpr: 1, x0: 0, y0: 0, ew: 0, eh: 0 };
  const pickers = new Set<(p: SkyPick) => void>();
  let disposed = false;
  let frame = 0;

  const store = new OverlayStore(() => schedule(false));

  const toCanvas = (x: number, y: number): [number, number] => {
    const { x0, y0, ew, eh } = dims;
    return [x0 + ((x + 2 * SQRT2) / (4 * SQRT2)) * ew, y0 + ((SQRT2 - y) / (2 * SQRT2)) * eh];
  };
  const fromCanvas = (px: number, py: number): [number, number] => {
    const { x0, y0, ew, eh } = dims;
    return [((px - x0) / ew) * 4 * SQRT2 - 2 * SQRT2, SQRT2 - ((py - y0) / eh) * 2 * SQRT2];
  };
  const project = (v: V3): [number, number] => {
    const [l, b] = vecToLb(v);
    const [x, y] = mollweideForward(l, b);
    return toCanvas(x, y);
  };

  function measure(): void {
    const w = Math.max(200, root.clientWidth || el.clientWidth || 600);
    const hh = opts.height ?? Math.round(w / 2);
    const dpr = Math.min(2, window.devicePixelRatio || 1);
    const m = 6;
    const ew = Math.min(w - 2 * m, 2 * (hh - 2 * m));
    const eh = ew / 2;
    dims = { w, h: hh, dpr, x0: (w - ew) / 2, y0: (hh - eh) / 2, ew, eh };
    canvas.style.width = `${w}px`;
    canvas.style.height = `${hh}px`;
    canvas.width = Math.round(w * dpr);
    canvas.height = Math.round(hh * dpr);
    root.style.height = `${hh}px`;
  }

  function renderMapLayer(): void {
    const { w, h: hh, dpr } = dims;
    const W = Math.round(w * dpr);
    const H = Math.round(hh * dpr);
    mapLayer = document.createElement('canvas');
    mapLayer.width = W;
    mapLayer.height = H;
    const ctx = mapLayer.getContext('2d');
    if (!ctx) return;
    const img = ctx.createImageData(W, H);
    const col = input ? makeColorizer(input, theme) : null;
    const nside = input ? nsideFromNpix(input.map.length) : 0;
    const blank = theme === 'dark' ? [20, 29, 46] : [236, 240, 246];
    for (let j = 0; j < H; j++) {
      for (let i = 0; i < W; i++) {
        const [x, y] = fromCanvas((i + 0.5) / dpr, (j + 0.5) / dpr);
        const lb = mollweideInverse(x, y);
        if (!lb) continue;
        const o = 4 * (j * W + i);
        if (col) {
          if (!col.color(lb2pix(nside, lb[0], lb[1]), img.data, o)) hatchMasked(img.data, o, i, j, theme, Math.round(6 * dpr));
        } else img.data.set(blank, o);
        img.data[o + 3] = 255;
      }
    }
    ctx.putImageData(img, 0, 0);
  }

  function polyline(ctx: CanvasRenderingContext2D, pts: V3[]): void {
    let prev: [number, number] | null = null;
    ctx.beginPath();
    for (const v of pts) {
      const p = project(v);
      if (!prev || Math.abs(p[0] - prev[0]) > dims.ew / 3) ctx.moveTo(p[0], p[1]);
      else ctx.lineTo(p[0], p[1]);
      prev = p;
    }
    ctx.stroke();
  }

  function marker(ctx: CanvasRenderingContext2D, p: [number, number], color: string, label?: string): void {
    ctx.beginPath();
    ctx.arc(p[0], p[1], 4.5, 0, 2 * Math.PI);
    ctx.fillStyle = color;
    ctx.fill();
    ctx.lineWidth = 2;
    ctx.strokeStyle = cssVar('--surface');
    ctx.stroke();
    if (label) {
      ctx.font = `600 10.5px ${cssVar('--font-mono')}`;
      ctx.fillStyle = cssVar('--text-1');
      ctx.textBaseline = 'middle';
      ctx.lineWidth = 3;
      ctx.strokeStyle = cssVar('--surface');
      ctx.strokeText(label, p[0] + 8, p[1]);
      ctx.fillText(label, p[0] + 8, p[1]);
    }
  }

  function draw(): void {
    const ctx = canvas.getContext('2d');
    if (!ctx) return;
    const { w, h: hh, dpr } = dims;
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    if (mapLayer) ctx.drawImage(mapLayer, 0, 0);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.lineJoin = 'round';
    ctx.lineCap = 'round';
    // graticule
    if (opts.grid !== false) {
      const g = graticule(opts.gridStepDeg ?? 30);
      ctx.strokeStyle = cssVar('--text-3');
      ctx.globalAlpha = 0.35;
      ctx.lineWidth = 1;
      g.lines.forEach((ln) => polyline(ctx, ln));
      ctx.globalAlpha = 0.7;
      polyline(ctx, g.equator);
      ctx.globalAlpha = 1;
    }
    // outline
    ctx.beginPath();
    ctx.ellipse(w / 2, hh / 2, dims.ew / 2, dims.eh / 2, 0, 0, 2 * Math.PI);
    ctx.strokeStyle = cssVar('--line-2');
    ctx.lineWidth = 1;
    ctx.stroke();
    // overlays
    for (const o of store.items) {
      const color = resolveCssColor(o.color);
      ctx.strokeStyle = color;
      ctx.lineWidth = 2;
      if (o.kind === 'circle') {
        polyline(ctx, smallCirclePoints(o.center, o.radiusDeg));
        if (o.label) marker(ctx, project(lbToVec(o.center[0], o.center[1])), color, o.label);
      } else if (o.kind === 'polygon') {
        const pts = polygonPoints(o.vertices);
        if (o.fill) {
          ctx.beginPath();
          pts.forEach((v, i) => {
            const p = project(v);
            if (i) ctx.lineTo(p[0], p[1]);
            else ctx.moveTo(p[0], p[1]);
          });
          ctx.closePath();
          ctx.globalAlpha = 0.18;
          ctx.fillStyle = color;
          ctx.fill();
          ctx.globalAlpha = 1;
        }
        polyline(ctx, pts);
      } else if (o.kind === 'chord') {
        // the chord passes through the observer; on a projection we join the two sky points
        const pa = project(lbToVec(o.a[0], o.a[1]));
        const pb = project(lbToVec(o.b[0], o.b[1]));
        ctx.globalAlpha = 0.8;
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.moveTo(pa[0], pa[1]);
        ctx.lineTo(pb[0], pb[1]);
        ctx.stroke();
        ctx.globalAlpha = 1;
        marker(ctx, pa, color, o.label);
        marker(ctx, pb, color);
      }
    }
    for (const o of store.items) if (o.kind === 'point') marker(ctx, project(lbToVec(o.lb[0], o.lb[1])), resolveCssColor(o.color), o.label);
  }

  function schedule(rebuildMap: boolean): void {
    if (disposed) return;
    if (rebuildMap) mapLayer = null;
    cancelAnimationFrame(frame);
    frame = requestAnimationFrame(() => {
      if (!root.isConnected) return view.dispose();
      if (!mapLayer) renderMapLayer();
      draw();
    });
  }

  const pickAt = (ev: MouseEvent): SkyPick | null => {
    const r = canvas.getBoundingClientRect();
    const [x, y] = fromCanvas(ev.clientX - r.left, ev.clientY - r.top);
    const lb = mollweideInverse(x, y);
    return lb ? pickValue(input, lb[0], lb[1]) : null;
  };
  canvas.addEventListener('click', (ev) => {
    const p = pickAt(ev);
    if (p) pickers.forEach((cb) => cb(p));
  });
  canvas.addEventListener('pointermove', (ev) => {
    const p = pickAt(ev);
    chrome.setHud(p ? formatPick(p, input?.unit) : '');
  });
  canvas.addEventListener('pointerleave', () => chrome.setHud(''));

  const ro = new ResizeObserver(() => {
    const before = dims.w;
    measure();
    schedule(before !== dims.w);
  });
  ro.observe(root);
  const offTheme = onThemeChange((t) => {
    theme = t;
    chrome.setColorbar(input, theme);
    schedule(true);
  });

  measure();
  schedule(true);

  const view: SkyView = {
    kind: 'mollweide',
    el: root,
    ...store.api(),
    setMap(next) {
      input = next;
      chrome.setColorbar(input, theme);
      schedule(true);
    },
    onPick(cb) {
      pickers.add(cb);
      return () => pickers.delete(cb);
    },
    resize() {
      measure();
      schedule(true);
    },
    dispose() {
      if (disposed) return;
      disposed = true;
      cancelAnimationFrame(frame);
      ro.disconnect();
      offTheme();
      pickers.clear();
    },
  };
  return view;
}
