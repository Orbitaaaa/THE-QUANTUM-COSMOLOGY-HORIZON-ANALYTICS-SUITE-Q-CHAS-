/**
 * Perceptual colormaps for maps and heatmaps, built in OKLab from the dashboard palette.
 *
 * * Diverging (CMB temperature, residuals): blue arm from the documented blue ramp, red arm with
 *   the same OKLCH lightness steps at the hue of the palette red, neutral gray midpoint. Each
 *   theme has its own steps (light: dark poles on a light midpoint; dark: bright poles on a dark
 *   midpoint).
 * * Sequential (magnitude): the blue ramp, lightest step near zero in light mode and darkest
 *   step near zero in dark mode, so "nothing" always recedes into the surface.
 */

export type Theme = 'dark' | 'light';
export type RGB = [number, number, number];

/** Documented sequential blue ramp, steps 100 -> 700. */
export const BLUE_RAMP = [
  '#cde2fb', '#b7d3f6', '#9ec5f4', '#86b6ef', '#6da7ec', '#5598e7', '#3987e5',
  '#2a78d6', '#256abf', '#1c5cab', '#184f95', '#104281', '#0d366b',
] as const;

const RED_POLE = { light: '#e34948', dark: '#e66767' } as const;
const MIDPOINT = { light: '#f0efec', dark: '#383835' } as const;
/** Blue-arm steps from the midpoint outwards (indices into BLUE_RAMP). */
const ARM_STEPS = { light: [2, 5, 7, 9, 12], dark: [10, 7, 5, 3, 1] } as const;

export function hexToRgb(hex: string): RGB {
  const h = hex.replace('#', '');
  const v = parseInt(h.length === 3 ? h.replace(/./g, (c) => c + c) : h, 16);
  return [(v >> 16) & 255, (v >> 8) & 255, v & 255];
}

export const rgbToHex = (c: RGB): string =>
  '#' + c.map((x) => Math.round(Math.min(255, Math.max(0, x))).toString(16).padStart(2, '0')).join('');

const toLinear = (c: number): number => {
  const x = c / 255;
  return x <= 0.04045 ? x / 12.92 : ((x + 0.055) / 1.055) ** 2.4;
};
const fromLinear = (x: number): number => 255 * (x <= 0.0031308 ? 12.92 * x : 1.055 * x ** (1 / 2.4) - 0.055);

export function rgbToOklab([r8, g8, b8]: RGB): [number, number, number] {
  const r = toLinear(r8);
  const g = toLinear(g8);
  const b = toLinear(b8);
  const l = Math.cbrt(0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b);
  const m = Math.cbrt(0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b);
  const s = Math.cbrt(0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b);
  return [
    0.2104542553 * l + 0.793617785 * m - 0.0040720468 * s,
    1.9779984951 * l - 2.428592205 * m + 0.4505937099 * s,
    0.0259040371 * l + 0.7827717662 * m - 0.808675766 * s,
  ];
}

/** OKLab -> linear sRGB (unclamped). */
function oklabToLinear([L, a, b]: [number, number, number]): [number, number, number] {
  const l = (L + 0.3963377774 * a + 0.2158037573 * b) ** 3;
  const m = (L - 0.1055613458 * a - 0.0638541728 * b) ** 3;
  const s = (L - 0.0894841775 * a - 1.291485548 * b) ** 3;
  return [
    4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s,
    -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s,
    -0.0041960863 * l - 0.7034186147 * m + 1.707614701 * s,
  ];
}

export const oklabToRgb = (lab: [number, number, number]): RGB =>
  oklabToLinear(lab).map((x) => fromLinear(Math.min(1, Math.max(0, x)))) as RGB;

const inGamut = (lab: [number, number, number]): boolean => oklabToLinear(lab).every((x) => x >= -1e-6 && x <= 1 + 1e-6);

/** Colour at OKLCH (L, C, h), reducing chroma until it fits in sRGB. */
function oklch(L: number, C: number, h: number): RGB {
  let c = C;
  for (let i = 0; i < 40; i++) {
    const lab: [number, number, number] = [L, c * Math.cos(h), c * Math.sin(h)];
    if (inGamut(lab)) return oklabToRgb(lab);
    c *= 0.95;
  }
  return oklabToRgb([L, 0, 0]);
}

const lch = (hex: string): { L: number; C: number; h: number } => {
  const [L, a, b] = rgbToOklab(hexToRgb(hex));
  return { L, C: Math.hypot(a, b), h: Math.atan2(b, a) };
};

/** Stops of the diverging map from -1 (blue pole) through 0 (midpoint) to +1 (red pole). */
export function divergingStops(theme: Theme): RGB[] {
  const blue = ARM_STEPS[theme].map((i) => BLUE_RAMP[i]!);
  const redHue = lch(RED_POLE[theme]).h;
  const redC = lch(RED_POLE[theme]).C;
  const red = blue.map((hex) => {
    const b = lch(hex);
    return rgbToHex(oklch(b.L, Math.max(b.C, 0.6 * redC), redHue));
  });
  const mid = MIDPOINT[theme];
  return [...[...blue].reverse(), mid, ...red].map(hexToRgb);
}

/** Interpolate a list of stops in OKLab into a lookup table of `n` RGB triplets. */
export function buildLut(stops: RGB[], n = 256): Uint8ClampedArray {
  const labs = stops.map(rgbToOklab);
  const out = new Uint8ClampedArray(n * 3);
  for (let i = 0; i < n; i++) {
    const u = (i / (n - 1)) * (labs.length - 1);
    const k = Math.min(labs.length - 2, Math.floor(u));
    const f = u - k;
    const a = labs[k]!;
    const b = labs[k + 1]!;
    const rgb = oklabToRgb([a[0] + f * (b[0] - a[0]), a[1] + f * (b[1] - a[1]), a[2] + f * (b[2] - a[2])]);
    out.set(rgb, 3 * i);
  }
  return out;
}

const lutCache = new Map<string, Uint8ClampedArray>();

export function divergingLut(theme: Theme): Uint8ClampedArray {
  const key = `div-${theme}`;
  let lut = lutCache.get(key);
  if (!lut) lutCache.set(key, (lut = buildLut(divergingStops(theme))));
  return lut;
}

export function sequentialLut(theme: Theme): Uint8ClampedArray {
  const key = `seq-${theme}`;
  let lut = lutCache.get(key);
  if (!lut) {
    const ramp = BLUE_RAMP.map(hexToRgb);
    lutCache.set(key, (lut = buildLut(theme === 'light' ? ramp : [...ramp].reverse())));
  }
  return lut;
}

/** Plotly colorscale ([[0, 'rgb(..)'], ..., [1, ...]]) sampled from a LUT. */
export function plotlyScale(lut: Uint8ClampedArray, samples = 11): [number, string][] {
  const n = lut.length / 3;
  return Array.from({ length: samples }, (_, i) => {
    const k = Math.round((i / (samples - 1)) * (n - 1));
    return [i / (samples - 1), `rgb(${lut[3 * k]},${lut[3 * k + 1]},${lut[3 * k + 2]})`] as [number, string];
  });
}

/** Index into a 256-entry LUT for value v on a symmetric scale [-limit, limit]. */
export const divergingIndex = (v: number, limit: number): number =>
  Math.max(0, Math.min(255, Math.round(((v / limit + 1) / 2) * 255)));
