import { describe, expect, it } from 'vitest';
import {
  BLUE_RAMP,
  divergingIndex,
  divergingLut,
  divergingStops,
  hexToRgb,
  oklabToRgb,
  plotlyScale,
  rgbToHex,
  rgbToOklab,
  sequentialLut,
} from '../src/lib/colormap';

const L = (rgb: [number, number, number]): number => rgbToOklab(rgb)[0];

describe('colormaps', () => {
  it('round-trips sRGB through OKLab', () => {
    for (const hex of ['#2a78d6', '#e34948', '#f0efec', '#383835', '#000000', '#ffffff']) {
      expect(rgbToHex(oklabToRgb(rgbToOklab(hexToRgb(hex))))).toBe(hex);
    }
  });

  it('diverging map: gray midpoint, blue and red poles with matched lightness per step', () => {
    for (const theme of ['light', 'dark'] as const) {
      const stops = divergingStops(theme);
      expect(stops).toHaveLength(11);
      const mid = rgbToOklab(stops[5]!);
      expect(Math.hypot(mid[1], mid[2])).toBeLessThan(0.01); // neutral, no hue
      for (let k = 1; k <= 5; k++) expect(L(stops[5 - k]!)).toBeCloseTo(L(stops[5 + k]!), 1);
      const blue = rgbToOklab(stops[0]!);
      const red = rgbToOklab(stops[10]!);
      expect(blue[2]).toBeLessThan(0); // b < 0: blue
      expect(red[1]).toBeGreaterThan(0.05); // a > 0: red
    }
    // light theme poles are darker than the midpoint, dark theme poles lighter
    expect(L(divergingStops('light')[0]!)).toBeLessThan(L(divergingStops('light')[5]!));
    expect(L(divergingStops('dark')[0]!)).toBeGreaterThan(L(divergingStops('dark')[5]!));
  });

  it('sequential map is monotone in lightness, receding into the surface at zero', () => {
    for (const theme of ['light', 'dark'] as const) {
      const lut = sequentialLut(theme);
      const ls = Array.from({ length: 256 }, (_, i) => L([lut[3 * i]!, lut[3 * i + 1]!, lut[3 * i + 2]!]));
      const dir = theme === 'light' ? -1 : 1;
      for (let i = 1; i < 256; i++) expect(dir * (ls[i]! - ls[i - 1]!)).toBeGreaterThanOrEqual(-1e-3);
    }
    expect(BLUE_RAMP).toHaveLength(13);
  });

  it('LUT helpers', () => {
    expect(divergingLut('dark')).toHaveLength(768);
    expect(divergingIndex(0, 100)).toBe(128);
    expect(divergingIndex(-1e9, 100)).toBe(0);
    expect(divergingIndex(1e9, 100)).toBe(255);
    const sc = plotlyScale(divergingLut('light'), 5);
    expect(sc[0]![0]).toBe(0);
    expect(sc[4]![0]).toBe(1);
    expect(sc[2]![1]).toMatch(/^rgb\(\d+,\d+,\d+\)$/);
  });
});
