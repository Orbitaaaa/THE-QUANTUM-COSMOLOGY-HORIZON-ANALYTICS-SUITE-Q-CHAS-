import { describe, expect, it } from 'vitest';
import {
  capacityCurve,
  capacityWindow,
  deSitterEntropy,
  minimalArea,
  schwarzschildEntropy,
} from '../src/lib/physics/capacity';
import {
  bhEntropy,
  curves,
  decodingInfidelityBound,
  haydenPreskillEmissions,
  hawkingEntropy,
  interiorVolume,
  islandMinEntropy,
  massFraction,
  pageCurve,
  pageCurveLower,
  pageMassFraction,
  pageTimeFraction,
  pageTransferSlope,
  totalQuanta,
} from '../src/lib/physics/page';

describe('Page law, [BH] Theorems C-F', () => {
  it('Page time 0.476 and mass fraction 0.806 ([BH] Fig. 1), as qchas.horizon.page', () => {
    expect(pageTimeFraction()).toBeCloseTo(0.47627090392600924, 14);
    expect(pageMassFraction()).toBeCloseTo(0.8060628405882714, 14);
    expect(pageTimeFraction().toFixed(3)).toBe('0.476');
    expect(pageMassFraction().toFixed(3)).toBe('0.806');
  });

  it('Hawking and capacity curves cross at the Page time; Page curve is their minimum', () => {
    const tp = pageTimeFraction();
    expect(hawkingEntropy(tp)).toBeCloseTo(bhEntropy(tp), 12);
    expect(islandMinEntropy(tp)).toBeCloseTo(0, 12);
    expect(hawkingEntropy(1)).toBeCloseTo(1.855, 12);
    expect(pageCurve(1)).toBe(0);
    for (const t of [0, 0.2, 0.5, 0.9]) expect(pageCurve(t)).toBe(Math.min(hawkingEntropy(t), bhEntropy(t)));
    expect(massFraction(0.875)).toBeCloseTo(0.5, 14);
    expect(pageTransferSlope(0.875)).toBeCloseTo(-4 / 3, 14);
  });

  it('2-design lower bound and interior volume agree with the Python', () => {
    expect(pageCurveLower(0.3, 100)).toBeCloseTo(39.25671272, 7);
    expect(pageCurveLower(0.476, 100)).toBeCloseTo(63.99413123, 7);
    expect(pageCurveLower(0.7, 100)).toBeCloseTo(44.81404747, 7);
    expect(interiorVolume(1)).toBeCloseTo((9 * Math.sqrt(3) * Math.PI) / 5, 12);
    expect(interiorVolume(1)).toBeCloseTo(9.794516566864775, 12);
  });

  it('curves() has the horizon.json.page layout', () => {
    const c = curves(5);
    expect(c.t).toEqual([0, 0.25, 0.5, 0.75, 1]);
    expect(c.interior_volume[4]).toBeCloseTo(1, 14);
    expect(c.beta).toBe(1.855);
    expect(Object.keys(c).sort()).toEqual(
      ['beta', 'bh_entropy', 'hawking', 'interior_volume', 'island_min', 'page', 'page_mass_fraction', 'page_time', 't'].sort(),
    );
  });

  it('quanta count and decoding bounds', () => {
    expect(totalQuanta(10)).toBeCloseTo(447, 10);
    expect(decodingInfidelityBound(2, 3)).toBe(1);
    expect(decodingInfidelityBound(2, 10)).toBe(2 ** -7);
    expect(haydenPreskillEmissions(3, 1 / 1024)).toBe(13);
  });
});

describe('capacity, [BH] Theorem G and [QG] Theorems C, U', () => {
  it('window capacity matches qchas.horizon.capacity.capacity_window', () => {
    const py: [number, number, number, number, number][] = [
      [10, 0.1, 19.605984708182902, 19.731648414326358, 25.2584049348618],
      [100, 1.0, 2505.444781258317, 2518.0111518726726, 2525.84049348619],
      [1.0, 0.01, -4.476875969323912, -4.475619332262475, 0.25258404934862],
      [1000, 5, 125653.5742168852, 125967.7334822449, 125977.86540895142],
    ];
    for (const [m, w, ln, g, ds] of py) {
      const c = capacityWindow(m, w);
      expect(c.ln_trace).toBeCloseTo(ln, 9);
      expect(c.theorem_G_v).toBeCloseTo(g, 8);
      expect(c.delta_S_BH).toBeCloseTo(ds, 8);
      expect(c.beta).toBeCloseTo(8 * Math.PI * m, 12);
    }
    const curve = capacityCurve(10, [0.1, 0.2]);
    expect(curve.width).toEqual([0.1, 0.2]);
    expect(schwarzschildEntropy(1)).toBeCloseTo(4 * Math.PI, 14);
  });

  it('de Sitter entropy and the area gap', () => {
    const ds = deSitterEntropy();
    expect(ds.log10_entropy).toBeCloseTo(122.52068405083706, 9);
    expect(ds.horizon_radius_Gpc).toBeCloseTo(5.381786087948357, 9);
    expect(minimalArea()).toBeCloseTo(5.965948868839189, 12);
    expect(minimalArea().toFixed(2)).toBe('5.97');
  });
});
