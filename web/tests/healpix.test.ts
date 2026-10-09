import { describe, expect, it } from 'vitest';
import {
  ang2pix,
  lb2pix,
  lbToVec,
  mollweideForward,
  mollweideInverse,
  npix,
  nsideFromNpix,
  pix2ang,
  pix2lb,
  vec2pix,
  vecToLb,
} from '../src/lib/healpix';

// reference values from healpy 1.x (pix2ang / ang2pix, RING)
const NSIDE1: [number, number, number][] = [
  [0, 0.84106867056793, 0.785398163397448], [1, 0.84106867056793, 2.356194490192345],
  [2, 0.84106867056793, 3.926990816987241], [3, 0.84106867056793, 5.497787143782138],
  [4, 1.570796326794897, 0.0], [5, 1.570796326794897, 1.570796326794897],
  [6, 1.570796326794897, 3.141592653589793], [7, 1.570796326794897, 4.71238898038469],
  [8, 2.300523983021863, 0.785398163397448], [9, 2.300523983021863, 2.356194490192345],
  [10, 2.300523983021863, 3.926990816987241], [11, 2.300523983021863, 5.497787143782138],
];
const NSIDE64: [number, number, number][] = [
  [0, 0.012757845597671, 0.785398163397448], [1, 0.012757845597671, 2.356194490192345],
  [3, 0.012757845597671, 5.497787143782138], [4, 0.025516210357419, 0.392699081698724],
  [100, 0.089334016244662, 3.702591341730828], [1000, 0.28160021642986, 5.462087227264072],
  [8064, 0.84106867056793, 0.012271846303085], [8069, 0.84106867056793, 0.134990309333936],
  [20000, 1.382179940619493, 3.939262663290326], [24576, 1.570796326794897, 3.153864499892878],
  [24577, 1.570796326794897, 3.178408192499048], [30000, 1.79132958771573, 4.319689898685965],
  [40000, 2.245927859731928, 4.724660826687775], [41087, 2.300523983021863, 6.270913460876501],
  [49051, 3.052258637345132, 2.580593965448759], [49147, 3.116076443232374, 5.890486225480862],
  [49151, 3.128834807992122, 5.497787143782138],
];
// random directions (theta, phi) and the healpy pixel for nside 1, 2, 64, 1024
const ANGLES: [number, number][] = [
  [1.547150874065, 2.071765476418], [0.448894531728, 4.95384364514], [2.362684197285, 1.905029296618],
  [0.457184468228, 2.849411276022], [1.956635874695, 0.842208822693], [1.724750880434, 2.532833593578],
  [0.856077355889, 1.278346978885], [1.753411344726, 1.648163326541], [1.471445590714, 4.714680286096],
  [2.808029646979, 1.761860188182], [1.039065193105, 3.048544801725], [1.494435512335, 6.162153563996],
];
const ANG2PIX: Record<number, number[]> = {
  1: [5, 3, 9, 1, 8, 9, 0, 5, 7, 9, 6, 4],
  2: [22, 3, 38, 1, 29, 31, 5, 30, 18, 45, 16, 27],
  64: [24020, 2490, 42143, 2443, 33698, 28391, 8372, 28867, 22336, 47777, 12028, 22907],
  1024: [6143302, 623371, 10769554, 645142, 8657445, 7253619, 2169665, 7433266, 5669889, 12236434, 3100611, 5814193],
};

describe('HEALPix RING', () => {
  it('nside 1: pixel 0 centre at theta = acos(2/3), phi = pi/4 and all twelve centres', () => {
    const [t0, p0] = pix2ang(1, 0);
    expect(t0).toBeCloseTo(Math.acos(2 / 3), 14);
    expect(p0).toBeCloseTo(Math.PI / 4, 14);
    for (const [pix, th, ph] of NSIDE1) {
      const [t, p] = pix2ang(1, pix);
      expect(t).toBeCloseTo(th, 12);
      expect(p).toBeCloseTo(ph, 12);
    }
  });

  it('nside 2: ring structure 4, 8, 8, 8, 8, 8, 4 with the HEALPix ring latitudes', () => {
    const rings = new Map<number, number[]>();
    for (let p = 0; p < npix(2); p++) {
      const [t, ph] = pix2ang(2, p);
      const key = Math.round(Math.cos(t) * 1e9);
      rings.set(key, [...(rings.get(key) ?? []), ph]);
    }
    const zs = [...rings.keys()].sort((a, b) => b - a).map((k) => k / 1e9);
    // z = 1 - i^2/(3 n^2) in the caps (i = 1), 4/3 - 2i/(3n) in the equatorial belt
    expect(zs).toEqual([11 / 12, 2 / 3, 1 / 3, 0, -1 / 3, -2 / 3, -11 / 12].map((z) => Math.round(z * 1e9) / 1e9));
    expect([...rings.keys()].sort((a, b) => b - a).map((k) => rings.get(k)!.length)).toEqual([4, 8, 8, 8, 8, 8, 4]);
    // shifted rings: the z = 1/3 ring starts at phi = 0, the z = 2/3 ring at pi/8
    expect(pix2ang(2, 12)[1]).toBeCloseTo(0, 14);
    expect(pix2ang(2, 4)[1]).toBeCloseTo(Math.PI / 8, 14);
  });

  it('nside 64 pixel centres match healpy', () => {
    for (const [pix, th, ph] of NSIDE64) {
      const [t, p] = pix2ang(64, pix);
      expect(t).toBeCloseTo(th, 12);
      expect(p).toBeCloseTo(ph, 12);
    }
  });

  it('ang2pix matches healpy for random directions at nside 1, 2, 64, 1024', () => {
    for (const [nside, expected] of Object.entries(ANG2PIX)) {
      ANGLES.forEach(([t, p], i) => expect(ang2pix(Number(nside), t, p)).toBe(expected[i]));
    }
  });

  it('round trips pix -> ang -> pix for every pixel of nside 1, 2, 4, 64 and samples of 8192', () => {
    for (const nside of [1, 2, 4, 64]) {
      for (let p = 0; p < npix(nside); p++) {
        const [t, ph] = pix2ang(nside, p);
        expect(ang2pix(nside, t, ph)).toBe(p);
      }
    }
    const big = 8192;
    for (const p of [0, 1, 7, 2 * big * (big - 1) - 1, 2 * big * (big - 1), 400_000_000, npix(big) - 2 * big * (big - 1), npix(big) - 1]) {
      const [t, ph] = pix2ang(big, p);
      expect(ang2pix(big, t, ph)).toBe(p);
    }
  });

  it('handles the poles, phi wrapping and invalid input', () => {
    expect(ang2pix(64, 0, 0)).toBeLessThan(4);
    expect(ang2pix(64, Math.PI, 0)).toBeGreaterThanOrEqual(npix(64) - 4);
    expect(ang2pix(64, 1.2, 0.3 + 2 * Math.PI)).toBe(ang2pix(64, 1.2, 0.3));
    expect(ang2pix(64, 1.2, 0.3 - 2 * Math.PI)).toBe(ang2pix(64, 1.2, 0.3));
    expect(() => ang2pix(3.5, 1, 1)).toThrow();
    expect(() => pix2ang(1, 12)).toThrow();
    expect(nsideFromNpix(49152)).toBe(64);
    expect(() => nsideFromNpix(100)).toThrow();
  });

  it('Galactic helpers follow the architecture convention', () => {
    const v = lbToVec(90, 0);
    expect(v[0]).toBeCloseTo(0, 15);
    expect(v[1]).toBeCloseTo(1, 15);
    expect(lbToVec(0, 90)[2]).toBeCloseTo(1, 15);
    const [l, b] = vecToLb(lbToVec(263.5, -37.25));
    expect(l).toBeCloseTo(263.5, 12);
    expect(b).toBeCloseTo(-37.25, 12);
    // theta = 90 - b, phi = l
    expect(lb2pix(64, 30, 20)).toBe(ang2pix(64, (70 * Math.PI) / 180, (30 * Math.PI) / 180));
    expect(vec2pix(64, lbToVec(30, 20))).toBe(lb2pix(64, 30, 20));
    const [l2, b2] = pix2lb(1, 0);
    expect(l2).toBeCloseTo(45, 12);
    expect(b2).toBeCloseTo(90 - (Math.acos(2 / 3) * 180) / Math.PI, 12);
  });
});

describe('Mollweide projection', () => {
  it('maps the centre, poles and the edge of the ellipse', () => {
    const [x0, y0] = mollweideForward(0, 0);
    expect(Math.abs(x0)).toBeLessThan(1e-15);
    expect(Math.abs(y0)).toBeLessThan(1e-15);
    expect(mollweideForward(0, 90)[1]).toBeCloseTo(Math.SQRT2, 12);
    const [xe] = mollweideForward(180 - 1e-9, 0);
    expect(Math.abs(xe)).toBeCloseTo(2 * Math.SQRT2, 6);
    // astronomical convention: positive l is drawn to the left
    expect(mollweideForward(30, 0)[0]).toBeLessThan(0);
    expect(mollweideForward(30, 0, { flipLon: false })[0]).toBeGreaterThan(0);
  });

  it('inverse undoes forward and rejects points outside the ellipse', () => {
    for (const [l, b] of [[0, 0], [10, 20], [123.4, -56.7], [270, 89], [359, -10], [200, 45]] as const) {
      const [x, y] = mollweideForward(l, b, { lon0: 0 });
      const back = mollweideInverse(x, y)!;
      expect(back[1]).toBeCloseTo(b, 9);
      expect(((back[0] - l + 540) % 360) - 180).toBeCloseTo(0, 8);
    }
    expect(mollweideInverse(2.9, 0)).toBeNull();
    expect(mollweideInverse(0, 1.5)).toBeNull();
  });
});
