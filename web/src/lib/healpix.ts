/**
 * HEALPix RING-scheme pixelization (Gorski et al. 2005, ApJ 622, 759) for any nside, Galactic
 * coordinate helpers and the Mollweide projection.
 *
 * Conventions follow docs/ARCHITECTURE.md: maps are RING ordered in Galactic coordinates, so
 * theta = 90 deg - b and phi = l (radians); the unit vector of (l, b) is
 * (cos b cos l, cos b sin l, sin b).
 */

const TWO_PI = 2 * Math.PI;
const HALF_PI = Math.PI / 2;
const DEG = Math.PI / 180;

export const npix = (nside: number): number => 12 * nside * nside;

export const nsideFromNpix = (n: number): number => {
  const nside = Math.round(Math.sqrt(n / 12));
  if (12 * nside * nside !== n) throw new Error(`${n} is not a valid HEALPix pixel count`);
  return nside;
};

function checkNside(nside: number): void {
  if (!Number.isInteger(nside) || nside < 1 || nside > 1 << 20) throw new Error(`invalid nside ${nside}`);
}

const isqrt = (x: number): number => {
  let r = Math.floor(Math.sqrt(x));
  while (r * r > x) r--;
  while ((r + 1) * (r + 1) <= x) r++;
  return r;
};

const mod = (a: number, n: number): number => ((a % n) + n) % n;

/** RING pixel index containing the direction (theta, phi) in radians (theta colatitude). */
export function ang2pix(nside: number, theta: number, phi: number): number {
  checkNside(nside);
  if (!(theta >= 0 && theta <= Math.PI)) throw new Error(`theta out of range: ${theta}`);
  const z = Math.cos(theta);
  const za = Math.abs(z);
  const tt = mod(phi, TWO_PI) / HALF_PI; // in [0, 4)
  if (za <= 2 / 3) {
    const temp1 = nside * (0.5 + tt);
    const temp2 = nside * z * 0.75;
    const jp = Math.floor(temp1 - temp2);
    const jm = Math.floor(temp1 + temp2);
    const ir = nside + 1 + jp - jm; // ring number counted from z = 2/3, in 1..2 nside + 1
    const kshift = 1 - (ir & 1);
    const ip = mod(Math.floor((jp + jm - nside + kshift + 1) / 2), 4 * nside);
    return 2 * nside * (nside - 1) + (ir - 1) * 4 * nside + ip;
  }
  const tp = tt - Math.floor(tt);
  // 1 - |z| = 2 sin^2(theta'/2), accurate near the poles
  const thetaP = z > 0 ? theta : Math.PI - theta;
  const s = Math.sin(thetaP / 2);
  const tmp = nside * Math.sqrt(6) * s;
  const jp = Math.floor(tp * tmp);
  const jm = Math.floor((1 - tp) * tmp);
  const ir = jp + jm + 1; // ring number counted from the closest pole
  const ip = mod(Math.floor(tt * ir), 4 * ir);
  return z > 0 ? 2 * ir * (ir - 1) + ip : npix(nside) - 2 * ir * (ir + 1) + ip;
}

/** Centre (theta, phi) in radians of a RING pixel. */
export function pix2ang(nside: number, pix: number): [number, number] {
  checkNside(nside);
  const np = npix(nside);
  if (!(Number.isInteger(pix) && pix >= 0 && pix < np)) throw new Error(`pixel ${pix} out of range`);
  const ncap = 2 * nside * (nside - 1);
  if (pix < ncap) {
    const iring = (1 + isqrt(1 + 2 * pix)) >> 1;
    const iphi = pix + 1 - 2 * iring * (iring - 1);
    const theta = 2 * Math.asin(iring / (Math.sqrt(6) * nside));
    return [theta, (iphi - 0.5) * (HALF_PI / iring)];
  }
  if (pix < np - ncap) {
    const ip = pix - ncap;
    const iring = Math.floor(ip / (4 * nside)) + nside;
    const iphi = (ip % (4 * nside)) + 1;
    const fodd = (iring + nside) & 1 ? 1 : 0.5;
    const z = ((2 * nside - iring) * 2) / (3 * nside);
    return [Math.acos(z), (iphi - fodd) * (Math.PI / (2 * nside))];
  }
  const ip = np - pix;
  const iring = (1 + isqrt(2 * ip - 1)) >> 1;
  const iphi = 4 * iring + 1 - (ip - 2 * iring * (iring - 1));
  const theta = Math.PI - 2 * Math.asin(iring / (Math.sqrt(6) * nside));
  return [theta, (iphi - 0.5) * (HALF_PI / iring)];
}

/** RING pixel of a Galactic direction (l, b) in degrees. */
export const lb2pix = (nside: number, l: number, b: number): number =>
  ang2pix(nside, Math.min(Math.PI, Math.max(0, (90 - b) * DEG)), l * DEG);

/** Galactic (l, b) in degrees of a RING pixel centre, l in [0, 360). */
export function pix2lb(nside: number, pix: number): [number, number] {
  const [theta, phi] = pix2ang(nside, pix);
  return [mod(phi / DEG, 360), 90 - theta / DEG];
}

/** Unit vector (cos b cos l, cos b sin l, sin b) of Galactic (l, b) in degrees. */
export function lbToVec(l: number, b: number): [number, number, number] {
  const lr = l * DEG;
  const br = b * DEG;
  return [Math.cos(br) * Math.cos(lr), Math.cos(br) * Math.sin(lr), Math.sin(br)];
}

/** Galactic (l, b) in degrees of any nonzero 3-vector, l in [0, 360). */
export function vecToLb(v: readonly number[]): [number, number] {
  const n = Math.hypot(v[0]!, v[1]!, v[2]!);
  const l = mod(Math.atan2(v[1]!, v[0]!) / DEG, 360);
  return [l, Math.asin(Math.max(-1, Math.min(1, v[2]! / n))) / DEG];
}

/** RING pixel containing a 3-vector direction. */
export function vec2pix(nside: number, v: readonly number[]): number {
  const [l, b] = vecToLb(v);
  return lb2pix(nside, l, b);
}

// --- Mollweide ------------------------------------------------------------------------

export interface MollweideOptions {
  /** central longitude in degrees (default 0) */
  lon0?: number;
  /** astronomical convention: longitude increases to the left (default true) */
  flipLon?: boolean;
}

/** Auxiliary angle t with 2t + sin 2t = pi sin(lat), by Newton iteration. */
function mollweideTheta(lat: number): number {
  if (Math.abs(Math.abs(lat) - HALF_PI) < 1e-12) return Math.sign(lat) * HALF_PI;
  const target = Math.PI * Math.sin(lat);
  let t = lat;
  for (let i = 0; i < 50; i++) {
    const f = 2 * t + Math.sin(2 * t) - target;
    const d = 2 + 2 * Math.cos(2 * t);
    if (d === 0) break;
    const step = f / d;
    t -= step;
    if (Math.abs(step) < 1e-14) break;
  }
  return t;
}

/**
 * Forward Mollweide projection of Galactic (l, b) degrees to (x, y) with x in [-2 sqrt 2, 2 sqrt 2],
 * y in [-sqrt 2, sqrt 2] (y up).
 */
export function mollweideForward(l: number, b: number, opts: MollweideOptions = {}): [number, number] {
  const { lon0 = 0, flipLon = true } = opts;
  let lam = mod(l - lon0 + 180, 360) - 180;
  lam *= DEG;
  const t = mollweideTheta(b * DEG);
  const x = ((2 * Math.SQRT2) / Math.PI) * lam * Math.cos(t);
  return [flipLon ? -x : x, Math.SQRT2 * Math.sin(t)];
}

/** Inverse Mollweide: (x, y) -> Galactic (l, b) degrees, or null outside the ellipse. */
export function mollweideInverse(x: number, y: number, opts: MollweideOptions = {}): [number, number] | null {
  const { lon0 = 0, flipLon = true } = opts;
  if ((x * x) / 8 + (y * y) / 2 > 1 + 1e-12) return null;
  const s = Math.max(-1, Math.min(1, y / Math.SQRT2));
  const t = Math.asin(s);
  const lat = Math.asin(Math.max(-1, Math.min(1, (2 * t + Math.sin(2 * t)) / Math.PI)));
  const ct = Math.cos(t);
  const xx = flipLon ? -x : x;
  const lam = ct < 1e-12 ? 0 : (Math.PI * xx) / (2 * Math.SQRT2 * ct);
  if (Math.abs(lam) > Math.PI + 1e-9) return null;
  return [mod(lam / DEG + lon0, 360), lat / DEG];
}
