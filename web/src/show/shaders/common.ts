/** GLSL chunks shared by the acts of the show. */

/** Galactic (x, y, z) <-> three.js (x, z, -y): Galactic north is "up". */
export const GALACTIC_GLSL = /* glsl */ `
const float PI = 3.141592653589793;
vec3 threeToGal(vec3 p) { return vec3(p.x, -p.z, p.y); }
vec2 galUv(vec3 g) {
  g = normalize(g);
  float l = atan(g.y, g.x);
  if (l < 0.0) l += 2.0 * PI;
  float b = asin(clamp(g.z, -1.0, 1.0));
  return vec2(l / (2.0 * PI), (b + 0.5 * PI) / PI);
}`;

/**
 * The show's CMB palette: deep-space midpoint, --cmb-cold #2f6bff and --cmb-hot #ff5a36 at one
 * display range, running to ice-white and gold-white at the extremes. Linear-light colours
 * above 1 bloom.
 */
export const CMB_COLOR_GLSL = /* glsl */ `
vec3 cmbColor(float v) {
  float s = clamp(abs(v), 0.0, 1.4);
  vec3 mid = vec3(0.006, 0.009, 0.024);
  vec3 cold = vec3(0.020, 0.110, 0.85);
  vec3 coldX = vec3(0.30, 0.70, 1.45);
  vec3 hot = vec3(0.95, 0.085, 0.025);
  vec3 hotX = vec3(1.50, 0.80, 0.25);
  vec3 a = v < 0.0 ? cold : hot;
  vec3 b = v < 0.0 ? coldX : hotX;
  vec3 c = mix(mid, a, smoothstep(0.0, 0.75, s));
  return mix(c, b, smoothstep(0.7, 1.35, s));
}`;

/** Hash and value noise. */
export const NOISE_GLSL = /* glsl */ `
float hash13(vec3 p) { p = fract(p * 0.1031); p += dot(p, p.zyx + 31.32); return fract((p.x + p.y) * p.z); }
vec3 hash33(vec3 p) {
  p = fract(p * vec3(0.1031, 0.1030, 0.0973));
  p += dot(p, p.yxz + 33.33);
  return fract((p.xxy + p.yxx) * p.zyx);
}
float vnoise(vec3 p) {
  vec3 i = floor(p); vec3 f = fract(p);
  f = f * f * (3.0 - 2.0 * f);
  float n000 = hash13(i), n100 = hash13(i + vec3(1,0,0)), n010 = hash13(i + vec3(0,1,0)), n110 = hash13(i + vec3(1,1,0));
  float n001 = hash13(i + vec3(0,0,1)), n101 = hash13(i + vec3(1,0,1)), n011 = hash13(i + vec3(0,1,1)), n111 = hash13(i + vec3(1,1,1));
  return mix(mix(mix(n000, n100, f.x), mix(n010, n110, f.x), f.y), mix(mix(n001, n101, f.x), mix(n011, n111, f.x), f.y), f.z);
}
float fbm(vec3 p) { float a = 0.5, s = 0.0; for (int i = 0; i < 5; i++) { s += a * vnoise(p); p *= 2.03; a *= 0.5; } return s; }`;

/**
 * A procedural deep sky for a direction d: stars from a hashed cell grid at two scales, plus a
 * dusty galactic band around the plane perpendicular to `bandAxis`.
 */
export const SKY_GLSL = /* glsl */ `
vec3 starLayer(vec3 d, float scale, float density, float bright, float sharp) {
  vec3 p = d * scale;
  vec3 cell = floor(p);
  vec3 col = vec3(0.0);
  for (int i = -1; i <= 1; i++) for (int j = -1; j <= 1; j++) for (int k = -1; k <= 1; k++) {
    vec3 c = cell + vec3(float(i), float(j), float(k));
    vec3 h = hash33(c);
    if (h.x > density) continue;
    vec3 sp = normalize(c + h);
    float dist = length(d - sp) * scale;
    float temp = hash13(c + 7.1);
    vec3 tint = mix(vec3(0.62, 0.74, 1.0), vec3(1.0, 0.82, 0.62), temp);
    float mag = pow(hash13(c + 3.7), 6.0);
    col += tint * bright * (0.25 + 4.0 * mag) * exp(-dist * dist * sharp);
  }
  return col;
}
vec3 deepSky(vec3 d, vec3 bandAxis) {
  d = normalize(d);
  float lat = dot(d, normalize(bandAxis));
  float band = exp(-lat * lat * 26.0);
  float dust = fbm(d * 3.1 + 4.0);
  float lanes = smoothstep(0.42, 0.75, fbm(d * 7.0 + 11.0));
  vec3 glow = vec3(0.040, 0.034, 0.055) * band * (0.25 + 0.9 * dust) * (1.0 - 0.8 * lanes * band);
  glow += vec3(0.020, 0.008, 0.036) * smoothstep(0.55, 0.9, fbm(d * 2.0 + 1.3));
  vec3 stars = starLayer(d, 200.0, 0.05 + 0.22 * band, 0.55, 110.0) + starLayer(d, 60.0, 0.05, 1.3, 420.0);
  return glow + stars;
}`;
