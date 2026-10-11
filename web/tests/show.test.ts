/**
 * The judge-facing show (src/show): routing between show and console, the exact physics the acts
 * draw, and the honesty rules of docs/SHOW.md applied to every line of copy.
 */

import { describe, expect, it } from 'vitest';
import { viewForHash } from '../src/routes';
import { bounceDensity } from '../src/show/acts/breakthrough-bounce';
import { ACTS, EVIDENCE, type Tag } from '../src/show/script';
import { hawkingEntropy, pageCurve, pageTimeFraction } from '../src/lib/physics/page';

describe('show routing', () => {
  it('lands on the show and keeps #/ routes for the console', () => {
    expect(viewForHash('')).toEqual({ mode: 'show', anchor: null });
    expect(viewForHash('#')).toEqual({ mode: 'show', anchor: null });
    expect(viewForHash('#/overview')).toEqual({ mode: 'console' });
    expect(viewForHash('#/bounce')).toEqual({ mode: 'console' });
    expect(viewForHash('#qes')).toEqual({ mode: 'console' });
  });

  it('maps plain anchors, aliases and timed QA anchors to acts', () => {
    expect(viewForHash('#discovery')).toEqual({ mode: 'show', anchor: 'discovery' });
    expect(viewForHash('#bounce')).toEqual({ mode: 'show', anchor: 'bounce' });
    expect(viewForHash('#blackhole')).toEqual({ mode: 'show', anchor: 'breakthrough' });
    expect(viewForHash('#experiment-6.5')).toEqual({ mode: 'show', anchor: 'experiment', at: 6.5 });
    expect(viewForHash('#nonsense')).toEqual({ mode: 'show', anchor: null });
  });
});

describe('exact bounce density', () => {
  const integrate = (T: number): number => {
    let s = 0;
    const n = 40000;
    const L = 6;
    for (let i = 1; i < n; i++) s += bounceDensity((i * L) / n, T) * (L / n);
    return s;
  };

  it('is normalised at every radiation time (unitary evolution)', () => {
    for (const T of [0, 0.7, Math.PI / 2 - 0.05, Math.PI / 2, 2.4]) expect(integrate(T)).toBeCloseTo(1, 4);
  });

  it('bounces at <a> = 0.165 l for a_max = 2 (bounce.json closed form)', () => {
    let m = 0;
    const n = 40000;
    for (let i = 1; i < n; i++) {
      const a = (i * 6) / n;
      m += a * bounceDensity(a, Math.PI / 2) * (6 / n);
    }
    expect(m).toBeCloseTo(0.16462, 4);
  });
});

describe('black-hole threads follow the Page curve', () => {
  it('re-routes no thread before the Page time and all of them at the end', () => {
    const rerouted = (t: number): number => hawkingEntropy(t) - pageCurve(t);
    expect(rerouted(0.3)).toBe(0);
    expect(rerouted(pageTimeFraction() - 1e-6)).toBeCloseTo(0, 6);
    expect(rerouted(1)).toBeCloseTo(hawkingEntropy(1), 12);
    expect(pageTimeFraction()).toBeCloseTo(0.476, 3);
  });
});

describe('show copy honours the binding facts', () => {
  const all = [
    ...ACTS.flatMap((a) => [a.title, ...a.beats.flatMap((b) => [b.title, b.body, ...(b.chips ?? []).map((c) => `${c.text} ${c.more ?? ''}`)]), ...a.truth.map((c) => `${c.text} ${c.more ?? ''}`)]),
    ...EVIDENCE.results.map((r) => `${r.k} ${r.v} ${r.d}`),
    ...EVIDENCE.limits,
  ].join('\n');

  it('never claims a detection', () => {
    expect(all).not.toMatch(/\b(we|Q-CHAS) (detect|discover)|torus (is|was) (detected|found)/i);
    expect(all).toMatch(/no torus/i);
    expect(all).toContain('p = 0.78');
    expect(all).toContain('0.79');
  });

  it('states the size bound and the orientation caveat of the 64x gain', () => {
    expect(all).toContain('0.985');
    expect(all).toMatch(/orientation/);
    expect(all).toMatch(/65×/);
    expect(all).toMatch(/multi-field primordial line mapped in 3D/);
  });

  it('tags every chip with a known honesty tag and labels the illustrative size', () => {
    const tags = new Set<Tag>(['ESTABLISHED', 'PROVED', 'COMPUTED', 'HYPOTHESIS', 'ILLUSTRATIVE', 'DATA']);
    for (const a of ACTS) for (const c of [...a.truth, ...a.beats.flatMap((b) => b.chips ?? [])]) expect(tags.has(c.tag)).toBe(true);
    const discovery = ACTS.find((a) => a.id === 'discovery')!;
    expect(discovery.beats.flatMap((b) => b.chips ?? []).some((c) => c.tag === 'ILLUSTRATIVE' && /excluded by data/.test(c.text))).toBe(true);
  });

  it('keeps the guided tour near one minute', () => {
    const total = ACTS.reduce((s, a) => s + a.duration, 0);
    expect(total).toBeGreaterThanOrEqual(45);
    expect(total).toBeLessThanOrEqual(70);
  });
});
