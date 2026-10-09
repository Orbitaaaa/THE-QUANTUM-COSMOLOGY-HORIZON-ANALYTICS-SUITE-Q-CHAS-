import { afterEach, describe, expect, it } from 'vitest';
import { BUNDLE_NAMES, clearCache, loadBinary, loadBundle, loadSkyMap, setFetcher } from '../src/lib/data';

type Route = { status?: number; body: string | ArrayBuffer; type?: string };

function mockFetch(routes: Record<string, Route>, calls: string[] = []): typeof fetch {
  return (async (input: RequestInfo | URL) => {
    const url = String(input);
    calls.push(url);
    const key = Object.keys(routes).find((k) => url.endsWith(k));
    if (!key) return new Response('<!doctype html><html></html>', { status: 200, headers: { 'content-type': 'text/html' } });
    const r = routes[key]!;
    return new Response(r.body, { status: r.status ?? 200, headers: { 'content-type': r.type ?? 'application/json' } });
  }) as typeof fetch;
}

afterEach(() => clearCache());

describe('bundle loaders', () => {
  it('lists every bundle of the data contract', () => {
    expect(BUNDLE_NAMES).toEqual([
      'manifest.json', 'topology.json', 'lattice.json', 'sky_map.json', 'cmb.json',
      'resonance.json', 'spherex.json', 'horizon.json', 'bounce.json', 'gw.json',
    ]);
  });

  it('reports a 404 as missing, never as data', async () => {
    setFetcher(mockFetch({ 'cmb.json': { status: 404, body: 'nope', type: 'text/plain' } }));
    const r = await loadBundle('cmb.json');
    expect(r.ok).toBe(false);
    if (!r.ok) expect(r.reason).toBe('missing');
  });

  it('treats a dev-server SPA fallback (index.html) as missing', async () => {
    setFetcher(mockFetch({}));
    const r = await loadBundle('gw.json');
    expect(r).toMatchObject({ ok: false, reason: 'missing' });
  });

  it('parses, validates and caches a bundle', async () => {
    const calls: string[] = [];
    const manifest = { generated_utc: '2026-10-09T12:00:00+00:00', qchas_version: '0.1.0', cmb_map_source: 'WMAP9-ILC', datasets: {}, bundles: [], warnings: [] };
    setFetcher(mockFetch({ 'manifest.json': { body: JSON.stringify(manifest) } }, calls));
    const a = await loadBundle('manifest.json');
    const b = await loadBundle('manifest.json');
    expect(a.ok && a.data.cmb_map_source).toBe('WMAP9-ILC');
    expect(b).toBe(a);
    expect(calls.filter((c) => c.endsWith('data/manifest.json'))).toHaveLength(1);
  });

  it('rejects malformed JSON and bundles lacking required keys', async () => {
    setFetcher(mockFetch({ 'topology.json': { body: '{"manifolds": [' }, 'horizon.json': { body: '{"qes": {}}' } }));
    expect(await loadBundle('topology.json')).toMatchObject({ ok: false, reason: 'invalid' });
    expect(await loadBundle('horizon.json')).toMatchObject({ ok: false, reason: 'invalid' });
  });

  it('decodes little-endian Float32 maps and checks their length', async () => {
    const nside = 1;
    const values = Float32Array.from({ length: 12 }, (_, i) => i - 5.5);
    const buf = new ArrayBuffer(48);
    const dv = new DataView(buf);
    values.forEach((v, i) => dv.setFloat32(4 * i, v, true));
    const mask = Uint8Array.from({ length: 12 }, (_, i) => (i % 2 ? 1 : 0));
    const meta = {
      nside, ordering: 'RING', frame: 'galactic', unit: 'uK', source: 'WMAP9-ILC', mask_file: 'sky_mask.bin',
      maps: [{ name: 'full', file: 'sky_map_full.bin', description: 'test', min: -5.5, max: 5.5 }],
    };
    setFetcher(
      mockFetch({
        'sky_map.json': { body: JSON.stringify(meta) },
        'sky_map_full.bin': { body: buf, type: 'application/octet-stream' },
        'sky_mask.bin': { body: mask.buffer, type: 'application/octet-stream' },
        'short.bin': { body: new ArrayBuffer(40), type: 'application/octet-stream' },
      }),
    );
    const r = await loadSkyMap();
    expect(r.ok).toBe(true);
    if (r.ok) {
      expect(Array.from(r.data.map)).toEqual(Array.from(values));
      expect(Array.from(r.data.mask!)).toEqual(Array.from(mask));
    }
    expect(await loadBinary('short.bin', 'f32', 12)).toMatchObject({ ok: false, reason: 'invalid' });
    expect(await loadBinary('absent.bin', 'f32')).toMatchObject({ ok: false, reason: 'missing' });
  });
});
