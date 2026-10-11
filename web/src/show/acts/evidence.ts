/**
 * Act IV · Evidence (docs/SHOW.md): a calm, information-dense landing in front of the mission
 * console. Live numbers from validation.json and manifest.json; doors into every console tool.
 */

import * as THREE from 'three';
import type { Manifest, ValidationBundle } from '../../lib/data';
import { h } from '../../lib/ui';
import { Act, type ActContext, type ActModule } from '../act';
import { skyDome } from '../gfx';
import { EVIDENCE } from '../script';

class Evidence extends Act {
  readonly id = 'evidence' as const;
  readonly duration = 5;
  private readonly sky: THREE.Mesh;

  constructor(ctx: ActContext, validation: ValidationBundle | null, manifest: Manifest | null) {
    super(ctx);
    this.camera.fov = 60;
    this.bloom = { strength: 0.5, radius: 0.4, threshold: 0.7 };
    this.sky = skyDome(1500, 0.35);
    this.scene.add(this.sky);
    ctx.hud.append(evidencePage(validation, manifest));
  }

  protected update(): void {
    const yaw = 0.2 + 0.01 * this.wall;
    this.camera.position.set(0, 0, 0);
    this.camera.lookAt(Math.sin(yaw), 0.15, Math.cos(yaw));
  }
}

function evidencePage(v: ValidationBundle | null, m: Manifest | null): HTMLElement {
  const match = v?.summary.match ?? 169;
  const total = v?.summary.total ?? 173;
  const groups = new Map<string, { n: number; ok: number }>();
  for (const r of v?.rows ?? []) {
    const g = groups.get(r.group) ?? { n: 0, ok: 0 };
    g.n++;
    if (r.verdict === 'match') g.ok++;
    groups.set(r.group, g);
  }
  const map = m?.cmb_map_source ?? 'WMAP9-ILC';
  const stat = (big: string, label: string, cls = ''): HTMLElement => h('div', { class: `ev-stat ${cls}`.trim() }, h('div', { class: 'ev-big' }, big), h('div', { class: 'ev-lab' }, label));
  return h(
    'div',
    { class: 'evidence show-ui' },
    h(
      'div',
      { class: 'ev-inner' },
      h('div', { class: 'cap-kicker' }, 'Act IV · Evidence'),
      h('h2', { class: 'ev-title' }, 'The evidence behind the show'),
      h('p', { class: 'ev-lead' }, EVIDENCE.lead),
      h(
        'div',
        { class: 'ev-stats' },
        stat(`${match} / ${total}`, 'paper numbers reproduced from scratch', 'hero'),
        stat(String(total - match), 'documented differences (grid bias in a printed table)'),
        stat('p = 0.79', `resonance search on ${map} (paper on Planck: 0.78)`),
        stat('0', 'detections claimed — the data show no torus'),
      ),
      h(
        'div',
        { class: 'ev-grid' },
        h('section', { class: 'ev-card' }, h('h3', null, 'Real data'), h('ul', { class: 'ev-list' }, ...EVIDENCE.data.map((d) => h('li', null, h('b', null, d.name), h('span', null, d.what), h('i', null, d.where))))),
        h('section', { class: 'ev-card' }, h('h3', null, 'Methods'), h('ul', { class: 'ev-list plain' }, ...EVIDENCE.methods.map((x) => h('li', null, x)))),
        h(
          'section',
          { class: 'ev-card wide' },
          h('h3', null, 'Results'),
          h('div', { class: 'ev-results' }, ...EVIDENCE.results.map((r) => h('div', { class: 'ev-res' }, h('div', { class: 'ev-rk' }, r.k), h('div', { class: 'ev-rv' }, r.v), h('div', { class: 'ev-rd' }, r.d)))),
        ),
        h(
          'section',
          { class: 'ev-card' },
          h('h3', null, 'Reproduction ledger'),
          h(
            'div',
            { class: 'ev-ledger' },
            ...[...groups.entries()].map(([g, c]) =>
              h('div', { class: 'ev-lrow' }, h('span', { class: 'ev-lg' }, g), h('span', { class: 'ev-lbar' }, h('span', { style: { width: `${(100 * c.ok) / Math.max(1, c.n)}%` } })), h('span', { class: 'mono' }, `${c.ok}/${c.n}`)),
            ),
          ),
        ),
        h('section', { class: 'ev-card' }, h('h3', null, 'Limits, stated plainly'), h('ul', { class: 'ev-list plain' }, ...EVIDENCE.limits.map((x) => h('li', null, x)))),
      ),
      h('h3', { class: 'ev-doors-title' }, 'Open the mission console'),
      h('div', { class: 'ev-doors' }, ...EVIDENCE.doors.map((d) => h('a', { class: 'ev-door', href: `#/${d.route}` }, h('span', { class: 'mono' }, d.code), h('span', null, d.label), h('span', { 'aria-hidden': 'true' }, '→')))),
      h('div', { class: 'ev-foot' }, 'Q-CHAS · NASA Space Apps Challenge · open source · data: NASA LAMBDA, NASA/IPAC IRSA, ESA PLA, GWOSC'),
    ),
  );
}

export const evidence: ActModule = {
  async mount(ctx: ActContext) {
    const [v, m] = await Promise.all([ctx.data.validation(), ctx.data.manifest()]);
    if (!ctx.stage) {
      // no WebGL: render the page alone
      ctx.hud.replaceChildren(evidencePage(v, m));
      return { id: 'evidence', duration: 0, t: 0, playing: false, onEnd: null } as never;
    }
    return new Evidence(ctx, v, m);
  },
};
