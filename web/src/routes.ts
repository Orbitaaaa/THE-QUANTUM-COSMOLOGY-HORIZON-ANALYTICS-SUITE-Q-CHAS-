/**
 * Route registry: the routes of docs/ARCHITECTURE.md, their navigation group, the bundles each
 * tool reads, and a lazy loader per panel (each panel is its own chunk).
 */

import type { NavGroupId, PanelModule, RouteId, RouteInfo } from './lib/panel';

export interface NavGroup {
  id: NavGroupId;
  code: string;
  title: string;
}

export const NAV_GROUPS: NavGroup[] = [
  { id: 'overview', code: 'Mission', title: 'Overview' },
  { id: 'm1', code: 'Module 1', title: 'Global Topology & All-Sky Survey Analyst' },
  { id: 'm2', code: 'Module 2', title: 'Dressed Horizon & Evaporation Simulator' },
  { id: 'ledger', code: 'Ledger', title: 'Claims & provenance' },
];

export const ROUTES: RouteInfo[] = [
  {
    id: 'overview',
    code: '0.0',
    label: 'Overview',
    group: 'overview',
    blurb: 'Mission status, the two analysis modules, source papers and the onboard physics self-test.',
    bundles: ['manifest.json', 'sky_map.json', 'lattice.json', 'resonance.json', 'cmb.json', 'horizon.json', 'bounce.json', 'gw.json'],
    live: true,
  },
  {
    id: 'manifolds',
    code: '1.1',
    label: 'Manifold Selector',
    group: 'm1',
    blurb: 'The 18 flat 3-manifolds against principles P1-P4: which survive, by which theorem.',
    bundles: ['topology.json', 'cmb.json'],
  },
  {
    id: 'lattice',
    code: '1.2',
    label: 'Lattice Shape Optimizer',
    group: 'm1',
    blurb: 'Standard-Model Casimir energy Z(2) over torus shapes; FCC/BCC minima, SC saddle, shape flow.',
    bundles: ['lattice.json', 'validation.json'],
    live: true,
  },
  {
    id: 'ghost-twin',
    code: '1.3',
    label: 'Ghost-Twin Tracker',
    group: 'm1',
    blurb: 'Matched circles, twin points n + g/χ* and ghost images on the real CMB map, with SPHEREx frame footprints.',
    bundles: ['sky_map.json', 'topology.json', 'cmb.json', 'spherex.json'],
    live: true,
  },
  {
    id: 'resonance',
    code: '1.4',
    label: 'Resonance Shell Inspector',
    group: 'm1',
    blurb: 'Topological resonance on dual-lattice shells: KL gain at known orientation, rank-N Woodbury test, pre-registered search.',
    bundles: ['resonance.json', 'cmb.json', 'validation.json'],
  },
  {
    id: 'sky-stats',
    code: '1.5',
    label: 'Sky Statistics',
    group: 'm1',
    blurb: 'Low-ℓ spectra, the quadrupole and S₁/₂ anomalies, aligned axes and the orientation-scanned torus likelihood.',
    bundles: ['cmb.json', 'validation.json'],
  },
  {
    id: 'bounce',
    code: '1.6',
    label: 'Quantum Big Bang',
    group: 'm1',
    blurb: 'Unitary Wheeler-DeWitt evolution through a = 0: the bounce scale, gravitons, mixmaster.',
    bundles: ['bounce.json'],
  },
  {
    id: 'capacity',
    code: '2.1',
    label: 'Crossed-Product Capacity',
    group: 'm2',
    blurb: 'Type III₁ → II∞ → II₁: the horizon capacity A/4G from an energy window, de Sitter entropy.',
    bundles: ['horizon.json', 'gw.json', 'validation.json'],
    live: true,
  },
  {
    id: 'page-curve',
    code: '2.2',
    label: 'Page Curve',
    group: 'm2',
    blurb: "Hawking's curve versus the Page transfer law S(R) ≤ min{S_H, A/4G}; island and interior volume.",
    bundles: ['horizon.json', 'gw.json', 'validation.json'],
    live: true,
  },
  {
    id: 'qes',
    code: '2.3',
    label: 'QES & Atmosphere',
    group: 'm2',
    blurb: 'The quantum extremal surface of an evaporating hole and the entangled quantum atmosphere.',
    bundles: ['horizon.json', 'validation.json'],
    live: true,
  },
  {
    id: 'hayden-preskill',
    code: '2.4',
    label: 'Hayden-Preskill Recovery',
    group: 'm2',
    blurb: 'Exactly unitary evaporation of qubits and decoding of a diary thrown in early or late.',
    bundles: ['horizon.json', 'validation.json'],
    live: true,
  },
  {
    id: 'gw-horizons',
    code: '2.5',
    label: 'GW Horizons (LIGO)',
    group: 'm2',
    blurb: 'LIGO strain and ringdowns of real mergers; remnant entropy, Hawking temperature, Page time.',
    bundles: ['gw.json', 'horizon.json', 'validation.json'],
  },
  {
    id: 'claims',
    code: 'L.1',
    label: 'Claims Ledger',
    group: 'ledger',
    blurb: 'Every statement with its status as the papers state it: proved, computed, measured or assumed.',
    bundles: ['topology.json', 'horizon.json', 'bounce.json', 'validation.json'],
  },
  {
    id: 'data',
    code: 'L.2',
    label: 'Data & Provenance',
    group: 'ledger',
    blurb: 'Datasets behind every number: archive URL, size, SHA-256, retrieval time and build warnings.',
    bundles: ['manifest.json'],
  },
];

export const DEFAULT_ROUTE: RouteId = 'overview';

export const PANEL_LOADERS: Record<RouteId, () => Promise<PanelModule>> = {
  overview: () => import('./panels/overview'),
  manifolds: () => import('./panels/manifolds'),
  lattice: () => import('./panels/lattice'),
  'ghost-twin': () => import('./panels/ghost-twin'),
  resonance: () => import('./panels/resonance'),
  'sky-stats': () => import('./panels/sky-stats'),
  bounce: () => import('./panels/bounce'),
  capacity: () => import('./panels/capacity'),
  'page-curve': () => import('./panels/page-curve'),
  qes: () => import('./panels/qes'),
  'hayden-preskill': () => import('./panels/hayden-preskill'),
  'gw-horizons': () => import('./panels/gw-horizons'),
  claims: () => import('./panels/claims'),
  data: () => import('./panels/data'),
};

export const routeById = (id: string): RouteInfo | undefined => ROUTES.find((r) => r.id === id);

// --- show vs console -----------------------------------------------------------------------

/**
 * Plain anchors owned by the judge-facing show (docs/SHOW.md "Routing"). Artifact hosts pass
 * only plain anchors, so the show never uses `#/`. `breakthrough` is Act 3A (black hole) and
 * `bounce` Act 3B (Big Bang); `blackhole` and `breakthrough-bounce` are accepted aliases.
 */
export const SHOW_ANCHORS = ['discovery', 'experiment', 'breakthrough', 'bounce', 'evidence'] as const;
export type ShowAnchor = (typeof SHOW_ANCHORS)[number];

const SHOW_ALIASES: Record<string, ShowAnchor> = {
  blackhole: 'breakthrough',
  'black-hole': 'breakthrough',
  'breakthrough-blackhole': 'breakthrough',
  'breakthrough-bounce': 'bounce',
};

export type AppView = { mode: 'show'; anchor: ShowAnchor | null; at?: number } | { mode: 'console' };

/**
 * Which top-level view a hash selects. `#/...` is always the mission console; an empty hash is
 * the show's landing; a show anchor opens that act (`#discovery-12` opens it paused at 12 s, for
 * QA screenshots); a bare console id (`#qes`, tolerated by `parseHash`) still opens the console;
 * anything else falls back to the show's landing.
 */
export function viewForHash(hash: string): AppView {
  const raw = hash.replace(/^#/, '');
  if (raw.startsWith('/')) return { mode: 'console' };
  let key = (raw.split(/[/?&]/)[0] ?? '').toLowerCase();
  if (!key) return { mode: 'show', anchor: null };
  let at: number | undefined;
  const timed = /^([a-z-]+?)-(\d+(?:\.\d+)?)$/.exec(key);
  if (timed && (SHOW_ANCHORS.includes(timed[1] as ShowAnchor) || SHOW_ALIASES[timed[1]!])) {
    key = timed[1]!;
    at = Number(timed[2]);
  }
  const anchor = (SHOW_ANCHORS as readonly string[]).includes(key) ? (key as ShowAnchor) : SHOW_ALIASES[key];
  if (anchor) return at === undefined ? { mode: 'show', anchor } : { mode: 'show', anchor, at };
  if (routeById(key)) return { mode: 'console' };
  return { mode: 'show', anchor: null };
}

/** Parse `#/route` (also tolerates `#route` and trailing segments such as `#/data/warnings`). */
export function parseHash(hash: string): { id: RouteId; rest: string[] } {
  const parts = hash.replace(/^#\/?/, '').split('/').filter(Boolean);
  const r = routeById(parts[0] ?? '');
  return { id: r ? r.id : DEFAULT_ROUTE, rest: r ? parts.slice(1) : [] };
}
