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
    bundles: ['manifest.json'],
    live: true,
  },
  {
    id: 'manifolds',
    code: '1.1',
    label: 'Manifold Selector',
    group: 'm1',
    blurb: 'The 18 flat 3-manifolds against principles P1-P4: which survive, by which theorem.',
    bundles: ['topology.json'],
  },
  {
    id: 'lattice',
    code: '1.2',
    label: 'Lattice Shape Optimizer',
    group: 'm1',
    blurb: 'Standard-Model Casimir energy Z(2) over torus shapes; FCC/BCC minima, SC saddle, shape flow.',
    bundles: ['lattice.json'],
    live: true,
  },
  {
    id: 'ghost-twin',
    code: '1.3',
    label: 'Ghost-Twin Tracker',
    group: 'm1',
    blurb: 'Matched circles, twin points n + g/χ* and ghost images of a source on the real sky map.',
    bundles: ['sky_map.json', 'cmb.json'],
    live: true,
  },
  {
    id: 'resonance',
    code: '1.4',
    label: 'Resonance Shell Inspector',
    group: 'm1',
    blurb: 'Topological resonance on dual-lattice shells: KL gain, rank-N Woodbury test, Planck search.',
    bundles: ['resonance.json'],
  },
  {
    id: 'sky-stats',
    code: '1.5',
    label: 'Sky Statistics',
    group: 'm1',
    blurb: 'Low-ℓ spectra, S1/2, KL information, the law-of-limitation ceiling and survey footprints.',
    bundles: ['cmb.json', 'sky_map.json', 'spherex.json'],
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
    bundles: ['horizon.json'],
    live: true,
  },
  {
    id: 'page-curve',
    code: '2.2',
    label: 'Page Curve',
    group: 'm2',
    blurb: "Hawking's curve versus the Page transfer law S(R) ≤ min{S_H, A/4G}; island and interior volume.",
    bundles: ['horizon.json'],
    live: true,
  },
  {
    id: 'qes',
    code: '2.3',
    label: 'QES & Atmosphere',
    group: 'm2',
    blurb: 'The quantum extremal surface of an evaporating hole and the entangled quantum atmosphere.',
    bundles: ['horizon.json'],
    live: true,
  },
  {
    id: 'hayden-preskill',
    code: '2.4',
    label: 'Hayden-Preskill Recovery',
    group: 'm2',
    blurb: 'Exactly unitary evaporation of qubits and decoding of a diary thrown in early or late.',
    bundles: ['horizon.json'],
    live: true,
  },
  {
    id: 'gw-horizons',
    code: '2.5',
    label: 'GW Horizons (LIGO)',
    group: 'm2',
    blurb: 'LIGO strain and ringdowns of real mergers; remnant entropy, Hawking temperature, Page time.',
    bundles: ['gw.json'],
  },
  {
    id: 'claims',
    code: 'L.1',
    label: 'Claims Ledger',
    group: 'ledger',
    blurb: 'Every statement with its status as the papers state it: proved, computed, measured or assumed.',
    bundles: ['topology.json', 'horizon.json'],
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

/** Parse `#/route` (also tolerates `#route` and trailing segments such as `#/data/warnings`). */
export function parseHash(hash: string): { id: RouteId; rest: string[] } {
  const parts = hash.replace(/^#\/?/, '').split('/').filter(Boolean);
  const r = routeById(parts[0] ?? '');
  return { id: r ? r.id : DEFAULT_ROUTE, rest: r ? parts.slice(1) : [] };
}
