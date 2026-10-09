/**
 * The panel contract shared by the router (main.ts) and every tool in src/panels/.
 *
 * A panel module exports `title` and `mount(el, ctx)`. `mount` fills `el` (the page body under
 * the page header) and may be async; whatever it registers through `ctx.onCleanup` /
 * `ctx.onThemeChange` is released when the user navigates away. See web/README.md.
 */

import type { BundleMap, BundleName, LoadResult, Manifest } from './data';
import type { Theme } from './theme';

export type RouteId =
  | 'overview'
  | 'manifolds'
  | 'lattice'
  | 'ghost-twin'
  | 'resonance'
  | 'sky-stats'
  | 'bounce'
  | 'capacity'
  | 'page-curve'
  | 'qes'
  | 'hayden-preskill'
  | 'gw-horizons'
  | 'claims'
  | 'data';

export type NavGroupId = 'overview' | 'm1' | 'm2' | 'ledger';

export interface RouteInfo {
  id: RouteId;
  /** short code shown in the navigation, e.g. "1.3" */
  code: string;
  /** navigation label */
  label: string;
  group: NavGroupId;
  /** one-line description (overview cards, placeholders) */
  blurb: string;
  /** bundles the tool reads; drives the readiness dot in the navigation */
  bundles: BundleName[];
  /** the tool also computes closed-form physics live in the browser */
  live?: boolean;
}

export interface PanelContext {
  route: RouteInfo;
  /** manifest.json, preloaded at start-up (null until `qchas build` has run) */
  manifest: Manifest | null;
  /** theme at mount time; subscribe with `onThemeChange` for updates */
  theme: Theme;
  /** cached bundle loader (never throws; check `.ok`) */
  load<K extends BundleName>(name: K): Promise<LoadResult<BundleMap[K]>>;
  /** run `fn` when the panel is unmounted (dispose charts, timers, WebGL) */
  onCleanup(fn: () => void): void;
  /** theme subscription released automatically on unmount */
  onThemeChange(cb: (theme: Theme) => void): void;
  /** replace the page-header subtitle */
  setSubtitle(text: string): void;
  /** put badges or small controls at the right of the page header */
  setHeaderExtras(nodes: Node[]): void;
  navigate(id: RouteId): void;
  /** true once the user has navigated away (check after awaits before touching the DOM) */
  readonly unmounted: boolean;
}

export interface PanelModule {
  title: string;
  /** page-header subtitle (defaults to the route blurb) */
  subtitle?: string;
  mount(el: HTMLElement, ctx: PanelContext): void | Promise<void>;
}
