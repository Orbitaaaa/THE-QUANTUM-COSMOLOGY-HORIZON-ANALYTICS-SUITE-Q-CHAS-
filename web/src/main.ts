/**
 * Q-CHAS mission console: application shell, hash router, navigation, status bar, theme.
 *
 * The shell preloads manifest.json once, then mounts one panel per route into the page body.
 * Leaving a route runs every cleanup the panel registered (charts, WebGL, timers).
 */

import './styles.css';
import { BUNDLE_NAMES, type Manifest, loadBundle } from './lib/data';
import type { PanelContext, RouteId, RouteInfo } from './lib/panel';
import { applyTheme, currentTheme, onThemeChange, storedTheme } from './lib/theme';
import { clear, fmtUtc, h, icon } from './lib/ui';
import { NAV_GROUPS, PANEL_LOADERS, ROUTES, parseHash, viewForHash } from './routes';
import type { ShowController } from './show/show';

applyTheme(storedTheme(), false);

const BUILT_BUNDLES = BUNDLE_NAMES.filter((n) => n !== 'manifest.json');

// --- shell ---------------------------------------------------------------------------------

function logo(): SVGSVGElement {
  const ns = 'http://www.w3.org/2000/svg';
  const svg = document.createElementNS(ns, 'svg');
  svg.setAttribute('viewBox', '0 0 32 32');
  svg.setAttribute('class', 'brand-mark');
  svg.setAttribute('aria-hidden', 'true');
  const parts: [string, Record<string, string>][] = [
    ['circle', { cx: '16', cy: '16', r: '13', fill: 'none', stroke: 'currentColor', 'stroke-width': '1.5' }],
    ['ellipse', { cx: '16', cy: '16', rx: '13', ry: '5', fill: 'none', stroke: 'currentColor', 'stroke-width': '1.2', opacity: '0.7' }],
    ['ellipse', { cx: '16', cy: '16', rx: '5', ry: '13', fill: 'none', stroke: 'currentColor', 'stroke-width': '1.2', opacity: '0.7' }],
    ['circle', { cx: '16', cy: '16', r: '2.4', fill: 'currentColor' }],
    ['circle', { cx: '25.5', cy: '9', r: '1.8', fill: 'currentColor' }],
  ];
  for (const [tag, attrs] of parts) {
    const el = document.createElementNS(ns, tag);
    for (const [k, v] of Object.entries(attrs)) el.setAttribute(k, v);
    svg.append(el);
  }
  return svg;
}

function statusField(key: string, extraClass = ''): { root: HTMLElement; value: HTMLElement } {
  const value = h('span', { class: 'v dim' }, '…');
  return { root: h('div', { class: `status-field ${extraClass}`.trim() }, h('span', { class: 'k' }, key), value), value };
}

const app = h('div', { class: 'app' });
const fCmb = statusField('CMB map');
const fBuild = statusField('Pipeline build', 'optional-narrow');
const fData = statusField('Datasets', 'optional');
const fBundles = statusField('Bundles', 'optional');
const fVersion = statusField('Version', 'optional-wide');
const warnValue = h('span', { class: 'v' }, '…');
const fWarn = h('a', { class: 'status-field', href: '#/data/warnings', title: 'Build warnings' }, h('span', { class: 'k' }, 'Warnings'), warnValue);
const clock = statusField('UTC');
const themeBtn = h('button', { type: 'button', class: 'icon-btn', 'aria-label': 'Toggle light or dark theme' });
const menuBtn = h('button', { type: 'button', class: 'icon-btn menu-btn', 'aria-label': 'Open navigation' }, icon('menu'));

const showLink = h('a', { class: 'show-link', href: '#discovery', title: 'Back to the show' }, '◂ Show');

const topbar = h(
  'header',
  { class: 'topbar' },
  menuBtn,
  showLink,
  h(
    'a',
    { class: 'brand', href: '#/overview', 'aria-label': 'Q-CHAS overview' },
    logo(),
    h('span', null, h('span', { class: 'brand-name' }, 'Q-CHAS'), h('span', { class: 'brand-sub' }, 'Mission console')),
  ),
  h(
    'div',
    { class: 'status-strip', role: 'status', 'aria-label': 'Pipeline status' },
    fCmb.root,
    fBuild.root,
    fData.root,
    fBundles.root,
    fVersion.root,
    fWarn,
    h('div', { class: 'status-right' }, clock.root, themeBtn),
  ),
);

const navLinks = new Map<RouteId, HTMLAnchorElement>();
const navDots = new Map<RouteId, HTMLElement>();
const sidebar = h('nav', { class: 'sidebar', 'aria-label': 'Tools' });
for (const g of NAV_GROUPS) {
  const group = h('div', { class: 'nav-group' });
  if (g.id !== 'overview')
    group.append(h('div', { class: 'nav-group-head' }, h('div', { class: 'nav-group-code' }, g.code), h('div', { class: 'nav-group-title' }, g.title)));
  for (const r of ROUTES.filter((x) => x.group === g.id)) {
    const dot = h('span', { class: 'nav-dot', 'aria-hidden': 'true' });
    const a = h(
      'a',
      { class: 'nav-link', href: `#/${r.id}` },
      h('span', { class: 'nav-code' }, r.code),
      h('span', { class: 'nav-label' }, r.label),
      dot,
    );
    navLinks.set(r.id, a);
    navDots.set(r.id, dot);
    group.append(a);
  }
  sidebar.append(group);
}
sidebar.append(
  h(
    'div',
    { class: 'sidebar-foot' },
    h('div', { class: 'legend-row' }, h('span', { class: 'nav-dot ready' }), 'inputs built'),
    h('div', { class: 'legend-row' }, h('span', { class: 'nav-dot partial' }), 'some inputs built'),
    h('div', { class: 'legend-row' }, h('span', { class: 'nav-dot' }), 'awaiting qchas build'),
    h('div', { style: { marginTop: '10px' } }, 'NASA Space Apps · Q-CHAS'),
  ),
);

const eyebrow = h('div', { class: 'page-eyebrow' });
const pageTitle = h('h1', { class: 'page-title', tabindex: '-1' });
const pageSub = h('p', { class: 'page-sub' });
const headExtras = h('div', { class: 'page-head-right' });
const pageBody = h('div', { class: 'page-body' });
const main = h(
  'main',
  { class: 'main', id: 'main' },
  h('div', { class: 'page' }, h('div', { class: 'page-head' }, h('div', null, eyebrow, pageTitle, pageSub), headExtras), pageBody),
);

app.append(topbar, sidebar, main);
app.hidden = true;
document.body.append(app);

// --- theme and clock -----------------------------------------------------------------------

function paintThemeButton(): void {
  const dark = currentTheme() === 'dark';
  themeBtn.replaceChildren(icon(dark ? 'sun' : 'moon'), dark ? 'Light' : 'Dark');
  themeBtn.setAttribute('aria-pressed', String(!dark));
}
themeBtn.addEventListener('click', () => {
  applyTheme(currentTheme() === 'dark' ? 'light' : 'dark');
  paintThemeButton();
});
paintThemeButton();

menuBtn.addEventListener('click', () => app.classList.toggle('nav-open'));

function tick(): void {
  clock.value.textContent = new Date().toISOString().slice(0, 19).replace('T', ' ');
  clock.value.classList.remove('dim');
}
tick();
setInterval(tick, 1000);

// --- status bar from the manifest -------------------------------------------------------------

function paintStatus(m: Manifest | null): void {
  for (const f of [fCmb, fBuild, fData, fBundles, fVersion]) f.value.classList.toggle('dim', !m);
  if (!m) {
    fCmb.value.textContent = 'not built';
    fBuild.value.textContent = 'run qchas build';
    fData.value.textContent = '—';
    fBundles.value.textContent = `0 / ${BUILT_BUNDLES.length}`;
    fVersion.value.textContent = '—';
    warnValue.replaceChildren(h('span', { class: 'indicator warning' }, 'no build'));
    for (const r of ROUTES) setDot(r, new Set());
    return;
  }
  fCmb.value.textContent = m.cmb_map_source ?? 'none';
  fBuild.value.textContent = fmtUtc(m.generated_utc);
  fData.value.textContent = String(Object.keys(m.datasets ?? {}).length);
  const built = new Set(m.bundles ?? []);
  fBundles.value.textContent = `${BUILT_BUNDLES.filter((b) => built.has(b)).length} / ${BUILT_BUNDLES.length}`;
  fVersion.value.textContent = m.qchas_version ?? '—';
  const n = m.warnings?.length ?? 0;
  const nErr = m.errors?.length ?? 0;
  const label = [nErr ? `${nErr} failed` : '', n ? `${n} active` : ''].filter(Boolean).join(' · ') || 'none';
  warnValue.replaceChildren(h('span', { class: `indicator ${nErr ? 'serious' : n ? 'warning' : 'good'}` }, label));
  fWarn.title = nErr ? `Build errors: ${m.errors!.map((e) => e.bundle).join(', ')}` : 'Build warnings';
  built.add('manifest.json');
  for (const r of ROUTES) setDot(r, built);
}

function setDot(r: RouteInfo, built: Set<string>): void {
  const dot = navDots.get(r.id);
  if (!dot) return;
  const have = r.bundles.filter((b) => built.has(b));
  const state = have.length === r.bundles.length ? 'ready' : have.length ? 'partial' : '';
  dot.className = `nav-dot ${state}`.trim();
  const missing = r.bundles.filter((b) => !built.has(b));
  navLinks.get(r.id)?.setAttribute('title', missing.length ? `Awaiting ${missing.join(', ')}` : 'Inputs built');
}

// --- router ------------------------------------------------------------------------------------

let manifest: Manifest | null = null;
let token = 0;
let cleanups: (() => void)[] = [];

function teardown(): void {
  for (const fn of cleanups.splice(0)) {
    try {
      fn();
    } catch (err) {
      console.error(err);
    }
  }
}

async function route(): Promise<void> {
  const { id, rest } = parseHash(location.hash);
  if (location.hash === '#/') history.replaceState(null, '', `#/${id}`);
  const info = ROUTES.find((r) => r.id === id)!;
  const my = ++token;
  teardown();
  app.classList.remove('nav-open');

  for (const [rid, a] of navLinks) {
    if (rid === id) a.setAttribute('aria-current', 'page');
    else a.removeAttribute('aria-current');
  }
  const group = NAV_GROUPS.find((g) => g.id === info.group)!;
  eyebrow.textContent = info.group === 'overview' ? 'Q-CHAS · Quantum Cosmology & Horizon Analytics Suite' : `${group.code} · ${group.title}`;
  pageTitle.textContent = info.label;
  pageSub.textContent = info.blurb;
  clear(headExtras);
  pageBody.replaceChildren(h('div', { class: 'page-loading' }, 'Loading tool…'));
  main.scrollTop = 0;

  let mod: Awaited<ReturnType<(typeof PANEL_LOADERS)[typeof id]>>;
  try {
    mod = await PANEL_LOADERS[id]();
  } catch (err) {
    // e.g. a stale chunk after a redeploy: say so instead of spinning forever
    console.error(err);
    if (my === token)
      pageBody.replaceChildren(
        h('div', { class: 'empty-state', role: 'alert' }, icon('warning'), h('div', null, h('div', { class: 'empty-title' }, 'This tool could not be loaded'), h('div', { class: 'empty-text' }, `${String(err)}. Reload the page to fetch the current build.`))),
      );
    return;
  }
  if (my !== token) return;
  pageTitle.textContent = mod.title;
  pageSub.textContent = mod.subtitle ?? info.blurb;
  document.title = `${mod.title} · Q-CHAS`;
  clear(pageBody);

  const ctx: PanelContext = {
    route: info,
    manifest,
    theme: currentTheme(),
    load: loadBundle,
    // after navigation away, release immediately instead of leaking into the next route
    onCleanup: (fn) => (my === token ? cleanups.push(fn) : fn()),
    onThemeChange: (cb) => {
      if (my === token) cleanups.push(onThemeChange(cb));
    },
    setSubtitle: (t) => (pageSub.textContent = t),
    setHeaderExtras: (nodes) => headExtras.replaceChildren(...nodes),
    navigate: (rid) => (location.hash = `#/${rid}`),
    get unmounted() {
      return my !== token;
    },
  };
  try {
    await mod.mount(pageBody, ctx);
  } catch (err) {
    console.error(err);
    if (my === token)
      pageBody.replaceChildren(
        h('div', { class: 'empty-state', role: 'alert' }, icon('warning'), h('div', null, h('div', { class: 'empty-title' }, 'This tool failed to load'), h('div', { class: 'empty-text' }, String(err)))),
      );
  }
  if (my === token && rest[0]) document.getElementById(rest[0])?.scrollIntoView({ block: 'start' });
}

// --- show (landing, plain anchors) vs console (`#/route`) ------------------------------------------

let show: ShowController | null = null;
let showLoading: Promise<ShowController> | null = null;
let manifestReady = false;

async function dispatch(): Promise<void> {
  const view = viewForHash(location.hash);
  if (view.mode === 'show') {
    if (!app.hidden) {
      token++;
      teardown();
      app.hidden = true;
    }
    document.documentElement.classList.add('show-mode');
    showLoading ??= import('./show/show').then((m) => m.mountShow(document.body));
    show = await showLoading;
    if (viewForHash(location.hash).mode === 'show') show.go(view.anchor, view.at);
    return;
  }
  show?.suspend();
  document.documentElement.classList.remove('show-mode');
  app.hidden = false;
  if (manifestReady) await route();
}

window.addEventListener('hashchange', () => void dispatch());
void dispatch();

void loadBundle('manifest.json').then((r) => {
  manifest = r.ok ? r.data : null;
  manifestReady = true;
  paintStatus(manifest);
  if (!app.hidden) void route();
});
