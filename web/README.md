# Q-CHAS dashboard (`web/`)

Vite + TypeScript (strict), no UI framework. Three.js for the globe, Plotly (cartesian bundle)
for charts, Vitest for tests. The site is static and loads the Python pipeline's bundles at
runtime from `./data/` (`web/public/data/`, written by `qchas build`; see the data contract in
`docs/ARCHITECTURE.md`).

```bash
npm install
npm run dev        # http://localhost:5173/#/overview
npm test           # vitest: HEALPix, physics ports, loaders, routes
npm run build      # tsc --noEmit && vite build -> dist/ (base './', works from a Pages subpath)
npm run preview    # serve dist/
```

Without bundles every view renders an explicit **"Data bundle not built yet: run qchas build"**
state. Never put placeholder or simulated numbers in a panel to fill that gap.

## Layout

```
src/main.ts            shell: status bar, sidebar, hash router, theme toggle, manifest preload
src/routes.ts          route registry (id, nav code, group, bundles read, lazy panel loader)
src/styles.css         design tokens (dark default + light), layout, panels, controls, badges
src/lib/panel.ts       the panel contract (PanelModule, PanelContext)
src/lib/data.ts        typed interfaces for every bundle + cached loaders
src/lib/ui.ts          h(), panel(), statTile(), slider(), select(), toggle(), segmented(), table(), badge(), note(), missingBundle()
src/lib/charts.ts      lineChart, scatterChart, barChart, heatmapChart (Plotly, lazy-loaded)
src/lib/colormap.ts    OKLab diverging (CMB) and sequential maps per theme
src/lib/healpix.ts     RING ang2pix/pix2ang (any nside), Galactic helpers, Mollweide
src/lib/sky.ts         sky views (globe or Mollweide) with one overlay API; sky-globe.ts is the Three.js part
src/lib/physics/       live closed-form physics, ports of the validated Python
src/panels/<route>.ts  one module per tool
tests/                 vitest specs
```

## Writing a panel

Each file in `src/panels/` exports `title` and `mount(el, ctx)`. Stubs use `mountStub(...)`;
replace that with your own `mount`. The route, its navigation entry and readiness dot already
exist in `src/routes.ts` (update `bundles` there if your tool reads different ones).

```ts
import { lineChart } from '../lib/charts';
import type { PanelContext } from '../lib/panel';
import { badge, grid, loadFailure, panel, slider } from '../lib/ui';

export const title = 'Page Curve';
export const subtitle = "Hawking's curve against the Page transfer law"; // optional

export async function mount(el: HTMLElement, ctx: PanelContext): Promise<void> {
  const res = await ctx.load('horizon.json');   // cached, never throws
  if (ctx.unmounted) return;                     // user navigated away during the await
  const p = panel({ title: 'Page law', subtitle: 'units of S_BH(M0)', badges: [badge('computed', '[BH] Fig. 1')], span: 8 });
  el.append(grid(p.root));
  if (!res.ok) return void p.body.append(loadFailure(res));   // "not built yet" or error state
  const { t, page, hawking } = res.data.page;
  const chart = lineChart(p.body, [
    { name: 'Page curve', x: t, y: page },        // slot 1
    { name: "Hawking's curve", x: t, y: hawking },// slot 2
  ], { xTitle: 'time / lifetime', yTitle: 'entropy / S_BH(M0)' });
  ctx.onCleanup(() => chart.destroy());
}
```

`PanelContext`:

| member | use |
|---|---|
| `route` | the `RouteInfo` (code, label, blurb, bundles) |
| `manifest` | `manifest.json`, preloaded (null before `qchas build`) |
| `load(name)` | typed, cached bundle loader returning `{ ok: true, data } \| { ok: false, reason, message }` |
| `onCleanup(fn)` | runs on navigation away: dispose charts, sky views, timers |
| `onThemeChange(cb)` | theme subscription, released automatically |
| `setSubtitle(text)`, `setHeaderExtras(nodes)` | page header |
| `navigate(id)` | go to another route |
| `unmounted` | check after every `await` before touching the DOM |

Binary maps: `loadSkyMap(name?)` in `lib/data.ts` returns `{ meta, entry, map: Float32Array, mask }`
(little-endian Float32, RING, length checked against `12 nside^2`).

## UI conventions

* **Provenance on every number.** Use `badge(kind, source)` with kind `proved`, `computed`,
  `measured`, `assumed`, `open`, `quoted` or `live`, e.g. `badge('quoted', '[TOP] Sec. 9')`.
  When a bundle carries `paper` values, show them next to ours. `classifyStatus(status)` maps a
  paper's printed status string to a kind (the Claims Ledger uses it).
* **Panels** sit in a 12-column `grid()`; `span` 3 to 12, collapsing to full width below 1240 px.
  `note(items, 'What the paper proves')` is the collapsible theorem/assumption box.
* **Strings are text.** `h()` and `table()` insert strings with `textContent`; do not use
  `innerHTML` with bundle content.
* **Status colours** (`statusPill('good' | 'warning' | 'serious' | 'critical', label)`) are
  reserved for state and always carry an icon and a label.

## Charts (`lib/charts.ts`)

`lineChart`, `scatterChart` (error bars, symmetric or asymmetric), `barChart`, `heatmapChart`
(`scale: 'sequential' | 'diverging'`, symmetric about `zmid`). They follow the dashboard's
data-viz rules, so panels need not restyle them:

* series colours come from the validated categorical palette in **fixed slot order**
  (`color: 1..8`; default = series index + 1; `'context'` for a de-emphasized gray). Colour
  follows the entity: keep a series on the same slot across charts and filters. Past eight
  series, fold into "other" or use small multiples. Never generate new hues.
* one y-axis per chart (no dual axes); two measures of different scale -> two charts;
* a legend appears automatically for two or more series; `refLines` draw labelled reference
  lines (Page time, certified values);
* every chart has a **Table** toggle with the same numbers, re-renders on theme change, resizes
  with its container and disposes itself when removed (`handle.destroy()` also works).

## Sky views (`lib/sky.ts`)

```ts
const { createSkyView } = await import('../lib/sky');           // three.js loads only here
const sky = await createSkyView(host, { mode: 'globe', view: 'outer', height: 420 });
sky.setMap({ map, mask, limit: 300, unit: 'μK' });            // diverging, symmetric limits
sky.addSmallCircle([l, b], radiusDeg, '--series-1', { label: 'C1' });
sky.addPoint([l, b], '--series-2', 'twin');
sky.addChordThroughOrigin(lbA, lbB, '--series-3');              // point and twin through the observer
sky.addPolygon(cornersGal, '--series-4', { fill: true });       // survey footprint
const off = sky.onPick(({ l, b, value, masked }) => { ... });   // click picking
(sky as GlobeView).setView('inner');                            // observer at the centre
ctx.onCleanup(() => sky.dispose());
```

Galactic `(l, b)` in degrees everywhere. Colours may be palette tokens (`'--series-1'`), which
follow the theme. Without WebGL `createSkyView` returns the Mollweide view (same API) and says
so on screen. In the outer view a chord makes the sphere semi-transparent so the line through
the observer is visible; `setOpacity` overrides it.

## Live physics (`lib/physics/`)

Pure functions, unit-tested against the Python and the papers:

| module | mirrors | main exports |
|---|---|---|
| `page.ts` | `qchas/horizon/page.py` ([BH] Theorems C-F) | `hawkingEntropy`, `pageCurve`, `pageCurveLower`, `islandMinEntropy`, `interiorVolume`, `pageTimeFraction`, `curves` |
| `qes.ts` | `qchas/horizon/qes.py` ([BH] Sec. 9) | `solveQes(M, c, rb_over_rs)` (Levenberg-Marquardt on the analytic log-space gradient), `PAPER_TABLE`, `depthFit`, `lagFit` |
| `lattice.ts` | `qchas/topology/lattice.py` ([TOP] Prop. E.2) | `epsteinZ2(B)` (Ewald), `rhombohedralBasis`, `cubicBasis`, `lllReduce`, `latticeVectors`, `shortestVector` |
| `circles.ts` | `qchas/topology/circles.py` ([TOP] Corollary F.1) | `circlePairs(basis, d/D, R)`, `matchedPoint`, `circlePoints`, `ghostImages` |
| `capacity.ts` | `qchas/horizon/capacity.py` ([BH] Theorem G, [QG] C, U) | `capacityWindow`, `capacityCurve`, `deSitterEntropy`, `minimalArea`, `ALGEBRA_LADDER` |
| `constants.ts` | `qchas/constants.py` | quoted numbers as `{ value, source }` |

Heavy or data-driven results (likelihoods, flows, transfer functions, strain fits) come from the
bundles, not from the browser.
