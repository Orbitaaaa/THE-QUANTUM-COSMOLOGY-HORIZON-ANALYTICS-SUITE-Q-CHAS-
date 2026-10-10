/**
 * Data & provenance: manifest.json rendered in full, every dataset with its archive URL, size,
 * SHA-256 and retrieval time, the bundles present, and the build warnings.
 */

import { BUILD_COMMAND, BUNDLE_NAMES, type DatasetRecord } from '../lib/data';
import type { PanelContext } from '../lib/panel';
import {
  copyButton,
  fmtBytes,
  fmtUtc,
  grid,
  h,
  icon,
  missingBundle,
  note,
  panel,
  statTile,
  statusPill,
  table,
  tiles,
} from '../lib/ui';

export const title = 'Data & Provenance';
export const subtitle =
  'The real datasets behind every number: where each file came from, its size and SHA-256, when it was retrieved, and what the build could not reach.';

const isUrl = (s: string): boolean => /^https?:\/\//i.test(s);

export async function mount(el: HTMLElement, ctx: PanelContext): Promise<void> {
  const res = await ctx.load('manifest.json');
  if (ctx.unmounted) return;

  if (!res.ok) {
    const p = panel({ title: 'Provenance record', subtitle: 'manifest.json' });
    p.body.append(
      res.reason === 'missing'
        ? missingBundle('manifest.json', `The provenance record is written by qchas fetch (URL, size and SHA-256 of every file in data/raw/provenance.json) and exported by ${BUILD_COMMAND}.`)
        : h('div', { class: 'empty-state' }, icon('warning'), res.message),
    );
    el.append(grid(p.root), grid(accessNote().root));
    return;
  }

  const m = res.data;
  const datasets = Object.entries(m.datasets ?? {}).map(([key, d]) => ({ key, ...d }));
  const totalBytes = datasets.reduce((s, d) => s + (d.bytes || 0), 0);
  const built = new Set(m.bundles ?? []);
  const warningsList = m.warnings ?? [];
  const errors = m.errors ?? [];
  // builder names in manifest.errors are bundle stems ("cmb"), except the sky map builder
  const failed = new Set(errors.map((e) => (e.bundle === 'sky' ? 'sky_map.json' : `${e.bundle}.json`)));

  const summary = panel({ title: 'Build summary' });
  summary.body.append(
    tiles(
      statTile({ label: 'Generated', value: fmtUtc(m.generated_utc, true), mono: true }),
      statTile({ label: 'Pipeline version', value: m.qchas_version, mono: true }),
      statTile({ label: 'CMB map source', value: m.cmb_map_source ?? 'none', mono: true }),
      statTile({ label: 'Datasets', value: datasets.length, sub: fmtBytes(totalBytes) }),
      statTile({ label: 'Bundles', value: `${BUNDLE_NAMES.filter((b) => b !== 'manifest.json' && built.has(b)).length} / ${BUNDLE_NAMES.length - 1}` }),
    ),
  );

  const warnings = panel({
    title: 'Build warnings',
    subtitle: 'Inputs the build host could not reach, what was used instead, and builders that failed',
    id: 'warnings',
    span: 6,
    stretch: true,
  });
  if (errors.length)
    warnings.body.append(
      h(
        'ul',
        { class: 'warning-list' },
        ...errors.map((e) =>
          h('li', null, icon('critical', 'build error'), h('span', null, h('strong', null, `${e.bundle}: `), e.error)),
        ),
      ),
    );
  if (warningsList.length)
    warnings.body.append(
      h('ul', { class: 'warning-list' }, ...warningsList.map((w) => h('li', null, icon('warning', 'warning'), h('span', null, w)))),
    );
  else if (!errors.length) warnings.body.append(statusPill('good', 'No warnings: every input was reachable'));

  const bundles = panel({ title: 'Bundles', subtitle: 'web/public/data/', span: 6, stretch: true, flush: true });
  bundles.body.append(
    table(
      [
        { key: 'name', label: 'Bundle', format: (r) => h('code', null, r.name) },
        {
          key: 'state',
          label: 'State',
          format: (r) =>
            r.ok ? statusPill('good', 'built') : failed.has(r.name) ? statusPill('serious', 'build failed') : statusPill('warning', 'not built'),
        },
      ],
      // the contract's bundles, then the extras the manifest lists (validation.json, binary maps)
      [
        ...BUNDLE_NAMES.filter((n) => n !== 'manifest.json'),
        ...[...built].filter((n) => n !== 'manifest.json' && !(BUNDLE_NAMES as string[]).includes(n)).sort(),
      ].map((name) => ({ name, ok: built.has(name) })),
    ),
  );

  const table1 = panel({ title: 'Datasets', subtitle: `${datasets.length} files with provenance`, flush: true });
  table1.body.append(
    table<DatasetRecord & { key: string }>(
      [
        {
          key: 'key',
          label: 'Dataset',
          className: 'wrap',
          format: (d) => h('div', null, h('div', { style: { fontWeight: '600' } }, d.mission || '—'), h('code', null, d.key)),
        },
        {
          key: 'description',
          label: 'Description',
          className: 'wrap wide',
          format: (d) => h('div', null, d.description, h('div', { class: 'small muted', style: { marginTop: '3px' } }, d.citation)),
        },
        { key: 'filename', label: 'File', className: 'break mono small', format: (d) => d.filename },
        { key: 'bytes', label: 'Size', align: 'num', format: (d) => fmtBytes(d.bytes) },
        {
          key: 'sha256',
          label: 'SHA-256',
          format: (d) =>
            d.sha256
              ? h('span', { class: 'hash', title: d.sha256 }, `${d.sha256.slice(0, 12)}…`, copyButton(d.sha256))
              : h('span', { class: 'muted' }, '—'),
        },
        { key: 'retrieved_utc', label: 'Retrieved', className: 'mono small nowrap', format: (d) => fmtUtc(d.retrieved_utc) },
        {
          key: 'url',
          label: 'Source',
          className: 'break small',
          format: (d) =>
            isUrl(d.url)
              ? h('a', { href: d.url, target: '_blank', rel: 'noopener noreferrer', title: d.url }, hostOf(d.url), ' ', icon('external'))
              : h('span', { class: 'muted' }, d.url || '—'),
        },
      ],
      datasets,
      'The manifest lists no datasets',
    ),
  );

  el.append(grid(summary.root), grid(warnings.root, bundles.root), grid(table1.root), grid(accessNote().root));
  const n = warningsList.length;
  ctx.setHeaderExtras([
    ...(errors.length ? [statusPill('serious', `${errors.length} builder${errors.length === 1 ? '' : 's'} failed`)] : []),
    statusPill(n ? 'warning' : 'good', n ? `${n} warning${n === 1 ? '' : 's'}` : 'no warnings'),
  ]);
}

function hostOf(url: string): string {
  try {
    const u = new URL(url);
    const path = u.pathname.split('/').filter(Boolean);
    return `${u.host}/…/${path[path.length - 1] ?? ''}`;
  } catch {
    return url;
  }
}

function accessNote() {
  const p = panel({ title: 'Archive access', subtitle: 'How the build obtains real data' });
  p.body.append(
    note(
      [
        'qchas fetch tries every mirror of every dataset and records URL, size and SHA-256 in data/raw/provenance.json; this page shows that record as exported to manifest.json.',
        'The Planck PR3 archive, the GWOSC API and IRSA are fetched on GitHub Actions runners. When one is unreachable the build falls back to data that is present (for example WMAP instead of Planck) and lists a warning here.',
        'No simulated or placeholder values ever stand in for an observation: a missing input produces a null field and a warning.',
      ],
      'Provenance policy',
      true,
    ),
  );
  return p;
}
