/**
 * Placeholder for a tool that has not been implemented yet: names the tool, lists what it will
 * show, the bundles it reads (with their build state) and the live physics modules available
 * to it. Replace the stub's `mount` with the real tool (see web/README.md).
 */

import { BUILD_COMMAND } from '../lib/data';
import type { PanelContext } from '../lib/panel';
import { badge, grid, h, icon, note, panel, statusPill, table } from '../lib/ui';

export interface StubSpec {
  /** what the finished tool will show, one item per line */
  plan: string[];
  /** live TypeScript physics modules the tool can use */
  physics?: { module: string; provides: string }[];
  /** the paper sections the tool implements */
  sources: string[];
}

export function mountStub(spec: StubSpec) {
  return async (el: HTMLElement, ctx: PanelContext): Promise<void> => {
    const status = await Promise.all(ctx.route.bundles.map(async (b) => ({ b, r: await ctx.load(b) })));
    if (ctx.unmounted) return;

    const intro = panel({
      title: `${ctx.route.code} · ${ctx.route.label}`,
      subtitle: 'Tool under construction: this console slot is reserved and wired to its data bundles.',
      badges: spec.sources.map((s) => badge('quoted', s)),
      span: 7,
      stretch: true,
    });
    intro.body.append(
      h(
        'div',
        { class: 'empty-state', role: 'status' },
        icon('wrench'),
        h(
          'div',
          null,
          h('div', { class: 'empty-title' }, `${ctx.route.label} is being implemented`),
          h('div', { class: 'empty-text' }, ctx.route.blurb),
        ),
      ),
      h('div', { style: { marginTop: '14px' } }, note(spec.plan, 'What this tool will show', true)),
    );

    const inputs = panel({
      title: 'Inputs',
      subtitle: `Bundles read from ./data/ (produced by ${BUILD_COMMAND})`,
      span: 5,
      stretch: true,
      flush: true,
    });
    inputs.body.append(
      table(
        [
          { key: 'b', label: 'Bundle', format: (r) => h('code', null, r.b) },
          {
            key: 'r',
            label: 'State',
            format: (r) =>
              r.r.ok
                ? statusPill('good', 'built')
                : r.r.reason === 'missing'
                  ? statusPill('warning', 'not built yet')
                  : statusPill('critical', r.r.reason),
          },
        ],
        status,
      ),
    );
    if (spec.physics?.length) {
      inputs.body.append(
        h('div', { style: { padding: '12px 14px 4px' }, class: 'muted small' }, 'Live in-browser physics available to this tool'),
        table(
          [
            { key: 'module', label: 'Module', format: (r) => h('code', null, r.module) },
            { key: 'provides', label: 'Provides', className: 'wrap' },
          ],
          spec.physics,
        ),
      );
    }
    el.append(grid(intro.root, inputs.root));
  };
}
