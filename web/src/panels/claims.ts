/**
 * Claims ledger: every statement of topology.json.claims and horizon.json.claims (plus
 * bounce.json.claims when present) with its status exactly as the papers' status tables give
 * it ([TOP] Sec. 13, [QG] Sec. 20, [BH] Sec. 16), filterable by provenance, paper and text.
 */

import type { Claim, LoadFailure } from '../lib/data';
import type { PanelContext } from '../lib/panel';
import {
  type ProvenanceKind,
  badge,
  classifyStatus,
  grid,
  h,
  loadFailure,
  missingBundle,
  note,
  panel,
  provenanceLabel,
  table,
} from '../lib/ui';

export const title = 'Claims Ledger';
export const subtitle =
  "Every statement Q-CHAS shows, with its status as the source paper's own status table gives it: theorem, computation, data analysis or hypothesis.";

interface Row extends Claim {
  kind: ProvenanceKind;
  paper: string;
  origin: string;
}

const KINDS: ProvenanceKind[] = ['proved', 'computed', 'measured', 'assumed', 'open', 'quoted'];

const paperKey = (source: string): string => /\[(TOP|QG|BH)\]/.exec(source)?.[1] ?? 'other';

export async function mount(el: HTMLElement, ctx: PanelContext): Promise<void> {
  const [topo, hor, bounce] = await Promise.all([ctx.load('topology.json'), ctx.load('horizon.json'), ctx.load('bounce.json')]);
  if (ctx.unmounted) return;

  const rows: Row[] = [];
  const add = (claims: Claim[] | undefined, origin: string): void => {
    for (const c of claims ?? [])
      rows.push({ ...c, kind: classifyStatus(c.status), paper: paperKey(c.source), origin });
  };
  if (topo.ok) add(topo.data.claims, 'topology.json');
  if (hor.ok) add(hor.data.claims, 'horizon.json');
  if (bounce.ok) add(bounce.data.claims, 'bounce.json');

  const failures = [topo, hor].filter((r): r is LoadFailure => !r.ok);
  const ledger = panel({
    title: 'Ledger',
    subtitle: rows.length
      ? `${rows.length} statements from ${[topo, hor, bounce].filter((r) => r.ok).map((r) => r.name).join(', ')}`
      : 'No claims loaded',
    flush: true,
  });

  if (!rows.length) {
    ledger.body.style.padding = '14px';
    ledger.body.append(
      failures.every((f) => f.reason === 'missing')
        ? missingBundle(['topology.json', 'horizon.json'], 'The ledger lists the status rows the Python pipeline exports from each paper; it is empty until they are built.')
        : h('div', { class: 'stack' }, ...failures.map(loadFailure)),
    );
    el.append(grid(ledger.root), grid(legendPanel().root));
    return;
  }

  // --- filters ---------------------------------------------------------------------------------
  const state = { kinds: new Set<ProvenanceKind>(), papers: new Set<string>(), q: '' };
  const counts = (k: ProvenanceKind): number => rows.filter((r) => r.kind === k).length;
  const kindChips = KINDS.filter((k) => counts(k) > 0).map((k) => {
    const chip = h('button', { type: 'button', class: 'chip', 'aria-pressed': 'false' }, provenanceLabel(k), h('span', { class: 'count' }, String(counts(k))));
    chip.addEventListener('click', () => {
      if (state.kinds.has(k)) state.kinds.delete(k);
      else state.kinds.add(k);
      chip.setAttribute('aria-pressed', String(state.kinds.has(k)));
      render();
    });
    return chip;
  });
  const papers = [...new Set(rows.map((r) => r.paper))].sort();
  const paperChips = papers.map((p) => {
    const chip = h('button', { type: 'button', class: 'chip', 'aria-pressed': 'false' }, p === 'other' ? 'Other' : `[${p}]`);
    chip.addEventListener('click', () => {
      if (state.papers.has(p)) state.papers.delete(p);
      else state.papers.add(p);
      chip.setAttribute('aria-pressed', String(state.papers.has(p)));
      render();
    });
    return chip;
  });
  const search = h('input', { type: 'search', placeholder: 'Search statements', 'aria-label': 'Search statements' });
  search.addEventListener('input', () => {
    state.q = search.value.trim().toLowerCase();
    render();
  });
  const shown = h('span', { class: 'muted small mono' });
  const filterBar = h('div', { class: 'filter-bar', style: { padding: '12px 14px 0' } }, ...kindChips, h('span', { class: 'spacer' }), ...paperChips, search);
  const tableHost = h('div');
  ledger.body.append(filterBar, h('div', { style: { padding: '0 14px 10px' } }, shown), tableHost);

  function render(): void {
    const q = state.q;
    const visible = rows.filter(
      (r) =>
        (!state.kinds.size || state.kinds.has(r.kind)) &&
        (!state.papers.size || state.papers.has(r.paper)) &&
        (!q || `${r.statement} ${r.status} ${r.source}`.toLowerCase().includes(q)),
    );
    shown.textContent = `${visible.length} of ${rows.length} shown`;
    tableHost.replaceChildren(
      table(
        [
          { key: 'n', label: '#', align: 'num', format: (r) => String(rows.indexOf(r) + 1) },
          { key: 'statement', label: 'Statement', className: 'wrap' },
          {
            key: 'status',
            label: 'Status (as printed)',
            className: 'wrap',
            format: (r) => h('div', null, badge(r.kind), h('div', { class: 'small text-2', style: { marginTop: '4px' } }, r.status)),
          },
          { key: 'source', label: 'Source', className: 'nowrap', format: (r) => h('code', null, r.source) },
        ],
        visible,
        'No statement matches these filters',
      ),
    );
  }
  render();
  ctx.setHeaderExtras([badge('quoted', '[TOP] Sec. 13 · [QG] Sec. 20 · [BH] Sec. 16')]);

  const parts: HTMLElement[] = [ledger.root];
  if (failures.length) {
    const miss = panel({ title: 'Not yet in the ledger' });
    miss.body.append(h('div', { class: 'stack' }, ...failures.map(loadFailure)));
    parts.push(miss.root);
  }
  el.append(grid(...parts), grid(legendPanel().root));
}

function legendPanel() {
  const p = panel({ title: 'How statuses are classified', subtitle: 'The badge summarizes the leading clause of the printed status; the full text is always shown.' });
  p.body.append(
    note(
      [
        'Proved: a theorem, proposition or lemma of the paper (possibly conditional on stated hypotheses, which the printed status names).',
        'Computed: a numerical result (linear theory, simulation, lattice computation) with the precision the paper states.',
        'Measured: a data analysis of public observations (Planck, WMAP, LIGO), calibrated as the paper describes.',
        'Assumed: a hypothesis, conditional statement, estimate or path-integral derivation that the papers do not prove as a theorem.',
        'Open: a quantity the papers say must still be measured.',
      ],
      'Badge legend',
      true,
    ),
  );
  return p;
}
