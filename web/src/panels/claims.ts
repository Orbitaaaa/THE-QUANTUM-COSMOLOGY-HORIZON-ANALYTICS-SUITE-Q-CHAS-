/**
 * Claims ledger: every statement of topology.json.claims and horizon.json.claims (plus
 * bounce.json.claims when present) with its status exactly as the papers' status tables give
 * it ([TOP] Sec. 13, [QG] Sec. 20, [BH] Sec. 16), filterable by provenance, paper and text.
 */

import './claims.css';
import type { Claim, LoadFailure, LoadResult, TopologyBundle, ValidationBundle, ValidationRow } from '../lib/data';
import type { PanelContext } from '../lib/panel';
import {
  type ProvenanceKind,
  append,
  badge,
  classifyStatus,
  fmtNum,
  grid,
  h,
  loadFailure,
  missingBundle,
  note,
  panel,
  provenanceLabel,
  statTile,
  statusPill,
  table,
  tiles,
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
  const [topo, hor, bounce, validation] = await Promise.all([
    ctx.load('topology.json'),
    ctx.load('horizon.json'),
    ctx.load('bounce.json'),
    ctx.load('validation.json'),
  ]);
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
    el.append(grid(ledger.root), grid(reproductionPanel(validation, topo.ok ? topo.data : null).root), grid(legendPanel().root));
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
  const jump = h('button', { type: 'button', class: 'btn' }, validation.ok ? `Reproduction ledger: ${validation.data.summary?.match ?? '—'} of ${validation.data.summary?.total ?? validation.data.rows.length} match` : 'Reproduction ledger');
  jump.addEventListener('click', () => document.getElementById('reproduction')?.scrollIntoView({ behavior: 'smooth', block: 'start' }));
  ctx.setHeaderExtras([badge('quoted', '[TOP] Sec. 13 · [QG] Sec. 20 · [BH] Sec. 16'), jump]);

  const parts: HTMLElement[] = [ledger.root];
  if (failures.length) {
    const miss = panel({ title: 'Not yet in the ledger' });
    miss.body.append(h('div', { class: 'stack' }, ...failures.map(loadFailure)));
    parts.push(miss.root);
  }
  el.append(grid(...parts), grid(reproductionPanel(validation, topo.ok ? topo.data : null).root), grid(legendPanel().root));
}

// --- reproduction ledger (validation.json) -----------------------------------------------------------

/** Ours and printed values: plain decimals in the normal range, 4 significant digits outside it. */
const fmtVal = (v: number): string => (Math.abs(v) >= 1e6 || (Math.abs(v) < 1e-3 && v !== 0) ? fmtNum(v, 4) : fmtNum(v, 6));

/**
 * Every number Q-CHAS recomputes from scratch against the value printed in the paper, with the
 * discrepancies we found stated next to both values.
 */
function reproductionPanel(res: LoadResult<ValidationBundle>, topo: TopologyBundle | null) {
  const p = panel({
    title: 'Reproduction ledger',
    subtitle: 'Every number Q-CHAS recomputes from scratch, next to the value printed in the paper',
    badges: [badge('computed', 'validation.json'), badge('quoted', 'printed values')],
    span: 12,
  });
  p.root.id = 'reproduction';
  if (!res.ok) {
    p.body.append(
      res.reason === 'missing'
        ? missingBundle('validation.json', 'qchas build writes validation.json, which compares each recomputed number with the printed one.')
        : loadFailure(res),
    );
    return p;
  }
  const V = res.data;
  const rows = V.rows;
  const total = V.summary?.total ?? rows.length;
  const nMatch = V.summary?.match ?? rows.filter((r) => r.verdict === 'match').length;
  const nDiff = V.summary?.differs ?? rows.filter((r) => r.verdict !== 'match').length;

  // groups in bundle order
  const groups: { name: string; n: number; match: number }[] = [];
  for (const r of rows) {
    let g = groups.find((x) => x.name === r.group);
    if (!g) groups.push((g = { name: r.group, n: 0, match: 0 }));
    g.n += 1;
    if (r.verdict === 'match') g.match += 1;
  }

  // the discrepancies, stated with both values
  const bcc = rows.find((r) => /Hessian eigenvalue at BCC/i.test(r.quantity) && r.note);
  const lambdaRows = rows.filter((r) => r.verdict !== 'match' && /Lambda|Λ/.test(r.group));
  const onlyOne = (topo?.manifolds ?? [])
    .filter((m) => Object.values(m.principles ?? {}).filter((v) => v !== 'allowed').length === 1)
    .map((m) => ({ label: m.label, fails: Object.entries(m.principles).find(([, v]) => v !== 'allowed')?.[0] ?? '' }));
  const diffs: HTMLElement[] = [];
  if (bcc)
    diffs.push(
      h(
        'li',
        null,
        h('strong', null, '[TOP] Prop. E.2, BCC Hessian. '),
        `The paper prints ${fmtVal(bcc.paper)} for the triple eigenvalue; we obtain ${fmtNum(bcc.ours, 6)} from both the analytic and a finite-difference Hessian. The printed value breaks the exact identity tr H = 2Z/3 by more than its rounding allows; ours satisfies it to 10⁻¹⁴. It is listed as a match (relative difference ${fmtNum(bcc.rel_diff, 2)}), and no conclusion of Prop. E.2 changes.`,
      ),
    );
  if (lambdaRows.length)
    diffs.push(
      h(
        'li',
        null,
        h('strong', null, '[QG] Theorem Λ table. '),
        `Our rerun on the paper's own finite-difference grid reproduces every printed number. The exact solution gives larger peaks of ⟨a⁻²⟩ (${lambdaRows.map((r) => `${fmtNum(r.ours, 4)} vs ${fmtNum(r.paper, 4)}`).join('; ')}): the grid sum misses a boundary term at a = 0. Exactly, ⟨a⟩_min tends to 0.1629 ℓ = 0.260 l_P. These are the ${lambdaRows.length} “differs” rows; the paper's qualitative claims hold with either set.`,
      ),
    );
  if (onlyOne.length)
    diffs.push(
      h(
        'li',
        null,
        h('strong', null, '[TOP] Secs. 4 and 5.6. '),
        `The sentence “every manifold except E1 fails at least two principles” does not hold for ${onlyOne.map((m) => m.label).join(', ')}, which by the paper's own Sec. 4 table fail only ${[...new Set(onlyOne.map((m) => m.fails))].join(', ')}. The Main Theorem is unaffected: both selecting pairs contain P2.`,
      ),
    );

  const state = { differsOnly: false, notedOnly: false, q: '' };
  const chip = (label: string, key: 'differsOnly' | 'notedOnly'): HTMLButtonElement => {
    const c = h('button', { type: 'button', class: 'chip', 'aria-pressed': 'false' }, label);
    c.addEventListener('click', () => {
      state[key] = !state[key];
      c.setAttribute('aria-pressed', String(state[key]));
      render();
    });
    return c;
  };
  const search = h('input', { type: 'search', placeholder: 'Search quantities', 'aria-label': 'Search reproduced quantities' });
  search.addEventListener('input', () => {
    state.q = search.value.trim().toLowerCase();
    render();
  });
  const shown = h('span', { class: 'muted small mono' });
  const host = h('div', { class: 'claims-repro-scroll' });

  function render(): void {
    const q = state.q;
    const visible = rows.filter(
      (r) =>
        (!state.differsOnly || r.verdict !== 'match') &&
        (!state.notedOnly || !!r.note) &&
        (!q || `${r.group} ${r.quantity} ${r.source} ${r.note}`.toLowerCase().includes(q)),
    );
    shown.textContent = `${visible.length} of ${rows.length} rows shown`;
    host.replaceChildren(
      table<ValidationRow>(
        [
          { key: 'group', label: 'Result', className: 'wrap claims-repro-group' },
          { key: 'quantity', label: 'Quantity', className: 'wrap' },
          { key: 'ours', label: 'Q-CHAS', align: 'num', format: (r) => fmtVal(r.ours) },
          { key: 'paper', label: 'Paper', align: 'num', format: (r) => fmtVal(r.paper) },
          { key: 'rel_diff', label: 'Rel. diff.', align: 'num', format: (r) => fmtNum(r.rel_diff, 2) },
          {
            key: 'verdict',
            label: 'Verdict',
            className: 'wrap',
            format: (r) => [statusPill(r.verdict === 'match' ? 'good' : 'warning', r.verdict), r.note ? h('div', { class: 'small text-2 claims-repro-note' }, r.note) : null],
          },
          { key: 'source', label: 'Source', className: 'nowrap', format: (r) => h('code', null, r.source) },
        ],
        visible,
        'No row matches these filters',
      ),
    );
  }
  render();

  append(
    p.body,
    tiles(
      statTile({ label: 'Numbers reproduced', value: String(total), sub: [badge('computed', 'validation.json'), 'recomputed from scratch'] }),
      statTile({ label: 'Match the printed value', value: String(nMatch), sub: [statusPill('good', 'match'), 'within tolerance or to the printed digits'] }),
      statTile({ label: 'Differ from the printed value', value: String(nDiff), sub: [statusPill(nDiff ? 'warning' : 'good', nDiff ? 'differs' : 'none'), 'each row says why'] }),
    ),
    h('p', { class: 'small text-2 claims-repro-lede' }, V.note),
    table<(typeof groups)[number]>(
      [
        { key: 'name', label: 'Result reproduced', className: 'wrap' },
        { key: 'n', label: 'Rows', align: 'num', format: (g) => String(g.n) },
        { key: 'match', label: 'Match', align: 'num', format: (g) => String(g.match) },
        { key: 'differ', label: 'Differ', align: 'num', format: (g) => (g.n - g.match ? h('strong', null, String(g.n - g.match)) : '0') },
      ],
      groups,
    ),
    diffs.length
      ? h(
          'div',
          { class: 'callout claims-repro-callout' },
          h('div', { class: 'claims-repro-callout-title' }, 'Where we differ from the printed papers'),
          h('ul', null, ...diffs),
          h('div', { class: 'row small' }, badge('computed', 'validation.json, topology.json'), badge('quoted', '[TOP] Prop. E.2, Secs. 4–5 · [QG] Theorem Λ')),
        )
      : null,
    h('div', { class: 'filter-bar claims-repro-filters' }, chip('Differs only', 'differsOnly'), chip('With a note', 'notedOnly'), h('span', { class: 'spacer' }), shown, search),
    host,
    V.skipped?.length ? h('div', { class: 'chart-caption' }, `Not compared in this build: ${V.skipped.join('; ')}.`) : null,
  );
  return p;
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
