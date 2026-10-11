/** Small DOM pieces of the show: honesty chips, buttons, meters. */

import { h } from '../lib/ui';
import { type Chip, TAG_INFO } from './script';

export function chip(c: Chip, extraClass = ''): HTMLElement {
  const el = h(
    'div',
    { class: `chip tag-${c.tag.toLowerCase()} ${extraClass}`.trim(), tabindex: c.more ? '0' : null, title: c.more ?? TAG_INFO[c.tag] },
    h('span', { class: 'chip-tag' }, c.tag === 'DATA' ? 'DATA' : c.tag),
    h('span', { class: 'chip-text' }, c.text),
  );
  if (c.more) el.append(h('span', { class: 'chip-more', role: 'note' }, c.more));
  return el;
}

export function button(label: string, onClick: () => void, cls = 'show-btn', aria?: string): HTMLButtonElement {
  const b = h('button', { type: 'button', class: cls, 'aria-label': aria ?? null }, label);
  b.addEventListener('click', onClick);
  return b;
}

export const fmt = (x: number, digits = 2): string =>
  !Number.isFinite(x) ? '—' : Math.abs(x) >= 100 ? x.toFixed(0) : Math.abs(x) >= 10 ? x.toFixed(1) : x.toFixed(digits);

/** The tag legend shown in the intro overlay. */
export function tagLegend(): HTMLElement {
  return h(
    'div',
    { class: 'tag-legend' },
    ...(['ESTABLISHED', 'PROVED', 'COMPUTED', 'HYPOTHESIS', 'ILLUSTRATIVE', 'DATA'] as const).map((t) =>
      h('span', { class: `chip mini tag-${t.toLowerCase()}`, title: TAG_INFO[t] }, h('span', { class: 'chip-tag' }, t)),
    ),
  );
}

/** Keep the filled part of a range input's track in sync with its value. */
export function syncFill(input: HTMLInputElement): void {
  const min = Number(input.min || 0);
  const max = Number(input.max || 100);
  const f = max > min ? ((Number(input.value) - min) / (max - min)) * 100 : 0;
  input.style.setProperty('--fill', `${f.toFixed(2)}%`);
}
