/**
 * Theme state: dark by default, persisted in localStorage, broadcast as a 'qchas:theme' event
 * so charts, sky views and colormaps can re-render with the new tokens.
 */

import type { Theme } from './colormap';

export type { Theme };

const STORAGE_KEY = 'qchas.theme';
const EVENT = 'qchas:theme';

export function storedTheme(): Theme {
  try {
    const v = localStorage.getItem(STORAGE_KEY);
    if (v === 'light' || v === 'dark') return v;
  } catch {
    /* storage unavailable (private mode, sandboxed preview) */
  }
  return 'dark';
}

export function currentTheme(): Theme {
  return document.documentElement.dataset.theme === 'light' ? 'light' : 'dark';
}

export function applyTheme(theme: Theme, persist = true): void {
  document.documentElement.dataset.theme = theme;
  if (persist) {
    try {
      localStorage.setItem(STORAGE_KEY, theme);
    } catch {
      /* ignore */
    }
  }
  document.dispatchEvent(new CustomEvent<Theme>(EVENT, { detail: theme }));
}

/** Subscribe to theme changes; returns the unsubscribe function. */
export function onThemeChange(cb: (theme: Theme) => void): () => void {
  const fn = (e: Event): void => cb((e as CustomEvent<Theme>).detail);
  document.addEventListener(EVENT, fn);
  return () => document.removeEventListener(EVENT, fn);
}

/** Resolved value of a CSS custom property on <html>, e.g. cssVar('--series-1'). */
export function cssVar(name: string): string {
  return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
}

/** The eight categorical series colours of the current theme, in fixed slot order. */
export function seriesColors(): string[] {
  return Array.from({ length: 8 }, (_, i) => cssVar(`--series-${i + 1}`));
}
