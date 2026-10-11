/**
 * The show controller (docs/SHOW.md): full-bleed stage, act navigation, captions and honesty
 * chips, the act transport (play/pause/scrub), the 60-second guided tour, the intro overlay,
 * keyboard control and the QA hook `window.__qchasShow`.
 */

import './show.css';
import { loadBundle, loadSkyMap } from '../lib/data';
import { h } from '../lib/ui';
import type { ShowAnchor } from '../routes';
import type { ActContext, ActHandle, ActModule, ShowData } from './act';
import { ACTS, type ActCopy, type ActId, BRAND, type Beat, actCopy } from './script';
import { Stage, webglSupported } from './stage';
import { button, chip, syncFill, tagLegend } from './ui';

const LOADERS: Record<ActId, () => Promise<ActModule>> = {
  discovery: () => import('./acts/discovery').then((m) => m.discovery),
  experiment: () => import('./acts/experiment').then((m) => m.experiment),
  breakthrough: () => import('./acts/breakthrough-blackhole').then((m) => m.blackhole),
  bounce: () => import('./acts/breakthrough-bounce').then((m) => m.bounce),
  evidence: () => import('./acts/evidence').then((m) => m.evidence),
};

const ORDER: ActId[] = ACTS.map((a) => a.id);
const TOUR_TOTAL = ACTS.reduce((s, a) => s + a.duration, 0);

export interface ShowController {
  go(anchor: ShowAnchor | null, at?: number): void;
  suspend(): void;
}

declare global {
  interface Window {
    __qchasShow?: {
      act(id: ActId): Promise<void>;
      seek(id: ActId, t: number): Promise<void>;
      play(): void;
      pause(): void;
      tour(): void;
      state(): { act: ActId | null; t: number; playing: boolean; touring: boolean };
    };
  }
}

function memo<T>(fn: () => Promise<T>): () => Promise<T> {
  let p: Promise<T> | null = null;
  return () => (p ??= fn());
}

function makeData(): ShowData {
  const bundle = <K extends Parameters<typeof loadBundle>[0]>(name: K) => memo(() => loadBundle(name).then((r) => (r.ok ? r.data : null)));
  return {
    sky: memo(() => loadSkyMap('lowl').then((r) => (r.ok ? r.data : null))),
    cmb: bundle('cmb.json'),
    resonance: bundle('resonance.json'),
    horizon: bundle('horizon.json'),
    bounce: bundle('bounce.json'),
    topology: bundle('topology.json'),
    validation: bundle('validation.json'),
    manifest: bundle('manifest.json'),
  };
}

const sleep = (ms: number): Promise<void> => new Promise((r) => setTimeout(r, ms));

export function mountShow(parent: HTMLElement): ShowController {
  const reducedMotion = window.matchMedia?.('(prefers-reduced-motion: reduce)').matches ?? false;
  const gl = webglSupported();
  const data = makeData();

  // --- DOM ---------------------------------------------------------------------------------
  const stageEl = h('div', { class: 'show-stage' });
  const hud = h('div', { class: 'show-hud' });
  const navLinks = new Map<ActId, HTMLAnchorElement>();
  const nav = h('nav', { class: 'show-nav', 'aria-label': 'Acts' });
  for (const a of ACTS) {
    const link = h('a', { href: `#${a.id}`, class: 'show-nav-link' }, h('span', { class: 'nav-num' }, a.number), h('span', { class: 'nav-name' }, a.nav), h('span', { class: 'nav-bar' }, h('span', { class: 'nav-fill' })));
    navLinks.set(a.id, link);
    nav.append(link);
  }
  const top = h(
    'header',
    { class: 'show-top show-ui' },
    h('a', { class: 'show-brand', href: '#', 'aria-label': 'Q-CHAS home' }, h('span', { class: 'brand-dot' }), h('span', { class: 'brand-word' }, BRAND.name), h('span', { class: 'brand-sub' }, 'Space Apps')),
    nav,
    h('a', { class: 'show-console', href: '#/overview' }, 'Mission console', h('span', { 'aria-hidden': 'true' }, ' →')),
  );
  const truth = h('div', { class: 'truth' });
  const beatChips = h('div', { class: 'beat-chips' });
  const chips = h('aside', { class: 'show-chips show-ui', 'aria-label': 'What is established, computed or assumed' }, truth, beatChips);
  const capKicker = h('div', { class: 'cap-kicker' });
  const capTitle = h('h2', { class: 'cap-title' });
  const capBody = h('p', { class: 'cap-body' });
  const caption = h('section', { class: 'show-caption', 'aria-live': 'polite' }, capKicker, capTitle, capBody);
  const statement = h('div', { class: 'show-statement', 'aria-live': 'polite' });

  // transport (free mode) ----------------------------------------------------------------------
  const playBtn = button('❚❚', () => togglePlay(), 'show-icon-btn', 'Play or pause');
  const scrub = h('input', { type: 'range', min: '0', max: '1', step: '0.001', value: '0', class: 'show-scrub', 'aria-label': 'Act timeline' }) as HTMLInputElement;
  const timeLabel = h('span', { class: 'show-time mono' }, '0.0 s');
  const tourBtn = button('▶ 60-second tour', () => startTour(), 'show-btn ghost');
  const transport = h('div', { class: 'show-transport' }, playBtn, scrub, timeLabel, tourBtn);

  // tour bar ---------------------------------------------------------------------------------------
  const segs = ACTS.map((a) => h('div', { class: 'tour-seg', style: { flexGrow: String(a.duration) }, title: a.nav }, h('div', { class: 'tour-fill' })));
  const tourPause = button('❚❚ Pause', () => (tourPaused ? resumeTour() : pauseTour()), 'show-btn');
  const tourbar = h(
    'div',
    { class: 'show-tourbar' },
    h('div', { class: 'tour-label' }, 'Guided tour'),
    h('div', { class: 'tour-track' }, ...segs),
    tourPause,
    button('Skip ⏭', () => skipChapter(), 'show-btn ghost'),
    button('Explore from here', () => exitTour(true), 'show-btn ghost'),
    button('✕', () => exitTour(false), 'show-icon-btn', 'Exit tour'),
  );
  const bottom = h('div', { class: 'show-bottom show-ui' }, transport, tourbar);

  // intro --------------------------------------------------------------------------------------------
  const intro = h(
    'div',
    { class: 'show-intro', role: 'dialog', 'aria-label': 'Welcome' },
    h(
      'div',
      { class: 'intro-inner' },
      h('div', { class: 'intro-eyebrow' }, `${BRAND.event} · ${BRAND.name}`),
      h('h1', { class: 'intro-title' }, 'The shape of space.', h('br'), 'The fate of information.', h('br'), 'The first instant.'),
      h('p', { class: 'intro-sub' }, `${BRAND.long}: three new results in mathematical physics, confronted with real NASA data — and honest about what the data can and cannot prove.`),
      h(
        'div',
        { class: 'intro-actions' },
        button('▶  Watch the 60-second tour', () => startTour(), 'show-btn primary big'),
        button('Explore the acts', () => explore(), 'show-btn big'),
        h('a', { class: 'show-btn big ghost', href: '#/overview' }, 'Mission console →'),
      ),
      h('div', { class: 'intro-legend' }, h('span', { class: 'legend-k' }, 'Every claim is tagged'), tagLegend()),
      h('div', { class: 'intro-data' }, 'Real data · NASA WMAP · Planck PR3 · SPHEREx QR3 · LIGO-Virgo strain'),
    ),
  );

  const nogl = h('div', { class: 'show-nogl' });
  const root = h('div', { class: 'show', 'data-act': '' }, stageEl, hud, statement, caption, chips, top, bottom, intro, nogl);
  parent.append(root);

  const stage = gl ? new Stage(stageEl) : null;

  // --- state -----------------------------------------------------------------------------------
  let current: { id: ActId; handle: ActHandle } | null = null;
  let token = 0;
  let touring = false;
  let tourPaused = false;
  let introOpen = false;
  let beatIndex = -1;
  let suspended = false;

  const ctx = (): ActContext => ({
    stage: stage!,
    hud,
    data,
    reducedMotion,
    userInput,
  });

  function userInput(): void {
    if (touring && !tourPaused) pauseTour();
  }

  async function switchTo(id: ActId, opts: { autoplay?: boolean; at?: number } = {}): Promise<void> {
    const my = ++token;
    root.dataset.act = id;
    paintNav(id);
    if (!gl || !stage) {
      showStatic(id);
      return;
    }
    if (current) {
      stage.fadeTo(0);
      await sleep(300);
      if (my !== token) return;
      current.handle.dispose();
      current = null;
      stage.setScene(null);
    }
    hud.replaceChildren();
    root.classList.add('loading');
    let handle: ActHandle;
    try {
      const mod = await LOADERS[id]();
      if (my !== token) return;
      handle = await mod.mount(ctx());
    } catch (err) {
      console.error(err);
      root.classList.remove('loading');
      showStatic(id, String(err));
      return;
    }
    if (my !== token) {
      handle.dispose();
      return;
    }
    root.classList.remove('loading');
    nogl.replaceChildren();
    current = { id, handle };
    handle.onEnd = () => onActEnd(id);
    stage.setScene(handle);
    beatIndex = -1;
    paintTruth(actCopy(id));
    if (opts.at !== undefined) {
      handle.seek(opts.at);
      handle.pause();
    } else if (opts.autoplay) handle.play();
    else if (reducedMotion) handle.seek(handle.duration);
    stage.fadeTo(1);
    stage.start();
    stage.renderNow();
  }

  function onActEnd(id: ActId): void {
    if (!touring || tourPaused) return;
    const next = ORDER[ORDER.indexOf(id) + 1];
    if (next) void goAct(next, { autoplay: true });
    else exitTour(true);
  }

  async function goAct(id: ActId, opts: { autoplay?: boolean; at?: number } = {}): Promise<void> {
    const hash = `#${id}`;
    if (location.hash !== hash && opts.at === undefined) history.replaceState(null, '', hash);
    await switchTo(id, opts);
  }

  // --- tour ------------------------------------------------------------------------------------------
  /** play an act from its start, reusing it when it is already mounted */
  function playFromStart(id: ActId): void {
    if (current?.id === id) {
      history.replaceState(null, '', `#${id}`);
      current.handle.seek(0);
      current.handle.play();
      beatIndex = -1;
    } else void goAct(id, { autoplay: true });
  }

  function startTour(): void {
    touring = true;
    tourPaused = false;
    closeIntro();
    root.classList.add('touring');
    root.classList.remove('tour-paused');
    tourPause.textContent = '❚❚ Pause';
    playFromStart('discovery');
  }

  function pauseTour(): void {
    tourPaused = true;
    current?.handle.pause();
    tourPause.textContent = '▶ Resume';
    root.classList.add('tour-paused');
  }

  function resumeTour(): void {
    tourPaused = false;
    root.classList.remove('tour-paused');
    tourPause.textContent = '❚❚ Pause';
    const hnd = current?.handle;
    if (!hnd) return;
    if (hnd.t >= hnd.duration - 1e-3) onActEnd(current!.id);
    else hnd.play();
  }

  function skipChapter(): void {
    if (!current) return;
    const next = ORDER[ORDER.indexOf(current.id) + 1];
    if (next) {
      tourPaused = false;
      root.classList.remove('tour-paused');
      tourPause.textContent = '❚❚ Pause';
      void goAct(next, { autoplay: true });
    } else exitTour(true);
  }

  function exitTour(keepPlaying: boolean): void {
    touring = false;
    tourPaused = false;
    root.classList.remove('touring', 'tour-paused');
    if (keepPlaying && current && current.handle.t < current.handle.duration) current.handle.play();
  }

  // --- intro -----------------------------------------------------------------------------------------
  function openIntro(): void {
    introOpen = true;
    root.classList.add('intro-open');
  }

  function closeIntro(): void {
    introOpen = false;
    root.classList.remove('intro-open');
  }

  function explore(): void {
    closeIntro();
    if (reducedMotion) void goAct('discovery', {});
    else playFromStart('discovery');
  }

  // --- chrome -----------------------------------------------------------------------------------------
  function paintNav(id: ActId): void {
    for (const [aid, a] of navLinks) {
      if (aid === id) a.setAttribute('aria-current', 'step');
      else a.removeAttribute('aria-current');
    }
  }

  function paintTruth(copy: ActCopy): void {
    truth.replaceChildren(...copy.truth.map((c) => chip(c)));
    capKicker.textContent = `Act ${copy.number} · ${copy.nav}`;
  }

  function paintBeat(copy: ActCopy, beat: Beat | undefined): void {
    beatChips.replaceChildren(...(beat?.chips ?? []).map((c) => chip(c, 'beat')));
    if (!beat) return;
    if (beat.statement) {
      statement.replaceChildren(h('div', { class: 'statement-kicker' }, beat.title), h('p', { class: 'statement-text' }, beat.body));
      root.classList.add('statement-on');
      caption.classList.add('hidden');
    } else {
      root.classList.remove('statement-on');
      caption.classList.remove('hidden');
      capTitle.textContent = beat.title;
      capBody.textContent = beat.body;
      caption.classList.remove('enter');
      void caption.offsetWidth; // restart the entrance animation
      caption.classList.add('enter');
    }
    capKicker.textContent = `Act ${copy.number} · ${copy.nav}`;
  }

  function afterFrame(): void {
    if (!current) return;
    const { id, handle } = current;
    const copy = actCopy(id);
    let bi = 0;
    copy.beats.forEach((b, i) => {
      if (handle.t >= b.at - 1e-6) bi = i;
    });
    if (bi !== beatIndex) {
      beatIndex = bi;
      paintBeat(copy, copy.beats[bi]);
    }
    const f = handle.duration > 0 ? handle.t / handle.duration : 0;
    if (document.activeElement !== scrub) scrub.value = String(f);
    syncFill(scrub);
    timeLabel.textContent = `${handle.t.toFixed(1)} / ${handle.duration.toFixed(0)} s`;
    playBtn.textContent = handle.playing ? '❚❚' : '▶';
    for (const [aid, a] of navLinks) {
      const fill = a.querySelector<HTMLElement>('.nav-fill');
      if (fill) fill.style.transform = `scaleX(${aid === id ? f : ORDER.indexOf(aid) < ORDER.indexOf(id) && touring ? 1 : 0})`;
    }
    if (touring) {
      const idx = ORDER.indexOf(id);
      segs.forEach((s, i) => {
        const fill = s.firstElementChild as HTMLElement;
        fill.style.transform = `scaleX(${i < idx ? 1 : i === idx ? f : 0})`;
      });
    }
  }
  if (stage) stage.afterFrame = afterFrame;

  function togglePlay(): void {
    const hnd = current?.handle;
    if (!hnd) return;
    if (touring) {
      if (tourPaused) resumeTour();
      else pauseTour();
      return;
    }
    if (hnd.playing) hnd.pause();
    else hnd.play();
  }

  scrub.addEventListener('input', () => {
    const hnd = current?.handle;
    if (!hnd) return;
    userInput();
    hnd.pause();
    hnd.seek(Number(scrub.value) * hnd.duration);
  });

  // --- static fallback (no WebGL) ------------------------------------------------------------------
  function showStatic(id: ActId, error?: string): void {
    const copy = actCopy(id);
    paintTruth(copy);
    if (id === 'evidence' && !error) {
      void LOADERS.evidence().then((m) => m.mount({ stage: null as unknown as Stage, hud, data, reducedMotion, userInput }));
      return;
    }
    nogl.replaceChildren(
      h(
        'div',
        { class: 'nogl-card' },
        h('div', { class: 'cap-kicker' }, `Act ${copy.number} · ${copy.nav}`),
        h('h2', { class: 'cap-title' }, copy.title),
        ...copy.beats.map((b) => h('div', { class: 'nogl-beat' }, h('h3', null, b.title), h('p', null, b.body))),
        h('p', { class: 'nogl-note' }, error ? `This act could not start (${error}).` : 'This browser has no WebGL, so the animation is replaced by its script.', ' ', h('a', { href: '#evidence' }, 'See the evidence →')),
      ),
    );
  }

  // --- keyboard ---------------------------------------------------------------------------------------
  function onKey(e: KeyboardEvent): void {
    if (suspended) return;
    const target = e.target as HTMLElement;
    if (target.closest('input, textarea, select') && e.key !== 'Escape') return;
    const idx = current ? ORDER.indexOf(current.id) : -1;
    if (e.key === ' ' && !target.closest('button, a')) {
      e.preventDefault();
      togglePlay();
    } else if (e.key === 'ArrowRight' && idx < ORDER.length - 1) {
      userInput();
      void goAct(ORDER[idx + 1]!, { autoplay: !reducedMotion || touring });
    } else if (e.key === 'ArrowLeft' && idx > 0) {
      userInput();
      void goAct(ORDER[idx - 1]!, { autoplay: !reducedMotion || touring });
    } else if (e.key === 'Escape') {
      if (touring) exitTour(false);
      else if (introOpen) explore();
    }
  }
  window.addEventListener('keydown', onKey);

  // --- QA hook ------------------------------------------------------------------------------------------
  window.__qchasShow = {
    act: (id) => goAct(id, { autoplay: false }),
    seek: async (id, t) => {
      if (current?.id !== id) await switchTo(id, { at: t });
      else {
        current.handle.pause();
        current.handle.seek(t);
      }
      stage?.renderNow();
    },
    play: () => current?.handle.play(),
    pause: () => current?.handle.pause(),
    tour: () => startTour(),
    state: () => ({ act: current?.id ?? null, t: current?.handle.t ?? 0, playing: current?.handle.playing ?? false, touring }),
  };

  return {
    go(anchor, at) {
      suspended = false;
      root.hidden = false;
      document.title = 'Q-CHAS · The shape of space, the fate of information';
      stage?.start();
      if (anchor === null) {
        if (touring) return;
        openIntro();
        if (current?.id !== 'discovery') void switchTo('discovery', {});
        return;
      }
      closeIntro();
      if (current?.id === anchor && at === undefined) return;
      if (at !== undefined) void switchTo(anchor, { at });
      else void switchTo(anchor, { autoplay: !reducedMotion || touring });
    },
    suspend() {
      suspended = true;
      if (touring) exitTour(false);
      current?.handle.pause();
      stage?.stop();
      root.hidden = true;
    },
  };
}

export const TOUR_SECONDS = TOUR_TOTAL;
