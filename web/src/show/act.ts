/**
 * The act contract of the show (docs/SHOW.md "Architecture"): every act is a pure function of
 * its act time t (so the tour and QA can seek to any moment), plus free-play state once the
 * judge takes over. `Act` implements the clock; subclasses implement `update`.
 */

import * as THREE from 'three';
import type { BounceBundle, CmbBundle, HorizonBundle, Manifest, ResonanceBundle, SkyMapData, TopologyBundle, ValidationBundle } from '../lib/data';
import type { ActId } from './script';
import { type BloomSettings, type Stage, type StageScene, disposeObject } from './stage';

/** Bundles the acts read, loaded once on first use. */
export interface ShowData {
  sky(): Promise<SkyMapData | null>;
  cmb(): Promise<CmbBundle | null>;
  resonance(): Promise<ResonanceBundle | null>;
  horizon(): Promise<HorizonBundle | null>;
  bounce(): Promise<BounceBundle | null>;
  topology(): Promise<TopologyBundle | null>;
  validation(): Promise<ValidationBundle | null>;
  manifest(): Promise<Manifest | null>;
}

export interface ActContext {
  stage: Stage;
  /** overlay layer for this act's HUD and controls (cleared on dispose) */
  hud: HTMLElement;
  data: ShowData;
  reducedMotion: boolean;
  /** report a judge interaction (pauses the guided tour) */
  userInput(): void;
}

export interface ActHandle extends StageScene {
  readonly id: ActId;
  /** scripted length in seconds */
  readonly duration: number;
  /** current act time in seconds */
  readonly t: number;
  readonly playing: boolean;
  seek(t: number): void;
  play(): void;
  pause(): void;
  /** called when the scripted part reaches its end while playing */
  onEnd: (() => void) | null;
  dispose(): void;
}

export interface ActModule {
  mount(ctx: ActContext): Promise<ActHandle>;
}

export abstract class Act implements ActHandle {
  abstract readonly id: ActId;
  abstract readonly duration: number;
  readonly scene = new THREE.Scene();
  readonly camera = new THREE.PerspectiveCamera(50, 1, 0.05, 4000);
  bloom: BloomSettings = { strength: 1.15, radius: 0.62, threshold: 0.1 };
  exposure = 1;
  onEnd: (() => void) | null = null;
  protected time = 0;
  protected isPlaying = false;
  /** wall-clock seconds since mount (for idle motion that is not part of the script) */
  protected wall = 0;
  protected width = 1;
  protected height = 1;
  protected pixelRatio = 1;
  private readonly disposers: (() => void)[] = [];

  constructor(protected readonly ctx: ActContext) {
    this.scene.background = new THREE.Color(0x03050b);
  }

  get t(): number {
    return this.time;
  }

  get playing(): boolean {
    return this.isPlaying;
  }

  frame(dt: number): void {
    this.wall += dt;
    if (this.isPlaying) {
      this.time = Math.min(this.duration, this.time + dt);
      if (this.time >= this.duration) {
        this.isPlaying = false;
        this.onEnd?.();
      }
    }
    this.update(dt);
  }

  seek(t: number): void {
    this.time = Math.min(this.duration, Math.max(0, t));
    this.update(0);
  }

  play(): void {
    if (this.time >= this.duration - 1e-6) this.time = 0;
    this.isPlaying = true;
  }

  pause(): void {
    this.isPlaying = false;
  }

  resize(width: number, height: number, pixelRatio: number): void {
    this.width = width;
    this.height = height;
    this.pixelRatio = pixelRatio;
  }

  /** register a cleanup that runs on dispose */
  protected own(fn: () => void): void {
    this.disposers.push(fn);
  }

  /** listen on an element for the lifetime of the act */
  protected listen<K extends keyof HTMLElementEventMap>(el: HTMLElement | Window, type: K, fn: (e: HTMLElementEventMap[K]) => void, opts?: AddEventListenerOptions): void {
    el.addEventListener(type, fn as EventListener, opts);
    this.own(() => el.removeEventListener(type, fn as EventListener, opts));
  }

  /** update the scene for the current time; dt is wall time (0 after a seek) */
  protected abstract update(dt: number): void;

  dispose(): void {
    for (const fn of this.disposers.splice(0)) {
      try {
        fn();
      } catch (err) {
        console.error(err);
      }
    }
    disposeObject(this.scene);
    this.ctx.hud.replaceChildren();
  }
}

/**
 * Drag-to-orbit for free play: pointer drags add yaw/pitch offsets with inertia. Returns the
 * current offsets; `active` is true while a drag is in progress.
 */
export class DragOrbit {
  yaw = 0;
  pitch = 0;
  private vy = 0;
  private vp = 0;
  private dragging = false;
  private lastX = 0;
  private lastY = 0;
  private readonly cleanup: () => void;

  constructor(el: HTMLElement, onInput: () => void, private readonly pitchLimit = 1.2) {
    const down = (e: PointerEvent): void => {
      if ((e.target as HTMLElement).closest('button, input, a, label, .show-ui')) return;
      this.dragging = true;
      this.lastX = e.clientX;
      this.lastY = e.clientY;
      onInput();
    };
    const move = (e: PointerEvent): void => {
      if (!this.dragging) return;
      const dx = (e.clientX - this.lastX) / Math.max(1, el.clientWidth);
      const dy = (e.clientY - this.lastY) / Math.max(1, el.clientHeight);
      this.lastX = e.clientX;
      this.lastY = e.clientY;
      this.vy = -dx * 4;
      this.vp = dy * 3;
      this.yaw += this.vy;
      this.pitch = Math.max(-this.pitchLimit, Math.min(this.pitchLimit, this.pitch + this.vp));
    };
    const up = (): void => {
      this.dragging = false;
    };
    el.addEventListener('pointerdown', down);
    window.addEventListener('pointermove', move);
    window.addEventListener('pointerup', up);
    this.cleanup = () => {
      el.removeEventListener('pointerdown', down);
      window.removeEventListener('pointermove', move);
      window.removeEventListener('pointerup', up);
    };
  }

  get active(): boolean {
    return this.dragging;
  }

  step(dt: number): void {
    if (this.dragging) return;
    const k = Math.exp(-dt * 3);
    this.vy *= k;
    this.vp *= k;
    this.yaw += this.vy * dt * 10;
    this.pitch = Math.max(-this.pitchLimit, Math.min(this.pitchLimit, this.pitch + this.vp * dt * 10));
  }

  dispose(): void {
    this.cleanup();
  }
}
