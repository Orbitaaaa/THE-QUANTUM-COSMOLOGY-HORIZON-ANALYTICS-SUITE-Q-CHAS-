/**
 * The shared WebGL stage of the show: one renderer, one post-processing chain
 * (RenderPass -> UnrealBloomPass at half resolution -> OutputPass (ACES, sRGB) -> finish pass
 * with vignette, film grain and a hint of chromatic aberration), a frame loop that drives the
 * active act, DPR/quality tiers, resize and tab-visibility handling.
 */

import * as THREE from 'three';
import { EffectComposer } from 'three/examples/jsm/postprocessing/EffectComposer.js';
import { OutputPass } from 'three/examples/jsm/postprocessing/OutputPass.js';
import { RenderPass } from 'three/examples/jsm/postprocessing/RenderPass.js';
import { ShaderPass } from 'three/examples/jsm/postprocessing/ShaderPass.js';
import { UnrealBloomPass } from 'three/examples/jsm/postprocessing/UnrealBloomPass.js';

export interface BloomSettings {
  strength: number;
  radius: number;
  threshold: number;
}

/** What an act hands the stage: a scene, a camera and a per-frame callback. */
export interface StageScene {
  scene: THREE.Scene;
  camera: THREE.PerspectiveCamera;
  bloom?: BloomSettings;
  /** exposure for the ACES tone mapping (default 1) */
  exposure?: number;
  /** advance by dt seconds of wall time and update the scene */
  frame(dt: number): void;
  resize?(width: number, height: number, pixelRatio: number): void;
}

const FINISH_SHADER = {
  uniforms: {
    tDiffuse: { value: null as THREE.Texture | null },
    uTime: { value: 0 },
    uVignette: { value: 0.9 },
    uGrain: { value: 0.045 },
    uAberration: { value: 0.0007 },
    uFade: { value: 1 },
  },
  vertexShader: /* glsl */ `
    varying vec2 vUv;
    void main() { vUv = uv; gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0); }`,
  fragmentShader: /* glsl */ `
    uniform sampler2D tDiffuse;
    uniform float uTime, uVignette, uGrain, uAberration, uFade;
    varying vec2 vUv;
    float hash(vec2 p) { p = fract(p * vec2(123.34, 456.21)); p += dot(p, p + 45.32); return fract(p.x * p.y); }
    void main() {
      vec2 c = vUv - 0.5;
      float r2 = dot(c, c);
      vec2 off = c * uAberration * (1.0 + 4.0 * r2);
      vec3 col;
      col.r = texture2D(tDiffuse, vUv + off).r;
      col.g = texture2D(tDiffuse, vUv).g;
      col.b = texture2D(tDiffuse, vUv - off).b;
      float vig = smoothstep(0.95, 0.18, r2 * uVignette * 2.2);
      col *= mix(1.0, vig, 0.85);
      float g = hash(vUv * 1000.0 + fract(uTime * 7.13)) - 0.5;
      col += g * uGrain * (0.35 + 0.65 * (1.0 - dot(col, vec3(0.3333))));
      gl_FragColor = vec4(col * uFade, 1.0);
    }`,
};

export function webglSupported(): boolean {
  try {
    const c = document.createElement('canvas');
    return !!(c.getContext('webgl2') || c.getContext('webgl'));
  } catch {
    return false;
  }
}

export class Stage {
  readonly canvas: HTMLCanvasElement;
  readonly renderer: THREE.WebGLRenderer;
  private readonly composer: EffectComposer;
  private readonly renderPass: RenderPass;
  private readonly bloom: UnrealBloomPass;
  private readonly finish: ShaderPass;
  private readonly ro: ResizeObserver;
  private current: StageScene | null = null;
  private raf = 0;
  private last = 0;
  private running = false;
  private elapsed = 0;
  private maxDpr = 1.5;
  private slowFrames = 0;
  private width = 1;
  private height = 1;
  /** crossfade of the whole frame (0 black, 1 full) */
  fade = 1;
  private fadeTarget = 1;
  /** called after every frame (the show updates captions and progress here) */
  afterFrame: (() => void) | null = null;

  constructor(private readonly container: HTMLElement) {
    this.renderer = new THREE.WebGLRenderer({ antialias: true, alpha: false, powerPreference: 'high-performance', preserveDrawingBuffer: true });
    this.renderer.setClearColor(0x03050b, 1);
    this.renderer.toneMapping = THREE.ACESFilmicToneMapping;
    this.renderer.toneMappingExposure = 1;
    this.renderer.outputColorSpace = THREE.SRGBColorSpace;
    this.canvas = this.renderer.domElement;
    this.canvas.className = 'show-canvas';
    this.canvas.setAttribute('aria-hidden', 'true');
    container.append(this.canvas);

    const empty = new THREE.Scene();
    const cam = new THREE.PerspectiveCamera();
    this.composer = new EffectComposer(this.renderer);
    this.renderPass = new RenderPass(empty, cam);
    this.bloom = new UnrealBloomPass(new THREE.Vector2(256, 256), 1.1, 0.6, 0.12);
    // bloom at half resolution: the blur mip chain starts at a quarter of the frame
    const bloomSetSize = this.bloom.setSize.bind(this.bloom);
    this.bloom.setSize = (w: number, h: number) => bloomSetSize(Math.max(2, Math.round(w / 2)), Math.max(2, Math.round(h / 2)));
    this.finish = new ShaderPass(FINISH_SHADER);
    this.composer.addPass(this.renderPass);
    this.composer.addPass(this.bloom);
    this.composer.addPass(new OutputPass());
    this.composer.addPass(this.finish);

    this.ro = new ResizeObserver(() => this.resize());
    this.ro.observe(container);
    this.resize();
    document.addEventListener('visibilitychange', this.onVisibility);
  }

  get size(): { width: number; height: number } {
    return { width: this.width, height: this.height };
  }

  get pixelRatio(): number {
    return this.renderer.getPixelRatio();
  }

  setScene(s: StageScene | null): void {
    this.current = s;
    if (s) {
      this.renderPass.scene = s.scene;
      this.renderPass.camera = s.camera;
      this.resize();
    }
  }

  /** fade the whole frame towards black (0) or full (1) over ~0.35 s */
  fadeTo(target: number): void {
    this.fadeTarget = target;
  }

  start(): void {
    if (this.running) return;
    this.running = true;
    this.last = performance.now();
    this.raf = requestAnimationFrame(this.loop);
  }

  stop(): void {
    this.running = false;
    cancelAnimationFrame(this.raf);
  }

  /** render one frame immediately (used after a seek while paused) */
  renderNow(): void {
    this.draw(0);
  }

  private readonly onVisibility = (): void => {
    if (document.hidden) cancelAnimationFrame(this.raf);
    else if (this.running) {
      this.last = performance.now();
      this.raf = requestAnimationFrame(this.loop);
    }
  };

  private readonly loop = (now: number): void => {
    if (!this.running) return;
    const dt = Math.min(0.1, Math.max(0, (now - this.last) / 1000));
    this.last = now;
    this.watchPerformance(dt);
    this.draw(dt);
    this.raf = requestAnimationFrame(this.loop);
  };

  private draw(dt: number): void {
    const s = this.current;
    this.elapsed += dt;
    this.fade += (this.fadeTarget - this.fade) * Math.min(1, dt * 7);
    if (Math.abs(this.fadeTarget - this.fade) < 0.002) this.fade = this.fadeTarget;
    if (!s) {
      this.renderer.setRenderTarget(null);
      this.renderer.clear();
      this.afterFrame?.();
      return;
    }
    s.frame(dt);
    this.afterFrame?.();
    const b = s.bloom ?? { strength: 1.1, radius: 0.6, threshold: 0.12 };
    this.bloom.strength = b.strength;
    this.bloom.radius = b.radius;
    this.bloom.threshold = b.threshold;
    this.renderer.toneMappingExposure = s.exposure ?? 1;
    const u = this.finish.uniforms as typeof FINISH_SHADER.uniforms;
    u.uTime.value = this.elapsed;
    u.uFade.value = this.fade;
    this.composer.render(dt);
  }

  /** drop to DPR 1 when frames take longer than 22 ms for about a second */
  private watchPerformance(dt: number): void {
    if (this.maxDpr <= 1 || dt === 0) return;
    this.slowFrames = dt > 0.022 && dt < 0.1 ? this.slowFrames + 1 : Math.max(0, this.slowFrames - 2);
    if (this.slowFrames > 50) {
      this.maxDpr = 1;
      this.slowFrames = 0;
      this.resize();
    }
  }

  resize(): void {
    const w = Math.max(1, this.container.clientWidth);
    const h = Math.max(1, this.container.clientHeight);
    const pr = Math.min(window.devicePixelRatio || 1, this.maxDpr);
    this.width = w;
    this.height = h;
    this.renderer.setPixelRatio(pr);
    this.renderer.setSize(w, h, false);
    this.composer.setPixelRatio(pr);
    this.composer.setSize(w, h);
    const s = this.current;
    if (s) {
      s.camera.aspect = w / h;
      s.camera.updateProjectionMatrix();
      s.resize?.(w, h, pr);
    }
  }

  dispose(): void {
    this.stop();
    this.ro.disconnect();
    document.removeEventListener('visibilitychange', this.onVisibility);
    this.composer.dispose();
    this.bloom.dispose();
    this.renderer.dispose();
    this.canvas.remove();
  }
}

/** Dispose every geometry, material and texture under an object. */
export function disposeObject(root: THREE.Object3D): void {
  root.traverse((o) => {
    const m = o as THREE.Mesh;
    m.geometry?.dispose();
    const mats = Array.isArray(m.material) ? m.material : m.material ? [m.material] : [];
    for (const mat of mats) {
      for (const v of Object.values(mat)) if (v instanceof THREE.Texture) v.dispose();
      const uniforms = (mat as THREE.ShaderMaterial).uniforms;
      if (uniforms) for (const u of Object.values(uniforms)) if (u?.value instanceof THREE.Texture) u.value.dispose();
      mat.dispose();
    }
  });
}
