# Q-CHAS · The Show (judge-facing experience)

The mission console stays as the **Evidence** layer. In front of it sits a cinematic,
interactive **show** in four acts that a judge understands in 60 seconds without the
mathematics. Every frame is driven by the real bundles in `web/public/data/` or by exact
closed-form physics; nothing is decorative fiction, and every claim carries an honesty tag.

## Audience, job, rules

* Audience: NASA Space Apps judges, on a laptop, often on a projector, 1–3 minutes each.
* Job: make them *feel* the idea (space may be a finite, repeating crystal; black holes give
  information back; the Big Bang bounces), then *trust* it (real NASA data, reproduced paper
  numbers, honest limits).
* Honesty tags on every scientific statement in the show, rendered as chips:
  `ESTABLISHED` (textbook / measured), `PROVED` (theorem in the papers), `COMPUTED`
  (we computed it; reproduced paper numbers), `HYPOTHESIS` (principle or assumption the
  result rests on), `ILLUSTRATIVE` (a visual exaggeration of scale; say which).
* Binding facts (never contradict): no torus evidence and no resonance in Planck (p = 0.78) or in
  our WMAP re-analysis (p = 0.79); data require d ≥ 0.985 D (circle searches; FCC ≥ 1.04 D in the
  paper's Planck likelihood); the 64× resonance gain (ours 65.0×) is **at known orientation**; a
  torus as large as 3 D needs a **multi-field primordial line mapped in 3D**. Any visual at
  d < 0.985 D is labelled "illustrative size, already excluded by data".

## Acts and beats

All copy lives in `web/src/show/script.ts` (single source; the science reviewer owns it).
Durations are for the guided tour; every act is also freely explorable.

### Act 1 · Discovery — "Space may be a crystal" (tour ≈ 18 s)
1. *Inside the oldest light.* Camera at the centre of a sphere textured with the real WMAP 9-yr
   CMB map (`sky_map_lowl.bin`, 2 ≤ ℓ ≤ 30, diverging CMB colormap, gentle rotation, faint
   starfield beyond). Caption: the sky NASA's WMAP measured; the light is 13.8 Gyr old.
2. *Light wraps around.* Camera pulls back out through the sphere wall (sphere fades to a
   translucent glowing shell). Copies of the same shell appear at FCC lattice translations
   (12 nearest neighbours first, then the next shells, staggered "ignite" animation). Each copy
   carries the **same** CMB texture (it is the same sky seen again). Where two shells intersect,
   the **exact** intersection circle ignites in seam-cyan: plane x·ĝ = |g|/2, radius
   √(χ*² − |g|²/4). Use the illustrative size d = 0.8 D (tag ILLUSTRATIVE).
3. *One shape survives.* The fundamental domain appears at the centre: the Voronoi cell of the
   lattice (FCC → rhombic dodecahedron; BCC → truncated octahedron; cubic → cube) as a glowing
   wireframe; a quick counter "18 flat shapes → 1" with the four principles as chips
   (P1 chirality ESTABLISHED, P2 finite thermodynamics HYPOTHESIS, P3 local Copernican HYPOTHESIS,
   P4 supersymmetry HYPOTHESIS·optional; Main Theorem PROVED). Slow orbit of the crystal.
4. *The seams on our sky.* Fly back to a single globe seen from outside. Matched circle pairs
   appear on it (physics/circles.ts), twin points pulse and are joined by luminous chords
   through the observer. Then d animates 0.8 D → 1.05 D: the circles shrink and vanish at d = D.
   Discovery statement (one sentence, large): *"Q-CHAS turns a theorem about the shape of space
   into an exact map of where its seams would lie on NASA's sky — and shows what the data can
   and cannot prove."* Reality chip: no seams in the data; d ≥ 0.985 D.

### Act 2 · Experiment — "Hunt for the seams" (tour ≈ 12 s, then free play)
* Left: the real-sky globe with live circles. Right: the **k-space crystal**: the dual-lattice
  wavevectors inside |κ|χ* ≤ ~14 as glowing points (size ∝ primordial power); a slow rotation.
* Controls (judge-sized, keyboard accessible): universe size d/D 0.6 → 1.5 (excluded range
  hatched), shape cubic / BCC / FCC, and a two-state switch **Standard search** (smooth spectrum,
  unknown orientation) ↔ **Resonance shell** (a narrow primordial line on one dual-lattice shell:
  6 vectors along cube axes / 8 along ⟨111⟩ / 12 along ⟨110⟩, known orientation).
* On "Resonance": the chosen shell's points flare into beams through the origin; the globe
  shows (blended over the real map, tagged ILLUSTRATIVE · hypothetical) the plane-wave pattern
  Σ cos(κ·χ* n̂ + φκ) that such a line would imprint; nothing like it is in the data.
* Live meter (fair comparison, same torus, same size): information K (nats, log axis) for the
  smooth spectrum from `cmb.kl` (interpolated in d) and with the line from
  `resonance.kl_table` (cubic, f = 0.5, |m|² = 16 shell; interpolate within 1.05–1.3 D; outside,
  show "not computed"), the KL = 1 detectability line, and the gain (×3 … ×65). Third bar,
  greyed: unknown orientation → most of the gain is lost (search p = 0.78 Planck / 0.79 WMAP).
* One-line takeaway: *"A narrow primordial line would make a torus up to 65× more visible — if
  we knew its orientation. Neither Planck nor WMAP shows one."*

### Act 3 · Breakthrough (the visual climax)
**A. "A black hole gives its information back"** (tour ≈ 18 s, then scrubbable)
* Scene: an evaporating Schwarzschild black hole against a procedural deep starfield and galactic
  band, with **screen-space gravitational lensing** (ray-deflection shader: Einstein ring,
  secondary images) and a bright photon ring; no accretion disk (an evaporating hole has none).
* The horizon is tiled with glowing hexagonal **Planck cells** (one quarter nat each: capacity
  A/4G), shader-based so tiles flicker gold; ILLUSTRATIVE cell size.
* Hawking quanta: particles born just outside the horizon stream outward (blue-white); each has
  an **infalling partner** inside (dim red) joined by a thin gold **entanglement thread**.
* Timeline t ∈ [0, 1] of the lifetime (auto-plays, judge can scrub): radius ∝ m = (1 − t)^{1/3};
  horizon cells go dark as area shrinks; the radiation cloud grows. HUD draws the curves live
  from physics/page.ts: Hawking's curve (field theory, keeps rising to 1.855 S_BH) and the Page
  curve min{S_H, A/4G}, turning over at **t = 0.476** (80.6 % of the mass left).
* At the Page time an **island** ignites: a magenta shell just inside the horizon
  (ILLUSTRATIVE depth; real depth 0.0399 c/S_BH of r_s), and threads re-route: partners of the
  earliest radiation migrate from the interior into the radiation cloud, so the radiation
  becomes entangled with itself (thread colour gold → cyan). At t → 1 the hole vanishes in a
  flash; every thread now lives inside the radiation: S(R) → 0.
* Interaction: **"Throw in a diary"**: a glowing cube falls in; its glyph re-emerges in the
  radiation according to the exact 20-qubit run (`horizon.json → evaporation`): thrown in after
  5 emissions → recovered at step 13 (after the Page time); after 14 → recovered at step 17
  (3 emissions later: Hayden–Preskill). Map the 21 emission steps onto the timeline.
* Truth panel (three chips, always visible): ESTABLISHED: black holes radiate; field theory on a
  fixed background forces the radiation entropy to keep rising ([BH] Theorem B). COMPUTED: an
  exactly unitary 20-qubit evaporation follows the Page curve (10.28 bits at the Page time, back
  to 1 bit; diaries recovered at steps 13 and 17), reproduced. HYPOTHESIS: the horizon has finite
  capacity A/4G and no independent interior (H*); the capacity is derived from the crossed
  product ([BH] Theorem G) and H* is proved for eternal black holes ([BH] Theorem K).

**B. "The Big Bang bounces"** (tour ≈ 8 s)
* A luminous 3D surface of |Ψ(a, T)|² over (a, T) for a_max = 2, computed live from the exact
  solution (odd combination of two coherent states of the half-line oscillator; see
  qchas/horizon/bounce.py), glowing ridge, interference fringes at the bounce.
* Beside it a translucent "universe bubble" whose radius is ⟨a⟩(T): it collapses, **bounces at
  ⟨a⟩_min = 0.163 ℓ ≈ 0.26 ℓ_P** and re-expands; a red ghost shows the classical a_max|cos T|
  crashing to a = 0. Graviton ripples (integer-frequency oscillators) run across the bubble
  unchanged through the bounce.
* Chips: PROVED unitary evolution through the singularity ([QG] Theorem Q); COMPUTED bounce
  scale (exact; the paper's printed table is grid-biased, we show both in Evidence);
  HYPOTHESIS closed radiation-filled universe with the radiation as clock; ESTABLISHED Planck
  shows no bounce feature (consistent: the bounce is transparent).

### Act 4 · Evidence (the console)
* A calm, information-dense landing for the console: what data (WMAP 9-yr, SPHEREx QR3 240
  frames, LIGO strain; Planck PR3 on the GitHub build), what method, **169 of 173 paper numbers
  reproduced** (the 4 differences explained), limitations, and a grid of doors into every
  console tool (#/overview, #/ghost-twin, … #/claims, #/data).

### Guided tour (≈ 60 s)
A "▶ Watch the 60-second tour" button (and auto-offer on first load) plays Act 1 → 2 → 3A → 3B →
4 with captions, chapter progress bar, pause/resume, skip, and "explore from here". Any user
input pauses the tour. Keyboard: Space pause, ←/→ chapters, Esc exit.

## Visual system
* Ground: deep space #03050b → #060a12. Palette tokens (dark-first; light theme only for
  Evidence, the acts stay dark): `--cmb-cold #2f6bff`, `--cmb-hot #ff5a36`, `--seam #59f3ff`
  (topology seams), `--info #ffc35a` (information / Planck cells / threads), `--island #ff4fd8`,
  `--text #e9eef8`, `--text-dim #8b97ad`.
* Type: display **Unbounded** (300–600) for act titles and the discovery statement; body
  **IBM Plex Sans**; data **IBM Plex Mono** (Google Fonts, with system fallbacks).
* Rendering: Three.js + `EffectComposer` (RenderPass → UnrealBloomPass (half res) → OutputPass),
  additive-blended GPU points, custom ShaderMaterials, ACES tone mapping, subtle vignette and
  film grain. Motion: slow, eased, cinematic (GSAP timelines), one orchestrated move at a time.
* Composition: full-bleed canvas, captions in the lower third with a soft gradient scrim,
  chips top-left, act navigation top-centre, "Mission console →" top-right.

## Architecture (`web/src/show/`)
| file | role |
|---|---|
| `show.ts` | landing controller: act nav, lazy-mounts acts, tour button, console link |
| `stage.ts` | shared WebGL stage: renderer, composer + bloom, DPR/quality tiers, resize, visibility pause, dispose, reduced-motion & no-WebGL fallbacks |
| `tour.ts` | timeline engine (chapters, beats, captions, progress UI, pause/skip, input-pauses) |
| `script.ts` | **all copy**: captions, statements, chip texts and tags, sources |
| `ui.ts`, `show.css` | captions, chips, controls, act nav; design tokens above |
| `acts/discovery.ts` | Act 1 |
| `acts/experiment.ts` | Act 2 |
| `acts/breakthrough-blackhole.ts` | Act 3A |
| `acts/breakthrough-bounce.ts` | Act 3B |
| `acts/evidence.ts` | Act 4 |
| `shaders/*.ts` | GLSL strings (lensing, Planck cells, CMB sphere, glow lines, points) |

Each act exports `{ id, title, mount(container, stage, data) -> ActHandle }` with
`ActHandle = { seek(t: number): void; play(): void; pause(): void; setInteractive(on): void;
dispose(): void; duration: number }` so the tour can drive it and QA can screenshot any moment.

Routing: the show is the landing page (empty hash) and owns plain anchors `#discovery`,
`#experiment`, `#breakthrough`, `#evidence` (artifact hosts pass only plain anchors). The console
keeps its `#/route` hashes. QA hook: `window.__qchasShow = { act(id), seek(id, t), play(), pause() }`.

## Performance and robustness
* 60 fps on an integrated laptop GPU at 1440×900: DPR ≤ 1.5 (quality tier lowers to 1.0 and
  halves particles if frame time > 22 ms for 1 s), bloom at half resolution, ≤ 60 k particles.
* Render only the visible act; pause when the tab is hidden; dispose geometries, materials,
  textures and composers on act change (no leaks across 20 act switches).
* `prefers-reduced-motion`: no auto-play, static keyframes, crossfades instead of flights.
* No WebGL: a designed static fallback per act (pre-rendered SVG/CSS) plus the Evidence console.
* Lazy-load: the show shell ships light; each act's code and three.js load on demand.
