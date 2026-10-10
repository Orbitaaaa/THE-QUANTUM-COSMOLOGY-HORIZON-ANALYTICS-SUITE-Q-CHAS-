/**
 * 2.5 GW Horizons: real LIGO strain and the remnant horizons ([QG] Secs. 17-19, [BH] Sec. 17).
 *
 * Event selector over gw.json events (LIGO tutorial data from GWOSC); whitened, band-passed H1
 * and L1 strain with L1 shifted by the measured lag and sign; the single-mode ringdown fit with
 * its false-alarm probability; remnant telemetry (Kerr thermodynamics with stated assumptions);
 * the quoted [QG] Sec. 19 nine-event echo table with our recomputed f220 and delay; and the quoted
 * Kerr spin-down table of [QG] Sec. 17. Bundle fields: gw.events[] (all), gw.paper_events,
 * gw.combined_R90, gw.paper_echo_summary, gw.paper_ringdown, gw.echo_bound_check, gw.gwtc;
 * horizon.kerr_spin_down; validation group "Kerr ringdown and echoes ([QG] Sec. 19)".
 */

import './m2.css';
import type { Layout } from 'plotly.js';
import { type ChartHandle, type LineSeries, lineChart, scatterChart } from '../lib/charts';
import type { GwBundle, GwEvent, GwPaperEvent, KerrSpinDownRow } from '../lib/data';
import type { PanelContext } from '../lib/panel';
import { type Child, append, badge, clear, fmtNum, fmtSci, grid, h, loadFailure, panel, segmented, statTile, statusPill, table } from '../lib/ui';
import { alignL1, fmtPow10, projectDampedSinusoid } from './horizon-calc';
import { keep, kv, proofList, tilesN, validationBlock } from './m2-ui';

export const title = 'GW Horizons (real LIGO data)';
export const subtitle =
  'Whitened LIGO strain of real binary black-hole mergers, the ringdown of the remnant horizon, and what the crossed-product framework says about it: no echoes, a Kerr spectrum, and a hole that spins down and evaporates.';

/** The paper calls the LIGO tutorial's LVT151012 by its catalogue name. */
const PAPER_NAME: Record<string, string> = { LVT151012: 'GW151012' };
const paperName = (n: string): string => PAPER_NAME[n] ?? n;

const pm = (v: number | null | undefined, e: number | null | undefined, d = 1): string =>
  v === null || v === undefined ? '—' : e === null || e === undefined ? v.toFixed(d) : `${v.toFixed(d)} ± ${e.toFixed(d)}`;

export async function mount(el: HTMLElement, ctx: PanelContext): Promise<void> {
  ctx.setHeaderExtras([badge('measured', 'LIGO strain (GWOSC)'), badge('quoted', '[QG] Secs. 17–19'), badge('computed', 'qchas')]);
  const [gwRes, hor, val] = await Promise.all([ctx.load('gw.json'), ctx.load('horizon.json'), ctx.load('validation.json')]);
  if (ctx.unmounted) return;
  if (!gwRes.ok) {
    const p = panel({ title: 'LIGO events', subtitle: 'gw.json' });
    p.body.append(loadFailure(gwRes));
    el.append(grid(p.root));
    return;
  }
  const gw: GwBundle = gwRes.data;
  const events = gw.events;

  // --- event selector -------------------------------------------------------------------------------------
  const selP = panel({
    title: 'Event',
    subtitle: `${events.length} binary black-hole mergers with H1 and L1 strain · source: ${events[0]?.source ?? 'GWOSC'}`,
    badges: [badge('measured', 'LIGO H1, L1')],
  });
  let current = events[0]?.name ?? '';
  const seg = segmented({
    label: 'LIGO event',
    options: events.map((e) => ({ value: e.name, label: e.name })),
    value: current,
    onChange: (v) => {
      current = v;
      renderEvent();
    },
  });
  const selStatus = h('div', { class: 'm2-pill-row', style: { marginBottom: '0' } });
  selP.body.append(h('div', { class: 'm2-pill-row' }, seg.root), selStatus);
  el.append(grid(selP.root));

  // --- strain + proves ------------------------------------------------------------------------------------
  const strainP = panel({
    title: 'Whitened strain around the merger',
    subtitle: 'H1 and L1, L1 shifted by the measured inter-detector lag and sign; band-passed for display',
    badges: [badge('measured', 'LIGO'), badge('computed', 'whitening, lag')],
    span: 8,
  });
  const zoomSeg = segmented({
    label: 'Time window',
    options: [
      { value: 'merger', label: 'Merger −140…+60 ms' },
      { value: 'full', label: 'Full ±200 ms' },
    ],
    value: 'merger',
    onChange: (v) => {
      zoom = v as 'merger' | 'full';
      drawStrain();
    },
  });
  let zoom: 'merger' | 'full' = 'merger';
  const strainHost = h('div');
  const strainCaption = h('div', { class: 'chart-caption' });
  strainP.body.append(h('div', { class: 'm2-pill-row' }, zoomSeg.root), strainHost, strainCaption);

  const provesP = panel({ title: 'What the paper proves', subtitle: '[QG] Secs. 17–20; what this page measures', span: 4, stretch: true });
  provesP.body.append(
    proofList(
      [
        { kind: 'computed', src: '[QG] Sec. 18 (2)', text: 'Prediction: no gravitational-wave echoes at any reachable amplitude. The horizon’s microstructure is a finite completion of dimension e^S, so a reflected wave has power |R|² ≤ e^(−S); the ringdown is the Kerr quasi-normal spectrum.' },
        { kind: 'measured', src: '[QG] Sec. 19.1', text: 'The paper’s own echo search in nine LIGO–Virgo mergers finds no echo: R < 0.28 combined (90 %, CLs), p = 0.78. Quoted below; Q-CHAS does not rerun that search.' },
        { kind: 'measured', src: '[QG] Sec. 19.2', text: 'Ringdowns of GW150914 and GW250114 agree with Kerr from late enough start times.' },
        { kind: 'proved', src: '[QG] Theorem K′', text: 'A Kerr hole spins down as it evaporates (photons, gravitons at every spin); scalar radiation alone drives it to a* = 0.555. Teukolsky computation, quoted below.' },
        { kind: 'measured', src: 'Q-CHAS', text: 'Here: the strain of four LIGO tutorial events, whitened by us, and a single-mode ringdown fit with a false-alarm probability from off-source noise. Only GW150914 gives a significant ringdown.' },
      ],
      true,
    ),
  );
  el.append(grid(strainP.root, provesP.root));

  // --- ringdown + telemetry ------------------------------------------------------------------------------------
  const rdP = panel({
    title: 'Ringdown fit',
    subtitle: 'Single damped sinusoid fitted to H1 after the peak, 20–1000 Hz, started 1–6 ms after the peak',
    badges: [badge('measured', 'H1'), badge('computed', 'fit, off-source p')],
    span: 6,
    stretch: true,
  });
  const telP = panel({
    title: 'Remnant telemetry',
    subtitle: 'The final black hole: Kerr numbers from the tutorial’s final mass and spin',
    badges: [badge('computed', 'Kerr formulas')],
    span: 6,
    stretch: true,
  });
  el.append(grid(rdP.root, telP.root));
  let strainChart: ChartHandle | null = null;
  let rdChart: ChartHandle | null = null;

  const ev = (): GwEvent => events.find((e) => e.name === current) ?? events[0]!;

  function drawStrain(): void {
    const e = ev();
    strainChart?.destroy();
    clear(strainHost);
    const s = e.strain;
    if (!s) {
      strainHost.append(h('div', { class: 'empty-state' }, `No strain exported for ${e.name}.`));
      strainCaption.textContent = '';
      return;
    }
    const fit = e.ringdown_fit;
    const tevent = e.tevent ?? 0;
    let tPeak = fit?.t_peak_gps !== undefined && e.tevent !== undefined ? fit.t_peak_gps - tevent : NaN;
    if (!Number.isFinite(tPeak)) {
      let iMax = 0;
      s.H1.forEach((v, i) => Math.abs(v) > Math.abs(s.H1[iMax]!) && (iMax = i));
      tPeak = s.t[iMax]!;
    }
    const ms = (t: number): number => (t - tPeak) * 1e3;
    const l1 = alignL1(s.t, s.L1, s.lag_ms ?? 0, s.lag_sign ?? 1);
    const series: LineSeries[] = [
      { name: 'H1 (Hanford)', x: s.t.map(ms), y: s.H1, color: 1, width: 1.5 },
      { name: `L1 (Livingston), shifted ${fmtNum(-(s.lag_ms ?? 0), 3)} ms${(s.lag_sign ?? 1) < 0 ? ', inverted' : ''}`, x: l1.t.map(ms), y: l1.y, color: 2, width: 1.5 },
    ];
    let proj: ReturnType<typeof projectDampedSinusoid> = null;
    if (fit && Number.isFinite(fit.f_Hz) && Number.isFinite(fit.tau_ms)) {
      const t0 = tPeak + fit.t0_ms_after_peak / 1e3;
      proj = projectDampedSinusoid(s.t, s.H1, fit.f_Hz, fit.tau_ms / 1e3, t0, (fit.window_ms ?? 20) / 1e3);
      if (proj)
        series.push({
          name: fit.significant ? `ringdown fit f = ${fmtNum(fit.f_Hz, 4)} Hz, τ = ${fmtNum(fit.tau_ms, 3)} ms` : `ringdown fit (not significant, p = ${fmtNum(fit.p_noise ?? NaN, 2)})`,
          x: proj.t.map(ms),
          y: proj.y,
          color: fit.significant ? 4 : 'context',
          dash: fit.significant ? 'solid' : 'dash',
          width: 2.5,
        });
    }
    strainChart = keep(
      ctx,
      lineChart(strainHost, series, {
        xTitle: `time from the H1 peak (ms)`,
        yTitle: 'whitened strain (noise σ before band-pass)',
        height: 360,
        xRange: zoom === 'merger' ? [-140, 60] : undefined,
        refLines: [{ axis: 'x', value: 0, label: 'peak' }],
        format: (v) => fmtNum(v, 4),
        ariaLabel: `Whitened H1 and L1 strain of ${e.name}`,
      }),
    );
    strainCaption.textContent =
      `Strain ${s.units ?? 'whitened'}, band-passed ${s.band[0]}–${s.band[1]} Hz for display, decimated to ${fmtNum(s.fs_export ?? s.fs, 4)} Hz. ` +
      `L1 lag ${fmtNum(s.lag_ms ?? NaN, 3)} ms (${s.lag_convention ?? 'positive: L1 later'}), correlation ${fmtNum(s.lag_corr ?? NaN, 2)}` +
      (Math.abs(s.lag_corr ?? 0) < 0.5 ? ': a weak H1–L1 correlation, so the lag and sign are not reliable for this quiet event. ' : '. ') +
      (proj && fit
        ? `Overlay: a damped sinusoid with the fitted f and τ (fit made on 20–1000 Hz data), started ${fmtNum(fit.t0_ms_after_peak, 3)} ms after the peak; its amplitude and phase are matched live to the displayed H1 trace by least squares.`
        : '');
  }

  function drawRingdown(): void {
    const e = ev();
    const fit = e.ringdown_fit;
    rdChart?.destroy();
    clear(rdP.body);
    if (!fit) {
      rdP.body.append(h('div', { class: 'empty-state' }, `No ringdown fit for ${e.name}.`));
      return;
    }
    const kf = fit.kerr_f220_Hz ?? e.f220_Hz;
    const kt = fit.kerr_tau220_ms ?? e.tau220_ms;
    append(
      rdP.body,
      h(
        'div',
        { class: 'm2-pill-row' },
        fit.significant
          ? statusPill('good', `significant ringdown · p_noise = ${fmtNum(fit.p_noise ?? NaN, 2)}`)
          : statusPill('warning', `not significant · p_noise = ${fmtNum(fit.p_noise ?? NaN, 2)}: indistinguishable from noise`),
        h('span', { class: 'small muted' }, `mode SNR ${fmtNum(fit.snr ?? NaN, 3)}, ${fit.n_noise ?? 0} off-source noise trials`),
      ),
      tilesN(
        2,
        statTile({ label: 'Fitted frequency', value: pm(fit.f_Hz, fit.f_err_Hz, 1), unit: 'Hz', badge: badge('measured'), sub: fit.f_pull !== null && fit.f_pull !== undefined ? `pull vs Kerr ${fmtNum(fit.f_pull, 2)}σ` : null }),
        statTile({ label: 'Fitted damping time', value: pm(fit.tau_ms, fit.tau_err_ms, 2), unit: 'ms', badge: badge('measured'), sub: fit.tau_pull !== null && fit.tau_pull !== undefined ? `pull vs Kerr ${fmtNum(fit.tau_pull, 2)}σ` : null }),
        statTile({ label: 'Kerr f₂₂₀ (prediction)', value: fmtNum(kf, 4), unit: 'Hz', badge: badge('computed', 'BCW fit') }),
        statTile({ label: 'Kerr τ₂₂₀ (prediction)', value: fmtNum(kt, 3), unit: 'ms', badge: badge('computed', 'BCW fit') }),
      ),
    );
    const host = h('div', { style: { marginTop: '12px' } });
    rdP.body.append(host);
    const ps = fit.per_start;
    if (ps && ps.t0_ms.length) {
      rdChart = keep(
        ctx,
        lineChart(
          host,
          [
            { name: 'frequency f / f_Kerr − 1', x: ps.t0_ms, y: ps.f_Hz.map((f) => f / kf - 1), color: 1, markers: true },
            { name: 'damping time τ / τ_Kerr − 1', x: ps.t0_ms, y: ps.tau_ms.map((t) => t / kt - 1), color: 2, markers: true, dash: 'dash' },
          ],
          {
            xTitle: 'fit start after the peak (ms)',
            yTitle: 'free mode / Kerr − 1',
            height: 250,
            refLines: [{ axis: 'y', value: 0, label: 'Kerr prediction' }],
            format: (v) => fmtNum(v, 3),
            ariaLabel: 'Fitted frequency and damping time relative to Kerr against start time',
            caption: fit.significant
              ? 'Like [QG] Fig. 8 (right) for GW250114: the free mode against the Kerr prediction from the inspiral, as the start moves later. The quoted values are medians over start times; the spread is the systematic error.'
              : 'Without a significant ringdown the fitted mode wanders with the start time: this is noise, not a measurement of the remnant.',
          },
        ),
      );
    }
    // the paper's ringdown table (quoted) for this event, when it has one
    const pr = gw.paper_ringdown?.rows?.find((r) => r.event === paperName(e.name));
    if (pr) {
      const pair = (v: unknown, d = 1): string => (Array.isArray(v) && v.length === 2 ? `${Number(v[0]).toFixed(d)} ± ${Number(v[1]).toFixed(d)}` : '—');
      rdP.body.append(
        h('div', { class: 'm2-section-label', style: { marginTop: '12px' } }, `The paper’s ringdown of ${pr.event} · [QG] Sec. 19.2 (quoted)`),
        table(
          [
            { key: 'q', label: 'Quantity' },
            { key: 'paper', label: 'Paper', align: 'num' },
            { key: 'ours', label: 'Q-CHAS', align: 'num' },
          ],
          [
            { q: 'free mode f (Hz)', paper: pair(pr.free_f, 0), ours: pm(fit.f_Hz, fit.f_err_Hz, 0) },
            { q: 'free mode τ (ms)', paper: pair(pr.free_tau_ms, 2), ours: pm(fit.tau_ms, fit.tau_err_ms, 2) },
            { q: 'Kerr f₂₂₀ from IMR (Hz)', paper: pair(pr.kerr_f220, 0), ours: fmtNum(kf, 4) },
            { q: 'Kerr τ₂₂₀ from IMR (ms)', paper: pair(pr.kerr_tau_ms, 2), ours: fmtNum(kt, 3) },
            { q: 'IMR remnant M_f (M☉, det.)', paper: pair(pr.imr_mf, 1), ours: fmtNum(e.mf, 4) },
            { q: 'IMR remnant a_f', paper: pair(pr.imr_af, 2), ours: fmtNum(e.af, 3) },
          ],
        ),
        h(
          'div',
          { class: 'chart-caption' },
          `Different pipelines: the paper fits its own IMRPhenomD remnant, starts the free mode +6 to +20 M_f after the peak and uses a measured noise covariance; Q-CHAS uses a remnant fit to the tutorial’s ${e.approximant ?? 'template'} parameters and starts 1–6 ms after the peak. The free-mode values agree within errors.`,
        ),
      );
    }
  }

  function drawTelemetry(): void {
    const e = ev();
    clear(telP.body);
    const frame = e.thermo_frame ?? e.mass_frame ?? 'detector';
    telP.body.append(
      tilesN(
        3,
        statTile({ label: 'Final mass M_f', value: fmtNum(e.mf, 4), unit: 'M☉', badge: badge('computed', 'remnant fit'), sub: `${e.mass_frame ?? 'detector'} frame; from m₁ = ${fmtNum(e.m1, 4)}, m₂ = ${fmtNum(e.m2, 4)} (tutorial)` }),
        statTile({ label: 'Final spin a_f', value: fmtNum(e.af, 3), badge: badge('computed', 'remnant fit'), sub: `aligned-spin fit from a₁ = ${fmtNum(e.a1 ?? NaN, 3)}, a₂ = ${fmtNum(e.a2 ?? NaN, 3)}` }),
        statTile({ label: 'Kerr f₂₂₀', value: fmtNum(e.f220_Hz, 4), unit: 'Hz', badge: badge('computed') }),
        statTile({ label: 'Kerr τ₂₂₀', value: fmtNum(e.tau220_ms, 3), unit: 'ms', badge: badge('computed') }),
        statTile({ label: 'Echo delay (l_P surface)', value: fmtNum(e.echo_delay_ms, 4), unit: 'ms', badge: badge('computed', '[QG] Sec. 19'), sub: 'if a reflecting surface sat one l_P outside the horizon' }),
        statTile({ label: 'Horizon entropy log₁₀(S_BH / nats)', value: fmtNum(e.S_BH_log10, 5), badge: badge('computed'), sub: e.S_BH_bits_log10 ? `log₁₀(S_BH / bits) = ${fmtNum(e.S_BH_bits_log10, 5)}` : null }),
        statTile({
          label: 'Echo-power bound log₁₀|R|² ≤ −S/ln 10',
          value: fmtSci(e.echo_bound_log10, 3),
          badge: badge('computed', '[QG] Sec. 18'),
          sub: `|R|² ≤ e^(−S) = ${fmtPow10(e.echo_bound_log10)}: no echo any instrument can reach`,
        }),
        statTile({ label: 'Hawking temperature T_H', value: fmtSci(e.T_H_K, 3), unit: 'K', badge: badge('computed') }),
        statTile({ label: 'Evaporation time t_evap', value: fmtPow10(e.t_evap_log10_yr), unit: 'yr', badge: badge('assumed', 'scalar emission'), sub: `Page time ${fmtPow10(e.page_time_log10_yr)} yr` }),
      ),
      h(
        'div',
        { class: 'callout m2-callout small' },
        h('strong', null, 'Assumptions. '),
        `Component masses and spins are the LIGO tutorial’s ${e.approximant ?? ''} template parameters (${e.mass_frame ?? 'detector'} frame); M_f and a_f follow from aligned-spin remnant fits (Barausse–Morozova–Rezzolla, Rezzolla et al.; ~1 % in mass, ~0.02 in spin). Entropy, temperature and lifetimes are evaluated in the ${frame} frame${e.redshift ? '' : ' (redshift not removed, so the source-frame mass is smaller by 1 + z)'}. ` +
          `Lifetimes: ${e.lifetime_model ?? 'massless scalar, Schwarzschild'}, for an isolated hole; with T_H far below the 2.7 K microwave background the real remnant absorbs more than it emits today. ` +
          'S_BH = A/4G with the Kerr area; the echo delay assumes a reflecting surface one Planck length outside the horizon, the structure the framework says is absent.',
      ),
    );
    const chk = gw.echo_bound_check;
    if (chk)
      telP.body.append(
        kv([
          [
            `${chk.m_sun} M☉ check`,
            [
              h('span', { class: 'm2-inline-num' }, `log₁₀|R|² ours ${fmtSci(chk.log10_R2_ours ?? NaN, 4)}, paper ${fmtSci(chk.log10_R2, 3)} `),
              badge('quoted', chk.source ?? '[QG] Sec. 18'),
            ],
          ],
        ]),
      );
  }

  function renderEvent(): void {
    const e = ev();
    const fit = e.ringdown_fit;
    clear(selStatus);
    append(
      selStatus,
      h('span', { class: 'small text-2' }, `${e.name}: m₁ = ${fmtNum(e.m1, 4)}, m₂ = ${fmtNum(e.m2, 4)} M☉ (${e.mass_frame ?? 'detector'} frame) → M_f = ${fmtNum(e.mf, 4)} M☉, a_f = ${fmtNum(e.af, 3)}`),
      fit ? (fit.significant ? statusPill('good', 'ringdown significant') : statusPill('warning', 'ringdown not significant')) : null,
      e.name !== paperName(e.name) ? h('span', { class: 'small muted' }, `called ${paperName(e.name)} in [QG] Sec. 19`) : null,
    );
    drawStrain();
    drawRingdown();
    drawTelemetry();
  }
  renderEvent();

  // --- echo table + R90 chart --------------------------------------------------------------------------------------
  const echoP = panel({
    title: 'No echoes in nine mergers · [QG] Sec. 19.1',
    subtitle: 'The paper’s echo search in public LIGO–Virgo strain (quoted), with our recomputed f₂₂₀ and echo delay from its remnant parameters',
    badges: [badge('quoted', '[QG] Sec. 19.1'), badge('computed', 'f₂₂₀, Δt')],
    span: 8,
    flush: true,
  });
  echoP.body.classList.add('m2-dense');
  const ours = new Set(events.map((e) => paperName(e.name)));
  const echoRows = gw.paper_events;
  const num = (v: number | undefined, d = 0): string => (v === undefined ? '—' : v.toFixed(d));
  const both = (paper: string, ours: string): Child[] => [h('span', { class: 'm2-inline-num' }, paper), h('div', { class: 'small muted' }, `ours ${ours}`)];
  const echoTable = table<GwPaperEvent>(
    [
      { key: 'event', label: 'Event', format: (r) => [r.event, ours.has(r.event) ? h('span', { class: 'small muted', title: 'LIGO strain of this event is shown above' }, ' ●') : null] },
      { key: 'snr', label: 'SNR', align: 'num', format: (r) => num(r.snr, 1) },
      { key: 'mf', label: 'M_f (M☉)', align: 'num', format: (r) => num(r.mf, 1) },
      { key: 'af', label: 'a_f', align: 'num', format: (r) => num(r.af, 2) },
      { key: 'f220', label: 'f₂₂₀ (Hz)', align: 'num', format: (r) => both(num(r.f220), num(r.f220_ours, 1)) },
      { key: 'dt_ms', label: 'Δt (ms)', align: 'num', format: (r) => both(num(r.dt_ms), num(r.dt_ms_ours, 1)) },
      { key: 'gamma', label: 'Γ', align: 'num', format: (r) => num(r.gamma, 2), title: 'greybody factor Γ(ω₂₂₀) of the (2,2) mode' },
      { key: 'echo_snr', label: 'echo SNR', align: 'num', format: (r) => num(r.echo_snr, 2) },
      { key: 'p', label: 'p', align: 'num', format: (r) => num(r.p, 2) },
      { key: 'R90', label: 'R (90 % UL)', align: 'num', format: (r) => r.R90_label ?? (r.R90 === null ? '—' : r.R90.toFixed(2)) },
    ],
    echoRows,
  );
  const sum = gw.paper_echo_summary;
  const vb = validationBlock(val, ['Kerr ringdown and echoes ([QG] Sec. 19)'], { summaryOnly: true });
  append(
    echoP.body,
    h('div', { style: { padding: '12px 14px 0' } }, vb),
    echoTable,
    h(
      'div',
      { class: 'chart-caption', style: { padding: '0 14px 12px' } },
      `● strain of this event shown above (LVT151012 = GW151012). M_f in M☉, detector frame. R: reflectivity of a surface one Planck length outside the horizon excluded above this value at 90 % (CLs); "> 1": the event alone does not constrain R. p: fraction of off-source trials with an equal or louder echo statistic. ` +
        `Combined: Σ SNR² = ${sum?.sum_snr2 ?? '—'}, p = ${sum?.p_combined ?? '—'} (Fisher ${sum?.p_fisher ?? '—'}); R < ${gw.combined_R90}; the echo amplitude of Abedi et al. is excluded above A = ${sum?.A90_gamma09 ?? '—'} (γ = 0.9). ` +
        `Our f₂₂₀ and Δt are recomputed from the paper’s own remnant parameters (${sum?.masses ?? 'detector frame'}); GW150914 Δt ours ${sum?.dt_gw150914_ms?.ours ?? '—'} ms vs Abedi ${sum?.dt_gw150914_ms?.abedi ?? '—'} ms.`,
    ),
  );
  const r90P = panel({
    title: 'Reflectivity limits',
    subtitle: '90 % upper limits on R, after [QG] Fig. 8 (left)',
    badges: [badge('quoted', '[QG] Sec. 19.1')],
    span: 4,
    stretch: true,
  });
  const sorted = [...echoRows].sort((a, b) => a.snr - b.snr);
  const lim = sorted.map((r, i) => ({ r, i })).filter((x) => x.r.R90 !== null);
  const unc = sorted.map((r, i) => ({ r, i })).filter((x) => x.r.R90 === null);
  keep(
    ctx,
    scatterChart(
      r90P.body,
      [
        { name: '90 % upper limit', x: lim.map((x) => x.r.R90 as number), y: lim.map((x) => x.i), text: lim.map((x) => `${x.r.event} (SNR ${x.r.snr})`), color: 1, size: 10 },
        { name: 'unconstrained (> 1)', x: unc.map(() => 1.05), y: unc.map((x) => x.i), text: unc.map((x) => `${x.r.event} (SNR ${x.r.snr}): > 1`), color: 'context', symbol: 'triangle-right', size: 11 },
      ],
      {
        xTitle: 'reflectivity R excluded above (90 %, CLs)',
        yTitle: '',
        xRange: [0, 1.15],
        height: 380,
        refLines: [{ axis: 'x', value: gw.combined_R90, label: `all 9: R < ${gw.combined_R90}` }],
        layout: {
          yaxis: { tickvals: sorted.map((_, i) => i), ticktext: sorted.map((r) => `${r.event} (${Math.round(r.snr)})`), range: [-0.6, sorted.length - 0.4] },
          margin: { l: 120, r: 18, t: 36, b: 46, pad: 0 },
        } as Partial<Layout>,
        format: (v) => fmtNum(v, 3),
        ariaLabel: 'Echo reflectivity upper limits per event',
        caption: 'Events ordered by network SNR (in brackets). Triangles at 1.05 stand for "> 1": no constraint from that event alone. Quoted from the paper; not an analysis run here.',
      },
    ),
  );
  el.append(grid(echoP.root, r90P.root));

  // --- Kerr spin-down -------------------------------------------------------------------------------------------
  const kerrP = panel({
    title: 'A real spinning horizon · Kerr spin-down',
    subtitle: 'd ln a* / d ln M from exact Teukolsky greybody factors: > 0 means the hole spins down as it evaporates',
    badges: [badge('quoted', '[QG] Sec. 17, Theorem K′')],
    span: 6,
    stretch: true,
  });
  const kerrSide = panel({ title: 'Spin-down table', subtitle: 'Quoted per field, with our consistency check g/f − 2', span: 6, stretch: true, flush: true });
  el.append(grid(kerrP.root, kerrSide.root));
  if (!hor.ok) {
    kerrP.body.append(loadFailure(hor));
    kerrSide.body.append(loadFailure(hor));
  } else {
    const k = hor.data.kerr_spin_down;
    const fields = [...new Set(k.rows.map((r) => r.field))];
    const yz = { zoom: false };
    const host = h('div');
    const zSeg = segmented({
      label: 'y range',
      options: [
        { value: 'full', label: 'All fields' },
        { value: 'zoom', label: 'Zoom: scalar fixed point' },
      ],
      value: 'full',
      onChange: (v) => {
        yz.zoom = v === 'zoom';
        drawKerr();
      },
    });
    let kChart: ChartHandle | null = null;
    const drawKerr = (): void => {
      kChart?.destroy();
      clear(host);
      kChart = keep(
        ctx,
        scatterChart(
          host,
          fields.map((f, i) => {
            const rows = k.rows.filter((r) => r.field === f && r.dlna_dlnM !== null);
            return { name: f, x: rows.map((r) => r.a_star), y: rows.map((r) => r.dlna_dlnM), color: i + 1, mode: 'lines+markers' as const, size: 8 };
          }),
          {
            xTitle: 'spin a*',
            yTitle: 'd ln a* / d ln M  (> 0: spin-down)',
            height: 320,
            yRange: yz.zoom ? [-1, 1.5] : undefined,
            xRange: [0.1, 1.02],
            refLines: [
              { axis: 'y', value: 0 },
              { axis: 'x', value: k.fixed_point_scalar, label: `scalar fixed point a* = ${k.fixed_point_scalar}` },
            ],
            format: (v) => fmtNum(v, 4),
            ariaLabel: 'Spin-down rate against spin for scalar, photon and graviton emission',
          },
        ),
      );
    };
    kerrP.body.append(h('div', { class: 'm2-pill-row' }, zSeg.root), host);
    drawKerr();
    const sd = k.spin_down_check;
    const facts = (k.facts ?? {}) as Record<string, Record<string, unknown> | undefined>;
    const gains = facts.superradiant_gain_a0999 as Record<string, number> | undefined;
    const born = facts['born_at_0.99'] as Record<string, number | string> | undefined;
    const tiles: HTMLElement[] = [
      statTile({ label: 'Scalar fixed point', value: k.fixed_point_scalar, unit: 'a*', badge: badge('quoted', '[QG] Sec. 17'), sub: 'Chambers–Hiscock–Taylor, from first principles' }),
    ];
    if (sd)
      tiles.push(
        statTile({
          label: `Mass left when a* falls ${sd.a_start} → ${sd.a_end}`,
          value: fmtNum(100 * sd.mass_fraction, 3),
          unit: '%',
          badge: badge('computed'),
          sub: [badge('quoted', '[QG] Sec. 17'), ` paper ${fmtNum(100 * sd.paper, 3)} %`],
        }),
      );
    if (born && typeof born['a_0.1_mass_fraction'] === 'number')
      tiles.push(statTile({ label: 'Mass left at a* = 0.1', value: fmtNum(100 * born['a_0.1_mass_fraction'], 3), unit: '%', badge: badge('quoted', '[QG] Sec. 17'), sub: `born at 0.99, ${String(born.fields ?? 'photons + gravitons')}` }));
    kerrP.body.append(h('div', { style: { marginTop: '12px' } }, tilesN(3, ...tiles)));
    const notes: Child[] = [
      'With photons and gravitons d ln a*/d ln M > 0 at every spin: the hole spins down monotonically and loses its spin long before its mass. ',
      gains ? `Maximal superradiant amplification at a* = 0.999: ${fmtNum(100 * (gains.scalar_l1m1 ?? NaN), 3)} % (scalar), ${fmtNum(100 * (gains.photon_l1m1 ?? NaN), 3)} % (photon), ${fmtNum(100 * (gains.graviton_l2m2 ?? NaN), 3)} % (graviton). ` : '',
      '"area grows": at high spin the horizon area grows while the hole radiates, as the area theorem requires.',
    ];
    kerrP.body.append(h('div', { class: 'chart-caption' }, ...notes));
    kerrSide.body.classList.add('m2-dense');
    kerrSide.body.append(
      table<KerrSpinDownRow>(
        [
          { key: 'field', label: 'Field' },
          { key: 'a_star', label: 'a*', align: 'num', format: (r) => fmtNum(r.a_star, 3) },
          { key: 'f', label: 'f = −M²Ṁ', align: 'num', format: (r) => fmtSci(r.f, 3) },
          { key: 'g', label: 'g = −(M/a*)J̇', align: 'num', format: (r) => (r.g === null ? '—' : fmtSci(r.g, 3)) },
          { key: 'd', label: 'd ln a*/d ln M', align: 'num', format: (r) => (r.dlna_dlnM === null ? '—' : (r.dlna_dlnM > 0 ? '+' : '') + r.dlna_dlnM.toFixed(3)) },
          { key: 'chk', label: 'g/f − 2', align: 'num', format: (r) => (r.dlna_dlnM_from_f_g === null ? '—' : fmtNum(r.dlna_dlnM_from_f_g, 3)), title: 'computed: g/f − 2 from the quoted f and g, a consistency check of the quoted column' },
          { key: 'gsl', label: 'GSL', title: 'generalized second law: Ṡ_rad / (−Ṡ_BH)', align: 'num', format: (r) => (r.area_grows ? 'area ↑' : r.gsl_ratio === null ? '—' : fmtNum(r.gsl_ratio, 3)) },
        ],
        k.rows,
      ),
      h(
        'div',
        { class: 'chart-caption', style: { padding: '0 14px 12px' } },
        `Quoted from ${k.source ?? '[QG] Sec. 17'} (Teukolsky greybody factors, ℓ ≤ 6). g/f − 2 is our check of the quoted d ln a*/d ln M from the quoted f and g. GSL: Ṡ_rad/(−Ṡ_BH); "area ↑": the area grows while the hole radiates. Sign-law violations: ${k.rows.reduce((a, r) => a + (r.sign_violations ?? 0), 0)} in all rows.`,
      ),
    );
  }

  // --- GWTC catalogue ------------------------------------------------------------------------------------------
  const catP = panel({ title: 'GWOSC event catalogue', subtitle: 'Every confident GWTC merger with its remnant horizon, when the build host can reach GWOSC', badges: [badge('measured', 'GWTC')] });
  if (!gw.gwtc || !gw.gwtc.length) {
    catP.body.append(
      h(
        'div',
        { class: 'callout' },
        statusPill('warning', 'not in this build'),
        ' ',
        'The GWOSC event API was not reachable from the build host, so gw.json carries gwtc = null. The catalogue (source-frame masses, remnant spin estimates, f₂₂₀, echo delays and horizon entropies for every binary black hole, neutron-star candidates flagged) appears here when ',
        h('code', null, 'qchas build'),
        ' runs where GWOSC is reachable. No catalogue numbers are shown in its place.',
      ),
    );
  } else {
    const cols: { key: string; label: string; d?: number }[] = [
      { key: 'name', label: 'Event' },
      { key: 'catalog', label: 'Catalogue' },
      { key: 'm1_source', label: 'm₁ (M☉, src)', d: 1 },
      { key: 'm2_source', label: 'm₂ (M☉, src)', d: 1 },
      { key: 'mf_source', label: 'M_f (M☉, src)', d: 1 },
      { key: 'chi_eff', label: 'χ_eff', d: 2 },
      { key: 'redshift', label: 'z', d: 2 },
      { key: 'af_estimate', label: 'a_f est.', d: 2 },
      { key: 'f220_Hz', label: 'f₂₂₀ (Hz)', d: 0 },
      { key: 'echo_delay_ms', label: 'Δt (ms)', d: 0 },
      { key: 'S_BH_log10', label: 'log₁₀ S_BH', d: 2 },
    ];
    catP.body.append(
      table<Record<string, unknown>>(
        [
          ...cols.map((c) => ({
            key: c.key,
            label: c.label,
            align: c.d === undefined ? ('text' as const) : ('num' as const),
            format: (r: Record<string, unknown>) => {
              const v = r[c.key];
              return typeof v === 'number' ? v.toFixed(c.d ?? 0) : v === null || v === undefined ? '—' : String(v);
            },
          })),
          { key: 'ns', label: '', format: (r) => (r.possible_neutron_star ? statusPill('warning', 'possible NS') : null) },
        ],
        gw.gwtc.slice(0, 300),
      ),
    );
  }
  el.append(grid(catP.root));
}
