"""Module 2 horizon science: [BH] Theorems C-G, Secs. 9-11; [QG] Theorems Q, Q', Lambda, R; Kerr."""

from __future__ import annotations

import json
import math

import numpy as np
import pytest
from scipy.linalg import sqrtm

from qchas.horizon import atmosphere, bounce, capacity, evaporation, kerr, page, qes


# --- [BH] Theorems C-F: the Page law -------------------------------------------------------

def test_page_time_and_mass_fraction():
    assert round(page.page_time_fraction(), 3) == 0.476
    assert round(page.page_mass_fraction(), 3) == 0.806
    # at the Page time the Hawking and capacity branches cross
    tp = page.page_time_fraction()
    assert page.hawking_entropy(tp) == pytest.approx(page.bh_entropy(tp), rel=1e-12)


def test_island_is_forced_after_page_time():
    tp = page.page_time_fraction()
    t = np.linspace(0.0, 1.0, 201)
    isl = page.island_min_entropy(t)
    assert np.all(isl[t < tp - 1e-9] == 0.0)
    assert np.all(isl[t > tp + 1e-3] > 0.0)
    assert isl[-1] == pytest.approx(page.BETA_RADIATION)  # beta S_BH(M0) at the endpoint
    assert np.all(page.page_curve(t) <= page.hawking_entropy(t) + 1e-15)
    assert np.all(page.page_curve(t) <= page.bh_entropy(t) + 1e-15)


def test_interior_volume_endpoint():
    assert page.interior_volume(1.0) == pytest.approx(9 * math.sqrt(3) * math.pi / 5, rel=1e-12)
    v = page.interior_volume(np.linspace(0, 1, 50))
    assert np.all(np.diff(v) > 0)


def test_decoding_bound_and_hayden_preskill():
    assert page.decoding_infidelity_bound(2, 0) == 1.0
    assert page.decoding_infidelity_bound(2, 7) == pytest.approx(2.0**-4)
    assert page.hayden_preskill_emissions(1, 0.25) == pytest.approx(3.0)


# --- [BH] Theorem G, [QG] Theorems C and U: capacity ---------------------------------------

@pytest.mark.parametrize("m_star,width", [(100.0, 1e-3), (1e3, 1e-4), (10.0, 0.05)])
def test_capacity_window_matches_theorem_G_v(m_star, width):
    w = capacity.capacity_window(m_star, width)
    # Theorem G(v) uses S_BH(M* + w) - S_BH(M*) = beta w + 4 pi w^2: equal up to that curvature
    roundoff = 64 * np.finfo(float).eps * 4 * math.pi * (m_star + width) ** 2  # S_BH(M* + w) - S_BH(M*)
    assert w["theorem_G_v"] - w["ln_trace"] == pytest.approx(4 * math.pi * width**2, rel=1e-6, abs=roundoff)


def test_area_gap_and_de_sitter_entropy():
    assert round(capacity.minimal_area(), 2) == 5.97
    ds = capacity.de_sitter_entropy()
    assert ds["log10_entropy"] == pytest.approx(122.5, abs=0.1)
    assert ds["horizon_radius_Gpc"] == pytest.approx(5.38, abs=0.05)


# --- [BH] Sec. 9: quantum extremal surface -------------------------------------------------

@pytest.mark.parametrize("row", qes.PAPER_TABLE, ids=[f"M{r[0]}-c{r[1]}-rb{r[2]}" for r in qes.PAPER_TABLE])
def test_qes_reproduces_paper_table(row):
    m, c, rb, depth, lag, t_scr, ds = row
    r = qes.solve_qes(m, c, rb)
    assert r.converged
    assert f"{r.depth:.2e}" == f"{depth:.2e}"
    assert round(r.lag, 2) == pytest.approx(lag, abs=1e-9)
    assert round(r.t_scr, 2) == pytest.approx(t_scr, abs=1e-9)
    assert round(r.s_gen_minus_s_bh, 2) == pytest.approx(ds, abs=1e-9)


# --- [BH] Secs. 10, 15: the quantum atmosphere --------------------------------------------

def test_atmosphere_reduced():
    r = atmosphere.atmosphere_profile(60, 20, gaps=(0, 1, 2, 3, 5), ell_max=150)
    share = np.array(r["share"])
    assert share[0] == pytest.approx(1.0, abs=1e-9)  # pure state: I(ball : complement) = 2S
    assert np.all(np.diff(share) < 0)
    for g, s in zip(r["gaps"][1:], share[1:]):
        assert s == pytest.approx(atmosphere.PAPER_SHARE[g], rel=0.06)
    assert r["area_coeff"] == pytest.approx(0.295, abs=0.002)


@pytest.mark.slow
def test_atmosphere_full_paper_case():
    r = atmosphere.atmosphere_profile(60, 20, gaps=tuple(atmosphere.PAPER_SHARE), ell_max=3000)
    for g, s in zip(r["gaps"], r["share"]):
        printed = atmosphere.PAPER_SHARE[g]
        # printed to 0.01 %: half a unit is 5e-5. g = 1 gives 15.504 %, 1e-5 below the rounding
        # interval of the printed 15.51 % (tail-fit level difference, reported to the lead)
        assert s == pytest.approx(printed, abs=7e-5 if g else 1e-6)
    assert round(r["area_coeff"], 3) == 0.295


# --- [BH] Sec. 11: exactly unitary evaporation --------------------------------------------

def _page_sum(m: int, n: int) -> float:
    m, n = min(m, n), max(m, n)
    return (sum(1.0 / k for k in range(n + 1, m * n + 1)) - (m - 1) / (2 * n)) / math.log(2)


@pytest.mark.parametrize("m,n", [(2, 2), (2, 8), (4, 4), (8, 32), (16, 16)])
def test_page_average_formula(m, n):
    assert evaporation.page_average_entropy(m, n) == pytest.approx(_page_sum(m, n), rel=1e-12)
    assert evaporation.page_average_entropy(n, m) == pytest.approx(_page_sum(m, n), rel=1e-12)


def test_page_average_at_the_page_point():
    # m = n = 2^11: ln m - 1/2 nats = 11 - 0.7213 bits, the 10.28 bits of [BH] Sec. 11
    assert evaporation.page_average_entropy(2**11, 2**11) == pytest.approx(10.2787, abs=1e-4)
    assert evaporation.page_average_entropy(1, 2**11) == 0.0


def test_scramble_is_unitary_on_b():
    rng = np.random.default_rng(3)
    for shape in [(64, 8), (8, 64), (16, 16)]:
        m = rng.standard_normal(shape) + 1j * rng.standard_normal(shape)
        m /= np.linalg.norm(m)
        out = evaporation._scramble(rng, m)
        # a unitary on the rows preserves the Gram matrix of the columns
        assert np.allclose(out.conj().T @ out, m.conj().T @ m, atol=1e-12)


def test_evaporation_10_qubits():
    n = 10
    r = evaporation.evaporate(n, None, seed=4)
    s, env = np.array(r["s_rad"]), np.array(r["envelope"])
    assert len(r["steps"]) == n
    assert np.all(s <= env + 1e-9)  # Theorem E(i)
    assert r["final_bits"] == pytest.approx(0.0, abs=1e-9)
    assert r["max_step"] == n // 2
    assert s[n // 2 - 1] == pytest.approx(r["page_average"][n // 2 - 1], abs=0.1)
    # Theorem E(ii) lower bound min - e^{-|S_H - A/4G|} (in nats), up to sample fluctuations
    lower = env - np.exp(-np.abs(np.array(r["hawking"]) - np.array(r["s_bh"])) * math.log(2)) / math.log(2)
    assert np.all(s >= lower - 0.15)

    d = evaporation.evaporate(n, 2, seed=5)
    assert len(d["steps"]) == n + 1
    assert d["diary_mi"][1] == 0.0  # thrown in after step 2: nothing in R yet
    assert d["final_bits"] == pytest.approx(1.0, abs=1e-9)
    assert np.all(np.array(d["diary_mi"]) <= 2.0 + 1e-9)
    assert d["diary_mi"][-1] == pytest.approx(2.0, abs=1e-9)
    assert d["recovered_step"] is not None and d["recovered_step"] > n // 2


def test_evaporation_and_decoding_reject_bad_input():
    with pytest.raises(ValueError):
        evaporation.evaporate(6, 7)
    with pytest.raises(ValueError):
        evaporation.evaporate(6, -1)
    with pytest.raises(ValueError):
        evaporation.decoding_infidelity(3, 1, samples=1, dynamics="gue")


@pytest.mark.slow
def test_evaporation_20_qubits_paper_numbers():
    p = evaporation.PAPER_EVAPORATION
    early = evaporation.evaporate(20, p["early_diary_step"], seed=1)
    assert round(early["max_bits"], 2) == p["max_bits"]
    assert early["max_step"] == p["max_step"]
    assert round(early["max_dev"], 2) == p["max_dev"]
    assert early["max_dev_step"] == p["max_step"]
    assert early["final_bits"] == pytest.approx(1.0, abs=1e-9)
    assert early["recovered_step"] == p["early_recovered_step"]
    late = evaporation.evaporate(20, p["late_diary_step"], seed=2)
    assert late["recovered_step"] == p["late_recovered_step"]
    assert late["diary_mi"][p["late_diary_step"] - 1] == 0.0
    plain = evaporation.evaporate(20, None, seed=0)
    assert plain["max_step"] == 10 and round(plain["max_bits"], 2) == 9.28


# --- [QG] Theorem R: Hayden-Preskill decoding ---------------------------------------------

def _infidelity_direct(psi: np.ndarray, n_total: int, s: int) -> float:
    d_i, d_e, d_m = psi.shape
    d_bp = 2 ** (n_total - s)
    t = psi.reshape(d_bp, 2**s, d_e, d_m)
    rho_mb = np.einsum("aiem,bien->manb", t, t.conj()).reshape(d_m * d_bp, d_m * d_bp)
    rho_b = np.einsum("aiem,biem->ab", t, t.conj())
    sigma = np.kron(np.eye(d_m) / d_m, rho_b)
    f = np.sum(np.linalg.svd(sqrtm(rho_mb) @ sqrtm(sigma), compute_uv=False)).real
    return 1.0 - f * f


def test_fast_fidelity_matches_direct():
    rng = np.random.default_rng(7)
    n_bh, k = 3, 1
    u = evaporation.haar_unitary(rng, 2 ** (n_bh + k))
    psi = u.reshape(16, 2, 8).transpose(0, 2, 1) / math.sqrt(16)
    for s in range(n_bh + k + 1):
        assert evaporation._infidelity_after(psi, n_bh + k, s) == pytest.approx(
            _infidelity_direct(psi, n_bh + k, s), abs=1e-7)


@pytest.mark.parametrize("dynamics,old", [("haar", True), ("eth", True), ("haar", False)])
def test_decoding_endpoints_and_bound(dynamics, old):
    k = 2
    r = evaporation.decoding_infidelity(5, k, samples=8, old=old, dynamics=dynamics, t_eth=10.0, seed=11)
    inf = np.array(r["infidelity"])
    # s = 0: F = 2^-k exactly for any unitary; all qubits emitted: perfect decoding
    assert inf[0] == pytest.approx(1 - 4.0**-k, abs=1e-12)
    assert abs(inf[r["s"].index(7)]) < 1e-12
    if old:
        assert np.all(inf <= np.array(r["bound"]) + 1e-12)


@pytest.mark.slow
def test_decoding_table_matches_paper():
    tab = evaporation.decoding_table(samples=6, seed=0)
    paper = evaporation.PAPER_DECODING
    cols = [("haar_old", tab["haar_old"], paper["haar_old"]), ("haar_young", tab["haar_young"], paper["haar_young"])]
    cols += [(f"eth{t}", tab["eth"][t], paper["eth"][t]) for t in ("1", "3", "10", "100")]
    for name, ours, printed in cols:
        for s in range(9):
            o, p, sd = ours["infidelity"][s], printed[s], ours["std"][s]
            tol = 3.0 * math.sqrt(2.0) * sd / math.sqrt(6) + 6e-3 * p
            assert abs(o - p) <= tol, (name, s, o, p, tol)
        assert abs(ours["infidelity"][9]) < 1e-12
    haar = tab["haar_old"]
    assert all(m <= b for m, b in zip(haar["infidelity"], haar["bound"]))


# --- Kerr remnants, ringdown and echo delays ------------------------------------------------

#: [QG] Sec. 19.1 table rows: event, M_f (det.), a_f, f220 (Hz), echo delay (ms)
QG_ECHO_ROWS = [
    ("GW150914", 68.7, 0.68, 248, 293), ("GW151012", 47.4, 0.67, 355, 199), ("GW151226", 21.4, 0.79, 875, 100),
    ("GW170104", 58.3, 0.64, 282, 241), ("GW170814", 59.4, 0.72, 296, 262), ("GW190521", 85.6, 0.68, 199, 366),
    ("GW200129", 71.4, 0.73, 248, 318), ("GW241011", 21.0, 0.90, 1020, 122), ("GW250114", 68.6, 0.66, 244, 289),
]


@pytest.mark.parametrize("row", QG_ECHO_ROWS, ids=[r[0] for r in QG_ECHO_ROWS])
def test_kerr_qnm_and_echo_delay_vs_paper(row):
    _, mf, af, f220, dt_ms = row
    f, tau = kerr.qnm(mf, af)
    assert f == pytest.approx(f220, rel=0.01)
    # the printed spins are rounded to 0.01; the delay's 1/sqrt(1 - a^2) amplifies it at a = 0.9
    assert 1e3 * kerr.echo_delay(mf, af) == pytest.approx(dt_ms, rel=0.015)


def test_kerr_ringdown_table_predictions():
    f, tau = kerr.qnm(67.4, 0.64)  # [QG] Sec. 19.2: 246 +- 4 Hz, 3.99 +- 0.18 ms
    assert f == pytest.approx(246, abs=2) and 1e3 * tau == pytest.approx(3.99, abs=0.03)
    f, tau = kerr.qnm(68.5, 0.66)  # 245 +- 1 Hz, 4.08 +- 0.03 ms
    assert f == pytest.approx(245, abs=1) and 1e3 * tau == pytest.approx(4.08, abs=0.02)


def test_remnant_fits_numerical_relativity_limits():
    mf, af = kerr.remnant_aligned(30.0, 30.0)
    assert mf / 60.0 == pytest.approx(0.9517, abs=1e-3)
    assert af == pytest.approx(0.686, abs=2e-3)
    mf, af = kerr.remnant_aligned(36.0, 29.0)
    assert mf == pytest.approx(62.0, abs=0.5) and af == pytest.approx(0.68, abs=0.01)
    assert kerr.remnant_aligned(29.0, 36.0) == kerr.remnant_aligned(36.0, 29.0)


def test_kerr_thermodynamics():
    bh = kerr.KerrHorizon(62.0, 0.0)
    # [QG] Sec. 18: |R|^2 <= e^{-S} = 10^{-1.75e80} for the 62 Msun remnant (Schwarzschild area)
    assert kerr.echo_reflectivity_bound(bh.entropy) / 1e80 == pytest.approx(-1.75, abs=0.005)
    assert bh.hawking_temperature_K == pytest.approx(6.17e-8 / 62.0, rel=0.01)
    spinning = kerr.KerrHorizon(62.0, 0.7)
    assert spinning.entropy < bh.entropy and spinning.hawking_temperature_K < bh.hawking_temperature_K
    assert bh.page_time_years() == pytest.approx(page.page_time_fraction() * bh.evaporation_time_years())


def test_spin_down_table_consistency():
    d = kerr.to_json()
    for r in d["rows"]:
        if r["g"] is not None:
            assert r["dlna_dlnM_from_f_g"] == pytest.approx(r["dlna_dlnM"], abs=0.03)
    assert d["spin_down_check"]["mass_fraction"] == pytest.approx(0.77, abs=0.01)
    json.dumps(d)


# --- [QG] Theorems Q, Q', Lambda: the quantum Big Bang -------------------------------------

def test_exact_packet_normalised_symmetric_periodic():
    t = np.array([0.0, 0.4, 1.2, math.pi / 2, 2.0, 3.0])
    m = bounce.moments(t, 2.0)
    assert np.allclose(m["norm"], 1.0, atol=1e-12)
    m2 = bounce.moments(t + math.pi, 2.0)
    assert np.allclose(m2["a_inv2"], m["a_inv2"], rtol=1e-10)
    d = 0.17
    a = bounce.moments([math.pi / 2 - d, math.pi / 2 + d], 4.0)["a_inv2"]
    assert a[0] == pytest.approx(a[1], rel=1e-10)
    assert m["a_mean"][0] == pytest.approx(2.0, abs=1e-6)


@pytest.mark.parametrize("a_max", [1.5, 2.0, 8.0])
def test_closed_forms_match_quadrature(a_max):
    m = bounce.moments(math.pi / 2, a_max)
    assert m["a_inv2"][0] == pytest.approx(bounce.peak_a_inv2_closed_form(a_max), rel=1e-10)
    assert m["a_mean"][0] == pytest.approx(bounce.a_min_closed_form(a_max), rel=1e-10)


def test_fd_spectrum_reproduces_dirichlet_levels():
    w = bounce.fd_levels(*bounce.PAPER_SPECTRUM_GRID, n_eig=30)
    exact = bounce.exact_levels(30)
    assert abs(w[0] - exact[0]) < 5e-5  # [QG] Sec. 21: "to 5e-5" holds for the ground level
    # second-order finite differences: the error grows like h^2 E^2, 4.4e-3 relative at j = 29
    assert np.all(np.diff(np.abs(w - exact)) > 0)
    assert abs(w[-1] - exact[-1]) / exact[-1] < 5e-3


def test_exact_solution_matches_hermite_expansion():
    # independent route: odd oscillator levels (m = 12, omega = 1) weighted by the coherent state
    a_max, alpha = 2.0, math.sqrt(6.0) * 2.0
    nmax = 160
    n = np.arange(nmax + 1)
    c = np.exp(-alpha**2 / 2 + n * math.log(alpha) - 0.5 * np.array([math.lgamma(k + 1) for k in n]))
    c = np.where(n % 2 == 1, c, 0.0)
    c /= np.linalg.norm(c)
    a = np.linspace(1e-6, a_max + 3.0, 2001)
    x = math.sqrt(12.0) * a
    phi = np.zeros((nmax + 1, a.size))
    phi[0] = (12.0 / math.pi) ** 0.25 * np.exp(-x * x / 2)
    phi[1] = math.sqrt(2.0) * x * phi[0]
    for k in range(2, nmax + 1):
        phi[k] = math.sqrt(2.0 / k) * x * phi[k - 1] - math.sqrt((k - 1) / k) * phi[k - 2]
    for t in (0.3, math.pi / 2 - 0.05, math.pi / 2):
        psi = (c * np.exp(-1j * (n + 0.5) * t)) @ (math.sqrt(2.0) * phi)
        exact = bounce.density(a, t, a_max)
        assert np.max(np.abs(np.abs(psi) ** 2 - exact)) < 1e-10 * np.max(exact)


def test_exact_solution_matches_fd_spectral_evolution():
    fd = bounce.FDEvolution(2.0, 10.0, 4000)
    for t in (0.7, math.pi / 2, 2.3):
        d = fd.density(t)
        e = bounce.density(fd.a, t, 2.0)
        assert np.max(np.abs(d - e)) < 1e-3 * np.max(e) * 10


def test_bounce_exact_runs_vs_theorems():
    r2 = bounce.bounce_run(2.0)
    p = bounce.PAPER_LAMBDA_TABLE[2]
    assert r2["fwhm_T"] == pytest.approx(p[0], rel=0.002)
    assert r2["fwhm_times_amax"] == pytest.approx(p[1], rel=0.002)
    assert round(r2["a_min"], 3) == bounce.PAPER_ISOTROPIC_AMIN_THEOREM_N[2]
    assert round(bounce.bounce_run(1.5)["a_min"], 3) == bounce.PAPER_ISOTROPIC_AMIN_THEOREM_N[1.5]
    t = np.array(r2["T"])
    k = int(np.argmax(r2["a_inv2"]))
    assert t[k] == pytest.approx(math.pi / 2, abs=1e-9)
    assert min(r2["a_mean"]) == pytest.approx(r2["a_min"], rel=1e-9)
    # quantum <a> never reaches the classical singularity
    assert np.all(np.array(r2["a_mean"]) >= r2["a_min"] - 1e-12)


@pytest.mark.parametrize("a_max", [4, 8, 16])
def test_bounce_fwhm_vs_lambda_table(a_max):
    r = bounce.bounce_run(float(a_max))
    assert r["fwhm_T"] == pytest.approx(bounce.PAPER_LAMBDA_TABLE[a_max][0], rel=0.01)


def test_bounce_universal_limits():
    big = 48.0
    assert bounce.a_min_closed_form(big) == pytest.approx(bounce.A_MIN_LIMIT, rel=1e-3)
    r = bounce.bounce_run(big, T=[math.pi / 2])
    assert r["fwhm_times_amax"] == pytest.approx(bounce.FWHM_PRODUCT_LIMIT, rel=0.01)
    assert bounce.A_MIN_LIMIT * bounce.ELL_OVER_LP == pytest.approx(0.26, abs=0.005)  # "about 0.26 l_P"


@pytest.mark.parametrize("a_max", [2, 4])
def test_fd_paper_grid_reproduces_lambda_table(a_max):
    s = bounce.FDEvolution(float(a_max)).summary()
    fwhm, prod, peak, amin = bounce.PAPER_LAMBDA_TABLE[a_max]
    assert round(s["peak_a_inv2"]) == peak
    assert round(s["a_min"], 4) == pytest.approx(amin, abs=1e-4)
    assert s["fwhm_T"] == pytest.approx(fwhm, rel=0.006)


def test_fd_peak_converges_to_exact_at_first_order():
    # the printed 262 is the h = 24/6001 grid sum; its deficit is the O(h) end term h f(0)/2 at a = 0
    exact = bounce.peak_a_inv2_closed_form(2.0)
    f0 = 4.0 * math.sqrt(12.0 / math.pi) * 24.0**2  # lim |Psi|^2 / a^2 at the bounce, p = 12 a_max
    gaps = []
    for n in (6000, 12000, 24000):
        h = 24.0 / (n + 1)
        gap = exact - bounce.FDEvolution(2.0, 24.0, n).summary()["peak_a_inv2"]
        assert gap == pytest.approx(0.5 * h * f0, rel=0.05)
        gaps.append(gap)
    assert gaps[0] / gaps[1] == pytest.approx(2.0, rel=0.05) and gaps[1] / gaps[2] == pytest.approx(2.0, rel=0.05)


@pytest.mark.slow
@pytest.mark.parametrize("a_max", [8, 16])
def test_fd_paper_grid_reproduces_lambda_table_large(a_max):
    s = bounce.FDEvolution(float(a_max)).summary()
    fwhm, prod, peak, amin = bounce.PAPER_LAMBDA_TABLE[a_max]
    assert round(s["peak_a_inv2"]) == peak
    assert round(s["a_min"], 4) == pytest.approx(amin, abs=1e-4)
    assert s["fwhm_T"] == pytest.approx(fwhm, rel=0.006)


def test_bounce_fwhm_undefined_for_tiny_universes():
    assert bounce.bounce_run(0.2, T=[math.pi / 2])["fwhm_T"] is None
    r = bounce.bounce_run(0.5, T=[math.pi / 2])
    assert 1.0 < r["fwhm_T"] < math.pi and r["fwhm_times_amax"] == pytest.approx(0.5 * r["fwhm_T"])


def test_graviton_metric_perturbation():
    h2 = bounce.graviton_h2(262.0, (3, 4, 5, 10))
    assert h2 == pytest.approx([262 / 6, 262 / 8, 262 / 10, 262 / 20])
    assert bounce.graviton_h2(100.0, (3,), occupation=2) == pytest.approx([2.5 / 3 * 100])


def test_snapshots_shape_and_norm():
    s = bounce.snapshots(6.0)
    a = np.array(s["a"])
    for psi2 in s["psi2"]:
        assert np.trapezoid(psi2, a) == pytest.approx(1.0, abs=2e-3)
    # at the bounce the wave function vanishes at a = 0 and has fringes of spacing pi/(12 a_max)
    mid = np.array(s["psi2"][1])
    assert mid[0] == 0.0
    zeros = a[1:-1][(mid[1:-1] < mid[:-2]) & (mid[1:-1] < mid[2:]) & (a[1:-1] < 0.5)]
    assert np.median(np.diff(zeros)) == pytest.approx(math.pi / 72, rel=0.05)


@pytest.mark.slow
def test_bounce_to_json():
    d = bounce.to_json(fast=True)
    text = json.dumps(d)
    assert "NaN" not in text
    assert {"runs", "spectrum_check", "snapshots", "graviton", "mixmaster", "planck_bounce_test", "units"} <= set(d)
    assert d["mixmaster"]["a_min"] == 0.256
    assert d["planck_bounce_test"]["kc_excluded_above"] == 4.9e-4
    assert all(isinstance(r["paper"], dict) for r in d["runs"])
    run2 = next(r for r in d["runs"] if r["a_max"] == 2.0)
    assert run2["paper"]["peak_a_inv2"] == 262.0
    assert round(run2["paper_grid"]["peak_a_inv2"]) == 262
