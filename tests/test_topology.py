"""Tests for qchas.topology: manifolds, lattice sums, shape flow, circles, vacuum ([TOP])."""

from __future__ import annotations

import json
import math

import numpy as np
import pytest
from scipy.integrate import quad

from qchas.constants import FCC_PREREGISTERED_AXES_GAL, HESSIAN_CERTIFIED, Z2_CERTIFIED
from qchas.topology import circles as C
from qchas.topology import lattice as L
from qchas.topology import manifolds as M
from qchas.topology import shapeflow as F
from qchas.topology import vacuum as V

# --- manifolds ------------------------------------------------------------------------


def test_catalogue():
    assert [m.label for m in M.MANIFOLDS] == [f"E{i}" for i in range(1, 19)]
    assert sum(m.closed for m in M.MANIFOLDS) == 10
    assert [m.label for m in M.MANIFOLDS if m.translations_only] == ["E1", "E11", "E16", "E18"]
    assert [m.label for m in M.MANIFOLDS if not m.orientable] == ["E7", "E8", "E9", "E10", "E13",
                                                                  "E14", "E15", "E17"]


def test_principle_table_matches_paper():
    for key, row in M.PRINCIPLE_TABLE_PAPER.items():
        for lab in M.expand_label_range(key):
            st = M.BY_LABEL[lab].principles()
            assert tuple(st[p] for p in M.PRINCIPLE_KEYS) == row, lab


def test_main_theorem_selection():
    assert M.survivors(["P2", "P3"]) == ["E1"]
    assert M.survivors(["P2", "P4"]) == ["E1"]
    assert M.survivors(["P1", "P2", "P3", "P4"]) == ["E1"]
    assert M.survivors(["P3"]) == ["E1", "E11", "E16", "E18"]  # Main Theorem (i)
    assert M.survivors(["P2"]) == [f"E{i}" for i in range(1, 11)]
    assert M.minimal_selecting_sets() == [["P2", "P3"], ["P2", "P4"]]
    with pytest.raises(ValueError):
        M.survivors(["P5"])


def test_logical_independence():
    assert M.failed_principles("E18") == ["P2"]
    assert M.failed_principles("E2") == ["P3", "P4"]
    assert M.failed_principles("E7") == ["P1", "P3", "P4"]
    assert M.satisfies("E7", ["P2"]) and not M.satisfies("E7", ["P1"])


def test_failure_counts():
    """[TOP] Secs. 4/5.6 say every manifold except E1 fails >= 2 principles; its own table
    shows E11, E16, E18 fail only P2. We assert what the table implies."""
    fails = {m.label: M.failed_principles(m.label) for m in M.MANIFOLDS}
    assert fails["E1"] == []
    single = sorted((k for k, v in fails.items() if len(v) == 1), key=lambda s: int(s[1:]))
    assert single == ["E11", "E16", "E18"]
    assert all(fails[k] == ["P2"] for k in single)
    assert all(len(v) >= 2 for k, v in fails.items() if k not in single + ["E1"])
    assert all(len(v) >= 2 for k, v in fails.items() if M.BY_LABEL[k].closed and k != "E1")


@pytest.mark.parametrize("label,order", [("E1", 1), ("E2", 2), ("E3", 4), ("E4", 3), ("E5", 6), ("E6", 4)])
def test_deck_groups(label, order):
    g = M.deck_group(label, L=1.3)
    assert g.order == order == M.BY_LABEL[label].holonomy_order
    assert all(np.isclose(np.linalg.det(A), 1.0) for A in g.holonomy)  # orientable
    for A in g.holonomy:  # holonomy normalizes T
        c = np.linalg.solve(g.translations.T, A @ g.translations.T)
        assert np.allclose(c, np.round(c))
    words = g.words(4)
    assert all(g.contains(A, b) for A, b in words)  # every element is (A, b_A + t)
    assert not g.contains(g.reps[-1][0], g.reps[-1][1] + 0.37) or order == 1
    assert M.fixed_point_free(g, radius=3.0)
    rng = np.random.default_rng(7)
    assert M.min_displacement(g, rng.uniform(-2, 2, (400, 3)), radius=2.0) > 0.05
    assert M.displacement_is_constant(g, rng.uniform(-2, 2, (50, 3))) == (label == "E1")
    assert M.parallel_spinor_dimension(g) == (2 if label == "E1" else 0)


def test_hw_generators_square_to_translations():
    g = M.deck_group("E6", L=1.0)
    for s, e in zip(g.generators, np.eye(3)):
        A2, b2 = M.compose(s, s)
        assert np.allclose(A2, np.eye(3)) and np.allclose(b2, e)


def test_topology_json_contract():
    d = M.to_json()
    assert set(d) == {"manifolds", "principles", "curved", "weeks", "evidence_torus",
                      "exclusion_torus", "other_flat", "claims"}
    need = {"label", "name", "closed", "orientable", "holonomy", "translations_only",
            "periodic_spin_possible", "principles"}
    assert all(need <= set(m) for m in d["manifolds"]) and len(d["manifolds"]) == 18
    assert d["manifolds"][0]["periodic_spin_possible"] is True
    assert [p["theorem"] for p in d["principles"]] == ["Theorem A", "Theorem C'", "Theorem B-loc", "Theorem D"]
    assert len(d["evidence_torus"]) == 12 and len(d["exclusion_torus"]) == 17
    assert len(d["other_flat"]) == 9 and len(d["curved"]) == 3 and len(d["weeks"]) == 8
    assert len(d["claims"]) == 29
    assert {"shape", "d_over_D", "lnB_planck", "median_if_true", "q05_if_true", "p_le_planck",
            "sims", "verdict"} <= set(d["exclusion_torus"][0])
    json.dumps(d)


def test_weeks_table_consistency():
    """Prop. 9.4 table: loop/D = systole/(2 chi_*/Rc); V_Weeks/V_ball with the hyperbolic ball
    volume pi (sinh 2x - 2x) Rc^3 (the Euclidean 4 pi x^3/3 misses the 0.015 row by 5%)."""
    for r in M.WEEKS:
        x = r["chi_over_Rc"]
        assert r["shortest_loop_over_D"] == pytest.approx(M.WEEKS_SYSTOLE / (2 * x), abs=0.035)
        ball = math.pi * (math.sinh(2 * x) - 2 * x)
        assert r["volume_ratio"] == pytest.approx(M.WEEKS_VOLUME / ball, rel=0.0018 / x + 0.005, abs=0.05)
        assert (r["circles"].startswith("no")) == (r["omega_k"] < M.WEEKS_OMEGA_K_THRESHOLD)
    ok = np.array([r["omega_k"] for r in M.WEEKS])
    x = np.array([r["chi_over_Rc"] for r in M.WEEKS])
    assert np.interp(M.WEEKS_SYSTOLE / 2, x, ok) == pytest.approx(0.0080, abs=2e-4)


# --- lattice --------------------------------------------------------------------------

@pytest.mark.parametrize("x", [0.05, 0.4, 1.0, 2.3])
def test_ewald_kernels_vs_quadrature(x):
    q = lambda f: quad(f, 1, np.inf, epsabs=0, epsrel=1e-13, limit=200)[0]
    assert L.phi1(x) == pytest.approx(q(lambda t: t * np.exp(-np.pi * x * t)), rel=1e-11)
    assert L.phi2(x) == pytest.approx(q(lambda u: u**-1.5 * np.exp(-np.pi * x * u)), rel=1e-10)
    assert L.dphi1(x) == pytest.approx(q(lambda t: -np.pi * t**2 * np.exp(-np.pi * x * t)), rel=1e-11)
    assert L.dphi2(x) == pytest.approx(q(lambda u: -np.pi * u**-0.5 * np.exp(-np.pi * x * u)), rel=1e-11)
    assert L.d2phi1(x) == pytest.approx(q(lambda t: np.pi**2 * t**3 * np.exp(-np.pi * x * t)), rel=1e-11)
    assert L.d2phi2(x) == pytest.approx(q(lambda u: np.pi**2 * u**0.5 * np.exp(-np.pi * x * u)), rel=1e-11)


@pytest.mark.parametrize("shape", ["SC", "BCC", "FCC"])
def test_z2_certified(shape):
    assert L.epstein_z2(L.cubic_basis(shape)) == pytest.approx(Z2_CERTIFIED[shape], rel=1e-12, abs=0)


def test_z2_gaps_and_rhombohedral_points():
    z = {s: L.epstein_z2(L.cubic_basis(s)) for s in L.SHAPES}
    assert z["BCC"] - z["FCC"] == pytest.approx(L.GAPS_CERTIFIED["BCC-FCC"], abs=5e-14)
    assert z["SC"] - z["BCC"] == pytest.approx(L.GAPS_CERTIFIED["SC-BCC"], abs=5e-14)
    for a, s in ((60.0, "FCC"), (90.0, "SC"), (math.degrees(math.acos(-1 / 3)), "BCC")):
        assert L.epstein_z2(L.rhombohedral_basis(a)) == pytest.approx(Z2_CERTIFIED[s], rel=1e-12)


def test_z2_invariances_and_direct_sum():
    rng = np.random.default_rng(3)
    B = rng.normal(size=(3, 3))
    z = L.epstein_z2(B)
    U = np.array([[1, 2, 0], [0, 1, 0], [3, 5, 1]])
    Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    assert L.epstein_z2(U @ B) == pytest.approx(z, rel=1e-12)
    assert L.epstein_z2(B @ Q) == pytest.approx(z, rel=1e-12)
    assert L.epstein_z2(2.0 * B, unit_covolume=False) == pytest.approx(
        L.epstein_z2(B, unit_covolume=False) / 16.0, rel=1e-12)
    Bu = L.normalize(L.lll_reduce(B)[0], "covolume")
    R = 25.0
    assert L.epstein_z2_direct(Bu, R) + 4 * math.pi / R == pytest.approx(z, abs=2e-3)


def test_dual_reduction_shortest():
    B = np.array([[1.0, 0.0, 0.0], [17.3, 1.0, 0.0], [5.0, 41.2, 0.7]])
    assert np.allclose(B @ L.dual_basis(B).T, 2 * np.pi * np.eye(3))
    Br, U = L.minkowski_reduce(B)
    assert np.allclose(U @ B, Br) and abs(round(np.linalg.det(U))) == 1
    v, lam = L.shortest_vector(B)
    g = L.lattice_vectors(B, lam * 1.5)
    assert lam == pytest.approx(np.linalg.norm(g, axis=1).min(), rel=1e-12)
    assert len(L.shortest_vectors(L.cubic_basis("FCC"))) == 12
    assert len(L.shortest_vectors(L.cubic_basis("BCC"))) == 8
    assert L.normalize(L.cubic_basis("FCC"), "shortest") @ np.ones(3) == pytest.approx(np.full(3, math.sqrt(2)))


def test_gram_gradient():
    rng = np.random.default_rng(5)
    B = L.normalize(L.cubic_basis("BCC") + 0.1 * rng.normal(size=(3, 3)), "covolume")
    h = L.gram(B)
    z, Mg = L.epstein_z2_gram(h, grad=True)
    assert np.trace(Mg @ h) == pytest.approx(-2 * z, rel=1e-12)  # Z is homogeneous of degree -2
    for _ in range(3):
        E = rng.normal(size=(3, 3))
        E = E + E.T
        eps = 1e-6
        fd = (L.epstein_z2_gram(h + eps * E) - L.epstein_z2_gram(h - eps * E)) / (2 * eps)
        assert np.trace(Mg @ E) == pytest.approx(fd, rel=1e-7)


@pytest.mark.parametrize("shape", ["SC", "BCC", "FCC"])
def test_theorem_E_and_hessian(shape):
    B = L.cubic_basis(shape)
    assert np.linalg.norm(L.shape_gradient(B)) < 1e-12  # Theorem E
    H = L.shape_hessian(B)
    ev = np.linalg.eigvalsh(H)
    paper = np.sort(HESSIAN_CERTIFIED[shape])
    tol = np.full(5, 6e-5)
    if shape == "BCC":
        tol[2:] = 3e-4  # [TOP] prints 3.3953; we get 3.39504 (see test_hessian_trace_identity)
    assert np.all(np.abs(ev - paper) <= tol)
    assert np.allclose(H, L.shape_hessian_fd(B, 1e-3), atol=2e-5)


def test_hessian_trace_identity():
    """tr Hess_S Z = (2/3) Z, exactly, since Hess of |g|^-4 traces to (2/3)|g|^-4.

    Fixes the BCC triple at 3.39504 (paper 3.3953 violates the identity by 8e-4)."""
    rng = np.random.default_rng(11)
    for B in [L.cubic_basis(s) for s in L.SHAPES] + [rng.normal(size=(3, 3))]:
        assert np.trace(L.shape_hessian(B)) == pytest.approx(2 / 3 * L.epstein_z2(B), rel=1e-12)
    ev = L.hessian_eigenvalues(L.cubic_basis("BCC"))
    assert ev[2:] == pytest.approx(3.39504, abs=1e-5)
    paper = HESSIAN_CERTIFIED["BCC"]
    assert sum(paper) - 2 / 3 * Z2_CERTIFIED["BCC"] == pytest.approx(8.2e-4, abs=1e-4)


def test_sc_certified_shear():
    E = L.traceless_symmetric_basis()
    H = L.shape_hessian(np.eye(3))
    assert H[2, 2] == pytest.approx(-1.834819428065220, abs=1e-13)
    assert E.shape == (5, 3, 3)
    assert np.allclose(np.einsum("aij,bij->ab", E, E), np.eye(5))


def test_degenerate_inputs_raise():
    for B in (np.zeros((3, 3)), np.array([[1.0, 0, 0], [0, 1, 0], [1, 1, 0]]), np.full((3, 3), np.nan)):
        for f in (L.epstein_z2, L.lll_reduce, L.normalize, L.dual_basis):
            with pytest.raises(ValueError):
                f(B)
    with pytest.raises(ValueError):
        L.epstein_z2_gram(np.diag([1.0, 1.0, -1.0]))
    with pytest.raises(ValueError):  # too flat to enumerate: refused instead of allocating GBs
        L.epstein_z2(np.diag([1.0, 1.0, 1e-7]))


def test_skewed_gram_is_reduced():
    """Z and M = dZ/dh are GL(3,Z)-covariant: Z(U h U^T) = Z(h), M_h = U^T M_{UhU^T} U; and
    homogeneous: Z(c h) = Z(h)/c^2, M(c h) = M(h)/c^3. The skewed input takes the LLL path."""
    B = L.normalize(L.cubic_basis("BCC") + 0.07 * np.random.default_rng(8).normal(size=(3, 3)), "covolume")
    h = L.gram(B)
    U = np.array([[7, 3, 0], [14, 5, 4], [3, 1, 1]])
    hs = U @ h @ U.T
    assert L._box_size(L._box_for(L.EWALD_RADIUS, np.diag(np.linalg.inv(hs)))) > 50 * L._BOX_REDUCE
    z, M = L.epstein_z2_gram(h, grad=True)
    zs, Ms = L.epstein_z2_gram(hs, grad=True)
    assert zs == pytest.approx(z, rel=1e-10)
    assert np.allclose(U.T @ Ms @ U, M, rtol=0, atol=1e-10 * np.abs(M).max())
    z3, M3 = L.epstein_z2_gram(3.7 * h, grad=True)
    assert z3 == pytest.approx(z / 3.7**2, rel=1e-12) and np.allclose(M3, M / 3.7**3, rtol=1e-10)
    Bs = U @ B
    g = L.lattice_vectors(Bs, 1.6)
    g0 = L.lattice_vectors(B, 1.6)
    assert len(g) == len(g0) and np.allclose(np.sort(np.linalg.norm(g, axis=1)), np.linalg.norm(g0, axis=1))


def test_rhombohedral_curve():
    c = L.rhombohedral_curve(np.arange(40.0, 119.6, 0.5))
    a, z = np.array(c["alpha_deg"]), np.array(c["Z"])
    i60, i90 = np.argmin(abs(a - 60)), np.argmin(abs(a - 90))
    assert z[i60] == pytest.approx(15.9621, abs=1e-4) and z[i90] == pytest.approx(16.5323, abs=1e-4)
    assert np.argmin(z) == i60
    assert z[i90] == z[(a > 70) & (a < 105)].max()  # SC is a local max along the family (a saddle)
    j = np.argmin(np.where((a > 100), z, np.inf))
    assert abs(a[j] - 109.47) < 0.6 and z[j] == pytest.approx(15.9719, abs=2e-4)
    assert z[0] > 17.0 and z[-1] > 17.0


def test_lattice_json():
    d = L.to_json(n_alpha=17)
    assert set(d) >= {"bases", "rhombohedral", "critical"}
    assert set(d["critical"]["FCC"]) >= {"Z", "Z_paper", "hessian", "hessian_paper"}
    for s, B in d["bases"].items():
        assert L.shortest_vector(np.array(B))[1] == pytest.approx(1.0)
    json.dumps(d)


# --- shape flow -----------------------------------------------------------------------

def test_haar_sampler_siegel():
    """Siegel: E #{g != 0 : |g| <= r} = 4 pi r^3 / 3 for Haar-random unimodular lattices."""
    r = F.siegel_check(np.random.default_rng(1), n=1000, t_max=F.T_MAX, radii=(0.8, 1.0))
    for mean, err, ref in zip(r["mean_count"], r["stderr"], r["siegel"]):
        assert abs(mean - ref) < 4 * err + 0.01 * ref  # T_MAX = 10: truncation bias below 1 %


def test_kak_density_and_reduction():
    t = F.sample_kak_t(np.random.default_rng(2), 500, 2.0)
    assert np.allclose(t.sum(axis=1), 0) and np.all(np.linalg.norm(t, axis=1) <= 2.0 + 1e-12)
    hs = F.haar_random_grams(np.random.default_rng(3), 20)
    for h in hs:
        assert np.linalg.det(h) == pytest.approx(1.0, rel=1e-9)
        assert F.orthogonality_defect(h) < 1.5


def test_flow_fast():
    hs = F.haar_random_grams(np.random.default_rng(42), 6)
    ends = []
    for h in hs:
        r = F.run_flow(h, kappa=1.0)
        assert r.converged and r.end in ("FCC", "BCC")
        assert np.linalg.det(r.h_end) == pytest.approx(1.0, rel=1e-9)
        assert r.Z[-1] <= r.Z[0]
        ends.append(r.end)
        g = F.run_flow(h, mode="gradient")
        assert g.converged and np.all(np.diff(g.Z) <= 1e-9)  # gradient flow: Z is a Lyapunov function
    assert "FCC" in ends


def test_inertial_flow_lyapunov():
    """Theorem H.2: W = Z + tr u^2 / (2 kappa) is non-increasing and det h stays 1."""
    from scipy.integrate import solve_ivp
    kappa = 2.0
    h0 = F.haar_random_grams(np.random.default_rng(17), 1)[0]
    sol = solve_ivp(F._rhs_factory(kappa, "inertial"), (0.0, 5.0), np.concatenate([h0.ravel(), np.zeros(9)]),
                    method="DOP853", rtol=1e-9, atol=1e-12, dense_output=True)
    W = []
    for t in np.linspace(0.0, 5.0, 41):
        y = sol.sol(t)
        h, X = y[:9].reshape(3, 3), y[9:].reshape(3, 3)
        u = np.linalg.solve(h, X)
        W.append(F.energy(F._sym(h)) + np.trace(u @ u) / (2 * kappa))
        assert np.linalg.det(h) == pytest.approx(1.0, abs=1e-6)  # integrator drift only
    assert np.all(np.diff(W) <= 1e-9) and W[-1] < W[0]


def test_flow_leaves_sc_and_stays_at_minima():
    rng = np.random.default_rng(9)
    S = rng.normal(size=(3, 3)) * 1e-3
    S = S + S.T - 2 * np.trace(S) / 3 * np.eye(3)
    r = F.run_flow(L.gram(L.deform(np.eye(3), S)), kappa=1.0)
    assert r.end in ("FCC", "BCC") and r.Z_end < 16.0  # SC is unstable (Theorem H.2)
    for s in ("FCC", "BCC"):
        r = F.run_flow(L.gram(L.normalize(L.cubic_basis(s), "covolume")), kappa=1.0)
        assert r.end == s and r.converged and r.tau[-1] <= 4.0


def test_relaxation_numbers():
    d = F.relaxation_numbers()
    assert d["computed"] is True
    c = d["ours"]
    assert c["kappa_cE_coeff"] == pytest.approx(d["kappa_cE_coeff"], abs=0.005)
    assert c["fcc_mu_min"] == pytest.approx(d["kappa_cE_rate_fcc"], abs=5e-4)
    assert c["TL4_per_dtau"] == pytest.approx(d["TL4_per_dtau"], abs=0.005)
    assert c["TL_today"] == pytest.approx(d["TL_today"], rel=0.05)


@pytest.mark.slow
def test_basins_n120():
    """[TOP] Prop. H.3: gradient flow from Haar-random lattices ends at FCC in 190/240 cases,
    at BCC otherwise. Same-seed regression check, compared within 3 joint sigma."""
    s = F.basin_statistics(n=120)
    assert s["mode"] == "gradient" and s["kappa"] is None and s["t_max"] == F.T_MAX
    assert s["counts"].get("other", 0) == 0 and "SC" not in s["counts"]
    assert s["unconverged"] == 0
    assert s["counts"]["FCC"] + s["counts"]["BCC"] == 120
    p = 190 / 240
    assert abs(s["fcc_fraction"] - p) < 3 * math.sqrt(p * (1 - p) * (1 / 120 + 1 / 240))
    d = F.to_json(stats=s)
    assert set(d["flows"]) >= {"n", "counts", "paper_counts", "trajectories"}
    assert d["flows"]["computed"] is True and len(d["flows"]["trajectories"]) == 24
    json.dumps(d, allow_nan=False)


# --- circles --------------------------------------------------------------------------

@pytest.mark.parametrize("shape,pairs,angles", [("SC", 3, [90.0]), ("BCC", 4, [70.5]), ("FCC", 6, [60.0, 90.0])])
def test_circle_table_0985(shape, pairs, angles):
    for d in (0.985, 0.99, 0.999):
        t = C.circle_table(shape, d)
        assert t["n_pairs"] == pairs and t["angles_between_centres_deg"] == angles
        assert t["radius_deg"] == pytest.approx(math.degrees(math.acos(d)), abs=1e-9)
    assert round(C.circle_table(shape, 0.985)["radius_deg"], 1) == 9.9
    assert C.circle_table(shape, 1.0)["n_pairs"] == 0 and C.circle_table(shape, 1.2)["n_pairs"] == 0


def test_circle_table_08():
    t = {s: C.circle_table(s, 0.8) for s in ("SC", "BCC", "FCC")}
    assert {s: v["n_pairs"] for s, v in t.items()} == C.FIG4_PAPER["pairs"]
    assert all(round(v["radius_deg"], 1) == 36.9 for v in t.values())
    sh = t["BCC"]["shells"]
    assert [x["n_pairs"] for x in sh] == [4, 3]
    assert sh[1]["radius_deg"] == pytest.approx(math.degrees(math.acos(0.8 * 2 / math.sqrt(3))))


def test_matched_circles_corollary_F1():
    R = np.array(C.fcc_preregistered_rotation()["R"])
    for p in C.circle_pairs(L.cubic_basis("FCC"), 0.9, R):
        g = np.array(p["g_sky"])
        assert np.allclose(p["centre_minus"], -np.array(p["centre_plus"]))  # antipodal
        pts, partners = C.circle_points(g, 36)
        for n, n2 in zip(pts, partners):
            m = C.matched_point(n, g)
            assert np.allclose(m, n2, atol=1e-12)  # zero twist: same azimuth
            assert np.allclose(C.CHI_STAR * m - C.CHI_STAR * n, g)  # same physical point
        assert np.allclose(pts @ np.array(p["centre_minus"]), math.cos(math.radians(p["radius_deg"])))
    with pytest.raises(ValueError):
        C.matched_point([0, 0, 1], [0.1, 0, 0])


def test_ghost_images():
    B = L.cubic_basis("SC")
    im = C.ghost_images((0.0, 0.0), 0.0, B, 0.8, horizon_over_D=0.5)
    assert len(im) == 7 and im[0]["is_source"] and im[0]["inside_horizon"]
    assert all(not e["inside_horizon"] and e["distance_over_D"] == pytest.approx(0.8) for e in im[1:])
    R = np.array(C.fcc_preregistered_rotation()["R"])
    im = C.ghost_images((120.0, 30.0), 0.2, L.cubic_basis("FCC"), 0.55, R, max_distance_over_D=1.5)
    x = 0.2 * C.lb_to_vec(120.0, 30.0)
    Bs = C.scaled_basis(L.cubic_basis("FCC"), 0.55)
    for e in im:
        y = x + R @ (np.array(e["n"]) @ Bs)
        assert np.linalg.norm(y) == pytest.approx(e["distance_over_D"])
        if not e["is_source"]:
            assert np.allclose(C.lb_to_vec(e["l"], e["b"]), y / np.linalg.norm(y))
        assert e["inside_horizon"] == (e["distance_over_D"] <= 0.5)
    n = np.array(np.meshgrid(*[np.arange(-8, 9)] * 3)).reshape(3, -1).T
    r = np.linalg.norm(x + (n @ Bs) @ R.T, axis=1)
    assert len(im) == int(np.sum(r <= 1.5))
    assert sum(e["inside_horizon"] for e in im) == int(np.sum(r <= 0.5)) >= 2
    im0 = C.ghost_images((10.0, 20.0), 0.0, B, 0.8)
    assert im0[0]["l"] is None and im0[0]["b"] is None
    json.dumps(im0, allow_nan=False)
    with pytest.raises(ValueError):
        C.ghost_images((10.0, 20.0), 0.1, B, 0.0)


def test_fcc_preregistered_rotation():
    r = C.fcc_preregistered_rotation()
    R = np.array(r["R"])
    assert np.allclose(R @ R.T, np.eye(3)) and np.linalg.det(R) == pytest.approx(1.0)
    assert r["residual_rms_deg"] < 0.06 and r["residual_max_deg"] < 0.1
    mapped = C.fcc_axes_lattice_frame() @ R.T
    target = C.lb_to_vec(*np.array(FCC_PREREGISTERED_AXES_GAL).T)
    assert np.all(np.abs(target @ mapped.T).max(axis=1) > math.cos(math.radians(0.1)))
    json.dumps(C.to_json())


# --- vacuum ---------------------------------------------------------------------------

def test_vacuum_e1_constant():
    g = V.slice_grid("E1", n=21)
    ex = np.array(g["excess"])
    assert np.ptp(ex) <= 1e-10 * abs(ex.mean())


def test_vacuum_e2_brighter_on_axes():
    g = V.slice_grid("E2", n=41)
    ex = np.array(g["excess"])
    xs = np.array(g["x"])
    assert np.ptp(ex) > 0.05 * ex.max()
    from scipy.ndimage import maximum_filter
    loc = ex == maximum_filter(ex, size=3, mode="nearest")
    maxima = sorted((round(xs[j], 6), round(xs[i], 6)) for i, j in zip(*np.nonzero(loc)))
    assert maxima == sorted((round(a, 6), round(b, 6)) for a, b in V.rotation_axes_xy("E2"))
    on_axis = ex[np.isclose(xs, 0.5)[:, None] & np.isclose(xs, 0.5)[None, :]][0]
    assert on_axis == pytest.approx(ex.max(), rel=1e-9)


def test_vacuum_kernels():
    r = np.array([1e-4, 1e-3])
    assert np.allclose(V.w0(r, 1.0) * 4 * np.pi**2 * r**2, 1.0, rtol=1e-5)  # Hadamard 1/(4 pi^2 r^2)
    spl = V.smoothed_w0_table(1.0, 0.01, 3.0)
    rr = np.array([0.5, 1.0, 2.0])
    assert np.allclose(np.exp(spl(rr)), V.w0(rr, 1.0), rtol=2e-3)
    sm = V.slice_grid("E2", n=21, smooth_sigma=0.35)
    raw = V.slice_grid("E2", n=21)
    assert 0 < np.ptp(sm["excess"]) < np.ptp(raw["excess"])
