"""Round 6p: independent mathematical audit (light checks; the full ledger is produced by scripts/round6p.py)."""
import warnings

from babyuniverse import audit6p as A

warnings.simplefilter("ignore")


def test_symbolic_audit():
    rows = (A.check_hayward_extremality_and_temperature() + A.check_closed_dust_inflection() + A.check_hydro_acceleration_equation()
            + A.check_dl_radius() + A.check_headon_4d())
    assert all(r["passed"] for r in rows), [r for r in rows if not r["passed"]]


def test_curvature_formulas_by_direct_riemann():
    rows = A.check_kretschmann_formula(4) + A.check_frw_kretschmann(4)
    assert all(r["passed"] for r in rows)


def test_independent_numerics():
    rows = A.check_n_acc_by_time_integration() + A.check_curvature_bound_analytic() + A.check_generator_vs_sigma()
    assert all(r["passed"] for r in rows), [r for r in rows if not r["passed"]]
