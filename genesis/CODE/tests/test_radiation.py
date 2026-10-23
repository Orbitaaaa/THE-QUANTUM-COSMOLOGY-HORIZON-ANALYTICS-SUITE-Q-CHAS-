"""Round 7: Davies-Fulling flux of an explicit 1+1 field on the collapse branch (laboratory state)."""
import warnings

from babyuniverse import radiation as Rd
from babyuniverse.model import model_from_dimensionless, ThinWallModel

warnings.simplefilter("ignore")


def test_collapse_branch_flux_tends_to_hawking_with_delayed_switch_on():
    for m in (model_from_dimensionless(0.6, 2.4), ThinWallModel(1.0, 0.3, 0.1)):
        f = Rd.flux_lab(m, n=3000)
        assert 0.95 < f["late_ratio"] < 1.01                    # kappa^2 / 48 pi, kappa = 1/4GM
        assert 3.0 < f["switch_on_in_4GM"] < 4.0
        assert -2.6 < f["excess_in_thermal_units"] < -1.4      # deficit relative to immediate thermal emission (daughter, D1)
        assert f["F_at_u0_over_FH"] < 0                        # negative-energy precursor at the turning-point signal


def test_eternal_history_ray_map_reproduces_hawking():
    f = Rd.flux(model_from_dimensionless(0.6, 2.4), n=1500)
    assert abs(f["late_ratio"] - 1) < 0.01
