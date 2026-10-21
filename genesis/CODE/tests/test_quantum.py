import warnings

import numpy as np

from babyuniverse.model import model_from_dimensionless
from babyuniverse import quantum as Q

warnings.simplefilter("ignore")


def test_hamilton_equation_all_sectors():
    m = model_from_dimensionless(-0.5, 4.5)
    for R, sg in ((0.5007, 1), (0.5007, -1), (1.5, 1), (1.5, -1)):
        c = Q.hamilton_equation_check(m, R, sign=sg)
        assert abs(c["dHdP"] / c["dRdt"] - 1) < 1e-7
        assert abs(c["H_at_P"] - m.M) < 1e-10
        if "H_closed" in c:
            assert abs(c["H_closed"] - m.M) < 1e-12
    m2 = model_from_dimensionless(0.6, 1.2)
    c = Q.hamilton_equation_check(m2, 0.495)
    assert c["sector"] == (-1, -1)
    assert abs(c["dHdP"] / c["dRdt"] - 1) < 1e-7


def test_difference_operator_symbol_on_euclidean_orbit():
    m = model_from_dimensionless(-0.5, 4.5)
    assert Q.eikonal_check(m) < 1e-10


def test_naive_quantization_inequivalent():
    m = model_from_dimensionless(0.0, 6.0)
    from babyuniverse.tunnelling import B_canonical
    Bc = B_canonical(m)
    assert abs(Q.B_naive(m, "membrane") / Bc - 1) > 0.1
    assert abs(Q.B_naive(m, "unit") / Bc - 1) > 0.01
