"""Every symbolic identity of JUNCTION_CONDITIONS.md / QUANTIZATION.md."""
import pytest

from babyuniverse import symbolic as S


@pytest.mark.parametrize("fn", [S.extrinsic_curvature_static, S.israel_equations, S.squaring_and_signs,
                                S.dimensionless_form, S.mass_formula, S.reduced_hamiltonian,
                                S.euclidean_angle_identity, S.closed_forms_M0, S.flat_space_limit,
                                S.wall_equation_of_motion])
def test_identity(fn):
    fn()  # raises AssertionError on failure
