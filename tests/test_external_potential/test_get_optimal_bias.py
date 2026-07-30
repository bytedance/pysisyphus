import pytest

from byted_pysisyphus.helpers import geom_loader
from byted_pysisyphus.calculators import XTB
from byted_pysisyphus.drivers.opt import get_optimal_bias
from byted_pysisyphus.testing import using


@using("xtb")
def test_get_optimal_bias():
    geom = geom_loader("lib:uracil_dimer.xyz", coord_type="redund")

    def calc_getter():
        return XTB(gfn="ff")

    opt_key = "lbfgs"
    opt_kwargs = {
        "mu_reg": 0.1,
        "dump": True,
        "max_cycles": 750,
    }

    _, k_opt, valid_k = get_optimal_bias(
        geom, calc_getter, opt_key, opt_kwargs, k_max=0.1, rmsd_target=0.05
    )
    assert k_opt == pytest.approx(-0.0375)
    assert valid_k
