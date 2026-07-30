import pytest

from byted_pysisyphus.calculators import ORCA
from byted_pysisyphus.helpers import geom_loader
from byted_pysisyphus.init_logging import init_logging
from byted_pysisyphus.testing import using


@using("orca")
def test_orca_stable():
    init_logging()
    geom = geom_loader("lib:ozone.xyz")

    orca_kwargs = {
        "keywords": "UHF def2-SVP",
        "pal": 1,
        "mult": 1,
        "charge": 0,
    }
    calc = ORCA(**orca_kwargs)
    geom.set_calculator(calc)

    unstable_energy = calc.get_energy(geom.atoms, geom.coords)["energy"]
    assert unstable_energy == pytest.approx(-224.0547837492)

    calc.do_stable = True
    stable_energy = calc.get_energy(geom.atoms, geom.coords)["energy"]
    assert stable_energy == pytest.approx(-224.1476205193)
