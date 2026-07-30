import pytest

from byted_pysisyphus.calculators.PyXTB import PyXTB
from byted_pysisyphus.helpers import geom_loader
from byted_pysisyphus.testing import using


@using("xtb")
@using("pyxtb")
@pytest.mark.parametrize(
    "gfn, ref_energy", [
        (0, -4.36679493),
        (1, -5.76865989),
        (2, -5.07043090),
        ("ff", -0.327181905),
    ]
)
def test_pyxtb(gfn, ref_energy):
    geom = geom_loader("lib:h2o.xyz")
    calc = PyXTB(gfn=gfn)
    geom.set_calculator(calc)

    energy = geom.energy
    forces = geom.forces
    assert energy == pytest.approx(ref_energy)
