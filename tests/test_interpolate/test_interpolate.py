import pytest

# from byted_pysisyphus.xyzloader import write_geoms_to_trj
from byted_pysisyphus.helpers import geom_loader
from byted_pysisyphus.interpolate.Interpolator import Interpolator
from byted_pysisyphus.interpolate.LST import LST
from byted_pysisyphus.interpolate.IDPP import IDPP
from byted_pysisyphus.interpolate.Redund import Redund
from byted_pysisyphus.interpolate.Geodesic import Geodesic
from byted_pysisyphus.testing import using


def test_idpp():
    initial = geom_loader("lib:09_htransfer_product.xyz")
    final = geom_loader("lib:10_po_diss_product_xtbopt.xyz")

    geoms = (initial, final)
    idpp = IDPP(geoms, 18, align=True)
    geoms = idpp.interpolate_all()
    # idpp.all_geoms_to_trj("idpp_opt.trj")

    assert len(geoms) == 20


@pytest.mark.parametrize(
    "interpol_cls",
    [
        Interpolator,
        LST,
        IDPP,
        Redund,
        pytest.param(Geodesic, marks=using("geodesic")),
    ],
)
def test_ala_dipeptide_interpol(interpol_cls):
    initial = geom_loader("lib:dipeptide_init.xyz")
    final = geom_loader("lib:dipeptide_fin.xyz")

    geoms = (initial, final)
    interpolator = interpol_cls(geoms, 28, align=True)
    geoms = interpolator.interpolate_all()
    # interpolator.all_geoms_to_trj("interpolated.trj")

    assert len(geoms) == 30
