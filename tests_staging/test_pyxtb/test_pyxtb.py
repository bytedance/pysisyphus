#!/usr/bin/env python3

import numpy as np
import pytest

from byted_pysisyphus.helpers import geom_from_library
from byted_pysisyphus.calculators.PyXTB import PyXTB
from byted_pysisyphus.calculators.XTB import XTB


# module load intel/2018.0.33
# module load xtb/6.2

@pytest.mark.skip
def test_pyxtb():
    geom = geom_from_library("benzene.xyz")
    pyxtb = PyXTB()
    geom.set_calculator(pyxtb)

    forces = geom.forces
    energy = geom.energy

    # xtb = XTB()
    # geom.set_calculator(xtb)
    # ref_forces = geom.forces
    # ref_energy = geom.energy
    # np.testing.assert_allclose(energy, ref_energy) 
    # np.testing.assert_allclose(forces, ref_forces, rtol=2e-2)


def test_pypsi4():
    from byted_pysisyphus.calculators.PyPsi4 import PyPsi4
    from byted_pysisyphus.calculators.Psi4 import Psi4

    kwargs = {
        "method": "scf",
        "basis": "sto-3g",
    }

    geom = geom_from_library("benzene.xyz")
    pypsi4 = PyPsi4(**kwargs)
    geom.set_calculator(pypsi4)

    forces = geom.forces
    energy = geom.energy

    # psi4 = Psi4(**kwargs)
    # geom.set_calculator(psi4)
    # ref_forces = geom.forces
    # ref_energy = geom.energy
    # np.testing.assert_allclose(energy, ref_energy) 
    # np.testing.assert_allclose(forces, ref_forces, rtol=2e-2)


if __name__ == "__main__":
    test_pyxtb()
    test_pypsi4()
