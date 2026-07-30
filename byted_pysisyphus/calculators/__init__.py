import logging

__all__ = [
    "AFIR",
    "AtomAtomTransTorque",
    "Composite",
    "ConicalIntersection",
    "DFTBp",
    "DFTD4",
    "Dimer",
    "Dummy",
    "EnergyMin",
    "EGO",
    "ExternalPotential",
    "FakeASE",
    "Gaussian09",
    "Gaussian16",
    "HardSphere",
    "IPIServer",
    "LennardJones",
    "MOPAC",
    "MultiCalc",
    "OBabel",
    "ONIOM",
    "OpenMolcas",
    "ORCA",
    "ORCA5",
    "Psi4",
    "PyPsi4",
    "PyXTB",
    "Remote",
    "TIP3P",
    "Turbomole",
    "TransTorque",
    "XTB",
    "GXTB",
    "CFOUR",
]


from byted_pysisyphus.calculators.AFIR import AFIR
from byted_pysisyphus.calculators.AtomAtomTransTorque import AtomAtomTransTorque
from byted_pysisyphus.calculators.Composite import Composite
from byted_pysisyphus.calculators.ConicalIntersection import ConicalIntersection
from byted_pysisyphus.calculators.DFTBp import DFTBp
from byted_pysisyphus.calculators.DFTD4 import DFTD4
from byted_pysisyphus.calculators.Dimer import Dimer
from byted_pysisyphus.calculators.Dummy import Dummy
from byted_pysisyphus.calculators.EnergyMin import EnergyMin
from byted_pysisyphus.calculators.EGO import EGO
from byted_pysisyphus.calculators.ExternalPotential import ExternalPotential
from byted_pysisyphus.calculators.FakeASE import FakeASE
from byted_pysisyphus.calculators.Gaussian09 import Gaussian09
from byted_pysisyphus.calculators.Gaussian16 import Gaussian16
from byted_pysisyphus.calculators.IPIServer import IPIServer
from byted_pysisyphus.calculators.HardSphere import HardSphere, PWHardSphere
from byted_pysisyphus.calculators.LennardJones import LennardJones
from byted_pysisyphus.calculators.MultiCalc import MultiCalc
from byted_pysisyphus.calculators.MOPAC import MOPAC
from byted_pysisyphus.calculators.Psi4 import Psi4
from byted_pysisyphus.calculators.OBabel import OBabel
from byted_pysisyphus.calculators.ONIOMv2 import ONIOM
from byted_pysisyphus.calculators.OpenMolcas import OpenMolcas
from byted_pysisyphus.calculators.ORCA import ORCA
from byted_pysisyphus.calculators.ORCA5 import ORCA5
from byted_pysisyphus.calculators.PyPsi4 import PyPsi4
from byted_pysisyphus.calculators.PyXTB import PyXTB
from byted_pysisyphus.calculators.Remote import Remote
from byted_pysisyphus.calculators.TIP3P import TIP3P
from byted_pysisyphus.calculators.TransTorque import TransTorque
from byted_pysisyphus.calculators.Turbomole import Turbomole
from byted_pysisyphus.calculators.XTB import GXTB, XTB
from byted_pysisyphus.calculators.CFOUR import CFOUR


logger = logging.getLogger("dimer")
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler("dimer.log", mode="w", delay=True)
fmt_str = "%(message)s"
formatter = logging.Formatter(fmt_str)
handler.setFormatter(formatter)
logger.addHandler(handler)
