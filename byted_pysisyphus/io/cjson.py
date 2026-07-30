import json

import numpy as np

from byted_pysisyphus.constants import ANG2BOHR
from byted_pysisyphus.Geometry import Geometry
from byted_pysisyphus.elem_data import INV_ATOMIC_NUMBERS


def parse_cjson(fn):
    with open(fn) as handle:
        cml = json.load(handle)
    atms = cml["atoms"]
    coords = np.array(atms["coords"]["3d"], dtype=float) * ANG2BOHR
    atom_nums = atms["elements"]["number"]
    atoms = tuple([INV_ATOMIC_NUMBERS[num].capitalize() for num in atom_nums])

    return atoms, coords.flatten()


def geom_from_cjson(fn, **kwargs):
    atoms, coords = parse_cjson(fn)
    geom = Geometry(atoms, coords.flatten(), **kwargs)
    return geom
