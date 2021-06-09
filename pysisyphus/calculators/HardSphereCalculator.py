# [1] https://doi.org/10.1002/jcc.26495
#     Habershon, 2021

import itertools as it

import numpy as np

from pysisyphus.helpers_pure import get_molecular_radius


class HardSphereCalculator:
    def __init__(self, geom, frags, kappa=1.0, permutations=False, **kwargs):
        """Intra-Image Inter-Molecular Hard-Sphere force.

        See A.2. in [1], Eq. (A1).
        """
        self.frags = frags
        self.kappa = kappa

        self.frag_num = len(self.frags)
        self.frag_sizes = np.array([len(frag) for frag in self.frags])
        self.pair_inds = np.array(list(it.combinations(range(self.frag_num), 2)))
        it_func = it.permutations if permutations else it.combinations
        self.pair_inds = np.array(list(it_func(range(self.frag_num), 2)))
        self.frag_inds = np.array([m for m, n in self.pair_inds])

        c3d = geom.coords3d
        frag_c3ds = [c3d[frag] for frag in self.frags]
        self.frag_radii = [get_molecular_radius(frag_c3d) for frag_c3d in frag_c3ds]
        self.radii_sums = np.array(
            [self.frag_radii[i] + self.frag_radii[j] for i, j in self.pair_inds]
        )

    def get_forces(self, atoms, coords, kappa=None):
        if kappa is None:
            kappa = self.kappa
        c3d = coords.reshape(-1, 3)

        centroids = np.array([c3d[frag].mean(axis=0) for frag in self.frags])
        mm, nn = np.array(self.pair_inds).T
        gdiffs = centroids[mm] - centroids[nn]
        gnorms = np.linalg.norm(gdiffs, axis=1)
        H = (gnorms < self.radii_sums).astype(int)
        N = np.zeros_like(H)
        for frag_ind, H_ in zip(self.frag_inds, H):
            N[frag_ind] += H_
        nonzero = H > 0
        N *= 3 * self.frag_sizes[self.frag_inds]
        phi = np.where(N > 0, kappa / N * (gnorms - self.radii_sums), 0)
        frag_forces = (phi * H / gnorms)[:, None] * gdiffs
        forces = np.zeros_like(c3d)
        # Distribute forces onto fragments
        for frag_ind, ff in zip(self.frag_inds, frag_forces):
            forces[self.frags[frag_ind]] += ff
        # As far as I can tell so far this is actually more like a gradient.
        # Moving into direction "forces" just increases it, so we multiply with
        # -1 to get actual forces.
        forces = -forces.flatten()
        return {"energy": 1, "forces": forces}
