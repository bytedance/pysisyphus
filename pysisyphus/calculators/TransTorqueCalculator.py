# [1] https://doi.org/10.1002/jcc.26495
#     Habershon, 2021


import numpy as np


def get_trans_torque_forces(
    mfrag,
    a_coords3d,
    b_coords3d,
    a_mats,
    b_mats,
    m,
    frags,
    N_inv,
    weight_func=None,
    skip=True,
    kappa=1,
):
    mcoords3d = a_coords3d[mfrag]
    gm = mcoords3d.mean(axis=0)

    if weight_func is None:

        def weight_func(m, n, a, b):
            return 1

    trans_vec = np.zeros(3)
    rot_vec = np.zeros(3)
    for n, nfrag in enumerate(frags):
        if skip and (m == n):
            continue
        amn = a_mats[(m, n)]
        bnm = b_mats[(n, m)]
        for a in amn:
            for b in bnm:
                rd = b_coords3d[b] - a_coords3d[a]
                gd = a_coords3d[a] - gm
                weight = weight_func(a, b, m, n)

                trans_vec += weight * abs(rd.dot(gd)) * rd / np.linalg.norm(rd)
                if np.isnan(trans_vec).any():
                    import pdb

                    pdb.set_trace()
                    pass
                rot_vec += weight * np.cross(rd, gd)
    trans_vec *= N_inv
    rot_vec *= N_inv
    forces = kappa * (np.cross(-rot_vec, mcoords3d - gm) + trans_vec[None, :])
    return forces


class TransTorqueCalculator:
    def __init__(
        self,
        frags,
        iter_frags,
        b_coords3d,
        a_mats,
        b_mats,
        weight_func=None,
        skip=True,
        kappa=1.0,
    ):
        """Translational and torque forces.
        See A.4. [1], Eqs. (A3) - (A5).
        """
        self.frags = frags
        self.iter_frags = iter_frags
        self.b_coords3d = b_coords3d
        self.a_mats = a_mats
        self.b_mats = b_mats
        self.weight_func = weight_func
        self.kappa = kappa
        self.skip = skip

        self.set_N_invs()

    def set_N_invs(self):
        Ns = np.zeros(len(self.frags))
        for m, mfrag in enumerate(self.frags):
            for n, _ in enumerate(self.iter_frags):
                if self.skip and (m == n):
                    continue
                amn = self.a_mats[(m, n)]
                bnm = self.b_mats[(n, m)]
                Ns[m] += len(amn) * len(bnm)
            Ns[m] *= 3 * len(mfrag)
        self.N_invs = 1 / np.array(Ns)

    def get_forces(self, atoms, coords, kappa=None):
        if kappa is None:
            kappa = self.kappa
        c3d = coords.reshape(-1, 3)
        forces = np.zeros_like(c3d)

        for m, mfrag in enumerate(self.frags):
            N_inv = self.N_invs[m]
            tt_forces = get_trans_torque_forces(
                mfrag,
                c3d,
                self.b_coords3d,
                self.a_mats,
                self.b_mats,
                m,
                self.iter_frags,
                N_inv,
                weight_func=self.weight_func,
                skip=self.skip,
            )
            forces[mfrag] = tt_forces

        return {"energy": 1, "forces": forces.flatten()}
