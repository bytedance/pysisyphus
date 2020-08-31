# [1] https://aip.scitation.org/doi/abs/10.1063/1.3664901
#     Behn, 2011, Freezing string method
# [2] https://aip.scitation.org/doi/pdf/10.1063/1.4804162
#     Zimmerman, 2013, Growing string with interpolation and optimization
#                      in internal coordiantes

import numpy as np

from pysisyphus.constants import ANG2BOHR
from pysisyphus.helpers import procrustes
from pysisyphus.optimizers.closures import bfgs_multiply
from pysisyphus.optimizers.Optimizer import Optimizer
from pysisyphus.optimizers.restrict_step import scale_by_max_step


class StringOptimizer(Optimizer):

    def __init__(self, geometry, gamma=1.25, max_step=0.1, stop_in_when_full=-1,
                 keep_last=10, lbfgs_when_full=True, gamma_mult=False, **kwargs):
        super().__init__(geometry, max_step=max_step, **kwargs)

        assert self.is_cos, \
            "StringOptimizer is only intended to be used with COS objects."

        # gamma = 1.25 Hartree/Bohr² ~ 5 Hartree/Angstrom²
        self.gamma = float(gamma)
        self.stop_in_when_full = int(stop_in_when_full)
        self.keep_last = int(keep_last)
        assert self.keep_last >= 0
        self.lbfgs_when_full = lbfgs_when_full
        if self.lbfgs_when_full and (self.keep_last == 0):
            print("lbfgs_when_full is True, but keep_last is 0!")
        self.gamma_mult = gamma_mult

        # Add one as we later subtract 1 before we check if this value is 0.
        self.stop_in = self.stop_in_when_full + 1
        self.is_cart_opt = self.geometry.coord_type == "cart"
        self.s_list = list()
        self.y_list = list()
        self.inds = list()

    def prepare_opt(self):
        if self.align and self.is_cart_opt:
            procrustes(self.geometry)

    def reset(self):
        pass

    def restrict_step_components(self, steps):
        too_big = np.abs(steps) > self.max_step
        self.log(f"Found {np.sum(too_big)} big step components.")
        signs = np.sign(steps[too_big])
        steps[too_big] = signs * self.max_step
        return steps

    def check_convergence(self, *args, **kwargs):
        # Normal convergence check with gradients etc.
        converged = super().check_convergence(*args, **kwargs)

        if self.geometry.fully_grown:
            # We only start decrementing the counter after the string is fully grown.
            self.stop_in -= 1
            self.log(f"String is fully grown. Stopping in {self.stop_in} cycles.")

        full_stop = self.geometry.fully_grown and (self.stop_in == 0)
        # full_stop will take precedence when True
        return full_stop or converged

    def optimize(self):
        new_image_inds = self.geometry.new_image_inds
        string_size_changed = len(new_image_inds) > 0

        if self.align and string_size_changed and self.is_cart_opt:
            procrustes(self.geometry)
            self.log("Aligned string.")

        forces = self.geometry.forces
        self.energies.append(self.geometry.energy)
        self.forces.append(forces)

        cur_size = self.geometry.string_size
        add_to_list = (
            self.keep_last > 0     # Only add to s_list and y_list if we want to keep
            and self.cur_cycle > 0 # cycles and if we can actually calculate differences.
            and (not self.lbfgs_when_full # Add when LBFGS is allowed before fully grown.
                 or self.lbfgs_when_full and self.geometry.fully_grown
                 and not string_size_changed # Don't add when string has to be fully grown
                                             # but grew fully in this cycle.
            )
        )
        if add_to_list:
            inds = list(range(cur_size))
            try:
                y = self.forces[-2] - forces
                s = self.coords[-1] - self.coords[-2]
            # Will be raised when the string grew in the previous cycle.
            except ValueError:
                cur_forces = np.delete(forces.reshape(cur_size, -1),
                                       new_image_inds, axis=0).flatten()
                y = self.forces[-2] - cur_forces
                cur_coords = np.delete(self.coords[-1].reshape(cur_size, -1),
                                       new_image_inds, axis=0).flatten()
                s = self.coords[-2] - cur_coords
                inds = np.delete(inds, new_image_inds)

            self.s_list.append(s)
            self.y_list.append(y)
            self.inds.append(inds)
            # Drop oldest vectors
            self.s_list = self.s_list[-self.keep_last:]
            self.y_list = self.y_list[-self.keep_last:]
            self.inds = self.inds[-self.keep_last:]

        # Results in steepest descent step for empty s_list & y_list
        step = bfgs_multiply(self.s_list, self.y_list, forces, gamma_mult=self.gamma_mult,
                             inds=self.inds, cur_size=cur_size, logger=self.logger)

        # Conjugate gradient step
        if self.keep_last == 0 and self.cur_cycle > 0 and not string_size_changed:
            cur_norm = np.linalg.norm(forces)
            prev_norm = np.linalg.norm(self.forces[-2])
            quot = min(cur_norm**2 / prev_norm**2, 1)
            step = forces / self.gamma + quot*self.steps[-1]

        step = scale_by_max_step(step, self.max_step)

        return step
