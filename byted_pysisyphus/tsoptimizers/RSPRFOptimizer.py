# See [1] https://pubs.acs.org/doi/pdf/10.1021/j100247a015
#         Banerjee, 1985
#     [2] https://aip.scitation.org/doi/abs/10.1063/1.2104507
#         Heyden, 2005
#     [3] https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.540070402
#         Baker, 1985
#     [4] https://link.springer.com/article/10.1007/s002140050387
#         Besalu, 1998


import numpy as np

from byted_pysisyphus.tsoptimizers.TSHessianOptimizer import TSHessianOptimizer


class RSPRFOptimizer(TSHessianOptimizer):
    def optimize(self):
        energy, gradient, H, eigvals, eigvecs, resetted = self.housekeeping()
        self.update_ts_mode(eigvals, eigvecs)

        # Transform gradient to eigensystem of hessian
        gradient_trans = eigvecs.T.dot(gradient)
        # Minimize energy along all modes, except the TS-mode
        min_indices = [i for i in range(gradient_trans.size) if i not in self.roots]
        # Maximize energy along all requested modes.
        max_indices = [i for i in range(gradient_trans.size) if i in self.roots]
        # Get line search steps, if requested.
        ip_step_trans, gradient_trans = self.step_and_grad_from_line_search(
            energy,
            gradient_trans,
            eigvecs,
            min_indices,
            max_indices,
        )

        """In the RS-(P)RFO method we have to scale the matrices with alpha.
        Unscaled matrix (Eq. 8) in [1]:
            (H  g) (x)          (S 0) (x)
                       = lambda
            (g+ 0) (1)          (0 1) (1)
        with
            S = alpha * Identity matrix
        and multiplying from the left with the inverse of the scaling matrix
            (1/alpha 0)

            (0       1)
        we get
            (1/alpha 0) (H  g) (x)          (x)
                                   = lambda
            (0       1) (g+ 0) (1)          (1)
        eventually leading to the scaled matrix:
            (H/alpha  g/alpha) (x)          (x)
                                   = lambda     .
            (g+             0) (1)          (1)
        """

        alpha = self.alpha0
        for mu in range(self.max_micro_cycles):
            self.log(f"RS-PRFO micro cycle {mu:02d}, alpha={alpha:.6f}")

            # Maximize energy along the chosen TS mode. The matrix is hardcoded
            # as 2x2, so only first-order saddle point searches are supported.
            H_aug_max = self.get_augmented_hessian(
                eigvals[max_indices], gradient_trans[max_indices], alpha
            )
            step_max, eigval_max, nu_max, self.prev_eigvec_max = self.solve_rfo(
                H_aug_max, "max", prev_eigvec=self.prev_eigvec_max
            )

            # Minimize energy along all modes, but the TS mode.
            H_aug_min = self.get_augmented_hessian(
                eigvals[min_indices], gradient_trans[min_indices], alpha
            )
            step_min, eigval_min, nu_min, self.prev_eigvec_min = self.solve_rfo(
                H_aug_min, "min", prev_eigvec=self.prev_eigvec_min
            )

            # Calculate overlap between directions over the course of the micro cycles
            # if mu == 0:
            # TODO: convert back to original space
            # ref_step_max = step_max.copy()
            # ref_step_min = step_min.copy()
            min_norm = np.linalg.norm(step_min)
            max_norm = np.linalg.norm(step_max)
            self.log(f"norm(step_max)={max_norm:.6f}")
            self.log(f"norm(step_min)={min_norm:.6f}")
            self.log(f"norm(step_max)/norm(step_min)={max_norm/min_norm:.2%}")
            # Calculate overlaps with originally proposed step in mu == 0
            # TODO: convert back to original space
            # max_ovlp = ref_step_max @ step_max
            # min_ovlp = ref_step_min @ step_min

            # As of Eq. (8a) of [4] max_eigval and min_eigval also
            # correspond to:
            # max_eigval = -forces_trans[max_indices].dot(max_step)
            # min_eigval = -forces_trans[min_indices].dot(min_step)

            # Create the full PRFO step
            step = np.zeros_like(gradient_trans)
            step[max_indices] = step_max
            step[min_indices] = step_min
            step_norm = np.linalg.norm(step)
            self.log(f"norm(step)={step_norm:.6f}")

            inside_trust = step_norm <= self.trust_radius
            if inside_trust:
                self.log(
                    "Restricted step satisfies trust radius of "
                    f"{self.trust_radius:.6f}"
                )
                self.log(
                    f"Micro-cycles converged in cycle {mu:02d} with "
                    f"alpha={alpha:.6f}!"
                )
                break

            # Derivative of the squared step w.r.t. alpha. Each summand
            #   g_i**2 / (eigval_i - rfo_eigval * alpha)**3
            # is singular when an RFO eigenvalue meets a Hessian eigenvalue. Near
            # convergence the numerator g_i -> 0 as well, so the physical limit of
            # the 0/0 term is 0. Evaluate it that way (safe divide) instead of
            # letting it produce nan, which otherwise propagates into alpha and
            # either blows the step up or leaves it degenerate (a null-space step
            # that does not move the geometry, ending in ZeroStepLength).
            def _dstep2(grad_sub, eigvals_sub, rfo_eigval, quad):
                num = grad_sub ** 2
                den = (eigvals_sub - rfo_eigval * alpha) ** 3
                quot = np.divide(num, den, out=np.zeros_like(num), where=(den != 0.0))
                return 2 * rfo_eigval / (1 + quad) * np.sum(quot)

            # max subspace
            dstep2_dalpha_max = _dstep2(
                gradient_trans[max_indices], eigvals[max_indices], eigval_max,
                step_max.dot(step_max) ** 2 * alpha,
            )
            # min subspace
            dstep2_dalpha_min = _dstep2(
                gradient_trans[min_indices], eigvals[min_indices], eigval_min,
                step_min.dot(step_min) * alpha,
            )
            dstep2_dalpha = dstep2_dalpha_max + dstep2_dalpha_min
            # If the derivative still vanishes or is non-finite, alpha cannot be
            # refined further; stop and let the trust-radius scaling below bound
            # the step.
            if (not np.isfinite(dstep2_dalpha)) or (dstep2_dalpha == 0.0):
                self.log("Restricted-step derivative not usable; stopping micro-cycles.")
                break
            # Update alpha
            alpha_step = (
                2 * (self.trust_radius * step_norm - step_norm**2) / dstep2_dalpha
            )
            if not np.isfinite(alpha_step):
                self.log("Non-finite alpha step; stopping micro-cycles.")
                break
            alpha += alpha_step

        # Right now the step is still given in the Hessians eigensystem. We
        # transform it back now.
        step += ip_step_trans
        step = eigvecs.dot(step)
        step_norm = np.linalg.norm(step)

        # Fall back to naive scaling when the restricted-step procedure did not
        # yield a sane step. This covers max_micro_cycles == 1 (RS disabled), an
        # early break from the singularity guard above, and RFO/alpha singularities
        # that leave the step non-finite.
        if not np.all(np.isfinite(step)):
            self.log(
                "Non-finite RS-PRFO step; falling back to a trust-radius step "
                "along the gradient direction."
            )
            grad_norm = np.linalg.norm(gradient)
            if grad_norm == 0.0:
                step = np.zeros_like(step)
            else:
                step = -gradient / grad_norm * self.trust_radius
            step_norm = np.linalg.norm(step)
        elif step_norm > self.trust_radius:
            step = step / step_norm * self.trust_radius
            step_norm = np.linalg.norm(step)
        self.log(f"norm(step)={np.linalg.norm(step):.6f}")

        # Eq. (6) from [4] seems erronous ... the prediction is usually only ~50%
        # of the actual change ...
        # predicted_energy_change = 1/2 * (eigval_max / nu_max**2 + eigval_min / nu_min**2)
        # self.predicted_energy_changes.append(predicted_energy_change)

        self.predicted_energy_changes.append(self.rfo_model(gradient, self.H, step))

        self.log("")
        return step
