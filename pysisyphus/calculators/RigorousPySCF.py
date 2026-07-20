import os
import shutil

import numpy as np
try:
    import cupy as cp
except ImportError:
    cp = None

from pyscf import gto, lib, dft, scf

from pysisyphus.calculators.OverlapCalculator import OverlapCalculator
from pysisyphus.helpers import geom_loader


class RigorousPySCF(OverlapCalculator):
    def __init__(
        self,

        # charge, # Handled by parent Calculator class
        # spin, # Handled by parent Calculator class (mult)
        # threads, # Handled by parent Calculator class (pal)
        # max_memory, # Handled by parent Calculator class (mem)
        method_xc,
        method_nlc,
        method_disp,
        grids_atom_grid,
        grids_level,
        nlcgrids_atom_grid,
        nlcgrids_level,
        basis,
        ecp,
        verbose,
        scf_conv_tol,
        direct_scf_tol,
        retry_soscf,
        basis_linear_dependency_threshold,
        with_df,
        auxbasis,
        with_gpu,
        with_lowmem,
        unrestricted,
        guessmix,
        grad_grid_response,
        hess_grid_response,
        with_tddft,
        tddft_options_tda,
        tddft_options_nstates,
        tddft_options_singlet,
        tddft_options_roots_for_tdgrad,
        tddft_options_lr_pcm,
        tddft_options_conv_tol,
        with_solvent,
        solvent_method,
        solvent_eps,
        solvent_solvent,
        solvent_lebedev_order,
        uniform_external_electric_field,

        ### Below are the options not specified by single point implementation, but are available in original pysisyphus-pyscf interface, and Henry think it makes sense to keep them.

        keep_chk,
        pruning,
        **kwargs,
    ):
        # The inputs are assumed to be valid, i.e. the validation step should be performed before calling this constructor.

        super().__init__(**kwargs)

        assert "method" not in kwargs
        assert "xc" not in kwargs
        assert "pseudo" not in kwargs

        self.pyscf_configs = {
            "method_xc": method_xc,
            "method_nlc": method_nlc,
            "method_disp": method_disp,
            "grids_atom_grid": grids_atom_grid,
            "grids_level": grids_level,
            "nlcgrids_atom_grid": nlcgrids_atom_grid,
            "nlcgrids_level": nlcgrids_level,
            "basis": basis,
            "ecp": ecp,
            "verbose": verbose,
            "scf_conv_tol": scf_conv_tol,
            "direct_scf_tol": direct_scf_tol,
            "retry_soscf": retry_soscf,
            "basis_linear_dependency_threshold": basis_linear_dependency_threshold,
            "with_df": with_df,
            "auxbasis": auxbasis,
            "with_gpu": with_gpu,
            "with_lowmem": with_lowmem,
            "unrestricted": unrestricted,
            "guessmix": guessmix,
            "grad_grid_response": grad_grid_response,
            "hess_grid_response": hess_grid_response,
            "with_tddft": with_tddft,
            "tddft_options_tda": tddft_options_tda,
            "tddft_options_nstates": tddft_options_nstates,
            "tddft_options_singlet": tddft_options_singlet,
            "tddft_options_roots_for_tdgrad": tddft_options_roots_for_tdgrad,
            "tddft_options_lr_pcm": tddft_options_lr_pcm,
            "tddft_options_conv_tol": tddft_options_conv_tol,
            "with_solvent": with_solvent,
            "solvent_method": solvent_method,
            "solvent_eps": solvent_eps,
            "solvent_solvent": solvent_solvent,
            "solvent_lebedev_order": solvent_lebedev_order,
            "uniform_external_electric_field": uniform_external_electric_field,

            "keep_chk": keep_chk,
            "pruning": pruning,
        }

        inherited_configs = {
            "charge": self.charge,
            "spin": self.mult - 1,
            "threads": self.pal,
            "max_memory": self.mem * self.pal,
        }

        print(f"PySCF configs = {self.pyscf_configs}, inherited configs = {inherited_configs}")

        self.out_fn = "pyscf.out"

        lib.num_threads(self.pal)

    @staticmethod
    def geom_from_fn(fn, **kwargs):
        geom = geom_loader(fn)
        geom.set_calculator(RigorousPySCF(**kwargs))
        return geom

    def get_driver(self, mol):
        xc = self.pyscf_configs["method_xc"]
        assert type(xc) is str and len(xc) > 0
        unrestricted = self.pyscf_configs["unrestricted"]
        assert type(unrestricted) is bool
        with_lowmem = self.pyscf_configs["with_lowmem"]
        assert type(with_lowmem) is bool
        grids_atom_grid = self.pyscf_configs["grids_atom_grid"]
        grids_level = self.pyscf_configs["grids_level"]
        nlcgrids_atom_grid = self.pyscf_configs["nlcgrids_atom_grid"]
        nlcgrids_level = self.pyscf_configs["nlcgrids_level"]
        nlc = self.pyscf_configs["method_nlc"]
        disp = self.pyscf_configs["method_disp"]
        with_df = self.pyscf_configs["with_df"]
        assert type(with_df) is bool
        auxbasis = self.pyscf_configs["auxbasis"]
        with_gpu = self.pyscf_configs["with_gpu"]
        assert type(with_gpu) is bool
        with_solvent = self.pyscf_configs["with_solvent"]
        assert type(with_solvent) is bool
        solvent_method = self.pyscf_configs["solvent_method"]
        solvent_eps = self.pyscf_configs["solvent_eps"]
        solvent_solvent = self.pyscf_configs["solvent_solvent"]
        solvent_lebedev_order = self.pyscf_configs["solvent_lebedev_order"]
        uniform_external_electric_field = self.pyscf_configs["uniform_external_electric_field"]
        basis_linear_dependency_threshold = self.pyscf_configs["basis_linear_dependency_threshold"]
        assert type(basis_linear_dependency_threshold) is float and basis_linear_dependency_threshold >= 0
        direct_scf_tol = self.pyscf_configs["direct_scf_tol"]
        assert type(direct_scf_tol) is float and direct_scf_tol >= 0
        scf_conv_tol = self.pyscf_configs["scf_conv_tol"]
        assert type(scf_conv_tol) is float and scf_conv_tol > 0
        with_tddft = self.pyscf_configs["with_tddft"]
        assert type(with_tddft) is bool

        ### Read parameters done

        assert with_lowmem is False, "Pysisyphus does not support with_lowmem yet."

        if xc.lower() == 'hf':
            mf = scf.UHF(mol) if unrestricted else scf.RHF(mol)
        else:
            mf = dft.UKS(mol, xc=xc) if unrestricted else dft.RKS(mol, xc=xc)

            if grids_level is not None:
                assert type(grids_level) is int
                mf.grids.level = grids_level
            elif grids_atom_grid is not None:
                assert type(grids_atom_grid) is not str
                mf.grids.atom_grid = grids_atom_grid
            else:
                raise ValueError('Both grids_level and grids_atom_grid are None, when DFT calculation is performed.')

            if mf._numint.libxc.is_nlc(mf.xc) or nlc:
                if nlcgrids_level is not None:
                    assert type(nlcgrids_level) is int
                    mf.nlcgrids.level = nlcgrids_level
                elif nlcgrids_atom_grid is not None:
                    assert type(nlcgrids_atom_grid) is not str
                    mf.nlcgrids.atom_grid = nlcgrids_atom_grid
                else:
                    raise ValueError('Both nlcgrids_level and nlcgrids_atom_grid are None, when DFT with NLC is specified.')

        mf.nlc = nlc
        mf.disp = disp

        if with_df:
            assert type(auxbasis) is str and len(auxbasis) > 0

            mf = mf.density_fit(auxbasis=auxbasis)

        if with_gpu:
            cp.get_default_memory_pool().free_all_blocks()
            mf = mf.to_gpu()

        if with_solvent:
            assert type(solvent_method) is str and len(solvent_method) > 0
            assert type(solvent_lebedev_order) is int and solvent_lebedev_order > 0

            solvent_method = solvent_method.upper()
            if solvent_method.endswith('PCM'):
                assert type(solvent_eps) is float and solvent_eps > 0
                mf = mf.PCM()
                mf.with_solvent.lebedev_order = solvent_lebedev_order
                mf.with_solvent.method = solvent_method.replace('PCM','-PCM')
                mf.with_solvent.eps = solvent_eps
            elif solvent_method.endswith('SMD'):
                assert type(solvent_solvent) is str and len(solvent_solvent) > 0
                mf = mf.SMD()
                mf.with_solvent.lebedev_order = solvent_lebedev_order
                mf.with_solvent.method = 'SMD'
                mf.with_solvent.solvent = solvent_solvent
            else:
                raise ValueError(f"Unrecognized solvent_method = {solvent_method}.")

        if uniform_external_electric_field:
            uniform_external_electric_field = np.array(uniform_external_electric_field, dtype = np.float64)
            assert uniform_external_electric_field.shape == (3,), f"Incorrect specification of external electric field = {uniform_external_electric_field}."

            assert with_gpu, "uniform_external_electric_field only supported if with_gpu."

            from gpu4pyscf.qmmm.external_field import add_external_field
            mf = add_external_field(mf, electric_field = uniform_external_electric_field)

        if not with_gpu:
            assert basis_linear_dependency_threshold == 0, "basis_linear_dependency_threshold only supported if with_gpu."
        else:
            import gpu4pyscf
            gpu4pyscf.scf.hf.remove_overlap_zero_eigenvalue = True
            gpu4pyscf.scf.hf.overlap_zero_eigenvalue_threshold = basis_linear_dependency_threshold

        mf.direct_scf_tol = direct_scf_tol
        mf.conv_tol = scf_conv_tol

        return mf

    def prepare_mol(self, atoms, coords, build=True):
        basis = self.pyscf_configs["basis"]
        assert type(basis) is str and len(basis) > 0
        ecp = self.pyscf_configs["ecp"]
        verbose = self.pyscf_configs["verbose"]
        assert type(verbose) is int and verbose >= 0

        ### Read parameters done

        mol = gto.Mole()
        mol.atom = [(atom, c) for atom, c in zip(atoms, coords.reshape(-1, 3))]
        mol.basis = basis
        if ecp is not None:
            assert type(ecp) is str and len(ecp) > 0
            mol.ecp = ecp
        mol.unit = "Bohr"
        mol.charge = self.charge
        mol.spin = self.mult - 1
        mol.symmetry = False
        mol.verbose = verbose
        mol.output = self.make_fn(self.out_fn)
        mol.max_memory = self.mem * self.pal
        if build:
            mol.build(parse_arg=False)
        return mol

    def generate_initial_guess(self, mf):
        # Generate the optional HOMO-LUMO mixed initial guess.

        guessmix = self.pyscf_config['guessmix']
        assert type(guessmix) is int or type(guessmix) is float
        unrestricted = self.pyscf_config['unrestricted']
        assert type(unrestricted) is bool

        ### Read parameters done

        if guessmix == 0:
            return None
        if not unrestricted:
            raise ValueError("guessmix != 0 requires unrestricted = True.")

        q = np.deg2rad(float(guessmix))
        mf_copy = mf.copy()
        mf_copy.init_guess_breaksym = 0
        mf_copy.max_cycle = 1
        mf_copy.kernel()

        mo_coeff = mf_copy.mo_coeff
        mo_occ = mf_copy.mo_occ

        if type(mo_coeff) in (tuple, list):
            assert len(mo_coeff) == 2 and mo_coeff[0].shape == mo_coeff[1].shape
        else:
            assert mo_coeff.ndim == 3 and mo_coeff.shape[0] == 2

        def _mix_one_spin(coeff, occ, angle):
            if cp is not None and isinstance(occ, cp.ndarray):
                occ = occ.get()
            occ = np.asarray(occ)
            occupied = np.where(occ > 0)[0]
            virtual = np.where(occ == 0)[0]
            if len(occupied) == 0 or len(virtual) == 0:
                raise ValueError('Unable to apply guessmix, cannot locate HOMO/LUMO from default initial guess.')
            homo_idx = occupied[-1]
            lumo_idx = virtual[0]

            c = np.cos(angle)
            s = np.sin(angle)
            mixed = coeff.copy()
            psi_homo = coeff[:, homo_idx].copy()
            psi_lumo = coeff[:, lumo_idx].copy()
            mixed[:, homo_idx] = c * psi_homo + s * psi_lumo
            mixed[:, lumo_idx] = -s * psi_homo + c * psi_lumo
            return mixed

        coeff_alpha = mo_coeff[0]
        Ca = _mix_one_spin(coeff_alpha, mo_occ[0], q)
        Cb = mo_coeff[1].copy()
        return mf_copy.make_rdm1((Ca, Cb), (mo_occ[0], mo_occ[1]))

    def prepare_input(self, atoms, coords, build=True):
        mol = self.prepare_mol(atoms, coords, build=build)
        assert mol._built, "Please call mol.build(parse_arg=False)!"
        return mol

    def store_and_track(self, results, func, atoms, coords, **prepare_kwargs):
        if self.track:
            self.store_overlap_data(atoms, coords)
            if self.track_root():
                # Redo the calculation with the updated root
                results = func(atoms, coords, **prepare_kwargs)
        results["all_energies"] = self.parse_all_energies()
        return results

    def get_energy(self, atoms, coords, **prepare_kwargs):
        assert "point_charges" not in prepare_kwargs, "point_charges is not supported in RigorousPySCF"

        mol = self.prepare_input(atoms, coords)
        mf = self.run(mol)
        results = {
            "energy": mf.e_tot,
        }
        results = self.store_and_track(
            results, self.get_energy, atoms, coords, **prepare_kwargs
        )
        return results

    def get_forces(self, atoms, coords, **prepare_kwargs):
        assert "point_charges" not in prepare_kwargs, "point_charges is not supported in RigorousPySCF"

        with_tddft = self.pyscf_configs["with_tddft"]
        assert type(with_tddft) is bool
        tddft_options_roots_for_tdgrad = self.pyscf_configs["tddft_options_roots_for_tdgrad"]
        grad_grid_response = self.pyscf_configs["grad_grid_response"]
        assert type(grad_grid_response) is bool

        ### Read parameters done

        mol = self.prepare_input(atoms, coords)
        mf = self.run(mol)
        grad_driver = mf.Gradients()
        if with_tddft:
            assert type(tddft_options_roots_for_tdgrad) is int and tddft_options_roots_for_tdgrad >= 0
            grad_driver.state = tddft_options_roots_for_tdgrad
        with_df = getattr(mf, 'with_df', None)
        if with_df:
            grad_driver.auxbasis_response = 1
        grad_driver.grid_response = grad_grid_response
        gradient = grad_driver.kernel()
        self.log("Completed gradient step")

        try:
            e_tot = mf._scf.e_tot
        except AttributeError:
            e_tot = mf.e_tot

        results = {
            "energy": e_tot,
            "forces": -gradient.flatten(),
        }
        results = self.store_and_track(
            results, self.get_forces, atoms, coords, **prepare_kwargs
        )
        return results

    def get_hessian(self, atoms, coords, **prepare_kwargs):
        assert "point_charges" not in prepare_kwargs, "point_charges is not supported in RigorousPySCF"

        hess_grid_response = self.pyscf_configs["hess_grid_response"]
        assert type(hess_grid_response) is bool

        ### Read parameters done

        mol = self.prepare_input(atoms, coords)
        mf = self.run(mol)
        hessian_driver = mf.Hessian()
        with_df = getattr(mf, 'with_df', None)
        if with_df:
            hessian_driver.auxbasis_response = 2
        hessian_driver.grid_response = hess_grid_response
        hessian = hessian_driver.kernel()

        hessian = hessian.transpose([0,2,1,3]).reshape([3 * mol.natm, 3 * mol.natm])
        results = {
            "energy": mf.e_tot,
            "hessian": hessian,
        }
        # results = self.store_and_track(
        # results, self.get_hessian, atoms, coords, **prepare_kwargs
        # )
        return results

    def run_calculation(self, atoms, coords, **prepare_kwargs):
        return self.get_energy(atoms, coords, **prepare_kwargs)

    def run(self, mol):
        mf = self.get_driver(mol=mol)

        if self.chkfile:
            new_chkfile = self.make_fn("chkfile", return_str=True)
            shutil.copy(self.chkfile, new_chkfile)
            self.chkfile = new_chkfile
            mf.chkfile = self.chkfile
            mf.init_guess = "chkfile"
            self.log(
                f"Using '{self.chkfile}' as initial guess."
            )

        keep_chk = self.pyscf_configs["keep_chk"]
        assert type(keep_chk) is bool
        if keep_chk and (self.chkfile is None):
            self.chkfile = self.make_fn("chkfile", return_str=True)
            try:
                os.remove(self.chkfile)
            except FileNotFoundError:
                self.log(f"Tried to remove '{self.chkfile}'. It doesn't exist.")
            self.log(f"Created chkfile '{self.chkfile}'")
            mf.chkfile = self.chkfile

        dm0 = self.generate_initial_guess(mf)
        mf.kernel(dm0)

        retry_soscf = self.pyscf_configs["retry_soscf"]
        assert type(retry_soscf) is bool
        if not mf.converged:
            if not retry_soscf:
                raise RuntimeError('SCF fail to converge.')
            mf = mf.newton().run()
            if not mf.converged:
                raise RuntimeError('Both the SCF and the fallback SOSCF fail to converge.')
            mf = mf.undo_newton()

        self.mf = mf.reset()  # release integrals and other temporary intermediates.
        if self.use_gpu:
            cp.get_default_memory_pool().free_all_blocks()

        self.calc_counter += 1

        return mf

    def parse_all_energies(self, exc_mf=None):
        if exc_mf is None:
            exc_mf = self.mf

        try:
            gs_energy = exc_mf._scf.e_tot
            exc_energies = exc_mf.e_tot
            all_energies = np.zeros(exc_energies.size + 1)
            all_energies[0] = gs_energy
            all_energies[1:] = exc_energies
        except AttributeError:
            gs_energy = exc_mf.e_tot
            all_energies = np.array((gs_energy,))

        return all_energies

    def prepare_overlap_data(self, path):
        raise NotImplementedError("prepare_overlap_data() is not implemented in RigorousPySCF")
        # gs_mf = self.mf._scf
        # exc_mf = self.mf

        # C = gs_mf.mo_coeff

        # first_Y = exc_mf.xy[0][1]
        # # In TDA calculations Y is just the integer 0.
        # if isinstance(first_Y, int) and (first_Y == 0):
        #     X = np.array([state[0] for state in exc_mf.xy])
        #     Y = np.zeros_like(X)
        # # In TD-DFT calculations the Y vectors is also present
        # else:
        #     # Shape = (nstates, 2 (X,Y), occ, virt)
        #     ci_coeffs = np.array(exc_mf.xy)
        #     X = ci_coeffs[:, 0]
        #     Y = ci_coeffs[:, 1]

        # all_energies = self.parse_all_energies(exc_mf)
        # return C, X, Y, all_energies

    def parse_charges(self):
        raise NotImplementedError("parse_charges() is not implemented in RigorousPySCF")
        # results = self.mf.analyze(with_meta_lowdin=False)
        # # Mulliken charges
        # charges = results[0][1]
        # return charges

    def get_chkfiles(self):
        return {
            "chkfile": self.chkfile,
        }

    def set_chkfiles(self, chkfiles):
        try:
            chkfile = chkfiles["chkfile"]
            self.chkfile = chkfile
            self.log(f"Set chkfile '{chkfile}' as chkfile.")
        except KeyError:
            self.log("Found no chkfile information in chkfiles!")

    def __str__(self):
        return f"RigorousPySCF({self.name})"
