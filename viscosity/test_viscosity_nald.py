#!/usr/bin/env python3
import sys
import unittest
import tempfile
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import viscosity_nald as vn
import diagonalize_channels as dc


class TestNALDViscosity(unittest.TestCase):
    def test_zero_frequency_single_mode(self):
        lam = np.array([4.0])          # s^-2
        gamma = np.array([2.0])        # N^2 kg^-1
        eta = vn.zero_frequency_viscosity_si(lam, gamma, 5.0, 3.0)
        self.assertAlmostEqual(eta, 0.075, places=14)

    def test_negative_inm_contributes_to_zero_limit(self):
        lam = np.array([-4.0, 4.0])
        gamma = np.array([2.0, 2.0])
        eta = vn.zero_frequency_viscosity_si(lam, gamma, 5.0, 3.0)
        self.assertAlmostEqual(eta, 0.15, places=14)

    def test_low_frequency_ratio_approaches_direct_limit(self):
        lam = np.array([-7.0, 4.0, 11.0])
        gamma = np.array([0.5, 2.0, 1.5])
        volume = 5.0
        nu = 3.0
        eta0 = vn.zero_frequency_viscosity_si(lam, gamma, volume, nu)
        omega = np.array([1.0e-7])
        loss = vn.loss_modulus_si(omega, lam, gamma, volume, nu)
        self.assertAlmostEqual(float(loss[0] / omega[0]), eta0, places=12)

    def test_explicit_cutoff_is_symmetric_in_lambda_sign(self):
        lam = np.array([-9.0, -1.0, 0.0, 1.0, 9.0])
        gamma = np.ones(5)
        kept_lam, _, info = vn.select_modes(
            lam, gamma, zero_tol=1e-12, cutoff_frequency=2.0
        )
        np.testing.assert_allclose(kept_lam, [-9.0, 9.0])
        self.assertEqual(info["removed_zero"], 1)
        self.assertEqual(info["removed_cutoff"], 2)

    def test_lammps_real_force_squared_per_amu_conversion(self):
        self.assertAlmostEqual(vn.GAMMA_REAL_MW_TO_SI, 2906915.7802365357, places=8)

    def test_affine_reader_sorts_atom_ids(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "af.data"
            path.write_text(
                "2 16.0 4.0 5.0 6.0\n"
                "1 12.0 1.0 2.0 3.0\n",
                encoding="utf-8",
            )
            mass, xi = dc.read_affine(path, skip_header=0)
            np.testing.assert_allclose(mass, [12.0, 16.0])
            np.testing.assert_allclose(xi, [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])

    def test_hessian_reader_accepts_trailing_blank_line(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "hessian.data"
            path.write_text(
                "1 0 0\n"
                "0 2 0\n"
                "0 0 3\n\n",
                encoding="utf-8",
            )
            h = dc.read_hessian(path, natoms=1)
            np.testing.assert_allclose(h, np.diag([1.0, 2.0, 3.0]))


if __name__ == "__main__":
    unittest.main()
