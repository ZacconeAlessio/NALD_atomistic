#!/usr/bin/env python3
"""Compute NALD shear or bulk loss modulus and zero-frequency viscosity.

The core response is evaluated in SI units,

    M''(Omega) = (1/V) sum_p Gamma_p nu Omega /
                 [ (lambda_p - Omega^2)^2 + (nu Omega)^2 ]

where M is G for shear or K for bulk. The zero-frequency viscosity is

    eta_s or zeta = (1/V) sum_p Gamma_p nu / lambda_p^2.

Two input conventions are supported:

* ``si``: lambda [s^-2], Gamma [N^2 kg^-1], V [m^3], Omega and nu [s^-1].
* ``lammps-real-massweighted``: lambda [THz^2 in the repository convention],
  Gamma [(kcal mol^-1 Angstrom^-1)^2 amu^-1], V [Angstrom^3], nu [s^-1].
  This is the convention produced by ``diagonalize_channels.py`` when the
  affine force is mass weighted before projection.

Negative instantaneous-normal-mode eigenvalues are retained. A low-frequency
cutoff is optional and explicit; no empirical cutoff is silently applied.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np

AVOGADRO = 6.02214076e23
AMU_KG = 1.66053906660e-27
KCAL_J = 4184.0
ANGSTROM_M = 1.0e-10
THZ_SI = 1.0e12
FORCE_REAL_TO_N = (KCAL_J / AVOGADRO) / ANGSTROM_M
GAMMA_REAL_MW_TO_SI = FORCE_REAL_TO_N**2 / AMU_KG
ANGSTROM3_TO_M3 = 1.0e-30


def _load_vector(path: str) -> np.ndarray:
    a = np.loadtxt(path, dtype=float)
    return np.asarray(a, dtype=float).reshape(-1)


def load_modes(eigenvalues_path: str, gamma_path: str) -> tuple[np.ndarray, np.ndarray]:
    lam = _load_vector(eigenvalues_path)
    gamma = _load_vector(gamma_path)
    if lam.size != gamma.size:
        raise ValueError(f"eigenvalues ({lam.size}) and gamma ({gamma.size}) differ in length")
    if not (np.all(np.isfinite(lam)) and np.all(np.isfinite(gamma))):
        raise ValueError("mode data contain non-finite values")
    if np.any(gamma < -1e-12):
        raise ValueError("Gamma must be non-negative because Gamma=(e.Xi)^2")
    gamma = np.maximum(gamma, 0.0)
    return lam, gamma


def convert_to_si(
    lam: np.ndarray,
    gamma: np.ndarray,
    volume: float,
    input_units: str,
) -> tuple[np.ndarray, np.ndarray, float]:
    if input_units == "si":
        return lam, gamma, volume
    if input_units == "lammps-real-massweighted":
        return (
            lam * THZ_SI**2,
            gamma * GAMMA_REAL_MW_TO_SI,
            volume * ANGSTROM3_TO_M3,
        )
    raise ValueError(f"unknown input unit convention: {input_units}")


def select_modes(
    lam_native: np.ndarray,
    gamma: np.ndarray,
    zero_tol: float,
    cutoff_frequency: float,
) -> tuple[np.ndarray, np.ndarray, dict[str, int]]:
    # Cutoffs are applied in the native input convention before SI conversion.
    abs_lam = np.abs(lam_native)
    keep_zero = abs_lam > zero_tol
    if cutoff_frequency > 0.0:
        keep_cutoff = np.sqrt(abs_lam) >= cutoff_frequency
    else:
        keep_cutoff = np.ones_like(keep_zero, dtype=bool)
    keep = keep_zero & keep_cutoff
    info = {
        "total": int(lam_native.size),
        "removed_zero": int(np.count_nonzero(~keep_zero)),
        "removed_cutoff": int(np.count_nonzero(keep_zero & ~keep_cutoff)),
        "kept": int(np.count_nonzero(keep)),
    }
    return lam_native[keep], gamma[keep], info


def loss_modulus_si(
    omega: np.ndarray,
    lam_si: np.ndarray,
    gamma_si: np.ndarray,
    volume_m3: float,
    nu_s1: float,
) -> np.ndarray:
    out = np.empty_like(omega)
    for j, om in enumerate(omega):
        den = (lam_si - om**2) ** 2 + (nu_s1 * om) ** 2
        out[j] = np.sum(gamma_si * nu_s1 * om / den) / volume_m3
    return out


def zero_frequency_viscosity_si(
    lam_si: np.ndarray,
    gamma_si: np.ndarray,
    volume_m3: float,
    nu_s1: float,
) -> float:
    return float(np.sum(gamma_si * nu_s1 / (lam_si**2)) / volume_m3)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--eigenvalues", required=True, help="one-column Hessian eigenvalue file")
    p.add_argument("--gamma", required=True, help="one-column affine-force correlator file")
    p.add_argument("--channel", choices=("shear", "bulk"), required=True)
    p.add_argument("--volume", type=float, required=True)
    p.add_argument("--nu", type=float, required=True, help="Markovian damping rate [s^-1]")
    p.add_argument(
        "--input-units",
        choices=("si", "lammps-real-massweighted"),
        default="lammps-real-massweighted",
    )
    p.add_argument(
        "--zero-tol",
        type=float,
        default=1.0e-10,
        help="remove |lambda| <= this tolerance, in native eigenvalue units",
    )
    p.add_argument(
        "--cutoff-frequency",
        type=float,
        default=0.0,
        help=("optional finite-size cutoff on sqrt(|lambda|), in native frequency units; "
              "0 keeps all nonzero modes"),
    )
    p.add_argument("--wmin", type=float, default=1.0e-6, help="minimum external frequency")
    p.add_argument("--wmax", type=float, default=1.0e2, help="maximum external frequency")
    p.add_argument("--nfrequency", type=int, default=300)
    p.add_argument("--output", default=None)
    return p


def main() -> None:
    args = build_parser().parse_args()
    if args.volume <= 0 or args.nu <= 0 or args.wmin <= 0 or args.wmax <= args.wmin:
        raise ValueError("volume, nu and frequencies must be positive, with wmax > wmin")
    if args.nfrequency < 2:
        raise ValueError("nfrequency must be at least 2")

    lam_native, gamma = load_modes(args.eigenvalues, args.gamma)
    lam_native, gamma, info = select_modes(
        lam_native, gamma, args.zero_tol, args.cutoff_frequency
    )
    lam_si, gamma_si, volume_m3 = convert_to_si(
        lam_native, gamma, args.volume, args.input_units
    )

    if args.input_units == "si":
        omega_native = np.geomspace(args.wmin, args.wmax, args.nfrequency)
        omega_si = omega_native
        frequency_label = "omega_s^-1"
        modulus_scale = 1.0
        modulus_label = "loss_modulus_Pa"
    else:
        omega_native = np.geomspace(args.wmin, args.wmax, args.nfrequency)
        omega_si = omega_native * THZ_SI
        frequency_label = "omega_THz_repo_convention"
        modulus_scale = 1.0e-9
        modulus_label = "loss_modulus_GPa"

    loss_pa = loss_modulus_si(omega_si, lam_si, gamma_si, volume_m3, args.nu)
    eta_omega = loss_pa / omega_si
    eta_zero = zero_frequency_viscosity_si(lam_si, gamma_si, volume_m3, args.nu)

    viscosity_label = "eta_s_Pa_s" if args.channel == "shear" else "zeta_Pa_s"
    out = Path(args.output or f"{args.channel}_viscosity.data")
    header = (
        f"{frequency_label} {modulus_label} {viscosity_label}\n"
        f"channel={args.channel}; input_units={args.input_units}; nu_s^-1={args.nu:.16g}; "
        f"zero_tol={args.zero_tol:.16g}; cutoff_frequency={args.cutoff_frequency:.16g}; "
        f"modes_total={info['total']}; modes_kept={info['kept']}; "
        f"removed_zero={info['removed_zero']}; removed_cutoff={info['removed_cutoff']}; "
        f"zero_frequency_viscosity_Pa_s={eta_zero:.16g}"
    )
    np.savetxt(out, np.column_stack((omega_native, loss_pa * modulus_scale, eta_omega)), header=header)

    print(f"Wrote {out}")
    print(
        f"Modes: {info['kept']}/{info['total']} kept "
        f"({info['removed_zero']} zero, {info['removed_cutoff']} below cutoff removed)"
    )
    print(f"Zero-frequency {viscosity_label}: {eta_zero:.12g}")
    if np.any(lam_native < 0):
        print("Negative instantaneous-normal-mode eigenvalues retained.")


if __name__ == "__main__":
    main()
