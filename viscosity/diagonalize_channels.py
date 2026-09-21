#!/usr/bin/env python3
"""Diagonalize a LAMMPS dynamical matrix and project shear/bulk affine forces.

Unlike the historical ``diagonalization.py``, this script does not hard-code N,
temperature, run number, or a single deformation channel. By default it uses
the mass-weighted affine force required by a mass-normalized dynamical matrix,

    Xi_tilde_(i,alpha) = Xi_(i,alpha) / sqrt(m_i [amu]),

so Gamma_p has units (LAMMPS-real force)^2 / amu and can be consumed directly
by ``viscosity_nald.py --input-units lammps-real-massweighted``.

The LAMMPS affine-force dump can be either
  id mass Xi_x Xi_y Xi_z   (recommended new format), or
  mass Xi_x Xi_y Xi_z      (historical format).
"""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np


def parse_channel(spec: str) -> tuple[str, Path]:
    if "=" not in spec:
        raise argparse.ArgumentTypeError("channel must be NAME=PATH, e.g. shear=AF_shear.data")
    name, path = spec.split("=", 1)
    name = name.strip()
    if not name:
        raise argparse.ArgumentTypeError("empty channel name")
    return name, Path(path)


def read_affine(path: Path, skip_header: int) -> tuple[np.ndarray, np.ndarray]:
    a = np.genfromtxt(path, skip_header=skip_header)
    if a.ndim == 1:
        a = a[None, :]
    if a.shape[1] >= 5:
        # id, mass, Xi_x, Xi_y, Xi_z
        a = a[np.argsort(a[:, 0])]
        mass = a[:, 1]
        xi = a[:, 2:5]
    elif a.shape[1] == 4:
        # historical: mass, Xi_x, Xi_y, Xi_z
        mass = a[:, 0]
        xi = a[:, 1:4]
    else:
        raise ValueError(f"{path}: expected 4 or 5 numeric columns, found {a.shape[1]}")
    if np.any(mass <= 0):
        raise ValueError(f"{path}: all atomic masses must be positive")
    return mass.astype(float), xi.astype(float)


def read_hessian(path: Path, natoms: int) -> np.ndarray:
    dim = 3 * natoms
    h = np.empty((dim, dim), dtype=float)
    with path.open("r", encoding="utf-8") as fh:
        for row in range(dim):
            vals = np.empty(dim, dtype=float)
            k = 0
            for _ in range(natoms):
                line = fh.readline()
                if not line:
                    raise ValueError(f"{path}: ended early while reading row {row}")
                block = np.fromstring(line, sep=" ")
                if block.size != 3:
                    raise ValueError(f"{path}: expected 3 Hessian entries per line")
                vals[k:k+3] = block
                k += 3
            h[row] = vals
        extra = fh.readline()
        if extra:
            raise ValueError(f"{path}: contains extra lines after expected {dim*natoms} blocks")
    return 0.5 * (h + h.T)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hessian", type=Path, required=True)
    p.add_argument("--affine", action="append", type=parse_channel, required=True,
                   help="repeatable NAME=PATH, e.g. shear=AF_shear.data")
    p.add_argument("--skip-header", type=int, default=9)
    p.add_argument("--outdir", type=Path, default=Path("."))
    p.add_argument("--legacy-unweighted-af", action="store_true",
                   help="reproduce historical raw-force projection (not recommended for new work)")
    args = p.parse_args()

    channels: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    natoms = None
    masses_ref = None
    for name, path in args.affine:
        mass, xi = read_affine(path, args.skip_header)
        if natoms is None:
            natoms = len(mass)
            masses_ref = mass
        elif len(mass) != natoms:
            raise ValueError("all affine-force files must contain the same number of atoms")
        elif not np.allclose(mass, masses_ref, rtol=1e-12, atol=1e-12):
            raise ValueError("atomic masses differ between affine-force files")
        channels[name] = (mass, xi)

    assert natoms is not None
    print(f"Reading {3*natoms} x {3*natoms} Hessian for N={natoms} atoms...")
    hessian = read_hessian(args.hessian, natoms)
    print("Diagonalizing symmetric dynamical matrix...")
    eigenvalues, eigenvectors = np.linalg.eigh(hessian)

    args.outdir.mkdir(parents=True, exist_ok=True)
    np.savetxt(args.outdir / "eigenvalues.data", eigenvalues)

    for name, (mass, xi) in channels.items():
        if args.legacy_unweighted_af:
            xi_vec = xi.reshape(-1)
            convention = "legacy-unweighted"
        else:
            xi_vec = (xi / np.sqrt(mass)[:, None]).reshape(-1)
            convention = "massweighted"
        projection = eigenvectors.T @ xi_vec
        gamma = projection**2
        np.savetxt(args.outdir / f"gamma_{name}.data", gamma)
        np.savetxt(
            args.outdir / f"eigen_gamma_{name}.data",
            np.column_stack((eigenvalues, gamma)),
            header=f"eigenvalue Gamma_{name}; affine_projection={convention}",
        )
        print(f"Wrote channel '{name}' ({convention} affine-force projection).")


if __name__ == "__main__":
    main()
