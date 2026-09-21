#!/usr/bin/env python3
"""Finite-size shear cutoff used in the 2025 JCP NALD viscosity paper.

Singh et al., J. Chem. Phys. 162, 244504 (2025), use

    omega_min = (2*pi/L) * sqrt(G_s/rho)

to discard frequencies that a finite box cannot support as propagating shear
modes.  This helper evaluates that expression either in SI units or in a
self-consistent reduced-unit system.

The paper gives this prescription for the SHEAR channel.  It does not provide
a bulk/longitudinal cutoff formula, so this helper intentionally does not invent
one.
"""

from __future__ import annotations

import argparse
import math


def shear_cutoff(box_length: float, density: float, shear_modulus: float) -> float:
    if box_length <= 0 or density <= 0 or shear_modulus < 0:
        raise ValueError("L and density must be positive and G_s must be non-negative")
    return (2.0 * math.pi / box_length) * math.sqrt(shear_modulus / density)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--box-length", type=float, required=True,
                   help="L: simulation-box size")
    p.add_argument("--density", type=float, required=True,
                   help="rho: mass density")
    p.add_argument("--shear-modulus", type=float, required=True,
                   help="G_s: zero-frequency shear modulus")
    p.add_argument("--units", choices=("si", "reduced"), default="si")
    args = p.parse_args()

    omega = shear_cutoff(args.box_length, args.density, args.shear_modulus)
    if args.units == "si":
        print(f"omega_min = {omega:.12g} s^-1")
        print(f"omega_min = {omega/1.0e12:.12g} THz (repository frequency scale)")
        print(f"lambda_cut = (omega_min/1e12)^2 = {(omega/1.0e12)**2:.12g} THz^2")
    else:
        print(f"omega_min = {omega:.12g} (reduced angular-frequency units)")
        print(f"lambda_cut = omega_min^2 = {omega**2:.12g} (reduced eigenvalue units)")


if __name__ == "__main__":
    main()
