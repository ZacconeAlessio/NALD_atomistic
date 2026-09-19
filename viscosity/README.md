# Shear and bulk viscosity extension (development)

This directory adds the two deformation channels needed to obtain viscosity
from atomistic NALD without changing the Hessian eigenbasis.

## Definitions

For shear,
`Xi_G = -df/dgamma_xz`, `Gamma_G,p = (e_p . Xi_G)^2`, and
`eta_s = lim_(omega->0) G''(omega)/omega`.

For isotropic compression/dilation,
`epsilon_v = Delta V/V`,
`Xi_K = -df/d epsilon_v`, `Gamma_K,p = (e_p . Xi_K)^2`, and
`zeta = lim_(omega->0) K''(omega)/omega`.

The bulk LAMMPS input uses **volumetric strain**, not linear strain. This avoids
the factor-of-three ambiguity in Xi_K (and factor-of-nine ambiguity in Gamma_K).

## Files

* `in.AF_shear`: centered finite difference for the xz shear affine force.
* `in.AF_bulk`: centered finite difference for isotropic volumetric strain.
* `viscosity_nald.py`: common shear/bulk loss-modulus and viscosity postprocessor.

## Important validation before quantitative use

The old epoxy example contains a system-specific damping coefficient and an
empirical low-frequency cutoff. Neither is used automatically here.

1. Converge `delta` / `deltaV` by repeating the affine-force calculation.
2. Remove only translational/numerical zero modes using a documented tolerance.
   Genuine low-frequency and negative INMs are retained.
3. Supply a damping/memory-kernel parameter appropriate to the material,
   temperature and force field. Do not reuse the epoxy value by default.
4. Establish the SI/reduced-unit conversion represented by `--prefactor` and
   `--eigenvalue-scale` by reproducing a known NALD loss-modulus calculation.
5. Average independent instantaneous configurations for liquids.

## Next validation step

Before merging this into the main documented workflow, compare the shear output
against the existing `G_dp.f90` on the epoxy example (with identical damping,
unit conversion and mode selection). Then test bulk response on a small system
against a direct finite-frequency/Green-Kubo benchmark.
