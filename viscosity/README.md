# Shear and bulk viscosity from atomistic NALD

This directory extends the atomistic NALD workflow to two deformation channels
using the same instantaneous Hessian eigenbasis:

- shear: `Xi_G = -df/dgamma`, `Gamma_G,p = (e_p . Xi_G)^2`,
  `eta_s = lim_(Omega->0) G''(Omega)/Omega`;
- bulk: `epsilon_v = Delta V/V`, `Xi_K = -df/d epsilon_v`,
  `Gamma_K,p = (e_p . Xi_K)^2`,
  `zeta = lim_(Omega->0) K''(Omega)/Omega`.

Defining the bulk perturbation through the **volumetric strain** `epsilon_v`
removes the common factor-of-three ambiguity associated with using the linear
strain in each box direction.

## Files

- `in.AF_shear` -- centered finite-difference affine force for xz shear.
- `in.AF_bulk` -- centered finite-difference affine force for isotropic volume strain.
- `diagonalize_channels.py` -- diagonalizes the LAMMPS mass-normalized dynamical
  matrix and projects one or more affine-force channels.
- `viscosity_nald.py` -- computes the loss modulus and its zero-frequency
  viscosity for shear or bulk.
- `validate_epoxy_reference.py` -- backward-compatibility check against the
  historical T=300 K epoxy `G_dp.f90` calculation.
- `legacy_epoxy_zero_frequency.py` -- reproduces the separately supplied direct
  epoxy zero-frequency viscosity formula while exposing its soft-mode cutoff.
- `finite_size_cutoff.py` -- evaluates the finite-size shear cutoff used in
  Singh et al., J. Chem. Phys. 162, 244504 (2025).
- `kobayashi_reference.csv` -- published 2026 simple-liquid bulk/shear benchmark values.
- `KOBAYASHI_BENCHMARK.md` -- benchmark protocol, limitations, and proposed NALD test.
- `validate_kobayashi_reference.py` -- checks the published relation
  `zeta=(K_inf-K_0) tau_zeta` against the transcribed benchmark table.
- `test_viscosity_nald.py` -- small synthetic tests of the modal formulas.
- `VALIDATION.md` -- numerical results of the epoxy reference check.

## 1. Generate the two affine-force fields

Run the shear and bulk LAMMPS inputs from the same undeformed configuration used
for the Hessian. Do not minimize or relax after applying the perturbation.
Both inputs use a centered finite difference and write atom ID, mass, and the
three components of `Xi`.

The bulk input applies isotropic scale factors

` s_+ = (1 + deltaV/2)^(1/3) ` and ` s_- = (1 - deltaV/2)^(1/3) `,

so the finite-difference variable is exactly `Delta V/V` to first order and the
plus/minus states are symmetric in volumetric strain.

## 2. Diagonalize once and project both channels

Example:

```bash
python diagonalize_channels.py \
  --hessian Hessian_T300.dat \
  --affine shear=AF_shear.data \
  --affine bulk=AF_bulk.data \
  --outdir modes_T300
```

LAMMPS `dynamical_matrix ... eskm` produces a mass-normalized dynamical matrix.
For a heterogeneous atomistic system the corresponding generalized affine force
must therefore be mass weighted. The new projection script uses
`Xi_i/sqrt(m_i)` by default before projection onto the eigenvectors. The old
repository script projected the raw affine force; use `--legacy-unweighted-af`
only when reproducing historical output.

The resulting `gamma_shear.data` and `gamma_bulk.data` have the same mode order
as `eigenvalues.data`.

## 3. Compute loss modulus and viscosity

For Markovian damping, the implemented modal response is

```text
M''(Omega) = (1/V) sum_p Gamma_p nu Omega /
             [ (lambda_p - Omega^2)^2 + (nu Omega)^2 ],
```

where `M=G` for shear and `M=K` for bulk. The direct zero-frequency limit is

```text
eta_s or zeta = (1/V) sum_p Gamma_p nu / lambda_p^2.
```

For the LAMMPS-real workflow above:

```bash
python viscosity_nald.py \
  --eigenvalues modes_T300/eigenvalues.data \
  --gamma modes_T300/gamma_shear.data \
  --channel shear --volume <A^3> --nu <s^-1> \
  --input-units lammps-real-massweighted
```

and replace `gamma_shear.data --channel shear` with
`gamma_bulk.data --channel bulk` for the bulk viscosity.

The conversion from LAMMPS real units to SI is explicit in the script; there is
no fitted numerical prefactor. The damping rate `nu` remains material,
temperature, force-field and thermostat/memory-kernel dependent and is never
silently taken from the epoxy example.

## Low-frequency modes and INMs

Negative instantaneous-normal-mode eigenvalues are retained. No empirical
low-frequency cutoff is applied by default.

A finite-size cutoff can nevertheless be physically required. Singh et al.,
J. Chem. Phys. 162, 244504 (2025), explicitly discard modes below the minimum
shear-wave frequency supported by the finite box,
`omega_min = (2 pi/L) sqrt(G_s/rho)`, where `L` is the box size, `rho` the mass
density, and `G_s` the zero-frequency shear modulus. The same paper evaluates
`G_s = G_inf - (V/k_B T)(<sigma_xy^2>-<sigma_xy>^2)`.

`finite_size_cutoff.py` evaluates this prescription and `viscosity_nald.py`
accepts the resulting value through `--cutoff-frequency`. Thus the cutoff is a
physical finite-size input, not a hidden numerical constant.

The 2025 paper gives this prescription for the **shear** channel only. It does
not specify the corresponding bulk/longitudinal cutoff. We therefore do not
hard-code a bulk analogue: that choice remains part of the bulk-viscosity
validation and should be checked against longitudinal response / system-size
scaling rather than copied from shear.

Numerical translational/Goldstone modes can separately be removed with
`--zero-tol`.

## Required validation for a new material

1. Check convergence with respect to `delta` and `deltaV`.
2. Verify that the affine-force files and Hessian come from exactly the same
   instantaneous configuration and atom ordering.
3. Establish the appropriate damping or memory kernel for the target material.
4. Document the low-frequency finite-size prescription, if one is used.
5. Average over independent instantaneous configurations.
6. Benchmark shear against Green-Kubo/NEMD where feasible; benchmark bulk
   against the corresponding pressure/volume-response calculation.

The historical epoxy data are intentionally not committed to this repository.
`validate_epoxy_reference.py` can be run on a local copy of those files to check
backward compatibility. `legacy_epoxy_zero_frequency.py` additionally makes the
cutoff sensitivity of the separately supplied direct-viscosity script explicit;
see `VALIDATION.md`.
