# Validation against the historical atomistic epoxy calculation

A local T=300 K dataset from the Vaibhav--Sirk--Zaccone atomistic epoxy workflow
was used to validate the loss-modulus arithmetic. The numerical data themselves
are not stored in this repository.

The reference directory contains 30,222 modes (`3 x 10,074`), with
`V = 96003.3129555773 A^3`, and the historical files `eigenvalues.data`,
`gamma.data`, and `G_dprime_T300.data`.

Running

```bash
python validate_epoxy_reference.py /path/to/run1
```

reproduces the current `G_dp.f90` calculation including its historical choices:

- `nu = 5.0e13 s^-1` (the value actually present in `G_dp.f90` and required to
  reproduce this supplied output);
- `|lambda| < 2.57` removed;
- the literal real-unit conversion constant `2.906e6`;
- the final mode omitted, matching the Fortran loop `i=1,m-1`.

For the supplied reference files the check gives:

```text
modes removed by historical |lambda|<2.57 cutoff: 33
max frequency absolute error: 4.547474e-13
max G'' relative error: 1.018206e-08
historical zero-frequency shear viscosity: 0.033905666086 Pa s
PASS: historical epoxy loss-modulus curve reproduced.
```

This establishes backward numerical compatibility, but the historical cutoff and
raw-force projection are **not** automatically carried into the new shear/bulk
workflow. The new implementation keeps INMs, makes any finite-size cutoff
explicit, and uses mass-weighted affine forces with the mass-normalized LAMMPS
dynamical matrix.

The factor `2.906e6` in the old Fortran can also be identified: to numerical
precision it is the SI conversion `(1 kcal mol^-1 A^-1)^2 / (1 amu)`. The new
code evaluates that conversion from physical constants rather than retaining a
magic number.
