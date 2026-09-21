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


## Additional direct-viscosity script supplied with the epoxy data

A second supplied script evaluates the zero-frequency expression directly rather
than first constructing `G''(Omega)`. Its constants match the historical
epoxy implementation (`N=10074`, `nu=5.0e13 s^-1`, `tu=1.0e12`,
`conv_gamma=2.906e6`), but it uses a different soft-mode window:
`|lambda| < 1.0` instead of the `|lambda| < 2.57` window in `G_dp.f90`.

Using the same supplied T=300 K run1 data gives

```text
legacy direct formula, |lambda| < 1.0  : eta = 0.147762800487 Pa s
legacy G_dp.f90 limit, |lambda| < 2.57 : eta = 0.033905666086 Pa s
no finite cutoff (zero tolerance only) : eta ~= 7.44e3 Pa s
```

Thus the zero-frequency viscosity is extremely sensitive to the soft-mode
prescription for this finite atomistic configuration.  The two historical
cutoffs are not mathematically equivalent and there is no basis in these files
alone for choosing one over the other.

The supplied direct script also contains two mechanical issues as written:
`range(1, 1)` executes zero iterations, and the file paths depend on `samp`
rather than the loop variable `tej`.  The reference implementation
`legacy_epoxy_zero_frequency.py` keeps the supplied algebra but removes those
loop/path accidents and exposes the cutoff explicitly.

For new shear or bulk calculations the finite-size low-frequency treatment
should therefore be stated as part of the physical model and tested for
convergence, rather than inherited from either legacy script.
