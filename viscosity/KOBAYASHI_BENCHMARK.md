# Kobayashi et al. JCP 2026 benchmark

Reference: H. Kobayashi, Y. Ishii, and N. Ohtori,
"Structural origin of the strong effect of attraction on bulk viscosity in simple liquids",
J. Chem. Phys. 165, 104506 (2026).

The file `kobayashi_reference.csv` transcribes the bulk/shear viscosity,
`K_inf-K_0`, and pressure-relaxation-time values reported in Table I of the
paper and Table S1 of its Supplementary Material. The viscosity table unit is
`10^-1 mPa s`; multiply by `1e-4` to obtain Pa s.

## Why this is a useful NALD bulk benchmark

The paper writes the Green-Kubo bulk viscosity as

```text
zeta = (K_inf - K_0) * tau_zeta,
K_inf - K_0 = V <(delta P)^2> / (k_B T).
```

This gives two independent physical pieces against which a NALD bulk
calculation can be interrogated: a static bulk-response amplitude and a
relaxation time. At 124 K and rho=1124.9 kg/m^3, weakening the pair attraction
from s=1.00 to s=0.00 changes the reported quantities approximately as follows:

- bulk viscosity: 3.79 -> 0.344 in units of 10^-1 mPa s (about -91%);
- shear viscosity: 0.958 -> 0.790 (about -18%);
- K_inf-K_0: 0.673 -> 0.3908 GPa (about -42%);
- tau_zeta: 0.56 -> 0.088 ps (about -84%).

Thus the strongest attraction dependence is dynamical: the pressure-relaxation
time changes much more than the static modulus difference.

## Simulation protocol reported in the paper

The production calculations use N=1372 Ar atoms. The three state points are
(140 K, 968 kg/m^3), (124 K, 1124.9 kg/m^3), and (90 K, 1390 kg/m^3).
For N=1372 and m_Ar=39.948 u, the corresponding cubic-box lengths inferred from
the reported densities are about 4.547, 4.325, and 4.031 nm, respectively.

The pair model is the Tang-Toennies / JHBV argon potential, decomposed into
WCA-like repulsive and attractive pieces and studied with attraction scale
s=1.00, 0.95, 0.90, 0.80, 0.50, and 0.00. An Axilrod-Teller-Muto three-body
term is also studied with the full pair potential. The ATM coefficient reported
in the paper is 7.32e-108 J m^9.

After equilibration, the paper reports a 2.5e7-step NVE production trajectory,
normally with dt=6.45 fs; dt=3 fs is used for the full-pair+three-body and
purely repulsive models. The pair cutoff is 1.7 nm and the three-body cutoff is
one quarter of the cell length. The simulations use an in-house MD code.

## Important limitation for exact reproduction

The complete numerical parameter set of the Tang-Toennies/JHBV pair potential
is not printed in the paper or its Supplementary Material; the paper refers to
its Ref. 34 for those parameters. Neither atomic configurations nor the in-house
MD source code are included in the supplied files. The Data Availability
statement says that supporting data are available from the corresponding author
upon reasonable request.

Therefore `kobayashi_reference.csv` is presently a **target dataset**, not an
end-to-end reproducible input deck.

## Immediate NALD test

The first quantitative test should use the 124 K, 1124.9 kg/m^3 state point
because all attraction strengths are tabulated there. For each interaction
model/configuration:

1. construct the instantaneous mass-normalized Hessian;
2. compute both shear and volumetric affine-force fields on the same snapshot;
3. project to obtain Gamma_shear,p and Gamma_bulk,p;
4. determine the finite-size low-frequency prescription from system-size
   convergence rather than copying the old epoxy cutoff;
5. compute eta and zeta and compare with `kobayashi_reference.csv`.

The paper also reports a finite-size warning: for the U2_s1.00_plus_U3 model,
the bulk viscosity changes by about 8% when N is increased to 32000. A rigorous
benchmark should therefore not interpret the N=1372 value as a thermodynamic-
limit result.

## Structural question opened by the comparison

Kobayashi et al. attribute the attraction sensitivity of bulk viscosity to slow
pressure relaxation coupled to low-q, long-range density fluctuations, whereas
the shear channel is much less sensitive. In NALD language, a natural quantity
to compare is how attraction redistributes the modal coupling
`Gamma_bulk,p` relative to `Gamma_shear,p`, especially among the lowest
collective modes. This comparison is a scientific target; it is not assumed by
the present implementation.
