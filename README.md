# NALD_atomistic_epoxy

Atomistic Non-Affine Lattice Dynamics (NALD) framework for computing viscoelastic moduli and mechanical response of cross-linked epoxy polymers and amorphous polymer systems.

The methodology combines atomistic molecular dynamics simulations, Hessian matrix analysis, instantaneous normal mode (INM) analysis, affine force-field correlations, and non-affine response theory to predict viscoelastic properties across broad frequency and time scales.

---

## Scientific Background

This repository implements the Non-Affine Lattice Dynamics (NALD) framework for atomistic epoxy thermosets and amorphous polymer systems.

The methodology combines:

- atomistic molecular dynamics simulations
- Hessian matrix construction
- instantaneous normal mode (INM) analysis
- affine force-field correlations
- non-affine response theory

to compute:

- storage modulus $G'(\omega)$
- loss modulus $G''(\omega)$
- plateau elastic modulus
- viscoelastic response over broad frequency ranges

A major advantage of the NALD framework is its ability to bridge the large time-scale gap between molecular dynamics simulations (high-frequency regime) and experimental mechanical spectroscopy (low-frequency regime).

The atomistic epoxy implementation follows the methodology developed in:

V. Vaibhav, T. W. Sirk, and A. Zaccone,  
*Time-Scale Bridging in Atomistic Simulations of Epoxy Polymer Mechanics Using Nonaffine Deformation Theory*,  
*Macromolecules* (2024).  
DOI: https://doi.org/10.1021/acs.macromol.4c01360

---

## Overview

The NALD framework predicts the viscoelastic response of amorphous materials from microscopic structural information.

The approach requires:

- vibrational density of states (vDOS)
- affine force-field correlator
- friction kernel or damping parameter

Starting from equilibrated atomistic configurations, the Hessian matrix is constructed and diagonalized to obtain the vibrational modes and eigenvectors of the system. These quantities are then used within non-affine response theory to compute the complex viscoelastic modulus.

The framework naturally incorporates both:

- stable vibrational modes
- unstable instantaneous normal modes (INMs)

which are important for finite-temperature amorphous systems.

---

## Workflow

### Step 1 — Build and equilibrate the epoxy network

Prepare a fully cross-linked epoxy structure and equilibrate it at the target temperature using atomistic molecular dynamics simulations.

Typical systems include:

- DGEBA epoxy monomers
- amine curing agents
- cross-linked thermoset networks

### Step 2 — Generate the Hessian matrix

Construct the Hessian matrix from equilibrated instantaneous configurations.

For a system containing $N$ atoms, the Hessian matrix has dimensions:

```math
3N \times 3N
```

The Hessian is evaluated from the second derivatives of the potential energy.

### Step 3 — Compute affine force fields

Apply a small affine shear deformation and compute the affine force-field vectors acting on the atoms.

### Step 4 — Perform normal mode analysis

Diagonalize the Hessian matrix to obtain:

- eigenvalues
- eigenvectors
- vibrational density of states (vDOS)
- instantaneous normal modes (INMs)

At finite temperature, both positive and negative eigenvalues can appear:

- positive eigenvalues → stable vibrational modes
- negative eigenvalues → unstable instantaneous normal modes

### Step 5 — Compute affine force correlators

Calculate the affine force-field correlator:

```math
\Gamma(\omega)
```

which characterizes the coupling between vibrational modes and non-affine deformation.

### Step 6 — Compute viscoelastic modulus

The frequency-dependent complex modulus is computed within the NALD framework:

```math
G^*(\Omega) = G'(\Omega) + iG''(\Omega)
```

where:

- $G'(\Omega)$ is the storage modulus
- $G''(\Omega)$ is the loss modulus

The viscoelastic response is evaluated over a broad frequency range.

### Step 7 — Compare with MD and experiments

The NALD predictions can be validated against:

- oscillatory shear molecular dynamics simulations
- dynamic mechanical analysis (DMA) experiments

This enables direct time-scale bridging between atomistic simulations and experimentally accessible frequencies.

---

## Viscoelastic Response

The NALD framework predicts the viscoelastic modulus from microscopic vibrational information using:

- Hessian eigenmodes
- affine force correlators
- non-affine particle dynamics

The storage modulus $G'(\Omega)$ typically exhibits:

- resonant peaks associated with molecular vibrations
- low-frequency plateau behavior
- temperature dependence below the glass transition

The loss modulus $G''(\Omega)$ characterizes dissipative relaxation processes and viscous damping.

---

## Instantaneous Normal Modes (INMs)

At finite temperature, instantaneous configurations contain unstable modes corresponding to negative Hessian eigenvalues.

These modes are known as:

```math
\text{Instantaneous Normal Modes (INMs)}
```

INMs are associated with:

- internal stresses
- anharmonicity
- local instabilities in the energy landscape

Their contribution is important for accurately describing the viscoelastic response of amorphous materials and polymer glasses.

---

## Computational Notes

For atomistic epoxy systems, the Hessian matrix is dense due to the presence of long-range Coulomb interactions.

Diagonalization therefore becomes computationally demanding in both memory usage and CPU time.

Efficient diagonalization libraries such as:

- LAPACK
- Intel MKL

are strongly recommended.

Large atomistic systems may require HPC resources with substantial memory.

---

## Requirements

The following software/packages are required:

- LAMMPS
- Python 3
- NumPy
- SciPy
- LAPACK
- Intel MKL (recommended)
- Fortran compiler (`gfortran` recommended)

---

## Typical Input Files

Depending on the workflow, the calculations may require:

- equilibrated atomistic configurations
- Hessian matrices
- affine force-field data
- eigenvalue/eigenvector files
- thermodynamic data
- oscillatory shear simulation data

Users may need to modify hard-coded parameters such as:

- temperature
- damping coefficient
- deformation frequency range
- directory structure
- force-field parameters

directly within the scripts.

---

## Applications

This framework can be applied to:

- cross-linked epoxy thermosets
- polymer glasses
- amorphous polymers
- viscoelastic spectroscopy
- atomistic polymer mechanics
- temperature-dependent elastic response
- time-scale bridging between MD and experiments

---

## Developed By

**Zaccone Group**  
University of Milan  
Italy

---

## References

1. V. Vaibhav, T. W. Sirk, and A. Zaccone,  
   *Time-Scale Bridging in Atomistic Simulations of Epoxy Polymer Mechanics Using Nonaffine Deformation Theory*,  
   *Macromolecules* (2024).  
   DOI: https://doi.org/10.1021/acs.macromol.4c01360

2. V. V. Palyulin, C. Ness, R. Milkus, R. M. Elder, T. W. Sirk, and A. Zaccone,  
   *Parameter-free predictions of the viscoelastic response of glassy polymers from non-affine lattice dynamics*,  
   *Soft Matter* **14**, 8475–8482 (2018).  
   DOI: https://doi.org/10.1039/C8SM01468J

3. A. Singh, V. Vaibhav, T. W. Sirk, and A. Zaccone,  
   *Viscosity of polymer melts using non-affine theory based on vibrational modes*,  
   *The Journal of Chemical Physics* **162**, 244504 (2025).  
   DOI: https://doi.org/10.1063/5.0272171

4. A. Zaccone,  
   *General theory of the viscosity of liquids and solids from nonaffine particle motions*,  
   *Physical Review E* **108**, 044101 (2023).  
   DOI: https://doi.org/10.1103/PhysRevE.108.044101

---

## License

This code is provided for academic and research purposes only.

Users are free to use and modify the code with appropriate citation to the original publications.
