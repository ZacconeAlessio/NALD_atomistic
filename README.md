# NALD atomistic (epoxy example)

Atomistic Non-Affine Lattice Dynamics (NALD) framework for computing viscoelastic moduli and mechanical response of cross-linked epoxy polymers and amorphous polymer systems. The methodology is general and can be applied to virtually any material (organic, inorganic) for which atomistic force-fields and simulations of the undeformed material are available. The case study of epoxy polymers is intended to be a fully worked-out case study to illustrate the method.

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
- friction kernel or damping parameter $\nu$

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

## Friction Parameter

The NALD framework includes a microscopic friction parameter $\nu$
that controls viscous damping in the generalized Langevin equation.

For the atomistic epoxy system, the value of $\nu$ is determined
by matching the NALD prediction of the viscoelastic modulus with
nonequilibrium molecular dynamics (MD) oscillatory shear calculations
in the high-frequency regime.

The high-frequency plateau and resonant peaks of the modulus are used
to calibrate the damping parameter.

Following the methodology of Vaibhav et al. [1], the friction parameter $\nu$ was estimated (for the particular epoxy polymer glass studied therein) as:

```math
\nu = 5.6 \times 10^{13} \ \mathrm{kg \ s^{-1}}
```

The precise value of the friction parameter depends on the
atomistic model, force field, temperature, and undeformed simulation protocol.

In practice, $\nu$ is calibrated by matching the NALD calculation 
of $G'$ with nonequilibrium MD
oscillatory shear calculations in the high-frequency regime where nonequilibrium MD is trustworthy.

---

## Viscoelastic Response

The NALD framework predicts the viscoelastic modulus from microscopic vibrational information using:

- Hessian eigenmodes
- affine force correlators
- non-affine particle dynamics
- friction parameter $\nu$ estimated from matching with MD data in the high-frequency plateau

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

````markdown
## Input Data Format

The NALD atomistic workflow requires atomistic configurations generated from molecular dynamics simulations together with the Hessian matrix and affine force-field data computed from the same configuration.

The repository expects the following directory structure:

```text
run1/
├── Hessian_T300.dat
├── AF_T300.dat

run2/
├── Hessian_T300.dat
├── AF_T300.dat
```

where each `runX` corresponds to an independent configuration or replica.

---

### Hessian Matrix

**File:**

```text
Hessian_T300.dat
```

The Hessian matrix is generated from an equilibrated atomistic configuration using the LAMMPS command:

```lammps
dynamical_matrix all eskm 1e-8 file hessian.data
```

The Hessian is a dense matrix of dimensions

```math
3N \times 3N
```

where \(N\) is the number of atoms.

For the epoxy system studied in the associated publication:

```text
N = 9920 atoms
```

which corresponds to a Hessian matrix of dimensions

```math
29760 \times 29760
```

The file is stored in block format. Each line contains three Hessian elements corresponding to one Cartesian degree of freedom. The Python script reconstructs the full matrix from these blocks and subsequently symmetrizes it before diagonalization.

---

### Affine Force Field

**File:**

```text
AF_T300.dat
```

The affine force field is generated by applying a very small affine shear deformation to the atomistic configuration and measuring the resulting forces.

The file is expected to contain:

```text
9 header lines
```

followed by

```text
N rows
```

containing:

```text
atom_ID    Fx    Fy    Fz
```

The script reads columns corresponding to the Cartesian affine-force components

```math
(F_x,F_y,F_z)
```

for each atom.

These values are assembled into a vector of dimension

```math
3N
```

which is subsequently projected onto the Hessian eigenvectors.

---

### Outputs

After diagonalization the code generates the following files.

#### Eigenvalues

```text
eigenvalues_T300.data
```

containing the Hessian eigenvalues

```math
\lambda_i
```

from which vibrational frequencies are obtained.

#### Affine Force Correlator

```text
gamma_T300.data
```

containing

```math
\Gamma_i
=
\left(\mathbf{e}_i \cdot \mathbf{\Xi}\right)^2
```

where:

- \(\mathbf{e}_i\) is the eigenvector of mode \(i\)
- \(\mathbf{\Xi}\) is the affine force field

#### Combined File

```text
eigen_gamma_T300.data
```

containing two columns:

```text
Eigenvalue    Gamma
```

which are subsequently used for computing:

- vibrational density of states (vDOS)
- affine force-field correlator
- storage modulus \(G'(\omega)\)
- loss modulus \(G''(\omega)\)

---

### Notes

- The Hessian matrix and affine force field must be computed from the same equilibrated configuration.
- Multiple independent configurations are strongly recommended to improve statistical averaging.
- The current implementation assumes a fixed atom count (`N = 9920`) and temperature (`T = 300 K`) in `diagonalization.py`; users should modify these values as required for their systems.
- Due to long-range Coulomb interactions in atomistic epoxy models, the Hessian matrix is dense and diagonalization may require substantial memory resources and HPC facilities.
````

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
   
 5. R. M. Elder, A. Zaccone, and T. W. Sirk,  
   *Identifying Nonaffine Softening Modes in Glassy Polymer Networks: A Pathway to Chemical Design*,  
   *ACS Macro Letters* **8**, 1160–1165 (2019).  
   DOI: https://doi.org/10.1021/acsmacrolett.9b00505

---

## License

This code is provided for academic and research purposes only.

Users are free to use and modify the code with appropriate citation to the original publications.
