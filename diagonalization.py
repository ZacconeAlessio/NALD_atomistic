import numpy as np
from numpy import linalg as LA

N = 5000
dim = 3 * N

for p in range(1, 3):   # run0 to run9

    print(f"Processing run {p}")

    # ---- Reset for each run ----
    hessian = np.zeros((dim, dim))
    count = 0
    shf = 0
    row = 0

    # ---- Read Hessian ----
    with open(f'./run{p}/Hessian_T0.5.dat', 'r') as f:
        for line in f:
            hessian[row, shf*3:(shf+1)*3] = np.fromstring(line, sep=' ')
            count += 1
            shf += 1

            if count % N == 0:
                row += 1
                shf = 0

    print("Matrix read. Symmetrizing...")

    # ---- Symmetrize ----
    hessian = 0.5 * (hessian + hessian.T)

    print("Diagonalizing...")

    # ---- Diagonalize (symmetric matrix) ----
    eigenvalues, eigenvectors = LA.eigh(hessian)

    print("Reading AF...")

    # ---- Read AF ----
    AF = np.zeros(dim)
    dataAF = np.genfromtxt(f'./run{p}/AF_T0.5.dat',
                           skip_header=9, usecols=(1,2,3))

    for i in range(N):
        AF[3*i + 0] = dataAF[i, 0]
        AF[3*i + 1] = dataAF[i, 1]
        AF[3*i + 2] = dataAF[i, 2]

    print("Projecting...")

    # ---- Compute gamma ----
    gamma = np.dot(eigenvectors.T, AF)**2

    # ---- Save individual files ----
    np.savetxt(f'./run{p}/eigenvalues_T0.5.data', eigenvalues)
    np.savetxt(f'./run{p}/gamma_T0.5.data', gamma)

    # ---- Save combined file (Eigen_gamm) ----
    combined = np.column_stack((eigenvalues, gamma))
    np.savetxt(f'./run{p}/eigen_gamma_T0.5.data',
               combined,
               header="Eigenvalue    Gamma",
               comments='')

    print(f"Run {p} completed\n")
