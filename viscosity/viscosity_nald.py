#!/usr/bin/env python3
"""NALD shear and bulk viscosity from precomputed mode data.

This script deliberately does not hard-code an epoxy damping coefficient or an
empirical low-frequency cutoff. It consumes Hessian eigenvalues lambda_p and
mode couplings Gamma_p=(e_p.Xi)^2 for either shear or bulk deformation.

The numerical prefactor/unit conversion depends on the units used to construct
the Hessian and affine force field. For that reason the user supplies
--prefactor. The same prefactor must reproduce the corresponding NALD loss
modulus in the chosen unit system.

For Markovian damping nu, the mode kernel used is
  M''(omega) = prefactor/V * sum_p Gamma_p * nu*omega /
               ((lambda_p-omega^2)^2 + (nu*omega)^2)
after all quantities have been converted to a mutually consistent unit system.

Then eta_s(omega)=G''(omega)/omega and zeta(omega)=K''(omega)/omega.
"""
import argparse
import numpy as np

def load_modes(path, zero_tol):
    a=np.loadtxt(path)
    if a.ndim != 2 or a.shape[1] < 2:
        raise ValueError("mode file must contain at least: eigenvalue Gamma")
    lam=a[:,0].astype(float)
    gamma=a[:,1].astype(float)
    keep=np.abs(lam) > zero_tol
    return lam[keep], gamma[keep], np.count_nonzero(~keep)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("mode_file", help="two columns: eigenvalue Gamma")
    p.add_argument("--channel", choices=("shear","bulk"), required=True)
    p.add_argument("--volume", type=float, required=True)
    p.add_argument("--nu", type=float, required=True,
                   help="Markovian damping in units consistent with eigenvalues/frequency")
    p.add_argument("--prefactor", type=float, default=1.0,
                   help="unit/mass prefactor multiplying the NALD mode sum")
    p.add_argument("--eigenvalue-scale", type=float, default=1.0,
                   help="multiply input eigenvalues by this factor before use")
    p.add_argument("--zero-tol", type=float, default=1e-10,
                   help="remove only numerical/Goldstone zero modes")
    p.add_argument("--wmin", type=float, default=1e-6)
    p.add_argument("--wmax", type=float, default=1e2)
    p.add_argument("--nfrequency", type=int, default=300)
    p.add_argument("--output", default=None)
    args=p.parse_args()

    lam,gamma,nzero=load_modes(args.mode_file,args.zero_tol)
    lam=lam*args.eigenvalue_scale
    w=np.geomspace(args.wmin,args.wmax,args.nfrequency)
    loss=np.empty_like(w)
    for j,omega in enumerate(w):
        den=(lam-omega**2)**2+(args.nu*omega)**2
        loss[j]=(args.prefactor/args.volume)*np.sum(gamma*args.nu*omega/den)

    viscosity=loss/w
    label="eta_s" if args.channel=="shear" else "zeta"
    out=args.output or f"{args.channel}_viscosity.data"
    np.savetxt(out,np.c_[w,loss,viscosity],
               header=f"omega loss_modulus {label}; removed_zero_modes={nzero}")
    print(f"Wrote {out}; removed {nzero} zero/Goldstone modes.")
    print(f"Lowest-frequency estimate: {label} = {viscosity[0]:.12g}")
    if np.any(lam < 0):
        print("NOTE: negative instantaneous-normal-mode eigenvalues were retained.")

if __name__=="__main__":
    main()
