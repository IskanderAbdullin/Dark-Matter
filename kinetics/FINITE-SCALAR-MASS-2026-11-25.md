# Finite daughter mass in energy-closed kinetics, 2026-11-25

For a scalar mass m and massless fermion daughter with dimensionless energies E and r (Tref units), let K=(x^2-m^2)/2. The invariant parent on-shell condition requires E >= (K/r + m^2*r/K)/2. The allowed rectangle measure in dE dr follows by subtracting the positive-part integrals at the upper/lower scalar energy boundaries. For fixed E, r lies between K/(E+p) and K(E+p)/m^2, p=sqrt(E^2-m^2). On that interval the primitive of E - K/(2r) - m^2*r/(2K) is E*r - (K/2)*ln(r) - m^2*r^2/(4K). The decay-angle Jacobian introduces 1/eta, eta=1-(m_phi/M)^2; kappa refers to the actual channel rest width. The two-dimensional collision matrix is now asymmetric, so daughter-fermion gain/loss uses its transpose.

The energy closure includes the scalar trace integral int f*m^2/E d^3q. Scalar energy/trace use eight-point shell quadrature; explicit scale-factor energy derivative cancels the scalar trace in the temperature equation. Number phase volume remains that of momentum shells. No post-production abundance normalization.

## Checks and interpretation

At x0=.1, kappa=.1383, lambda=100 and 320 cells, m_phi/M=0 vs 1e-9 gives omega=.11984834 vs .11984833, hot fraction .43653426 vs .43653401. This validates the massless production limit for an electroweak-scale parent (86 eV / 86 GeV = 1e-9). Ratio .001 changes omega only about .04%; ratio .01 changes omega to .12415 and hot fraction to .41253. At .03 the 160-cell run yielded omega=.15043, hot fraction=.29619; a 320-cell run timed out and this point is not converged. The last two ratios are numerical stress tests, not cosmologically consistent claims (constant g=106.75 inappropriate at keV parent temperatures).

Run from kinetics:
- OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python massive_daughter_bath.py --n 160 --mass-ratios 0 1e-9 1e-5 .001 .01 .03 --output massive160.json
- OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python massive_daughter_bath.py --n 320 --mass-ratios 0 1e-9 .001 .01 --rtol 5e-7 --output massive320resolved.json

Independent allowed-area quadrature sampled at four x values, detailed balance, directional Jacobian, daughter number and continuity passed. Remaining finite-mass tasks: finite fermion daughter, thermal dispersion masses, realistic g(T), parent chemical and kinetic collision rates. The massless scalar production approximation at high parent scale is not currently the dominant uncertainty; earlier initial production temperature and actual replenishment rate are.
