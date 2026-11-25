# Energy-closed FRW production, 2026-11-25

## Model

Use a reference temperature Tref proportional to 1/a, x=M/Tref, q=p/Tref, theta=T_bath/Tref. The bath contains g=106.75 radiation degrees of freedom, excluding the explicitly tracked parent/daughters. Parent has kinetic-equilibrium Fermi occupation with chemical potential fixed by its comoving number. Chemical relaxation lambda=Gamma_chem/H is a parameter, not a derived microphysical reaction.

In units Tref^4, total energy R is R_bath(theta)+R_parent(x,N,theta)+R_daughters(f). FRW continuity implies dR/dln(a)=rho_parent-3P_parent in the massless-daughter limit. Friedmann H/H_ref=sqrt(R/R_bath(theta=1)). The scalar temperature equation follows by differentiating R and isolating the logarithmic temperature derivative; an independent ODE variable integrates the trace. This closes the total background energy budget without claiming that the assumed kinetic-equilibrium closure is exact. The discrete pair collision kernel is the same detailed-balanced approximation as before. Parent number and energy are distinct variables through the chemical potential.

## Reproduction

Run from kinetics/ with OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1:

- python energy_closed_bath.py --n 240 --relaxation 0 10 100 --output energy240.json
- python energy_closed_bath.py --n 400 --relaxation 10 100 --rtol 5e-7 --output energy400.json
- python energy_closed_bath.py --n 320 --kappa .1383 --relaxation 100 --keep-psd --output energy320matched.json
- python energy_closed_bath.py --n 480 --kappa .1383 --relaxation 100 --rtol 5e-7 --keep-psd --output energy480matched.json

## Result

At x0=.1, kappa=.1347, lambda=100, closing energy reduces omega from the earlier finite-bath ~.1201 to .1145. The 400-cell hot fraction p/T_bath>.1 is .45774. Retuning the physical width (not the PSD) to kappa=.1383 yields omega=.119849 (320) and .120147 (480), with hot fractions .43653/.43593. Theta at freeze-out is 1.00298. Final normalized continuity residual is -1.45e-9 at 480 cells. Daughter number mismatch -1.15e-16, Pauli ceiling .99429. Stable massless chi2 would give Delta Neff=.06977. At lambda=0 the initial parent number is conserved with each daughter number; at the original width omega=.18493, hot fraction .75193.

The hot fraction is a diagnostic cut, not a fitted cosmological hot species fraction. The remaining limitations are constant g, chemical source closure, instantaneous parent kinetic equilibrium, vacuum massless-daughter approximation and unknown parent/daughter Lagrangian. Parent abundance may have a tiny negative late-time numerical residue treated as empty. The energy residual tests the actual numerical integration rather than assuming exact machine conservation.

## Next

Export the freeze-out spectrum in q_actual=q_ref/theta_f (occupation unchanged); both number and radiation moments transform. Compare CLASS at matched abundance and independently check its intermediate-scale enhancement. Introduce finite daughter mass and constrain the allowed mass ratio and reheating initial condition.
