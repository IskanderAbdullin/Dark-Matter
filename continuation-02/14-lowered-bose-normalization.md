# Lowered Bose halo: correct quantum normalization and its effect on a fiducial particle mass

Follow-up derivation, 2026-11-21. Density distribution I(psi;eta) is defined in manuscript 11.
For f(p)=((exp((E-mu)/kT)-1)^-1 - const) restricted by escape E<E_t and q=p/sqrt(2 m kT),
\[
 n=\frac{1}{(2\pi\hbar)^3}4\pi(2m k_BT)^{3/2} I(\psi;\eta)
  =\frac{\sqrt2}{\pi^2}(mV_0/\hbar)^3 I(\psi;\eta),
\]
and with lambda_dB=sqrt(2pi)hbar/(mV0) the dimensionless central phase-space density is
\[
\nu_{\rm lowered}=n_0\lambda_{\rm dB}^3=\frac{4}{\sqrt\pi}I(\psi_0;\eta).
\]
The mass relation is **exactly** the same normalization identity as in the original source,
\[
 m=\left[\frac{\sqrt{\pi/2}\,\hbar^3}
 {G R^2 V_0\,\nu_{\rm lowered}}\right]^{1/4},
\]
but nu_lowered is NOT simply Li_{3/2}(exp(psi0-eta)). Removing unbound momentum modes and subtracting their cutoff occupancy changes the normalization.

At eta=4, and using only the SOURCE'S F583-4 scale R=0.6 kpc, V0=44.1 km/s as a hypothetical common reference (NOT a refit), arithmetic gives:

| psi0 | nu_lowered | mass corresponding to those R,V0 |
|---:|---:|---:|
| 1.0 | 0.007923 | 367.34 eV |
| 2.0 | 0.066950 | 215.45 eV |
| 3.0 | 0.314096 | 146.39 eV |
| 3.8 | 1.165564 | 105.48 eV |

In the formal psi0->eta=4 central critical limit, the integral I_c is finite (near q=0 the apparent 1/q² BE singularity is canceled by the q² measure). Direct integration yields I_c=1.0872575050, nu_lowered,c=2.4536774359 and the corresponding hypothetical mass ~87.57 eV for those same R,V0. These are **not observational estimates**, since refitting new halo shapes would change R,V0. At fixed V0,R, a stronger escape cutoff requires larger particle mass for a lower central phase-space density.

Physical mass scale if R=0.6 kpc and V0=44.1 km/s is M0=V0² R/G≈2.713×10^8 solar masses. For the four equilibrium solutions in manuscript 11, the illustrative (not fitted) boundary radii would be 3.58, 5.92, 9.39, 16.42 kpc and masses 0.586, 1.648, 3.098, 5.088 billion solar masses. It is invalid to identify these numbers with F583-4's actual virial mass.

One must approach the saturated limit carefully; all four numerical BVP solutions from manuscript 11 used psi0<eta. Their center is smooth. The formal equality psi0=eta admits the same potential nonanalyticity as the ideal Bose distribution, and its full boundary-value solution remains uncomputed.
