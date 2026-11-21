# Toy colder Bose spectra at unchanged 86-eV mass and abundance

**NOT a particle-production mechanism or physical phase-transition model.** This follows up the measured suppression of linear power for a standard early thermal real 86-eV scalar using CLASS. It isolates the importance of the initial *shape and occupied phase-space amplitude*, with all scientific limitations explicit.

Take a one-real-scalar reference thermal Bose occupation f_B(x)=1/(exp(x)-1), and introduce a positive artificial coldness parameter 0<alpha<=1:

\[
\boxed{f_\alpha(p)=\alpha^{-3}f_B[p/(\alpha T_b)]},
\]

where T_b is the *reference boson temperature* before/after decoupling under free redshift. In the CLASS numerical PSD file (q=p/(alpha T_b)), include the quantum-cell factor:

\[
f_{\rm CLASS,\alpha}(q)=\frac{\alpha^{-3}}{(2\pi)^3[\exp(q)-1]},
\qquad T_{\rm ncdm}=\alpha T_b.
\]

This is **not a thermal equilibrium Bose distribution** with chemical potential; its occupancy amplitude is artificially multiplied and its momentum scale reduced.

From a change of variables p=alpha*x*T_b, the number density is independent of alpha (the factor alpha^-3 cancels the shrinking momentum volume), so the present omega_B h² is the same ~0.1191294 for the same 86-eV particle, while the early relativistic density is proportional to alpha:

\[
\boxed{n_\alpha=n_1,\quad \rho_{\rm rel,\alpha}=\alpha\,\rho_{\rm rel,1},\quad
\langle p\rangle_\alpha=\alpha\langle p\rangle_1,\quad
\Delta N_{\rm eff,\alpha}=\alpha\Delta N_{\rm eff,1}.}
\]

The nonrelativistic rms momentum is proportional to alpha, but the full cosmological free-streaming length is NOT exactly alpha times the reference value because the time of nonrelativistic transition also changes.

| alpha | Relative relic number | Relative early radiation | Early DeltaNeff |
|---:|---:|---:|---:|
| 1 | 1 | 1 | 0.064815 |
| 0.5 | 1 | 0.5 | 0.032408 |
| 0.3 | 1 | 0.3 | 0.019445 |
| 0.1 | 1 | 0.1 | 0.006482 |

The integral identities and table were independently checked numerically in Python. The [CLASS scan code](cold_psd_scan.py) runs exactly these four PSDs as noncold species in otherwise matched linear cosmologies and compares each to a cold reference. The physically important question is what values of alpha still suppress galactic structures and whether a **specified** particle Lagrangian could generate the required alpha and high occupancy without violating collisional, condensate and cosmological constraints.

Do not claim that the Sakurai–Yin or other stimulated production models actually generate this particular f_alpha. It is a useful **shape-sensitivity toy**, not a viable-model demonstration.
