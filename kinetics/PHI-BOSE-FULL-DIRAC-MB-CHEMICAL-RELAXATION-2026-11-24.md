# Independent classical full-Dirac Yukawa chemical relaxation: one true conserved U1 and four damping modes (2026-11-24 UTC)

This independent fixed-temperature **MAXWELL–BOLTZMANN**, *not quantum FD/BE*, benchmark computes a physically normalized vacuum-tree four-channel chemical relaxation matrix for species χ1,antiχ1,χ2,antiχ2,real φ (masses .1/.02/86e-9 GeV and y=.001). The on-shell channels are separate Dirac conjugate 1↔2 χ1→χ2+φ and antiχ1→antiχ2+φ, plus χ1antiχ1→2φ via virtual χ2 and χ2antiχ2→2φ via virtual χ1. Exact S is defined in phi_bose_yukawa_charge_audit.py.

For an ideal MB thermal manifold at common T, n_i=g_i m_i² T K2(m_i/T)/(2π²) and χ_i=∂n_i/∂μ_i=n_i/T; g=2 separately for each Dirac particle/antiparticle and g=1 for φ. The equilibrium decay rate *per charge state* is γ_dec=nχ1 Γrest K1(m1/T)/K2(m1/T); independent vectorized t/u amplitudes and properly identical-normalized Gondolo–Gelmini s integrals give γ_pair1 and γ_pair2.

The linearized four reaction sources near μ=0 satisfy R_r=−(γ_r/T)(S_r⋅δμ). Hence the symmetric physical chemical relaxation matrix L=χ^−1/2 S^T diag(γ_r/T) S χ^−1/2 is positive semidefinite and preserves exactly the vector Dirac U1=(1,−1,1,−1,0). There is just **one zero chemical eigenmode** in the connected full four-channel network; old positive sums of 1↔2 are not true global Noether charges.

At T1GeV, g*rho=61.75 (illustrative external RD bath), the rates γ per physical volume and time are:

| reaction | MB γ (GeV^4) |
|---|---:|
| χ1→χ2φ | 2.74580594711e−11 |
| antiχ1→antiχ2φ | 2.74580594711e−11 |
| χ1antiχ1→2φ | 3.49120144421e−16 |
| χ2antiχ2→2φ | 2.08863494488e−16 |

With H_RD=1.068430987449e−18GeV, the physical MB linear chemical eigenvalues (in units H) are approximately **0 (U1), 462.9526, 8979.282, 2.539711×10^8, 7.612737×10^8**. The two pair processes are order y^4 and much weaker than the y^2 decays, but still can cause substantial chemical relaxation in this prescribed MB thermal background. At T=.1GeV, g*=17.25, γ_decay=1.67723410e−13, γ_pairheavy=5.60131396e−21, γ_pairlight=7.85972209e−21 GeV4, and the smallest positive chemical eigenvalue is about 1.931×10^3 H within the same MB assumptions. The single near-zero numeric eigenvalue fluctuates in sign at machine precision compared with enormous fast modes; its exact mathematical zero is from S U1=0, not the sign of a floating eigensolver.

**CRITICAL:** λ/H values here are NOT physical chemical freezeout temperatures, proof of Bose thermalization of arbitrary f(p), or an actual 86-eV dark relic. Assumed all five populations remain MB thermally shaped at fixed T with χ=n/T, and a supplied external H,g*. This deliberately excludes FD blocking, BE enhancement and possible condensate. With no momentum-changing elastic plasma, scalar IR modes can remain inaccessible despite a finite “chemical eigenvalue” on the assumed common-μ manifold. Rates are vacuum tree, no thermal mass/quartic/SM portal, energy/entropy exchange or self-consistent FRW. The separate physically normalized FIVE-species FD/BE quantum event model and its genuine quantum susceptibility eigenvalues are computed by phi_bose_full_yukawa_chemical_relaxation.py, and should be compared by actual CI before promotion.

Reproduction: kinetics/phi_bose_full_yukawa_mb_chemical_oracle.py and tests/test_phi_bose_full_yukawa_mb_chemical_oracle.py, source spin amplitude kinetics/phi_bose_pair_independent_vectorized.py (both initial-mass orderings). Numerical values reproduced by independent scipy integral with physical K1/K2 and tanh CM angle. Run python -m pytest -q tests/test_phi_bose_full_yukawa_mb_chemical_oracle.py and python -m kinetics.phi_bose_full_yukawa_mb_chemical_oracle.

NEXT_ACTION: compare full FD/BE **same physical event** chemical linear response and observe Bose enhancement and finite-P/IR convergence; then integrate all four channels in physically DOS-normalized positive finite-step five-field solver with independently validated quantum KL, include actual in-medium χ/φ and finite SM+FRW, absolute Ω and observed likelihood.
