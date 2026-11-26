# Momentum-resolved FD electron inverse portal at vector pole

Extends the previous finite-BW MB inversion diagnostic with an independent quantum electron occupation correction, **at the vacuum pole only** (narrow-width approximation). Empty dark final states, massless e, unchanged benchmark. This is not the full FD Breit-Wigner result, not a time-dependent collision solution.

For fixed two dark incoming energies E,E', s=M², lab pair momentum P=sqrt((E+E')²-s), beta=sqrt(1-4m²/s), cosψ=(E-E')/(P beta), where ψ is angle of dark incoming CM direction to boost. At measure-zero P=0 set arbitrary axis. The massless-electron vector spin kernel is A+B cos²θ with A=1+4m²/s, B=1-4m²/s. For outgoing electron CM cosine x to boost, azimuth average cos²θ=cos²ψ x²+(1-cos²ψ)(1-x²)/2. Lab electron energies E±=(E+E'±Px)/2. FD product divided by MB exp(-(E+E')/T) equals W=[1+exp(-E+/T)]^-1[1+exp(-E-/T)]^-1, equivalently reverse electron Pauli weights by detailed balance at zero chemical potential. Hence replace MB pole kernel by F=int_-1^1 (A+B<c²>)W dx/[2(A+B/3)], exactly 1/4<=F<=1. Resonance angular s integral becomes π/(MΓ) M²(M²+2m²) F restricted to exact partner energy range. Finite electron mass and all other visible channels omitted.

## Local results

| T GeV | pphi/T | J1/J2 MB pole | J1/J2 FD pole |
|---|---|---|---|
| 0.4 | 0.012 | 1.000175705757 | 1.000175706735 |
| 0.6882487043 | 0.001 | 1.000092398365 | 1.000092457352 |
| 0.6882487043 | 0.012 | 1.000162794623 | 1.000162861431 |
| 1 | 0.001 | 1.000139871422 | 1.000140240347 |
| 1 | 0.012 | 1.000147705659 | 1.000148065278 |
| 2 | 0.001 | 1.000091964820 | 1.000092495948 |
| 5 | 0.0001 | 1.000059796699 | 1.000063721246 |
| 10 | 0.0001 | 1.000037267049 | 1.000035181844 |

All eight selected positive MB pole inversion points remain positive under FD pole correction. This is not a claim about all momentum/temperature points or finite-width FD. 32→64 nested partner and angular Gauss orders produce ratio changes <=1e-15 on these samples. Rate check at T1,m=.02, MB NWA3.44711583713e-16 GeV^4 matches earlier independent spin-3 rate3.447115846e-16 to ~3e-9 relative. FD rate3.38456405819e-16 (~1.81% suppression); T5 rate8.53220745825e-12 vs MB1.29266918674e-11 (~34.0% suppression). This correction cannot be represented by a universally fixed factor. Positive inversion is tiny (3.5e-5 to 1.76e-4 on these samples); even modest differential finite-width or other channel corrections can matter. For instance at T=.4, alpha=.012 MB full BW ratio1.0022229216 vs pole1.0001757058, so pole approximation is inadequate there for precise inversion.

Reproduce with python -m kinetics.phi_bose_portal_fd_pole_gate. No new remote CI claim. Next: finite-width FD angular kernel, controlled pole/nonpole comparison and coupled depletion. No physical relic/BEC/CLASS.
