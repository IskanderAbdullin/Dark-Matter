# Energy-moment bound on repeated scalar multiplication, including fermion conversion and lepton elastic heating

## Scope
Initially empty homogeneous sector, fixed RD g*=10.75, TR 5–100 MeV, vacuum tree benchmark y=.001, m1=.1,m2=.02,mphi=86e-9 GeV. Include prompt or delayed 1↔2 decay/fusion, χ_i antiχ_i↔φφ, fermion pair up/down conversion via scalar exchange, number-preserving elastic redistribution in dark sector, lepton pair injection via vector, and χ_i lepton elastic scattering through same massive vector. All Bose factors on final φ can be arbitrary. Fermi occupations unpolarized, Pauli<=1. No scalar self number-changing (e.g. loop 2→4), other visible/hadronic sources, external coherent fields, additional A' processes/thermal widths, or changed entropy and FRW. The minimal Lagrangian also permits additional processes: this is a specified network, NOT a general proof for complete QFT.

## Exact number ledger and energy-weighted reaction ceilings
Only decay increases dark total particle number by +1; fusion decreases −1. Pair conversions preserve total count. Integrated net decays <= number of directly injected χ1+anti1 + two times gross reactions creating χ1 pairs from φφ or light fermions. Therefore Nφ<=Ntotal<=Ninitial_injected+N1_initial_injected+2Nrecycle; conservatively base yield <=4 times total injected pair event yield. Holds with arbitrary inverse cycles of the listed channels. No assumption of monotonic total count or thermal dark distribution.

The previous arbitrary-momentum vacuum forward bound for χ1 antiχ1→φφ is K1=y^4/(8π m2²). Detailed balance of cross sections gives reverse σv_CM=8 βχ1/βφ σv_forward_CM (spin multiplicities 4, factor2 for identical final forward). Hence reverse σv_lab<=Kinv=y^4/(π m2²)=7.9577471546e-10. Threshold E E'>=m1² since s<=4EE'. Identical initial φ event factor1/2, final Pauli<=1. Consequently Rφφ→1anti1 <= Kinv ρφ²/(2m1²), independent of arbitrarily large fφ. This uses energy density, not a momentum-shape or number-density cutoff.

Missing fermion endothermic conversion must also be counted: χ2 antiχ2→χ1 antiχ1 and same-sign χ2χ2→χ1χ1, anti analog. Distinguishable t-channel trace gives |M_t|²=y^4[((M+m)²−t)/(mphi²−t)]². Put M=m1,m=m2,B=(M+m)²−mphi², Δ=M²−m². For s>=4M², t_+=−Δ²/(s/2−M²−m²+s βiβf/2) <=−Δ²/s; βi>=sqrt(1−m²/M²). Since (1+B/z)²<=2(1+B²/z²), integrated t kernel <=2s(1+B²/Δ²). Exact σv_CM= y^4/(8π s² βi) integral dt (1+B/z)², giving uniform KD=y^4(1+B²/Δ²)/(16π M² β_min)=6.5989963646e-12. Same-sign t/u Cauchy and identical final factor give <=2 KD. Summing all three gross event rates with initial symmetry factors gives <=KD (ρχ2+ρanti2)²/M². A 1000-point exact integrated t-kernel scan from threshold to s=1e10 GeV² has sampled max2.135076e-12 below proof ceiling, not used to establish proof.

Because scalar and light fermion energies are disjoint, combined number increase rate from heavy-pair recreation <=max(Kinv,2KD) ρD²/m1², rather than summing separate maximal energies. This includes both mechanisms and all repeated cycles within specified network. No reliance on fφ being small.

## Energy injection and elastic heat ceiling
Define X=ρD/sSM. For positive kinetic pressure, dX/dt <= Q_pair/s + Q_heat/s, neglecting visible losses. With fixed vacuum matrix element MB forward source, Q_pair=T² dγ/dT exactly by Boltzmann weight. In invariant integral replace K1(u) factor by T[2K1(u)+u K0(u)]. Local six derivative checks max relative2.70e−8. FD forward and Pauli final cannot increase this positive MB energy source.

For distinct fermion lepton t-channel vector scattering, define S=s−mχ²−ml², U=u−mχ²−ml². |U|<=S, t<=0, spin average amplitude²=2C[S²+U²+2t(mχ²+ml²)]/(t−Mv²)²<=4C S²/Mv⁴. Integrate allowed t interval λ/s and multiply vM<=2: σv<=C S/(2π Mv⁴)<=2C Eχ El/(π Mv⁴). Positive transfer to χ per event <=El; sum e,mu,tau and charges, spin2 each, MB second-energy moment <=144 T⁵/π². Thus Q_heat<=κ ρD, κ=288 C T⁵/(π³ Mv⁴). Dark internal collisions do not inject energy. Gronwall: X(t)<=exp(J_R) X_pair(t), J_R=integral κ dt=63.6890223 TR³ for TR in GeV. At TR20MeV J=.0005095, i.e. <=0.051% multiplicative increase of energy envelope; this does not imply zero momentum redistribution. Bound covers charged leptons only, not hadronic heat.

## Result
Use 160 logT interpolated MB source, solve normalized cumulative X and integral of e^-z X², z=ln(TR/T), with analytic positive tail after z5. Source T<.5MeV omitted as exponentially negligible in these benchmarks; numerical interpolation not certified interval enclosure. ∫s dt=(2π²/45) sqrt(g*) Mpl TR/1.66. Multiply base + recycling number envelope by exp(2J_R), deliberately overcounting base too. Scalar Ω upper = mphi s0 Y_upper/rhoc even when not cold. No physical lower bound or actual prediction.

| TR MeV | Xpair final GeV | cumulative energy factor | Omega_phi upper no heat | J_R | upper incl lepton heat |
|---|---|---|---|---|---|
| 5 | 1.70777e-11 | 0.853727 | 2.816791e-8 | 7.9611e-6 | 2.816836e-8 |
| 8 | 4.32215e-10 | 0.803660 | 6.087488e-7 | 3.2609e-5 | 6.087885e-7 |
| 10 | 1.45645e-9 | 0.781139 | 2.073478e-6 | 6.3689e-5 | 2.073742e-6 |
| 15 | 1.00328e-8 | 0.748978 | 3.156208e-5 | 2.1495e-4 | 3.157565e-5 |
| 20 | 3.73945e-8 | 0.742201 | 4.382308e-4 | 5.0951e-4 | 4.386776e-4 |
| 25 | 1.12073e-7 | 0.750824 | 4.742767e-3 | 9.9514e-4 | 4.752216e-3 |
| 30 | 2.90290e-7 | 0.760389 | 3.830565e-2 | 1.7196e-3 | 3.843762e-2 |
| 50 | 3.96986e-6 | 0.758814 | 1.187723e+1 | 7.9611e-3 | 1.206785e+1 |
| 100 | 8.66828e-5 | 0.728741 | 1.087583e+4 | 6.3689e-2 | 1.235327e+4 |

At TR20MeV, even arbitrary stimulation, arbitrary dark momentum redistribution, repeated listed cycles and maximal lepton elastic energy gain, scalar Ωh²<=4.38678e-4, >273× below .120. At TR30 upper .03844, still below .120; above 50 upper loose and not exclusion. This is a substantially stronger and broader conditional low-TR gate than one-pass assumptions. It does NOT include scalar self number-changing or complete SM/hadronic/thermal vector processes, and it does not prove actual initial conditions. Next quantify missing injection/heat or alternative interactions; no CLASS/BEC/relic validation. Reproduce python -m kinetics.phi_bose_portal_recycling_energy_bound. Local analytic and SciPy numerical checks; remote CI not verified.
