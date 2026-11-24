# Low-reheat single-pass cold-particle count bound with arbitrary stimulation

## Scope
Initially empty sector, instantaneous visible reheating TR, constant g*S=g*rho=10.75 RD and conserved visible comoving entropy, e/mu/tau annihilation via benchmark finite vacuum vector. Fermion injection followed only by redshifting, single chi1->chi2 phi decay, and terminal chi_i anti-chi_i -> 2 phi (and possible removal); no kinetic reenergization/elastic momentum transfer, no inverse decays or recycling/multiplication, no hadronic/other injection. Allow arbitrary final Bose enhancement of these allowed one-pass events; count maximum possible cold scalar number rather than compute rate. Pauli blocking of injected fermions can only reduce forward source. This is NOT a full-model exclusion. 

Chosen rms cold diagnostic q_rms<=1.14255e-6 eV (match thermal 5.7keV FD rms, not actual observational criterion). Markov requires at least half the particle count below qcut=sqrt(2)*q_rms. With Tgamma0=2.34865e-4eV, gS0=3.91, aT=Tgamma0(gS0/gS)^(1/3), this is p<=alpha*T with alpha=.0096379014. At each later T<=TR the upper cold physical momentum decreases. At TR cold emission requires Echi1>=Edec_min(k) from exact massive decay, or total chi pair E>=2m cosh(acosh(m/mphi)-acosh(Ephi/mphi)) for Ephi<m (exact vacuum pair kinematics). For the terminal pair at least one incoming particle has E>=half Epair,min; m=.02GeV gives lowest threshold. At most one cold scalar per terminal pair because 2Ephi_cut<2m. A fermion cannot have gained energy since injection under the stated one-pass history (redshift and daughter energies <= mother). Thus every cold direct daughter requires a distinct injected chi1 above Edec_min(TR); every cold terminal daughter consumes at least one distinct injected fermion above Epair_min(TR)/2. Counting both channels can overcount, which is safe for upper bound. Both thresholds are minimal at TR in this cold range.

## Uniform analytic source tail envelope

For each charged lepton reverse process chi anti -> l anti-l, spin-averaged CM sigma vM=C/(6pi) s/D(s)*(1+2mchi²/s)*beta_l*(1+2ml²/s). Last product <=1; first bracket <=1.5. Lab vM<=CM value. Hence K_l=C/(4pi) max s/D=C/(4pi)*(sqrt(Mv^4+A²)+Mv²)/(2A²)=2.5927917121e-10 GeV^-2, A=MvΓ. This deliberately allows the resonant maximum even for all low-T tail particles. Detailed balance at empty dark final gives forward MB occupation source <= exp(-E/T)*K_l*nbar_MB <= exp(-E/T)*K_l*2T³/pi². Visible FD occupation <= MB; physical visible masses do not increase bound. Sum 3 charged lepton species. Integrating over one charge sector g=2 at E>=E_* and using p²dp=pE dE<=E²dE yields gamma_tail <=6 K_l T^6/pi^4 exp(-x)(x²+2x+2), x=E_*/T. Integrating through constant-g* RD, Y_tail <= COEF*TR exp(-x_R)(x_R+2), COEF=6 K_l Mpl/[pi^4(2pi²g*/45)1.66 sqrt(g*)]=7597338.405. Exact integral identity int_x^inf exp(-u)(1+2/u+2/u²)du=exp(-x)(1+2/x). This is a per charge and per species bound. Sum cold yield <= 2 Y_tail(Edec_min)+4 Y_tail(Epair_min/2). The cold mass density parameter upper is mphi s0 Ycold/rhoc; if all DM is scalar and chosen rms criterion met, it must be >=.060. No assumptions on actual stimulation amplitude or efficiency enter the count bound.

| TR MeV | cold pmax MeV | Edec_min GeV | one hard pair E_min GeV | log10 Omega_cold upper |
|---|---|---|---|---|
| 5 | 0.048190 | 49.8034 | 4.15030 | -351.01596 |
| 8 | 0.077103 | 31.1272 | 2.59396 | -131.54581 |
| 10 | 0.096379 | 24.9018 | 2.07519 | -80.94750 |
| 15 | 0.14457 | 16.6013 | 1.38350 | -31.05055 |
| 20 | 0.19276 | 12.4510 | 1.03767 | -13.64460 |
| 25 | 0.24095 | 9.96093 | 0.830177 | -5.62152 |
| 30 | 0.28914 | 8.30086 | 0.691858 | -1.28398 |
| 50 | 0.48190 | 4.98084 | 0.415269 | 4.96061 |
| 100 | 0.96379 | 2.49117 | 0.207996 | 7.56290 |

At TR20MeV, at least one hard fermion E>=1.03767GeV (51.88 TR) for a terminal pair cold product; direct decay parent E>=12.451GeV. Even the deliberately generous resonance-everywhere bound yields Omega_cold<=2.26674e-14, >2.6e12 below .060. At 30MeV upper=.0520026; the analytic upper first crosses .060 at 30.095246MeV under fixed inputs. Above this the bound is uninformative, not viability. Other earlier one-pass total number bounds can be stronger there.

Local exact antiderivative numerical checks at x=2,10,50 agree within 5e-16; physical thresholds and finite width maximum explicit. Source and JSON reproducible. Next actual redistribution, number multiplication, larger TR/alternative dark production and stable partner depletion. No CLASS/BEC/actual abundance claim.