# Empty-sector portal inversion and reduced positive-gain duration (2026-11-24 UTC)

## Physical question and exact identities

In the benchmark real neutral scalar (86 eV), same-charge Dirac chi1=100 MeV and chi2=20 MeV, y=.001, vector M=10 GeV, gD=.1, epsilon=1e-5, Gamma=.005305164626 GeV, both fermions are injected by the electron portal. For one decay/fusion event the quantum factor is exactly f1(1-f2)+(f1-f2) fphi. Therefore the stimulated coefficient is proportional to the integral of f1(E)-f2(E-Ephi) over the exact finite-mass parent interval [Emin,Emax]. At common FD temperature, mu1-mu2<=Ephi implies nonpositive gain. Neither that thermal statement nor the calculation below establishes the actual nonthermal occupation.

At empty dark sector, massless MB electrons, fixed bath T, vacuum finite Breit-Wigner and per charge/spin occupation, detailed balance gives source
S_m(E)=gbar C exp(-E/T)/(192 pi^3 E p) * integral_m^infty dE' exp(-E'/T) integral_{s-}^{s+} ds s(s+2m^2)/[(s-M^2)^2+(M Gamma)^2].
Here gbar=2, C=gD^2 epsilon^2 e^2, s±=2m²+2(EE'±pp'). The s integral uses a closed rational/log/atan expression near resonance and positive Gauss quadrature away. At E=m exactly the implementation returns zero rather than the limiting source; this single measure-zero endpoint is not used as a rest-mode result. Integrating over g p²dp/(2pi²) gives gamma per species charge. Independent previous integrated rate at T=.03: 1.56395543439e-31 vs 1.56395543266e-31 GeV^4; T=.1: 2.42043315166e-27 vs 2.42043419797e-27; T=1: 3.450407733e-16 vs 3.45040772412e-16. Last two relative mismatches 4.33e-7 and 2.6e-9.

## Duration bound (reduced, not full model)

Let J_i be energy-integrated initial injection source in the exact scalar kinematic interval (shift chi2 by Ephi), and C1=int S1(E)/gamma(E) dE, gamma(E)=Gamma_rest*m1/E. Under fixed injection, spontaneous parent decay, no chi2 depletion, no inverse and no Bose backreaction, integrated f1<=min(J1*t,C1), f2>=J2*t (decay feed only increases it). Set K=(g1/gphi)Gamma_rest*m1²/(2*kstar*p_phi*E_phi), with g1=4 and gphi=1. Positive gain exponent <=K*C1²*(J1-J2)_+/(2 J1 J2). If J1<=J2 no positive gain. This assumes those reduced occupation inequalities, not full nonlinear dynamics.

Stronger: at least half the positive source weight lies at E<=2<E> (Markov). Define tau=C1/J1, delta=(J1-J2)_+/J1. Using 1-(1-exp(-x))/x = integral_0^1(1-exp(-xu))du >= x/[2(1+x)], get integrated f1<=J1*t*[1-t/(4(2*tau+t))]. When delta<1/4 positive gain cannot survive past tcut=8*tau*delta/(1-4*delta). Then exponent <= K*(J1-J2)_+*tcut²/2. Take minimum of these bounds. All finite interval integrals numerically truncate the exponentially suppressed far tail at 60/80 T; this is a numerical diagnostic, not a certified interval enclosure. No assumption that stimulation by itself respects the reduced parent inequality once inverse reactions populate the parent.

## Results

42 sampled (T,p_phi/T) points, T=.1,.2,.4,.6882487043,1,2,5,10 GeV, alpha=.0001,.0003,.001,.003,.012,.03 (only alpha>=.003 for first two temperatures). Largest sharpened exponent per T:
| T GeV | alpha | J1/J2 | exponent upper in reduced dynamics |
|---|---|---|---|
| .1 | all sampled | <=1 | 0 |
| .2 | .03 | 1.0014323742 | 1.90989e-16 |
| .4 | .012 | 1.0022229216 | 3.61554e-13 |
| .6882487043 | .001 | 1.0000882959 | 1.32670e-11 |
| 1 | .001 | 1.0001397122 | 4.10170e-9 |
| 2 | .001 | 1.0000920030 | 1.16079e-7 |
| 5 | .0001 | 1.0000597798 | 3.46364e-5 |
| 10 | .0001 | 1.0000372254 | 1.72079e-5 |

At .6882487043 and alpha=.012, ratio=1.0002046816, exponent<=1.16255e-12. At 5 GeV alpha=.0001 the weaker triangle estimate is 151.4, hence cannot on its own reject substantial amplification. Rapid loss compared with the tiny inversion sharply changes this estimate. Maximum sharpened exponent on sampled points <3.47e-5. Higher quadrature order checked selected .688 and 5 GeV examples; not a proof over a continuous domain.

## Limitations and next

FD corrections can change the tiny positive difference. All non-electron visible sources, temperature-dependent widths, inverse and dark pair channels, feedback, Hubble evolution, continuous parameter extrema and actual fphi are unresolved. In particular do NOT infer no Bose burst for full model, physical BEC, actual relic or CLASS readiness. Next calculate genuine FD electron injection, compare its sign to MB, then include other channels and depletion self-consistently. Stable-heavy population and entropy gates remain separate necessary constraints.

Reproduce: python kinetics/phi_bose_portal_initial_gain_gate.py (full 42-point scan, nested SciPy integrals); selected outputs in JSON. The source was executed locally; no new remote CI claimed.
