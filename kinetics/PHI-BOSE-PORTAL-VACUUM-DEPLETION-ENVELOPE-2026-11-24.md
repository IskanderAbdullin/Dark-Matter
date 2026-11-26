# Arbitrary-momentum vacuum annihilation envelope and optimistic low-reheat depletion

Real scalar Yukawa, stable symmetric χ2 mass m=.02GeV, virtual χ1 mass M=.1GeV, y=.001, scalar μ=86e-9GeV. Exact spin-averaged squared single diagram numerator, q=p1−k1, t=q², d=m²+t−μ²:
Tr[(slash p2−m)(slash q+M)(slash p1+m)(slash q+M)]/4 = s(M²−t)/2−(d+2mM)²/2. This follows from p1·q=d/2, p2·q=−d/2, p1·p2=(s−2m²)/2. Thus single |M_t|²<=y⁴ s/[2(M²−t)]. Both t and u <=0 for physical equal incoming/equal outgoing masses: βiβf <=1−2(m²+μ²)/s from AM-GM. Cauchy |M_t+M_u|²<=2(|M_t|²+|M_u|²) yields combined <=y⁴ s[(M²−t)^−1+(M²−u)^−1]. For identical scalar final states, σv_CM=βf/(32πs) integral_-1^1 dcos |M|². Consequently exact vacuum envelope σv_CM<= y⁴/(8πM²)=3.9788735773e−12 GeV^-2. In any frame v_Moller/v_CM=s/(4E E')<=1, so same upper bound for arbitrary momentum distributions. No kinetic-equilibrium assumption, no claim saturation. Spin-independent unpolarized occupations; no stimulated final factors. If fphi<=F pointwise on accessible final modes, rate envelope grows at most by (1+F)². Other channels not included.

## Conditional injection history and deliberately maximal removal

Assume prescribed symmetric injection-only MB lepton history (same as low-reheat gate), all χ1/anti1 convert immediately into χ2/anti2; no prior χ1 annihilation, no inverse, no χ2→SM or additional number-changing removal. Constant g*=10.75 RD, no entropy change. For single-charge Y= n/s, source I(t)>=0, the actual vacuum scalar loss <=s K Y², K=bound above. Comparison ODE with loss exactly sKY² gives a lower bound on the surviving population within these assumptions. If total fixed injection Y_in all appears at reheating, minimum final Y_in/(1+Y_in K integral s dt), since earlier injection maximizes loss. Integral s dt=(2π²/45) sqrt(g*) Mpl TR/1.66. Independently use distributed injection I(z) interpolated on 160 logT nodes [.0005,.1]GeV and solve comparison ODE z=ln(TR/T); integrate to z5, then analytic pure-loss tail. Source below cutoff omitted; negligible for sampled history. This is a conditional comparison, not observed or quantum self-consistent yield. Constant g* above MeV illustrative.

| TR MeV | ΩH injection only | instantaneous lower survival | distributed lower survival | ΩH lower distributed | F required to exceed for 99% scalar mass (optimistic instantaneous bound) |
|---|---|---|---|---|---|
| 5 | 0.00327246 | 0.999937 | 0.999947 | 0.00327229 | 428320 |
| 8 | 0.0680109 | 0.997918 | 0.998349 | 0.0678986 | 74276.6 |
| 10 | 0.202833 | 0.992282 | 0.994087 | 0.201634 | 38469.0 |
| 15 | 1.05826 | 0.942618 | 0.958093 | 1.01391 | 13750.5 |
| 20 | 3.01785 | 0.812043 | 0.861265 | 2.59917 | 7051.26 |
| 25 | 6.97312 | 0.599330 | 0.684638 | 4.77406 | 4148.62 |
| 30 | 14.5121 | 0.374593 | 0.466635 | 6.77184 | 2624.83 |
| 50 | 120.998 | 0.0413211 | 0.0704023 | 8.51854 | 703.397 |
| 100 | 1482.07 | 0.00175637 | 0.00575123 | 8.52373 | 141.317 |

At TR20MeV, vacuum arbitrary-momentum optimistic distributed loss leaves >=86.1265% of injection (ΩH>=2.599 in restricted MB history). Necessary one-pass survival <=8.6869e−8 for φ>=99% by mass. If final stimulation bounded by (1+F)^2, even instantaneous injection needs F>7051 to have a chance at that survival fraction at 20MeV. A high occupation at inaccessible k does not imply this enhancement. Conversely actual stimulated rates may exceed vacuum bound, actual FD injection and other SM/number-multiplication reactions change this history. These results are not full exclusion or a validated relic.

Local exact trace checks at 35 physical angles/energies, positive envelope and analytic instant-bound comparison. Distributed injection differs from previous quadrature source by <1e−6 relative at sampled temperatures. Reproduce python -m kinetics.phi_bose_portal_vacuum_depletion_envelope. Next: cold kinematic access and number-changing, then actual multi-channel kinetic evolution. CLASS-ready=false.
