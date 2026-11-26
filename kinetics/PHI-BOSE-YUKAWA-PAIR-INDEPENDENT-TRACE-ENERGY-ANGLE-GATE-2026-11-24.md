# Independent vectorized physical Yukawa pair 2-to-2 invariant trace, MB event energy and angular-convergence gate (2026-11-24)

The existing physical χ2+antiχ2→2 real φ off-diagonal Dirac Yukawa source in phi_bose_pair_annihilation_mb.py includes correct spin averages, identical 2φ final factor and a tanh angular mapping. We independently reconstructed the spin trace as a vectorized Minkowski 2x2 interference contraction, avoiding the original nested per-event p/q four-vectors, and independently integrated the Maxwell–Boltzmann GG reaction density and its total physical event energy as minus the derivative in inverse temperature.

For sqrt(s)=.05,.1,.5,1,5,100GeV the new and old vacuum sigma agree at nangle128 within 5e-8 relative. The MB forward event density is computed with physical spin states gchi*ganti=4, identical outgoing scalar 1/2 embedded in sigma. The energy moment obeys the exact derivative identity for a fixed vacuum rate:
  Q_event=-d gamma/d(1/T) = T² d gamma/dT.
The integral energy kernel is T*[2*K1(z)+z*K0(z)], with z=sqrt(s)/T. It is a physical incoming-pair total energy weak moment, not yet separate 4-body lab-frame occupation and heat exchange.

Independently converged 128-node tanh angular integration yields:
| T (GeV) | gamma forward GeV^4 | mean total incoming energy GeV | Gamma forward per chi2/H with illustrative RD gstar |
|---:|---:|---:|---:|
| .1 | 7.85972180345e-21 | .4467886841 | 6936.017 (gstar17.25) |
| 1 | 2.08863491011e-16 | 4.4090619073 | 964.782 (gstar61.75) |
| 100 | 6.66192746689e-8 | 415.3259516 | 23.4023 (gstar106.75) |

Angular node counts 32/48/64/96/128 are saved in the independent code and tests. At T=.1 the 48-point result 7.86615035e-21 is **0.0818% above** 128-node, even though forward/reverse detailed balance of the SAME coarse event weights is exactly conserved. At T=1 the 48-point rate differs by ~0.0404%; convergence at 96 nodes agrees with 128 to <1e-7. Thus exact detailed balance/entropy is NOT evidence of angular integral convergence. The physical rates above are forward-only ideal MB, not full FD/BE quantum, thermal width, χ+SM/FRW or actual 86eV cold relic.

The earlier microscopic reciprocity and thermal-inverse calculation in phi_bose_pair_mb_detailed_balance.py derives sigma_inverse/sigma_forward=8 pchi²/pphi², with the *incoming identical real scalar pair 1/2!* cancelling spin/final factors exactly under common invariant equilibrium s weights. A single shared event weight is necessary to avoid spurious equilibrium flux when independent adaptive numerical reverse integrations do not have bitwise-identical results.

Source: kinetics/phi_bose_pair_independent_vectorized.py; independent tests: tests/test_phi_bose_pair_independent_vectorized.py; baseline source: kinetics/phi_bose_pair_annihilation_mb.py; physical reciprocal MB source: kinetics/phi_bose_pair_mb_detailed_balance.py. Run python -m pytest -q tests/test_phi_bose_pair_independent_vectorized.py and python -m kinetics.phi_bose_pair_independent_vectorized.

NEXT_ACTION: build a fully physical microreversible 2↔2 FD/BE lab-event operator preserving real Dirac U1, energy and quantum H, integrate with 1↔2 and actual thermal plasma; no BEC/Ω/CLASS-ready.
