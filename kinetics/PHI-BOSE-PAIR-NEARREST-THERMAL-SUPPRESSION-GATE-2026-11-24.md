# Exact 86-eV near-rest χ2+antiχ2→2φ Yukawa threshold and FD/BE thermal pointwise suppression (2026-11-24)

The [exact χ2 pair global massive scalar support theorem](PHI-BOSE-PAIR-MASSIVE-SCALAR-GLOBAL-IR-SUPPORT-2026-11-24.md) proves the χ2 pair can reach a scalar exactly at 86-eV rest with incoming χ2+antiχ2 momenta ≈4.651TeV EACH, versus the original χ1→χ2φ parent threshold55.814TeV. This opens lower kinematic IR but does NOT establish a substantial real near-rest rate.

The full 2→2 invariant pair threshold can be INVERTED analytically: for any one φ of energy Eφ≤K, with 86eV=mφ≤K<mχ and equal-mass χ and antiχ incoming, total incoming LAB energy must obey

    E1+E2 >= 2mχ cosh[arcosh(mχ/mφ)-arcosh(K/mφ)].

At K=mφ this gives χ2-pair E_total≥9302.32558GeV, each parent p*=4651.16279GeV. At K=86.1eV the χ2 pair still requires E_total≥**8864.4138775GeV**, and χ1 pair requires ≥**221610.34694GeV**. This is an exact attainable invariant-mass/boost lower bound, not an artifact of a sparse 5D Gauss grid. Proof follows by maximizing the outgoing daughter boost at fixed LAB total E and s≥4mχ², using the threshold aligned equal-parent configuration.

Now independently impose a **prescribed** ideal FD χ/anti bath at temperature Tχ, both chemical potentials0, and a **prescribed** physical BE real scalar bath Tφ and subcritical μφ<mφ. For ANY on-shell event with one φ≤K, the dimensionless quantum FORWARD gain occupation factor obeys:

    F=fχ(E1)fanti(E2)(1+bφ(E3))(1+bφ(E4))
     ≤ exp[-Etotal,min/Tχ] /
        [(1-exp[-(mφ-μφ)/Tφ])
         (1-exp[-(Etotal,min-K-μφ)/Tφ])].

The reverse factor R=(1-fχ)(1-fanti)b(E3)b(E4) admits a similar upper bound. We compute their logarithms directly using exact stable expm1 to avoid catastrophic underflow; these are **pointwise multiplicative factors**, NOT bounds on the full integrated rate (which also needs the W event measure, initial phase-space volume, thermal widths and model parameters).

For χ2 pair and K=86.1eV, common Tχ=Tφ:

| Prescribed T | log10 of pointwise forward factor upper bound |
|---:|---:|
| .1GeV | **−38491.5948** |
| 1GeV | **−3842.70053** |
| 100GeV | **−29.43216** |
| 1000GeV | **+6.21574** |

The 1GeV result: parent pair FD factor ≤exp(−8864.41), about 10^−3849.766; even a low-energy thermal Bose factor ≤T/mφ≈1.16×10^7 leaves the forward upper factor below 10^−3842.70 PER EVENT. For the heavier χ1 pair the corresponding 1GeV upper factor is below roughly 10^−96237. These are stronger conditions than simply testing kinematical reach: a reaction can be on shell yet exponentially unpopulated in a given fixed T bath. With near-critical but strictly subcritical μφ=mφ−1e−16GeV and T=1GeV the bound is still below 10^−3830.

At a hypothetical early T~1TeV the Boltzmann suppression becomes far weaker and Bose factors very large, so a cold near-rest daughter might participate; BUT **a large upper bound does not imply net production**. At common thermal and chemical equilibrium the full F−R cancels exactly by detailed balance. An actual nonthermal 86eV relic needs the full quantum inverse kinetics, in-medium χ/φ masses, all 1↔2 and 2↔2 channels, finite SM heat and FRW. The low-temperature bound cannot be used to exclude high-temperature Bose bursts, actual BEC or general 86-eV DM.

Reproduce source kinetics/phi_bose_pair_nearrest_thermal_gate.py, tests tests/test_phi_bose_pair_nearrest_thermal_gate.py. Test verifies kinematic threshold via independently evaluated global pair support, the full 86.1eV suppression logs, strictly subcritical BE and the high-T caveat. This stage is a **physical on-shell kinematic + conditional thermal occupation bound**, no imaginary high-P numerical rate or observational likelihood. Keep class_ready=false.

NEXT_ACTION: retain exact pair kmin(P) energy-cell boundaries in the five-field strong solver and test actual near-critical finite-step realizability/IR convergence, while identifying physical thermal χ/φ dispersion, scalar quartic/SM portal, full FRW, frozen f(q), Ω and existing CLASS observational comparison. Don't infer condensate from a source peak or claim 1GeV tail bounds at TeV.
