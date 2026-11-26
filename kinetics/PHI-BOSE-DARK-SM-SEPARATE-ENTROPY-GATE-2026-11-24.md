# Independent conditional dark/SM entropy split after electron vector decoupling: a factor-two temperature ambiguity (2026-11-24)

## Underlying physical problem

The minimal real-φ dark sector consists of two distinct Dirac fermions χ1,χ2 (particle and antiparticle degrees 4 each) and one real scalar φ (g=1), initially all relativistic at a hypothetical common dark+SM temperature Tdec near 0.7GeV. The exact event-level χ1↔χ2φ and χχbar↔φφ conserve dark U1 and dark total energy; the proposed gD=.1,ε=1e−5,mV10GeV EM-mixed vector has a physical electron transport-weighted MB collision-frequency proxy νe/H≈1 at **T=0.6882487043GeV** if an externally constant g*ρ=61.75 is assumed. This is NOT a physical kinetic drag eigenvalue or a proof of equilibrium before Tdec; e,mu,hadrons, mediators and plasma can modify the crossing. At x=mχ2/T=8, TSM=2.5MeV, the same electron νe/H≈1.5565565e−7 while exact Yukawa χ2pair annihilation forward/H≈4.969. The electron vector interaction cannot on its own maintain common Tdark=TSM to that epoch. Merely increasing epsilon to make this **frequency proxy** reach H would demand ε≈**0.02534647** at x8 (factor2534.6 vs nominal 10−5), and that is NOT a phenomenologically permitted coupling claim.

## Rigorous conditional separate-entropy identity

Suppose sectors actually decouple at Tdec and thereafter each **comoving entropy** is separately conserved, with ideal thermal internally equilibrated but otherwise disconnected baths. Since s_i=(2π²/45)g*S,i(T_i)T_i³, equality of a³s_i at the decoupling event and later implies

    ξ ≡ Tdark/TSM = [
       g*S,SM(TSM) / g*S,SM(Tdec)
       × g*S,dark(Tdec) / g*S,dark(Tdark)
    ]^(1/3).

This follows algebraically without inventing a rate or future dark density. For a relativistic 1-real-scalar+two-Dirac field content,
    g*S,dark(Tdec) = 1 + 7/8*(4+4) = **8**.

Take only an explicitly declared illustrative SM entropy change
    g*S,SM(Tdec)=61.75 → g*S,SM(MeV)=10.75.
Actual smooth QCD and finite-neutrino kinetic decoupling g*S(T) need empirical tables; these two numbers are NOT a computed SM history. Two extreme dark entropy outcomes are:

| Independently postulated final dark thermal dof | Conditional ξ=Tdark/TSM |
|---|---:|
| g*S,dark remains 8 (no internal dark rest-mass entropy transfer) | **0.5583722596874** |
| g*S,dark becomes 1 (ALL initially dark fermion entropy thermalizes into real φ) | **1.1167445193749** |

The latter is physically incompatible with a **substantial residual stable χ2**, so cannot be used as this model's actual late dark temperature; it is a deliberately extreme energy/entropy example. The two cases differ by exactly factor **2**, from the cubic dark field dof ratio 8/1, not numerical noise. Both assume internal equilibrium and separate comoving entropy conservation, which may not hold near real freezeout. Actual dark scalar IR modes can be collision-inaccessible at finite P, and chemical potentials can retain entropy/nonthermal structure. Thus the real Tdark/T_SM may be different, and these values are **NOT rigorous bounds** outside this restricted thermal two-bath scenario.

## Consequence and next scientific gate

The earlier constant-g* common-SM/dark-T χ2 symmetric MB yield Ωχ2+anti2 h²≈879.96 remains a **legitimate falsifiable thermal branch**, but **cannot** be promoted to the explicit A'10GeV cosmological outcome. The full conditional evolution must maintain the SM+dark heat/entropy ledger, dark U1, nonthermal χ1/χ2/φ FD/BE and reverse reactions through QCD, electron mass, neutrino decoupling, FRW redshift and q=a p. It must compute the total Ωcharged+Ωφ, and if physical y/eps/field content are constrained, evaluate halo/CLASS likelihood only after physical f(q) is numerically converged.

Reproduction: [exact algebraic source](phi_bose_dark_sm_entropy_split_gate.py) and [four numerical algebra tests](../tests/test_phi_bose_dark_sm_entropy_split_gate.py); [independent physical electron/Yukawa rate source](phi_bose_chi2_thermal_branch_consistency.py) and [reproduction](../tests/test_phi_bose_chi2_thermal_branch_consistency.py), [machine JSON](results/phi-bose-dark-SM-entropy-split-2026-11-24.json). Electron-crossing focused CI #45454728281 succeeded; the new entropy tests are in the next focused workflow (must check its result before citing green). No actual BEC/86eV Ω/CLASS-ready or galaxy observation claim.

## Finite Dirac mass entropy correction (actual thermal χ1 and χ2, 2026-11-24)
Do not confuse the **exact massless limit** g*S,dark=8 with a mass-independent actual thermal value at T=.6882487043GeV. Evaluating the physical zero-chemical-potential FD Fermi integral for each Dirac particle+antiparticle of g=4 (mass χ1=.1GeV and χ2=.02GeV):
  g*S,FD(m/T)=45g/(4π⁴) ∫_0^∞ dy y² [sqrt(y²+(m/T)²)+y²/(3sqrt(y²+(m/T)²))] / [exp(sqrt(y²+(m/T)²))+1]
gives **g*S,dark=7.991670794968991**, 0.1042% below the massless 8. This is a useful independent finite-mass normalization and does not justify either dark endpoint scenario. New source function physical_initial_dark_gstar_s and fifth numerical regression verify independent Gauss adaptive FD entropy; latest CI result should be checked separately.
