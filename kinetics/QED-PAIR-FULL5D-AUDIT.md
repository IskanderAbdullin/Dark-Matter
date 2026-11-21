# First full 5D finite-temperature quantum QED e−e+↔γγ reaction moments (2026-11-22)

## Claim boundary and pre-registered checks

This stage now integrates a genuine on-shell **five-dimensional**
tree-level QED electron/positron pair-annihilation and inverse
two-photon collision **thermal FD/BE moment**. Unlike the preceding
loss-only σv diagnostic, both photon Bose enhancement and inverse
γγ→e−e+ Pauli-blocked gain are included, with exact massive
electron and relativistic photon kinematics.

But this is a **restricted thermal occupation subfamily**, not yet
a callable arbitrary nonthermal C_e−(p), C_e+(p), Cγ(k), nor a
time integrator, coupled χ1–SM rate, neutrino plasma or cosmological
freeze-out. The electron/positron thermal distributions have common
Te, photons Tγ, and each sector uniform μ. No thermal radiative
corrections or real collisional thermalization rates beyond QED
tree-level 2↔2 are included.

Predeclared numerical checks:
1. With initially empty photons (bγ=0), full five-dimensional QED
   integration equals the **independently coded vacuum total σ(s)
   Møller-loss quadrature** from earlier to relative 2×10⁻⁸ at
   Te/me=.025 and .1, at nrad32/nangle24.
2. For equal Te=Tγ and μe−+μe+−2μγ=0, **net pair reaction exactly
   zero** despite nonzero forward and inverse amplitudes.
3. Positive chemical affinity gives pair annihilation and total
   entropy increase; negative affinity / hotter photon bath gives
   pair creation and positive entropy. Exact net electric charge
   and total energy; no free/reset thermal reservoir.
4. nrad24→32 changes the Te/me=.1, affinity 0.02 benchmark
   by <5×10⁻⁸; include cutoff/angles and temperature sensitivity,
   not a claim of universal grid convergence.

## On-shell CM-to-laboratory QED integral

E−=√(p²+me²), E+=√(k²+me²),
s=(E−+E+)²−|p+k|², β*=√(1−4me²/s),
and incoming relative μ=p̂·k̂. In CM the outgoing photons have
ω*₁=ω*₂=√s/2. The incoming electron's CM angle to the total
boost is cosα=[γ(p_parallel−βB E−)]/p*, where
p_parallel=p(p+kμ)/|p+k|, p*=√(s−4me²)/2.
At photon CM scattering cosθ*=c and azimuth φ,

    cosα' = cosα c + √(1−cos²α)√(1−c²)cosφ,
    ω₁=(E−+E+)/2 [1+βB cosα'],
    ω₂=E−+E+−ω₁.

This construction conserves four-momentum for every sampled
event without energy-cell center interpolation. With
D=1−β*²c² and H=(2−D)/D+
2(1−β*²)/D−2(1−β*²)²/D², the microscopic
spin/polarization-averaged angular cross section including
**identical photons 1/2** obeys

    d(σ v_M,lab)/dc = πα²/[2E− E+] H(β*,c).

The full collision reaction density (units energy⁴) is

    R_net = [4/(8π⁴)]∫dp p²∫dk k²∫dμ∫dc∫dφ/(2π)
       × πα²/[2E−E+] H × (F−B),

where F=f−(E−)f+(E+)(1+bγ(ω₁))(1+bγ(ω₂)) and
B=(1−f−)(1−f+)bγ(ω₁)bγ(ω₂). Electron/positron
spin factor 4, phase-space measures and photon identical
factor are all explicit; no extra γ polarisation factor
may be multiplied after amplitude spin/polarisation sum.

Thermal on-shell F/B =
exp{(μe−+μe+)/Te+(E−+E+)(1/Tγ−1/Te)}
for μγ=0. To avoid cancellation near equilibrium use
F−B=−F expm1[−ln(F/B)].

**Physical weak-moment ledgers at fixed volume:**

    Ndot_e−=Ndot_e+=−R_net, Ndot_γ=+2R_net,
    Qdot_e−=−∫E− dR_net, Qdot_e+=−∫E+ dR_net,
    Qdot_γ=+∫(ω₁+ω₂)dR_net,
    Qdot_total=0, Qdot_e−−Qdot_e+=0 (charge),
    Sdot=∫(F−B) ln(F/B)×[positive measure] ≥0.

These are **instantaneous integrated moments**, not state
evolution and not an assumed constant χ1 BGK relaxation.

## Independent numerical results and negative controls

Tensor Gauss nrad=32 on p/T∈[0,32], relative incoming μ and
outgoing c nangle=24, CM azimuth nφ=12, me=1 energy unit,
α=1/137.035999084. Discretization of phi is well behaved
in tested cases; some low-T configurations need radial
cutoff/angle scans before global accuracy claims.

| Te/me | Tγ/Te | (μe−+μe+)/Te | R_net/me⁴ | Qγ/me⁵ | Sdot/me⁴ |
|---|---:|---:|---:|---:|---:|
| .1 | 1 | 0 | 0 | 0 | 0 |
| .1 | 1 | +.02 | +1.51124370585e−19 | +3.51205647510e−19 | +3.02248741169e−21 |
| .2 | 1 | +.02 | +3.29059132436e−14 | +8.79183866571e−14 | +6.58118264872e−16 |
| 1 | 1 | +.04 | +6.35830164194e−8 | +3.59017590002e−7 | +2.54332065678e−9 |
| .1 | 1.1 | 0 | −5.53323540797e−17 | −1.30728742646e−16 | +1.18844311497e−16 |
| .1 | .9 | 0 | +6.90348970587e−18 | +1.60628729867e−17 | +1.78476366518e−17 |

At .1, thermal μe−=+.02, μe+=−.02 (electric
charge chemical potential) and zero photon μ, **R_net=0
exactly**, although R_forward=R_inverse≈7.48091e−18.
The equal-T +.02 affinity case satisfies
R_net/R_forward=1−e⁻^.02, and
Sdot=.02 R_net; this is a physical detailed-balance gate,
not a nonzero rate at equilibrium.

**Independent normalization:** sending Tγ→0 (empty initial
photons, no final Bose) reproduces previous 3D
electron-FD×QED σv total-loss calculation to
~1e−14 relative at Te/me=.025,.05,.1 with matching
p/T cutoff=32 and refined angles. There is no
posthoc rate symmetrization.

**Bose physics is quantitatively important:** at thermal
equilibrium, forward reaction (canceled exactly by inverse)
relative to empty-photon forward is ~1.000165545 at T/me=.1,
~1.01025133 at .2, and **1.34326026 at 1**. An
incoming-only Fermi loss model omits this photon final-state
enhancement, and reporting that forward rate as a net
chemical source is demonstrably false.

**Numerical convergence:** for Te/me=.1, μsum/Te=.02,
R_net nrad12→16→24→32→40:
1.51117968505e−19 →
1.51124354754e−19 →
1.51124370403e−19 →
1.51124370585e−19 →
1.51124370650e−19, with nangle24 and nφ12.
At nrad32, nangle10→16→24→36 shifts
1.51124370686e−19 to 1.51124370580e−19.
At nrad32, cutoff p/T=24→28→32→36 yields
1.51124053051e−19 →
1.51124361643e−19 →
1.51124370585e−19 →
1.51124370789e−19. The n32 cutoff=32 to36
shift is ~1.35e−9 relative, but cutoff=24 can
cause ~2e−6 residual. These are **one benchmark**
with thermal occupations and do not certify
uniform nonthermal collision error.

Conservation of Qγ+Qe−+Qe+ and net charge holds to
floating roundoff for all tabulated nonzero cases;
normalized total-energy residual ≤1.6e−15,
Sdot positive for both forward and inverse directions.
The exact local balance is ensured by energy-conserving
CM kinematics; no arbitrary balanced pair fitting.

## Reproduction and next required physical work

[Executable source](qed_pair_full5d.py),
[tests](../tests/test_qed_pair_full5d.py),
[compact data](results/qed-pair-full5d-2026-11-22.json).

    python -m pytest -q tests/test_qed_pair_full5d.py
    python -m kinetics.qed_pair_full5d --Te-over-me 0.1 --Tgamma-over-me 0.1 --pair-affinity-over-Te 0.02

The source uses α fixed to the vacuum fine-structure value,
tree-level QED m_e and photon dispersion; no finite-T QED
corrections, photon thermal mass, charge-neutral proton
species, weak neutrino processes, evolving photon entropy,
Friedmann H, or chemical freezeout time integration.
One must next construct genuine nonthermal electron,
positron and photon **mode-resolved strong operators**
with Pauli/Bose bounds, adapted phase-space conservative
quadrature and stiff dynamics; compare with these weak
moments, verify exact charge/total energy/entropy ledgers
and inverse-reaction relaxation. Then integrate ν and
gauge-consistent χ1–e/ν UV completion, physical χe
relaxation spectrum, χ1↔χ2+φ Bose production and
cosmological freezeout/observables.

References: [Geant4 reference QED annihilation](https://geant4.web.cern.ch/documentation/pipelines/master/prm_html/PhysicsReferenceManual/electromagnetic/electron_incident/annil.html)
and [Grohs et al. PRD 101, 063507 (2020)](https://doi.org/10.1103/PhysRevD.101.063507).
No 86-eV ΩDM, macroscopic BEC, frozen f(q), CLASS
or physical χ1 BGK γel can be claimed.
