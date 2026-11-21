# QED e−e+↔γγ microscopic vacuum-rate gate (2026-11-22)

## Status, hypothesis, falsification controls

This is a **scoped intermediate** microphysics result for a real
Standard-Model reaction omitted by the equilibrium massive e± heat
reference. It does **NOT** yet calculate a finite-temperature net
Boltzmann collision operator, chemical freeze-out, an evolving SM bath
or physical χ1–SM relaxation.

Hypothesis: the independently spin/polarization-averaged exact
tree-level e−e+→γγ amplitude yields the same total vacuum
annihilation σ(s) as Dirac/Heitler's electron-at-rest expression
(implemented in the Geant4 physics reference), *including the 1/2
identical final-photon factor*; its nonrelativistic σ v_M tends
to πα²/me². A positive FD-initial-state, **no-Bose loss-only**
rate-density can then be computed as a bounded *diagnostic*
with explicit finite-temperature mass and the proper Møller flux.
The quantum Bose/Pauli event ledger must obey detailed balance at
μe−+μe+−2μγ=0, even for nonzero charge-neutral μe−=−μe+.

Pre-declared gates: angular integration vs exact σv and Geant4
rest-frame expression agrees <1e−11 relative at β=.01,.1,.5,.8,.95;
the β→0 threshold within 2e−8 for β≤1e−4;
incoming loss quadrature n16→32 at T/me=.025 within 1e−8;
positive quantum entropy event at nonzero chemical affinity.
These gates test only the stated partial kernel and numerical factors.

## Independent kinematic and QED formulas

For me=m, s≥4m², β=√(1−4m²/s), c=cosθ* and
D=1−β²c², spin-averaged / final-photon-polarisation-summed

    |M|²_avg = 4e⁴ [(2−D)/D
                       +2(1−β²)/D−2(1−β²)²/D²].

After integrating the identical-photon 1/2,

    σ(s) v_M,CM = (2πα²/s)∫_{−1}¹dc H(β,c)
                 = (πα²/m²) (1−β²)/2 ×
                    {[(3−β⁴)/β] ln[(1+β)/(1−β)]
                     −2(2−β²)},

with v_M,CM=2β. Near threshold the independent Taylor
expansion of the ratio to πα²/m² is

    1+β²−(12/5)β⁴+(52/105)β⁶+(4/105)β⁸+O(β¹⁰),

avoiding catastrophic cancellation. The *independent*
Geant4/Heitler **rest-electron** total σ uses γ=s/(2m²)−1,

    σlab(γ)=πα²/[m²(γ+1)]×
       [(γ²+4γ+1)/(γ²−1) acosh(γ)
         −(γ+3)/√(γ²−1)].

The relation σlab(γ)×2β = σ(s)v_M,CM is frame-invariant
under the correct Møller flux convention. In an arbitrary lab,
v_M,lab=√[(p−·p+)²−m⁴]/(E−E+)
       =sβ/(2E−E+), hence

    σ(s)v_M,lab = [σ(s)v_M,CM]s/(4E−E+).

This converts Lorentz-invariant σ(s) to the actual fluid
collision rate without falsely using a special-relativistic
relative speed βrel=2β/(1+β²) as the Møller flux.

For incoming thermal FD distributions, with electron and
positron spin degeneracy 2 each, the **loss-only**
reaction density is

    R_loss = (4/8π⁴) ∫dp p² ∫dk k² ∫_{−1}¹dμ
               f−(p)f+(k) σ(s)v_M,lab.

The 4/(8π⁴) coefficient contains the two incoming spin
states, the pair of angular measures and the (2π)⁻⁶
normalisation; there is **no incoming 1/2** since e− and e+
are distinguishable. The photon identical-particle factor
was *already included* in σ. The computed R_loss
does **not** contain [1+bγ(k1)][1+bγ(k2)] and omits
the inverse γγ→e−e+ with Pauli blocking: it is **not
an equilibrium net collision or actual relaxation rate**.

For a single on-shell event, full quantum statistical factors
are F=f−f+(1+b1)(1+b2),
B=(1−f−)(1−f+)b1 b2. By exact shell conservation,

    F/B=exp[(μe−+μe+−2μγ)/T].

Thus even μe−≠0 is consistent with pair equilibrium when
μe+=−μe− and μγ=0; (F−B) ln(F/B)≥0.
Only this **local** check is implemented, not the
full quantum collision integral with photon phase space.

## Numerical checks (dimensionless and independently reproducible)

Gauss angular amplitude integration n160 vs analytic σv
and **independent** Geant4 Heitler σ agree to <8e−13
relative at β=.01,.1,.5,.8,.95. Threshold normalized
σv/(πα²/m²)→1 as β→0.

| T/me | R_loss/(n_e− n_e+ πα²/me²) | incoming Gauss order |
|---|---:|---|
| .025 | .99485197502049 | 32 |
| .05 | .98201716589639 | 32 |
| .1 | .94132778017644 | 32 |
| .2 | .83614414744417 | 32 |
| 1 | .29534354374317 | 32 |

This is finite-temperature **kinematic** modification of the
initial FD loss. The ratio is NOT a thermal chemical
relaxation coefficient: quantum photon stimulation, inverse,
Pauli final states and evolving e± numbers are absent. At
T/me=.025, orders 16,24,32 agree within ~1e−13 relative,
but at T/me=1 the order-16 value differs from order-32
by ~2.3e−4 relative; convergence is temperature-dependent.

**Nontrivial negative control:** setting R_net=R_loss
in common-temperature thermal equilibrium is WRONG.
The inverse reaction is exactly the same as loss once
photon Bose and electron final Pauli terms are included.
For on-shell energies (2,3→2.5,2.5) at T=1,
μe−=.2, μe+=−.2, μγ=0, F=B=0.0065937702566995;
shifting μe−=.21 yields F/B=exp(.01)
and (F−B) ln(F/B)=6.617415314×10⁻⁷.
A loss-only rate would spuriously deplete thermal pairs.

Reproducibility: [code](qed_pair_annihilation_control.py),
[regression](../tests/test_qed_pair_annihilation_control.py),
[compact data](results/qed-pair-annihilation-control-2026-11-22.json).

    python -m pytest -q tests/test_qed_pair_annihilation_control.py
    python -m kinetics.qed_pair_annihilation_control --T-over-me 0.025

References: [Geant4 Physics Reference Manual,
positron/electron annihilation](https://geant4.web.cern.ch/documentation/pipelines/master/prm_html/PhysicsReferenceManual/electromagnetic/electron_incident/annil.html)
for independent Heitler σ; [Grohs et al.,
Electron-positron annihilation freeze-out in the early universe,
Phys. Rev. D 101, 063507 (2020)](https://doi.org/10.1103/PhysRevD.101.063507).
Neither externally validates an unknown dark portal.

## Next P0 and forbidden extrapolations

Use the *differential* QED |M|² to build the actual
finite-temperature forward-minus-reverse e+e−↔γγ 5D
collision, spin/photon degeneracy and symmetry weights,
Bose/Pauli including finite moving-plasma photons,
with finite e−/e+/γ energy/charge/entropy and ν heat.
Independently confirm weak/strong moments, γγ inverse,
q/angle/time/cutoff convergence and EM chemical freezeout,
rather than inferring it from R_loss. Then match a truly
gauge-invariant χ1–SM e/ν portal and couple physical χ1↔χ2+φ
Bose production + FRW. No physical Γel, 86-eV ΩDM,
BEC, frozen f(q) or CLASS-ready result is established.
