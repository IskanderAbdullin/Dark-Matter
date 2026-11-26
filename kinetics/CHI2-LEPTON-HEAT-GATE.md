# Missing χ2–SM collision under equal-charge real-φ gauge completion (2026-11-22)

## Physical selection gate and hypothesis

The preceding [exact gauge audit](REAL-PHI-YUKAWA-GAUGE-CHARGE-AUDIT.md)
found that the original renormalizable Yukawa term
y φ barχ2 χ1 with one **real** U(1)′-neutral Z2-odd φ
forces qχ1=qχ2 for vectorlike Dirac charge eigenstates.
Thus a nonzero χ1 vector Z′ portal also couples χ2
to electrons and the active ν sectors; a χ1-only
elastic collision operator is not a complete version
of this simplest gauge model.

The next question is whether the χ2 scattering
channel can be dismissed on **mass kinematics alone**,
even before calculating its actual abundance.
We calculate a conditional full-Pauli 5D finite-vector
weak heat for two distinct masses with equal gauge
coupling and identical *assumed thermal* FD μ=0
occupations. This tests only the mass sensitivity
of the hypothetical daughter channel; **it is
not** its real decay-produced fχ2(q), number density,
physical Γ_el, relic or cosmic heat.

## Benchmark and independent angular reference

Mass units m_e=1 (physical m_e=0.51099895 MeV),
mχ1=.1 m_e, mχ2=.02 m_e,
mφ=86 eV=0.000168297801786 m_e.
The vacuum χ1→χ2+φ threshold is open.
The φ=86-eV benchmark here is **only a mass input**
to a threshold check, not an ΩDM calculation.

Lepton baths at Tℓ=.1m_e, both fermions at
Tχ1=Tχ2=.09m_e, G=.001m_e^−2,
mV=.3m_e, lepton degrees g_e±=4 and
g_(νe+νμ+anti)=4. Exact relativistic CM
boosts and finite-t vector amplitude, Pauli
forward/reverse, independent adaptive polar
scattering angle. The neutrino channel is
massless; the electron channel retains m_e.
A common photon/neutrino bath temperature is
**assumed** here: there is no actual ν decoupling,
plasma kinetic time or charge/energy backreaction.

The source calls the same verified 5D
adaptive kernel with M=mχ1 and M=mχ2
**separately**, as required for the two gauge
charge-equal elastic species. No rate averaging
or BGK source is introduced.

## Measured conditional thermal heat

| Incoming nrad=nμ | Qχ1←2ν | Qχ2←2ν | Qχ2/Qχ1 |
|---:|---:|---:|---:|
| 6 | +1.76508422593e−18 | +2.24918101293e−18 | 1.27426 |
| 8 | +1.83292071242e−18 | +2.34110208231e−18 | 1.27725 |
| 10 | +1.82051274037e−18 | +2.32968791034e−18 | 1.27969 |

All Q use common m_e^5-based benchmark units.
Both ν thermal heat signs are positive because
the bath is warmer than the hypothetical χ
distributions. The χ2 channel is **not
suppressed** at these specific masses if it
had a comparable thermal occupancy: its heat
is about 28% larger than χ1's.

For the massive e± bath at incoming nrad=6
and8, respectively:

    n6: Qχ1(e)=+1.03644783425e−21,
        Qχ2(e)=+1.11895983257e−21,
        daughter/mother=1.07961.

    n8: Qχ1(e)=+1.02380011920e−21,
        Qχ2(e)=+1.09739978213e−21,
        daughter/mother=1.07189.

The *absolute* neutrino Q changes by several
percent under 6→8 radial refinement,
and its 8→10 shift remains ~0.5%. The
**ratio** varies ~0.4% across n6–10, but
this is not a global cutoff/angular error
bound. At this T, e± heat is strongly
Boltzmann suppressed. No physical χ2
chemical potential, density or χ1 decay
history has been established. These ratios
could be altered substantially by the
genuine nonthermal daughter distribution,
by active ν thermodynamics, mediator bounds
and time-dependent occupation.

## Consequence and next physical test

The *structural* χ2 vector collision is mandatory
for the qχ1=qχ2≠0 real-φ gauge completion.
The conditional mass scan shows no generic
kinematic argument allowing its omission.
A complete continuum quantum Boltzmann system
must therefore evolve **both** χ1 and χ2
elastic SM 2→2 operators together with
χ1↔χ2+φ Bose/Pauli decay/inverse decay
and e±/ν/γ energy, charge, entropy and
Hubble expansion. This is not done yet.
If a different UV theory suppresses
χ2 coupling, its new fields and changed
Yukawa must be written explicitly and all
rates rederived.

Executable [source](chi2_lepton_heat_gate.py),
[unit checks](../tests/test_chi2_lepton_heat_gate.py),
[compact data](results/chi2-lepton-heat-2026-11-22.json).

    python -m pytest -q tests/test_chi2_lepton_heat_gate.py
    python -m kinetics.chi2_lepton_heat_gate --M2-me 0.02 --nr 6 --nmu 6

**Terminal gates remain open.** No physical
γ_el spectrum, 86-eV ΩDM, macroscopic BEC,
frozen phase-space distribution, CLASS/P(k)
or observational fit follows from this
conditional thermal mass comparison.
