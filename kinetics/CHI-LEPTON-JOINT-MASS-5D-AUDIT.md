# Conditional finite-vector 5D Pauli weak heat across BOTH χ masses and massless/massive lepton (2026-11-22)

## Physical obligation, hypothesis and gates

For a **neutral real φ** with Yukawa χ1↔χ2+φ and equal vectorlike χ charges, both χ1 and χ2 scatter through the portal. A 5D rate for χ1 alone cannot validate the true finite-SM heat/bath, especially after the χ2 population changes. The previous independently verified energy-guided momentum and forward-log-transfer CM quadratures are now applied to **two distinct χ masses** Mχ=.1,.02 and **two lepton masses** mℓ=0,.5 in the common synthetic reference energy unit T=1. These are *conditional channels*, not an actual χ1→χ2+φ daughter abundance, thermal electron/muon or validated gauge coupling.

Prespecified numerical controls: exact on-shell Eχ+Eℓ=Eχ′+Eℓ′, fully Pauli-blocked F−R and H=(F−R)ln(F/R)≥0 for every channel, independent common-temperature FD fixed point at **distinct** μχ=.2, μℓ=−.1 including massive lepton, plus a minimum two scheme levels (coarse/mid) for 2×2 masses and both smooth/oscillatory prescribed shapes. Three levels for select mass channels expose whether apparently stable heat is a robust continuum statement. Never combine simultaneous changes in cutoff, angular and radial grids into one claimed error estimate.

## 5D physical model and numerical conventions

The incoming Gauss radial p,k+cosμ and outgoing CM log transfer variable y=log1p[(2p*²/mV²)(1−cosθ)] + azimuth are those of the earlier independent full relativistic χℓ finite-vector 2→2 weak source, with mV=.1, G=.001. The vector-tree spin trace and full forward/reverse Pauli are unchanged. The historical gχgℓ=4×4 is a *numerical pilot convention* and CANNOT be used as the correct spin/flavor/helicity multiplicity for each active SM neutrino/electron/positron channel. All Q and entropy values are in arbitrary common units; no physical Γel/H can be deduced.

For each mass pair/shape, define:
- **coarse:** cutoff12, radial E-guided segments12×local4, nμ6, mapped ncos12, nφ12;
- **mid:** cutoff12, segments12×local5, nμ8, mapped ncos16, nφ16;
- **fine** (selected cases): cutoff16, segments16×local5, nμ12, mapped ncos32, nφ24.

These levels jointly change numerical variables, so their differences are **sensitivity** not isolated a posteriori convergence bounds. Earlier exact angular I0/I1 and radial ∫p²dp proofs provide the separate analytic normalization controls.

## Measured weak Qχ across channels

All values below are Qχ in units 10⁻¹⁴ at coarse/mid; full Pauli entropy is positive, eventwise pair energy defect is at floating precision in every tested case.

| Mχ | mℓ | Smooth coarse / mid | Oscillatory coarse / mid |
|---:|---:|---:|---:|
| .1 | 0 | +9.013006 / +8.995854 | +0.470172 / +0.466254 |
| .1 | .5 | +8.156784 / +8.144078 | +0.164349 / +0.174428 |
| .02 | 0 | +9.084849 / +9.066973 | +0.509656 / +0.500993 |
| .02 | .5 | +8.199324 / +8.185587 | +0.188347 / +0.193652 |

For the smooth synthetic state at the **mid grid**, mℓ=.5 suppresses Qχ by 9.47% at Mχ=.1 and by 9.72% at Mχ=.02 relative to massless mℓ=0. This does **not** establish an electron thermal correction in nature: the e± bath, actual T, chemical species, physical couplings/decays and χ daughter distribution were not independently evolved. Mχ=.02 versus .1 changes massless smooth Q by ~0.79% at the mid test, but that is not a universal χ2/χ1 coupling ratio.

Selected third-level comparisons:
- Mχ=.1,mℓ=0 smooth: Qcoarse=9.013006456×10⁻¹⁴, Qmid=8.995854126×10⁻¹⁴, Qfine=8.997892033×10⁻¹⁴;
- Mχ=.1,mℓ=0 oscillatory: Qcoarse=4.701720499×10⁻¹⁵, Qmid=4.662544974×10⁻¹⁵, Qfine=4.668530928×10⁻¹⁵;
- Mχ=.02,mℓ=.5 oscillatory: Qcoarse=1.883469591×10⁻¹⁵, Qmid=1.936516185×10⁻¹⁵, Qfine=1.970260932×10⁻¹⁵.

**Negative:** even the refined *selected* massive-lepton oscillatory case shifts Q by +1.74% mid→fine, and simultaneous cutoff/angular/radial changes prevent attributing this to one physical effect. Consequently these rates cannot calibrate γelastic or γdrag at a claimed 1% physical accuracy. At finite lepton mass, shared T and μχ=.2, μℓ=−.1 still give Q≈0 and entropy≈0 as required by microscopic energy conservation, but chemical fixed point does not certify correct absolute normalization/precision.

## Reproducibility, physical open gates and next step

[Joint-channel code](chi_lepton_joint_mass_gate.py), [focused tests](../tests/test_chi_lepton_joint_mass_gate.py), [machine-readable conditional values](results/chi-lepton-joint-mass-5d-2026-11-22.json), [permanent dated stage](../archive/stages/2026/11/2026-11-22-chi-lepton-joint-mass.md).

    OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_chi_lepton_joint_mass_gate.py
    OPENBLAS_NUM_THREADS=1 python -m kinetics.chi_lepton_joint_mass_gate --quick

Three local tests passed. The *full* optional CLI performs eight cases at three grid levels and is heavier; do not trigger unnecessary large GitHub Actions matrices. Previously compiled independent mV=.01/.1 angular data and previous sign/cutoff negatives remain part of this audit.

**Next P0:** derive a shared-lepton **paired strong nonlinear** Cχ1,Cχ2,Cℓ source with truly evolving χ1 and χ2 separate abundances, correct physical spin/flavor and gauge charges, and maintain total Nχi/Nℓ plus energy/Pauli/H. Couple genuine χ1↔χ2+real φ Bose/Pauli, on-shell mediator and complete e±γν weak thermodynamics in FRW/H before calculating physical γelastic, 86-eV f(q), ΩDM, BEC or CLASS-ready perturbations. Current result remains an instrument/scoped negative + positive control, NOT terminal research.
