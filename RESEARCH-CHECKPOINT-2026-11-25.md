# Research checkpoint — 2026-11-25

Canonical research goal is the physical chain from microphysics to production spectrum, linear transfer, nonlinear structure and galactic equilibrium, as in [the previous theoretical roadmap](THEORETICAL-STATUS-2026-11-21.md). Do not treat a successful numerical workflow as confirmation of the 86 eV particle or the halo model.

## Evidence and negative results

- Fixed-bath direct daughter-energy solver: integrated and balanced closures converge in moments; [evidence](kinetics/DAUGHTER-ENERGY-FV-2026-11-25.md).
- FRW kinetic evolution, finite chemical parent, energy-closed model: at matched abundance x0=.1, p/T>.1 fraction about 44%; [energy closure](kinetics/ENERGY-CLOSED-BATH-2026-11-25.md).
- Finite vacuum scalar mass at m/M=1e-9 negligible in tested regime; high-ratio run timed out, not converged. [Evidence](kinetics/FINITE-SCALAR-MASS-2026-11-25.md).
- Earlier x0=.005 on two meshes gives omega=.11997, hot fraction .07847, peak occupation ~1.5e17; x0=.002 exploratory only. [Evidence](kinetics/REHEATING-MATCHED-2026-11-25.md).
- CLASS automatic q quadrature failed, transformed quadrature passed moments. Newtonian two-tracer solution failed long-wave normalization even with no fluid closure; independent synchronous/hybrid controls resolved the symptom. [Failure](cosmology/CLASS-NORMALIZATION-NEGATIVE-2026-11-25.md), [control](cosmology/CLASS-NORMALIZATION-RECOVERED-2026-11-25.md).
- [Verified corrected fixed-bath CLASS run 45625224041](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/45625224041): half-power at z=5 ~.6984 Mpc^-1, ratio at k=10 ~.0371, mesh/partition drift below .41%. No observation likelihood. [Retained numerical summary](cosmology/results/kinetic-psd-class-synchronous-45625224041.json).
- [Two-mass geometry](kinetics/TWO-MASS-GEOMETRY-2026-11-25.md) verified reference integration but not evolution.
- Corrected energy-closed [run 45625224293](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/45625224293) succeeded; half-power k=.68988 Mpc^-1 at z=5 on 480 cells, with mesh/partition drift <.5% before half. Full hierarchy changes power at k=1–1.5 substantially: [corrected results and limitation](cosmology/ENERGY-CLOSED-CLASS-CORRECTED-2026-11-25.md).
- Historical failed earlier-onset [run 45625464842](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/45625464842) completed both kinetic and CLASS calculations but failed an incorrect final assertion requiring a half-crossing. Both reported crossings are null on .003–20 Mpc^-1; at z=5 the 480-cell P(k=10)/P_CDM=.597707, 320-cell .597024. Logs are retained; full transfer JSON was not written. Correct the right-censored guard and rerun; no final full-grid convergence claim yet.

## Next actions (in order)

1. Collect the rerun of earlier-onset CLASS with corrected right-censoring guard and test full-range mesh agreement. Collect targeted hierarchy precision run. Refine the perturbation hierarchy around the half crossing and at high k (the new fluid control shows 16–26% differences at k=1–1.5). If early spectrum has no crossing in tested range, extend the physical k range rather than assuming a result.
2. Implement finite fermion daughter in the conservative asymmetric energy collision kernel using the tested geometry; extend total stress-energy to both traces. Derive temperature-dependent dispersion/chemical rates for a physically specified Lagrangian, with the early burst's tiny IR scale as the decision criterion.
3. Only after conditional production/transfer checks, compare admissible scenarios to CMB/Lyman-alpha and dwarf-galaxy likelihoods and close the nonlinear halo/relaxation link.

A new research program is not necessary. Continue this repository, branch main, preserve previous evidence. Computer use is optional and restricted to explicitly authorized project location; no personal machine resources required for this checkpoint.

A [targeted full hierarchy convergence workflow](.github/workflows/class-hierarchy-precision.yml) is prepared to compare lmax=30/45 and 120/240 momentum nodes at k=.003–1.5 Mpc^-1, without the expensive and currently unnecessary high-k run. It saves raw comparison results even if the convergence threshold fails.


## Follow-up: completed workflows

[Verified follow-up](cosmology/KINETIC-CLASS-FOLLOWUP-2026-11-25.md): early x0=.005 has no half-power crossing up to 20 Mpc^-1 on two meshes (full-grid drift .251%); P/P_CDM=.597708 at k=10,z=5, conditional default fluid closure. Late x0=.1 full hierarchy yields k_half,P=.717581 Mpc^-1 (z=5), stable to .054% between two precision settings up to k=1.5. The default-fluid half power .68988 is superseded by this more accurate local calculation. Next: finite fermion mass/degeneracy and its physical cosmological fate; early full-hierarchy validation remains open.

## Two-mass local continuation

[Two-mass kinetic solver and scan](kinetics/TWO-MASS-BATH-2026-11-25.md): mass degeneracy sharply suppresses the hot component but makes stable daughter matter density untenable; .99 ratio integration failed. Next implement energy-closed chi2 bath relaxation with independent number-source accounting.

## Daughter bath interaction

[Energy-closed detailed-balanced chi2 relaxation](kinetics/DAUGHTER-BATH-2026-11-25.md): three validated exploratory cases return the massive fermion to the bath. Kappa=.05 overproduces the scalar. Match abundance by changing the channel width; verify independent meshes and physical rates. GitHub reference update encountered a server-side error; local commits in the authorized workspace preserve the work.


## Degeneracy + bath return + Yukawa-rate candidate

See kinetics/YUKAWA-RATE-AND-DEGENERACY-2026-11-25.md and kinetics/results/yukawa-degeneracy-2026-11-25.json (actual local outputs). For M=86 TeV, m2/M=0.9, kappa=0.00224, kappaSM1=4000, kappa2=0.03, two independent meshes n=240,320 give Omega_phi h^2=0.120113,0.120224 and hot cut p/T>0.1 fractions 0.00351646,0.00351514; scalar number source, Pauli and energy checks pass. Three-scale analytic PSD export passes; CLASS transfer workflow added, outcome pending. The parent is initially thermal and assumed kinetic-equilibrium, SM rates approximate, and the illustrative uncancelled seesaw contribution is too large. Next: empty-parent initialization, weaker parent couplings and actual CLASS result; no claim of microscopic viability.

Publication limitation: OAuth credential rejects workflow changes (missing workflow scope). The executable workflow is saved as a proposed configuration in cosmology/workflows/daughter-bath-class.proposed.yml; it is not installed or running. CLASS results remain uncomputed. Original local workflow commit c2b5790 preserved on checkpoint/unpublished-workflow-20261125.

## Empty-parent and thermal screen

See kinetics/EMPTY-PARENT-THERMAL-SCREEN-2026-11-25.md. Initial empty x0=.01 with kappa_SM1=47 gives Omega=.08354 at previous width; stronger parent couplings overproduce. Equilibrium simple-Yukawa thermal-mass/peak ratios reach 0.12-0.44 early vs 0.00354 late. Need genuine thermal dispersion and momentum-dependent parent production; do not interpret early results as predictions. Existing regression suite 16/16 after explicit halo initial integration step. CLASSpp 26 Windows wheel installed in isolated project Python; older 23.2 PSD-reader failure and new 26 alias incompatibility identified. Testing the corrected backend.

## Weak-coupling fit and actual local CLASSpp

See kinetics/WEAK-PARENT-FIT-2026-11-25.md, results/weak-parent-fit-2026-11-25.json. Empty-parent x0=.01, kappa_SM1=47, physical kappa=.0043950631; n=240/320 Omega=.121077/.121459, hot=.00363545/.00362122. Thermal benchmark mass/peak ~.26, not controlled. Bath-only parent-shape test TV 3.87%, energy mismatch 4.43% at same number for weak source. Actual CLASSpp 26 calculation of previous strong thermal-parent PSD completed to k=1.5, P/P_CDM near unity and partition drift .2023%; broad and hierarchy controls pending. Tests 20/20. Later weak empty-parent fit and full parent momentum/thermal physics are next priorities.
