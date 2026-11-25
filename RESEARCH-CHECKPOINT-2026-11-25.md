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
- Earlier-onset [run 45625464842](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/45625464842) was running at checkpoint creation; collect actual output, do not infer success.

## Next actions (in order)

1. Collect earlier-onset CLASS results and test normalization/half-power convergence. Refine the perturbation hierarchy around the half crossing and at high k (the new fluid control shows 16–26% differences at k=1–1.5). If early spectrum has no crossing in tested range, extend the physical k range rather than assuming a result.
2. Implement finite fermion daughter in the conservative asymmetric energy collision kernel using the tested geometry; extend total stress-energy to both traces. Derive temperature-dependent dispersion/chemical rates for a physically specified Lagrangian, with the early burst's tiny IR scale as the decision criterion.
3. Only after conditional production/transfer checks, compare admissible scenarios to CMB/Lyman-alpha and dwarf-galaxy likelihoods and close the nonlinear halo/relaxation link.

A new research program is not necessary. Continue this repository, branch main, preserve previous evidence. Computer use is optional and restricted to explicitly authorized project location; no personal machine resources required for this checkpoint.
