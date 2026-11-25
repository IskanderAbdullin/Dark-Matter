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

## Late weak-parent and actual CLASS checkpoint

# Late weak-parent candidate and local CLASS validation (2026-11-25)

## Late weak-parent kinetic fit
M1=86 TeV, m2/M1=0.9, m_phi=86 eV; x0=0.1, initially empty parent, kappa_SM1=47, kappa2=0.03, kappa=0.005278180162600994. Width was fit on 160 cells and held fixed on 240/320 cells. Omega h^2=0.1207153994 and 0.1209727907, hot fraction q/T_actual>0.1=0.003565679 and 0.003556194. Nominal uncancelled seesaw mass ~0.05024 eV. Equilibrium thermal scalar mass/peak ratio ~0.00359 (contrast early weak case ~0.26); this is only a simple-Yukawa thermal diagnostic. Independent bath-only shape screen: peak TV distance 0.01818, energy discrepancy 2.04%, number discrepancy <0.1%. It excludes dark depletion and cannot establish closure in the stimulated regime. Unprotected loop scale ~28 MeV versus 86 eV; protective symmetry remains essential.

## Local CLASS-community 26.0.0
For the earlier strong thermal-parent candidate (not the late weak candidate), actual synchronous linear transfer was calculated at z=0,5 and k=0.003--20/Mpc on both 240 and 320 cell PSDs. At z=5 k=10/Mpc the ratios P/P_CDM are 0.989236364 and 0.988734742. Fine-mesh minimum sampled ratio 0.983762068 at z=5; no half-power crossing in range. Independent three-tracer partition drift <=0.44%. These conditional collisionless linear results do not establish viability or observational likelihood. Small apparent large-scale enhancement is not physically validated.

Thermal Bose reference comparison with archived official CLASS 3.4.3: half-amplitude crossing z=5 at k=0.42961044 vs 0.43275181/Mpc (0.726% shift); abundance and early radiation agree. One deeply suppressed thermal mode at k=0.78978243 has invalid power ratio; recorded as null, not a physical value.

## Negative results and next action
Late weak candidate PSD export failed transformed quadrature tolerance: third tracer moment errors -0.00297, -0.02920, -0.000316. No CLASS spectrum for that candidate. Repartition/rescale tracer or refine quadrature and check moments before solver. Full hierarchy alternate run failed because CLASS-community 26 renamed tolerance input; adapter patched to source-confirmed tol_perturb_integration, not rerun. Prioritize late weak candidate with momentum-resolved parent including depletion and self-consistent finite-T self energies. No observational viability claim.

Raw numerical records: results/late-weak-and-class-2026-11-25.json. Local environment and public research branch preserve all completed calculations. Main not updated (server-side ref error); workflow draft not activated (missing workflow scope).

## Quadrature follow-up

Both late weak-parent meshes and both partitions fail 120/240-node moment checks. All four pass at 480 nodes (max relative error 2.28e-5). Added explicit --bins 480 option with 960-node alternate check; analytic export passed both meshes. Recorded full checks in kinetics/results/late-weak-quadrature-2026-11-25.json and kinetics/results/late-weak-export-2026-11-25.json. The earlier failed coarse check is preserved. Next run actual CLASS at these validated settings; transfer and full hierarchy remain untested for this candidate.

Final regression: 24 passed (16.74 s). Strong-candidate CLASS cross-mesh max relative transfer-power drift: 0.0008230023 at z=0 and 0.0008229460 at z=5. Full hierarchy and late weak CLASS transfer remain untested.

## Actual late weak-parent CLASS

See cosmology/LATE-WEAK-CLASS-VERIFIED-2026-11-25.md. Old 480-node run failed due to nonfinite state; nonzero subnormal tail causes q/f0 overflow. Slow negligible positive floor and q^5 f hot-scale selection resolve it; actual two-mesh CLASS at 120/240 nodes gives P/P_CDM(z=5,k=10)=.98518/.98453, minimum fine .97567 up to k=20/Mpc. Conditional only; full hierarchy in progress.

## First fully momentum-resolved parent

A fractional-stoichiometry finite-volume collision operator interpolates the parent energy between adjacent shells and uses geometric products of parent availability/blocking. At positive occupation, forward/backward log ratios enforce exact discrete detailed balance and particle/energy conservation. A tiny C1 continuation below parent occupancy 1e-15 is optional; the exact and regularized n=64 solutions through x=.3 differ in Omega proxy by less than 9e-7 relative. Both show a stimulated burst, but the Omega proxy at x=.3 is only .074999; this is an intermediate state, not a relic prediction. Pauli checks and FRW residuals pass; spontaneous vacuum normalization is discretization limited (up to 12% on this coarse mesh at x=20). Fine and long integrations are running, not yet validated. Data: kinetics/results/resolved-parent-early-2026-11-25.json. Code: kinetics/resolved_parent_bath.py.

## Published portal rates and momentum-resolved evolution

See kinetics/LAINE-AND-RESOLVED-PARENT-2026-11-25.md and kinetics/results/laine-and-resolved-parent-2026-11-25.json. Actual Ghisoiu-Laine high-T table changes g=4 fixed-width Omega=.115007 (320), retune .006 gives .119965 (240), while g=2 fixed-width underproduces .038792. Fully resolved n=80 reaches x=80 with .112865; strict n=120 budget-limited at x=7.758, late energy-matched hybrid to x=80 gives .115815. Early hybrid source mismatch 80% is a negative result; late handoff mismatch 0.39%. Fine and two-spin scan pending. These are conditional rates/kinetic models, not physical viability.

## Spin scan, late-handoff refinement and identical-pair interpretation

See kinetics/SPIN-AND-HANDOFF-2026-11-25.md. g=2 width scan up to .12 gives at most .076764 for tested single pair. g=4 exactly equals two independent identical Majorana pairs in the symmetric occupation subspace (numerical operator test passes), but flavor/neutrino mass and dark thermal self-energies remain unresolved. New table + explicit parent hybrid gives .098908 at old width; n=320 number closure at .006 gives .120218, and n=240 CLASS power ratio .984315 at k=10,z=5. Older-rate hybrid n=200 .117509; same-grid two late switches agree at 1.6e-5 relative. Adaptive two-check handoff and separate cross-kernel guard added; retuning pending. Existing and initial new tests 29 passed, extended targeted 6 passed.

## Final width/transfer checkpoint

[Final evidence and next priorities](kinetics/FINAL-AUTONOMOUS-RUN-2026-11-25.md). Published-rate closure CLASS converged on two meshes, but new resolved width trials rejected: positivity failure at 0.006 and 28.22% handoff mismatch at 0.012. No accepted full-kinetic abundance fit.

## Parent-shell aliasing and independent angular quadrature

See kinetics/PARENT-GRID-AND-ANGLE-2026-11-25.md and raw results. Found spurious legacy parent-shell loads up to 12543 times physical width; guarded BDF reproduces and rejects accepted negative parent step. Uniform parent grid reduces max load ratio to 3.60. Independent parent-angle kernel with separate heavy grids reaches x=80 quickly; at kappa=.006 n120/240 Omega=.116493/.115436, unlike old shared-grid hybrid .102325/.104607. Analytical spontaneous spectrum favors separated rectangle kernel (n480 peak bias .01565%) over angle (.3385%). Need full comparison at independent meshes and widths; no accepted physical fit yet.

## Matched-width independent kinetic kernels

See kinetics/CONVERGED-INDEPENDENT-KERNELS-2026-11-25.md. At fixed kappa=.00685, full angle n240/320 yields .1203498/.1202148, same-grid independent rectangle hybrid .1197704/.1198517. Strict exact dense n240 .1203496. CLASS conditional P/P_CDM at z=5,k=10/Mpc=.983907/.983511. This supersedes earlier shared-grid positivity failures as a numerical solution, not a physical viability result. The highest-priority next issue is early scalar self-interaction and thermal dispersion, including possible condensate and number-changing reactions. Full CLASS hierarchy remains unverified.

## Contact self-interaction: physical bottleneck

See kinetics/CONTACT-ELASTIC-AND-CONDENSATE-LIMITS-2026-11-25.md. Fixed-PSD quantum net elastic rate of relative q2 change is +5.1875/H even for sigma/m=1e-7 cm2/g; independent conservative scheme +5.2334/H. At .001 bare diagnostic +3.56e4/H, uncontrolled lambda*f~11. Collisionless assumption therefore requires new physical justification; this is not a solved interacting relic or exclusion. Ideal number-conserving condensation/free-streaming limit yields P/P_CDM(z5,k10)=.968895; ideal tight fluid sound horizon to equality .12828 Mpc. BEST (2026 exact kinematics and stiff integration) and CosmoLattice identified for next-stage implementation.

## Autonomous continuation: final verified state

Full regression after BOM correction: 40 passed in 48.67 s. All four new elastic diagnostics and five kernel tests included. Two independent scalar-collisionless kinetic discretizations at matched width converge near Omega h^2=.12; strong early contact interaction invalidates extrapolating the collisionless result to halo-thermalizing models without new physics. Even the weak contact screen at 1e-7 gives relative q2 change ~10.19/H during burst x=.2 and ~5.19/H at endpoint, independently checked ~5.23/H. Explicit relaxation trial exceeded 90 s; no integrated interacting evolution. Full CLASS hierarchy trial stopped without final result; no claim of observational viability. BEST v1.3.5 (c70875be714e7a6bc6d38578c34bb0dc9f5aefa9) checked out and condensate example inspected; mpi4py/vegas missing. Prioritize its controlled independent validation, then self-consistent interacting production, zero mode and number-changing reactions. Details and raw results: kinetics/CONTACT-ELASTIC-AND-CONDENSATE-LIMITS-2026-11-25.md and kinetics/results/contact-elastic-and-condensate-limits-2026-11-25.json.

## BEST benchmark and conservative stiff relaxation

See kinetics/BEST-AND-IMPLICIT-ELASTIC-FIRST-2026-11-25.md. BEST v1.3.5 actually runs in isolated Intel MPI environment and matches moderate direct moment benchmark to .265%, but high-occupancy actual PSD n160 gives spurious number and energy source +6.86/H, +109.98/H. Do not use that MC result for physical elastic evolution. Deterministic analytic Jacobian BDF completes fixed-mass fixed-geometry relaxation at sigma/m=1e-7; n64 over .1 Hubble time changes q2 by 3.3-3.5%, with machine-level number/energy conservation. Remesh errors remain material. Fine grid and coupled massless production/elastic trials in progress.
