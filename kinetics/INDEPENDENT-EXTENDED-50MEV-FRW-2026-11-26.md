# Longer controlled evolution from 50 MeV reheating

The earlier 45 MeV controls sample only a minority of the integrated visible
source. This check extends the same independent two-species calculation to
reference visible endpoints of 30 and 20 MeV, and includes a 15 MeV vacuum
source control. The visible plasma evolves with energy feedback and exterior
source cohorts. Final scale factors are set by the corresponding reference
temperature on the adiabatic visible trajectory; actual endpoints are lower
by the quoted thermal feedback. Retained comoving support is [.05,16] in
units of T_R unless noted.

| Nominal endpoint | 24-step scalar energy fraction (matched vacuum, quantum evolution) | Quantum-run visible temperature shift | Max number-source difference vs independent integral |
| ---: | ---: | ---: | ---: |
| 30 MeV | 0.28138597, 0.28138704 | -1.91e-5 | 3.5e-4 |
| 20 MeV | 0.25913742, 0.25913867 | -2.25e-5 | 1.01e-3 |
| 15 MeV | 0.24230079, not run | -2.37e-5 (vacuum) | 1.65e-3 (vacuum) |

The 20 MeV 48-step matched vacuum calculation gives 0.25906222, so
time refinement exceeds the matched source correction. At 20 MeV the
24-step vacuum run with a finite-volume first cell extending to zero
gives 0.25913177 and 1.37e-7 of scalar particles below q/T_R=0.01.
This agrees in scale with the early soft-cell result. It is not precision
convergence of a differential spectrum or proof of a singular mode near
zero. The same a⁴ energy and particle-number trace ledgers include
exterior source injection; the 48-step control still conserves represented
counts to 3.2e-13 relative.

The independent temperature integral uses the same vacuum portal and a
prescribed adiabatic background. It is a useful check of cumulative
production and the time measure, not the exact solution with
occupational blocking or visible energy feedback. The maximum complete
portal count and injected energy differences in the 20 MeV 24-step
quantum run are 1.01e-3 and 8.08e-4 respectively. Retained numbers are not
substituted for the complete source: the low/high daughter moments and
the initial pulse are explicitly included in this comparison. These are
not final cosmological abundances.

## Reproduce

Prepare the external `EOS.dat` with the [shared benchmark setup](INDEPENDENT-BENCHMARK-SETUP.md).

In `kinetics`, use the external central EOS:

```sh
export OPENBLAS_NUM_THREADS=1
COMMON=(EOS.dat --physical-amplitude --visible-source --visible-backreaction --trace-exterior --with-elastic --empty-pulse 1e-8 --support-low .05 --vacuum-scaled-residual)
for END in .6 .4; do
  for SRC in matched-vacuum-order quantum-visible-source; do
    python -W error joint_resonant_frw_split_reference.py "${COMMON[@]}" --"$SRC" --empty-stop "$END" --steps 12 24
  done
done
python -W error joint_resonant_frw_split_reference.py "${COMMON[@]}" --matched-vacuum-order --empty-stop .3 --steps 24
python -W error joint_resonant_frw_split_reference.py "${COMMON[@]}" --matched-vacuum-order --empty-stop .4 --steps 48
python -W error joint_resonant_frw_split_reference.py "${COMMON[@]}" --matched-vacuum-order --empty-stop .4 --steps 24 --support-low .005 --allow-empty-shells --include-zero-cell
python -W error compare_extended_frw_source_reference.py EOS.dat results/visible-production-cumulative-2026-11-26.json results/extended-frw-20-quantum-2026-11-26.json results/extended-frw-30-quantum-2026-11-26.json results/extended-frw-15-vacuum-2026-11-26.json
```
