# A matched early-expansion test with exactly the same inner state

Previous q16 vs q32 runs used independently prepared source archives.
This stricter comparison starts with the same populated embedded q32
archive in both runs. In the cropped control only its four original
outer shells are removed; every surviving entropy coefficient is
identical. The removed fractions are 2.2135e-5 scalar and 3.4822e-5
fermion by number. No artificial vacuum cell or positive floor is
added. The cropped and full runs conserve their own retained
counts and use their own dark density in H.

With the same exact dark operators and compact elastic panels, evolve
to b=1.1066819197 on matching time grids:

| Intervals | T_phi kept wide (MeV) | T_phi cropped (MeV) | Fraction in first scalar bin, wide / cropped |
| ---: | ---: | ---: | ---: |
| 4 | 9.928019 | 9.930221 | 1.98047e-8 / 1.97913e-8 |
| 8 | 9.919155 | 9.921529 | 1.98384e-8 / 1.98216e-8 |

The scalar effective-temperature response to retaining the exterior
is about -0.024% at eight intervals. The first-bin fraction changes
by about 0.085%, within numerical and model limitations; comparing
fractions also includes the small normalization change in particle
number. The temperature time-step shift from four to eight intervals
is several times larger than the cutoff response. These are controlled
short-expansion results, not evidence about an arbitrary late tail.
The sampled on-shell channel does not directly access this first bin.

## Initial-guess diagnostic for stiff intervals

An optional fallback solves successively larger collision intervals
as **initial guesses**, always targeting the same old moments.
Only the final root with the full interval is applied. Invariants and
the exact residual criterion are unchanged. In a manufactured
fixed-state check the directly solved and continuation-solved moments
agree to 3.8e-15 in the local count normalization, and the exact
centered residual is 1.3e-14. Saving the rounded sum of the reference
and perturbation alone and then subtracting loses enough precision
that recomputing the stiff residual gives 3.8e-9 on this test. The
optional diagnostic exposes the separately stored accepted reference,
perturbation and parent ratio for same-step verification. It does
not change the physical state used by the expansion solver.

## Reproduce

Prepare [both EOS inputs](INDEPENDENT-BENCHMARK-SETUP.md). From
`kinetics`, use the source-generated embedded q32 archive:

```sh
export OPENBLAS_NUM_THREADS=1
BASE=(EOS.dat results/wide32-embedded-vacuum-2026-11-26.json results/wide32-embedded-vacuum-shape-2026-11-26.json)
for N in 4 8; do
 python -W error late_archived_shape_continuation_reference.py "${BASE[@]}" --steps "$N" --stop 1.1066819197 --prompt-until 20 --radiation-table EOS.late.dat --dark-gravity --compact-elastic
 python -W error late_archived_shape_continuation_reference.py "${BASE[@]}" --steps "$N" --stop 1.1066819197 --prompt-until 20 --radiation-table EOS.late.dat --dark-gravity --compact-elastic --crop-edge 16
done
python -W error prompt_initial_guess_audit_reference.py EOS.dat empty.json
```
