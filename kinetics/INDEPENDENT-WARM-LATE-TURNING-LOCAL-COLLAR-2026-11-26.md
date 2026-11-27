# Dynamical momentum-split control for the late .68 GeV turning clock

At the same 7.394 MeV fixed background and same six-step backward-Euler
clock, split four momentum shells in the narrow threshold neighborhood
(41 to 45 cells). Keep all previous boundaries and the cut at
0.7148213914435242 GeV. The initial occupations are identical physical
functions: restrict the affine kinetic-entropy profile into both halves
of each selected shell. The relative Hubble normalization is unchanged
at recorded precision; the largest species moment restriction L1
difference is 1.12e-22 relative, audited shell by shell.

All six steps solve independently on the refined phase. The previous
coarse-grid endpoints supply only Newton guesses, never the target
moments; each target comes from the preceding accepted refined phase.
Species counts and total closed-sector energy agree with the common
initial phase to 1.2e-15 relative; the maximum refined discrete excess
ledger difference is <1.7e-32 in units of initial dark energy.
The signed excess flux stays positive through Hdt=7.0e-8 and is negative
at 8.5e-8 as on the original six-step grid. The largest difference in the
excess above the fixed cut is 0.034% of the original sampled peak.

This controls a targeted change of momentum discretization at fixed
time mesh only, not full resolution convergence, a continuum change
time, extension of the far momentum support, or non-prompt reactions.
Time discretization moves the bracket more in the
[12 and 24 step comparisons](INDEPENDENT-WARM-LATE-TAIL-TIME-REFINEMENT-2026-11-26.md)
than this local change does at six steps.

~~~sh
a=results/warm-corrected-anchor-signed-floor-2026-11-26.json
seed=results/warm-late-aligned-turning-seed-2026-11-26.json.gz
mkdir -p coarse fine
python -W error warm_late_cut_turning_clock_reference.py EOS.dat "$a" "$seed" coarse --substeps 6 >coarse/report.json
python -W error warm_late_turning_local_collar_replay_reference.py EOS.dat coarse fine >fine/report.json
python -W error warm_late_turning_local_collar_audit_reference.py EOS.dat coarse fine
~~~

[Audited comparison](results/warm-late-turning-local-collar-audit-2026-11-26.json)
and [hashed refined phases](results/warm-late-turning-local-collar-six-2026-11-26.json.gz).
