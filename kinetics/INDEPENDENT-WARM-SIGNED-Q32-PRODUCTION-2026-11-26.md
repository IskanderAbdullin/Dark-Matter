# Extend warm production support without changing the inner bins

Repeat the [48-step signed warm production benchmark](INDEPENDENT-50MEV-PAIR-FROM-EMPTY-2026-11-26.md)
starting from an empty dark sector at 50 MeV. Keep the entire
q/T_R ≤ 16 grid unchanged, including the first resolved cell,
and add five outer bins up to q/T_R = 32. Preserve FD forward and
inverse lepton source, prompt return, below-pole elastic scattering,
scalar pairs, and the visible energy/entropy feedback. Both
versions reach 11.999710 MeV at the same scale factor.

| At the 12 MeV checkpoint | Original support | Extended support |
| --- | ---: | ---: |
| Resolved momentum bins | 19 | 24 |
| Inclusive stable fixed-count Ωh² | 87.187259 | 87.182991 |
| Inclusive scalar fixed-count Ωh² | 0.00015009616 | 0.00015011451 |
| Exterior stable fixed-count Ωh² | 0.0129962 | 1.95e-8 |
| Archived exterior energy (comoving a⁴ρ) | 2.04e-12 | 5.15e-18 |

Enlarging the support captures nearly all those previously
unresolved particles and changes the inclusive stable result by
−4.90e-5 fractionally and the scalar result by +1.22e-4.
The gross inverse-to-forward pair ratio peaks at about 2.94
in both runs, so it would not be legitimate to remove the
inverse term. The normalized joint energy defects are
−5.96e-14 and −5.20e-14. The new shape preserves populated outer
bins for later use, rather than inferring their occupancies
from the old exterior energy ledger. At the new endpoint, the
resolved outer bins contain 2.42e-5 of retained scalar number
and 3.78e-5 of retained stable number. Their maximum occupations
are 6.35e-11 and 6.13e-11, respectively. These are measured in
the extended solution, not estimates obtained by putting the old
omitted count into those cells.

These are finite-resolution and fixed-count outputs at 12 MeV.
They are not a terminal relic, a change of the particle model,
or proof that all later exterior effects are negligible.
There is no change to the initial vacuum, density normalization,
or time mesh in this support-only control.

~~~sh
export OPENBLAS_NUM_THREADS=1
python -W error joint_resonant_frw_split_reference.py EOS.dat \
  --physical-amplitude --visible-source --visible-backreaction \
  --trace-exterior --with-elastic --empty-pulse 1e-8 --empty-stop .24 \
  --support-low .005 --vacuum-scaled-residual --allow-empty-shells \
  --include-zero-cell --pair-conversion --steps 48 --quantum-visible-source \
  --empty-support-high 32 --empty-outer-bins 5 --save-snapshots q32-shape.json >q32.json
python -W error matched_warm_outer_support_reference.py EOS.dat \
  results/50MeV-from-empty-pair-quantum-2026-11-26.json \
  results/50MeV-from-empty-pair-quantum-shape-2026-11-26.json q32.json q32-shape.json
~~~
