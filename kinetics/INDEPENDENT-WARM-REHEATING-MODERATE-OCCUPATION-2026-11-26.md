# Continue the early-seeded thermal comparison from .28 to moderate occupation

A wider control keeps the same early FD leptonic
sources and exact signed thermal affinities. It uses
MB event factors with FD/BE equilibrium moments, like
the [earlier seed test](INDEPENDENT-WARM-REHEATING-EARLY-SEED-SENSITIVITY-2026-11-26.md).
The numerical comparison now stops at 40 MeV or an
imposed scalar phase-factor ceiling of .65. It does
not assume that reaching .28 is a physical transition.

All five near-empty initializations reach 40 MeV in
this **reduced** model without meeting .65. The sampled
maximum scalar phase factor is about .639 around 62 MeV.
For the 99.98 MeV seed, the last retained light and scalar
equivalent yields are respectively 4.2863e−4 and
2.2127e−4; the last scalar phase factor is .61481.
An independent saturated thermal normal-branch
solution at the same light-equivalent number and energy
gives a scalar-invariant/capacity ratio .7521 at visible
63.35 MeV and .7467 at 40 MeV. That is not a
bound against any nonthermal condensate.

Interpolation of the recorded states at visible 60 MeV
gives S≈1.741e−4 and a dark temperature of about
14.95 MeV, compared with S≈5.87e−5 and a
different dark temperature in the independently
imposed 60 MeV seed. This difference is a useful check
on early population feedback in the chosen closure.
It does not show that the thermal model is valid
throughout the cooling epoch or determine the true
yield.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error warm_prompt_reheating_early_seed_sensitivity_reference.py EOS.late.dat --peak .65 > mids.json
python -W error warm_prompt_reheating_saturation_distance_reference.py EOS.late.dat mids.json
~~~

[Summary and links by initialization](results/warm-prompt-reheating-mid-occupation-summary-2026-11-26.json);
[saturation comparison](results/warm-prompt-reheating-mid-occupation-saturation-2026-11-26.json).
Raw arms:
[99.98](results/warm-prompt-reheating-mid-occupation-arm1-2026-11-26.json.gz),
[99.90](results/warm-prompt-reheating-mid-occupation-arm2-2026-11-26.json.gz),
[99.50](results/warm-prompt-reheating-mid-occupation-arm3-2026-11-26.json.gz),
[99.00](results/warm-prompt-reheating-mid-occupation-arm4-2026-11-26.json.gz),
[97.00](results/warm-prompt-reheating-mid-occupation-arm5-2026-11-26.json.gz).
