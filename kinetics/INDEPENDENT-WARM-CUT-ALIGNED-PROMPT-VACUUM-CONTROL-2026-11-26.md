# Resolve the excess-energy threshold and the near-vacuum implicit solve

At the late corrected q64 anchor, place the
[known exterior receipts](INDEPENDENT-WARM-LATE-RECORDED-TAIL-REACTIVATION-2026-11-26.md)
in the test packet distribution and pad the momentum grid to
1.05 GeV. Evolve the closed prompt channel on a fixed clock
without visible terms, pair or elastic scattering. A backward
Euler root always targets the exact **old** number and kinetic
moments. For active near-empty shells, only the Newton guess
uses a positive reconstruction on the scale of the predicted
incoming moment. Exactly zero old moments therefore remain
zero in the residual and in the initial conservation ledger;
no arbitrary floor is inserted in the physical reference.

## The threshold must be a shell edge

For a fixed energy cut C=.714821 GeV, split the momentum
mesh at both mass-specific thresholds sqrt(C²−m²).
Restrict the archived affine entropy coefficients to the
new subcells without changing their functions. Then
max(E−C,0) is affine in each species' kinetic coordinate on
every shell, so its integral is an exact linear combination
of the resolved number and kinetic moments. Without the
split, it is not: an early explicit difference of two
unaligned reconstructed spectra gave 2.52e-10 per Hubble
time, while the separate angular calculation gave 1.48e-10.

With the aligned shell mesh, the initial native flux in the
.68 GeV case is +1.477763634e-10 per retained-energy Hubble
unit and the independent angular fold gives +1.477763644e-10.
The first accepted interval to Hdt=1e-11 produces
1.4772722e-21 of initial retained energy above the cut in
the excess metric. This is a short step validation, not yet
a converged lifetime.

## Assemble the same direct Jacobian locally

Evaluate the event flux in direct physical occupations,
rather than dividing nearly empty channels by a thermal
reference profile. The optional faster Jacobian accumulates
exact derivatives only within each parent phase panel, then
projects to the same global moments. On the three archived
test states the rates are bit-for-bit identical, the largest
column-relative derivative discrepancy is below 9e-14,
and independent finite differences agree within 6e-10.
Assembling the tested full basis took about 5–6 times less
wall time in this cloud run. This is not a timing guarantee.

Replaying all three exact-vacuum implicit microsteps with
the local Jacobian gives exactly the same recovered shell
moments as the original direct solver. The count and energy
defects are below 6e-16. No root tolerance was relaxed.

The earlier [unbounded-exterior rate upper](INDEPENDENT-WARM-LATE-INVERSE-MOMENT-ENVELOPE-2026-11-26.md)
uses the original retained cell maxima. It must not be
applied automatically to these later populations: by
Hdt=1e-8, the .68 case already exceeds selected old maxima
even in the retained cells. The number above the cap remains
small (5.79e-10 of the original scalar count and 2.41e-10 of
the original light count), but their largest relative
increases occur precisely in the initially dilute cells.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
s=results/warm-late-cut-aligned-hard068-test-step0-2026-11-26.json.gz
t=results/warm-late-cut-aligned-hard068-test-step1-2026-11-26.json.gz
u=results/warm-late-cut-aligned-hard068-test-step4-2026-11-26.json.gz
v=results/warm-late-unaligned-hard080-test-step0-2026-11-26.json.gz
python -W error warm_compact_absolute_prompt_reference.py EOS.dat "$s" "$u" "$v"
python -W error warm_compact_absolute_prompt_accepted_control_reference.py EOS.dat "$s" "$t"
python -W error warm_late_retained_cap_after_reactivation_reference.py \
  results/warm-multishell-corrected-q64-fine-stage-48-2026-11-26.json.gz "$s" "$u"
~~~

[Jacobian controls](results/warm-compact-direct-prompt-jacobian-2026-11-26.json);
[accepted numerical control](results/warm-compact-direct-prompt-accepted-2026-11-26.json);
[original cap persistence test](results/warm-late-retained-cap-persistence-2026-11-26.json).

## First sign bracket on one accepted time ladder

The aligned .68 GeV choice has positive signed excess transfer
at Hdt=1e-8 (+9.155e-11 in the units above) and negative
transfer at Hdt=1e-7 (−1.433e-12). The independent angle
fold gives the same signs. Thus there is at least one
reversal between those samples in this selected fixed-clock
prompt-only run. It does not locate the first or unique
turning time. By the latter sample the total accumulated excess
is 2.2752e-18 of the initial represented dark energy. This is
the directly integrated state moment, not a rate multiplied by
a guessed duration. The values use three equal implicit
substeps per reporting interval and are not a time-refined
comparison.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error warm_corrected_late_prompt_cut_aligned_reference.py EOS.dat \
  results/warm-corrected-anchor-signed-floor-2026-11-26.json \
  results/warm-multishell-corrected-q64-fine-stage-48-2026-11-26.json.gz trial068 \
  --hard .68 --substeps 3 --predictor
python -W error warm_late_cut_aligned_prompt_clock_audit_reference.py trial068 --intervals 5
~~~

[Accepted bracket](results/warm-late-cut-aligned-hard068-first-sign-bracket-2026-11-26.json);
[raw state at Hdt=1e-7](results/warm-late-cut-aligned-hard068-test-step5-2026-11-26.json.gz).
