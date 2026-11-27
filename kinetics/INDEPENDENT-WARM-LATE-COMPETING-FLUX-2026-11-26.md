# Competing dark reactions across the late .68 GeV turning interval

Use the same reactivated phase from the 7.394 MeV archive and the
same physical energy cut, `C = 0.7148213914435242` GeV. Measure
`W = sum of (E − C)_+ over daughter particles`. The known stored packets were
reactivated before this interval, and subsequent exterior output in
the pair reaction returns no incoming particles.

## Frozen channel comparison

Compare four exact saved phases of the six-step prompt run at Hdt
(1.0,7.0,8.5,10.0) times 10⁻⁸. Use the same norm ρH,
with ρ the initial represented dark energy. In the closed
below-pole elastic operator the ordinary lab-energy cap
√s + 32T does not reach the observable: all outgoing and
incoming energies are below C. The zero in that calculation is a
kinematic truncation, not evidence of a vanishing physical rate.
Reopen the lab boost up to the ceiling set by the closed momentum
support of both species (lab energy at most 2.1002 GeV). Still
restrict the invariant mass to s < 0.82 M². Keep inverse elastic
events at the same lab energy and independently verify the event
W ledger against the affine moment projection.

With the greater lab support, the first elastic signed current is
(1.81–2.41)e-18; it is negative at each of the
three later samples under both quadratures. The signed pair current
starts at (8.33–8.81)e-19. Counting both
directions of W transfer, the summed elastic and pair gross
currents are at most 3.57e-6 of the absolute prompt
current across these samples and orders. The prompt current changes
from +2.454e-12 to -2.139e-12 at the
two samples bounding the turn. Pair output beyond the resolved
mesh contributes at most 1.81e-28 to the newly outgoing
W current. These are sampled rates, not a bound between samples.
In particular, a small direct tail current need not imply small
heat transfer or small rates in the cold bulk.

## One matched expanding block

On the six-step phase at Hdt = 7.0e−8, add
Δlog a = 1.5e−8 and compare three matched calculations.
Use the same first and second redshifts, fixed physical cut and
prescribed visible entropy background. Compute W by integrating
the reconstructed profile with the cut as a quadrature boundary,
not by reusing the original cut-aligned linear weights after redshift.
The first redshift agrees with the independent collisionless Liouville
calculation to 3.7e-14 relative in W. The pair
reaction uses two explicit midpoint half-kicks around the implicit
prompt and optional closed elastic collision. The prompt-only frozen
endpoints serve only as solver guesses. The largest measured count
and expansion-corrected energy discrepancies are below 2.0e-14.

| Selected reactions | Endpoint W / initial dark energy | Rel. change from expanding prompt |
| --- | ---: | ---: |
| Prompt | 2.3956624643e-18 | 0 |
| Prompt + pair | 2.3956624391e-18 | -1.0524e-8 |
| Prompt + pair + closed elastic | 2.3956624314e-18 | -1.3738e-8 |

Free redshift alone reaches 2.4277413884e-18 from
2.4277453379e-18. Thus the decline in the matched block
is principally reaction-driven. The prompt-plus-pair count exchange
is about +8.40e-9 of the initial scalar number and
-5.14e-9 of the initial light-fermion number, not a
change in their summed number when including pair output.

As an independent local clock control, start from the common
12-step archived phase at Hdt=7.0e-8. Advance the same
Δlog a in one block or in two half blocks, with the second
block starting from the actual accepted first endpoint, including
any positive exterior pair receipts. Both prompt-only and full
settings end with smaller W. Halving the block size increases the
endpoint by 0.2544% in both settings; the addition of pair and elastic
changes it by only −1.38e−8 and −1.57e−8
in the one- and two-block controls. The qualitative result survives
this halving, while the time-discretization shift is substantially
larger than the measured channel correction. The two trajectories
do not establish a continuum relative correction.

This is one local window in a known reactivated stress shape, not an
evolution of the unknown exterior or a cosmological prediction.
Elastic includes only the closed projected below-pole channel.
Prompt return is the stationary-parent projection and does not
evolve visible scattering, slow source, full heavy recycling,
scalar number-changing interactions or perturbations.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MALLOC_ARENA_MAX=2
a=results/warm-corrected-anchor-signed-floor-2026-11-26.json
seed=results/warm-late-aligned-turning-seed-2026-11-26.json.gz
mkdir -p turn6 turn12 competing
for n in 6 12; do
  python -W error warm_late_cut_turning_clock_reference.py EOS.dat "$a" "$seed" "turn$n" --substeps "$n" >"turn$n/report.json"
done
for k in elastic pair; do
  python -W error warm_late_reactivated_competing_flux_reference.py EOS.dat "competing/$k-low.json" "$k" turn6/step-{0,4,5,6}.json
  python -W error warm_late_reactivated_competing_flux_reference.py EOS.dat "competing/$k-high.json" "$k" turn6/step-{0,4,5,6}.json --order-mode 1
done
run() { python -W error warm_late_reactivated_expanding_competing_reference.py EOS.dat "$@"; }
run turn6/step-4.json turn6/step-5.json competing/expanding-prompt.json prompt
run turn6/step-4.json turn6/step-5.json competing/expanding-pair.json prompt-pair
run turn6/step-4.json turn6/step-5.json competing/expanding-all.json all --extent 285.
python -W error warm_late_competing_flux_audit_reference.py turn6/step-4.json competing
for mode in prompt all; do
  run turn12/step-8.json turn12/step-10.json "competing/compare-$mode-coarse.json" "$mode" --extent 285.
  run turn12/step-8.json turn12/step-9.json "competing/compare-$mode-mid.json" "$mode" --extent 285.
  run turn12/step-9.json turn12/step-10.json "competing/compare-$mode-fine.json" "$mode" --extent 285. --predecessor "competing/compare-$mode-mid.json"
done
python -W error warm_late_competing_clock_audit_reference.py turn12/step-8.json competing
~~~

[Audit table](results/warm-late-competing-flux-audit-2026-11-26.json),
[clock-halving audit](results/warm-late-competing-clock-audit-2026-11-26.json),
[frozen diagnostics and accepted states](results/warm-late-competing-flux-tapes-2026-11-26.json.gz),
[clock-halving states](results/warm-late-competing-clock-tapes-2026-11-26.json.gz).
