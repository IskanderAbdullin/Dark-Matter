# Warm continuation to 0.414 MeV and the finite-support pair limit

All runs start with the existing signed warm $T_R=50$ MeV
benchmark, not with a new high-reheat setup. The first two
stages use the [3.43 MeV](INDEPENDENT-WARM-COUPLED-B35-2026-11-26.md)
and [1.206 MeV](INDEPENDENT-WARM-COUPLED-B10-2026-11-26.md) restart
surfaces. As a higher time-resolution check, the 48+42-step branch
gives conditional stable (Omega h^2=57.3573) at 1.206 MeV,
0.157% below the 32+28 result of 57.4472. This does not test
independent momentum refinement or omitted physics.

Restart the latter branch at 1.206 MeV and expand by a further
factor of 3. At 0.413763 MeV the stable fixed-count proxy is
57.04316 at 20 steps and 57.04240 at 32 steps. The extra losses
are 0.70339% and 0.70472%, respectively, counting prior exterior
cohorts only once. The scalar fixed-count proxy is about
2.7972e-4; the maximum scalar occupation in the
resolved 32-step state is 7.29e-4. The cumulative
count audit checks each exported state and records the files'
hashes. There is no inverse flux from the retained scalar support
in either final-stage run.

## Exact support reason, and a conditional remaining-loss interval

Assume exterior populations do not re-enter or stimulate the represented
collision kernel. The scalar and stable distributions have the same
maximum represented momentum (P/b) at expansion (b) from the
1.206 MeV anchor, with (P=17.466) MeV. The scalar edge energy
$\sqrt{m_\phi^2+(P/b)^2}$ is less than $m_\chi$. Two retained
scalars cannot reach the stable-pair threshold, and at least one
scalar produced by every stable pair is above the edge.

For two massive incident particles with $p_i\le p=P/b$, the lower
lab-frame daughter energy satisfies
$$
E_{\phi,\min}\ge (E_1+E_2-|\mathbf p_1+\mathbf p_2|)/2
\ge\sqrt{m_\chi^2+p^2}-p.
$$
Compare that expression with the scalar edge. Equality occurs
at (p=)11.5470 MeV, so for (b>1.5125682) **both** forward
scalars are outside the retained support. This is a kinematic
property of the represented kernel, not proof that there is no
physical inverse or stimulated term from a broader spectrum.

In the [previous quadratic loss bound](INDEPENDENT-WARM-PAIR-LATE-SUPPORT-BOUND-2026-11-26.md),
the uniform Bose factor can therefore be replaced by
$1+F(1-b_{\rm crit}^{-3})$, where $F$ bounds only the retained
occupations. Here $1-b_{\rm crit}^{-3}=0.7110$.
For example, the conditional lower stable proxy rises to 38.85
at (F=1). At the newer 0.414 MeV anchor both outgoing states are
already outside support. The remaining loss integral is at most
0.00865, so the stable fixed-count proxy stays in **[56.55, 57.04]**
for this restricted continuation: no inverse, no active exterior,
no new channels or visible entropy production, and the same
photon expansion floor. These are not observational bounds.
Ordinary finite-support evolution tracks missing particles in count
and energy ledgers; it does not track their momentum distributions.

## Retain the new outgoing high-momentum scalars as a separate control

To test the artificial exterior cut, evolve only the signed pair
operator on the same prescribed radiation clock in two 20-step runs.
Start from the exact 1.206 MeV retained distribution. In one run use
the original momentum cutoff (q/T_R=16). In the other use 96,
adding initially empty shells. It captures essentially all new
daughters, so newly produced scalars can enter the inverse operator
and stimulate later forward reactions. This control omits the
elastic and prompt blocks and the dark contribution to (H),
and must be compared with its matched pair-only run.

| Pair-only support | Stable loss from this interval | Peak reverse/forward | New exterior scalar number / initial retained scalar number |
| --- | ---: | ---: | ---: |
| Original, up to 16 | 0.630569% | 0 | 0.0055891 |
| Padded, up to 96 | 0.630568% | 4.223e-6 | 6.0e-24 |

The two inclusive stable output counts differ fractionally by
about 1.23e-8 in this control. The padded run
does not restore any exterior population produced before 1.206 MeV.
This tests feedback from *newly* produced tails, not general
convergence in support, collision quadrature, or the cosmology.

## Reproduce

```sh
export OPENBLAS_NUM_THREADS=1
r=results/warm-coupled-pair
python -W error coupled_warm_pair_continuation_reference.py EOS.dat \
  "$r-b35-48-shape-2026-11-26.json" --stop 2.857142857142857 \
  --steps 42 --save-snapshots fine-1206-shape.json >fine-1206-rows.json
python -W error coupled_warm_pair_continuation_reference.py EOS.dat \
  "$r-b10-28-shape-2026-11-26.json" --stop 3 --steps 20 32 \
  --save-snapshots late-shapes.json >late-rows.json
# Exact single-grid archived inputs used for the stage and bound audit:
python -W error warm_coupled_multistage_audit_reference.py \
  results/50MeV-from-empty-pair-budget-2026-11-26.json \
  results/50MeV-from-empty-pair-quantum-shape-2026-11-26.json \
  --stage "$r-b35-32-2026-11-26.json" "$r-b35-32-shape-2026-11-26.json" \
  --stage "$r-b10-28-2026-11-26.json" "$r-b10-28-shape-2026-11-26.json" \
  --stage "$r-b30-32-2026-11-26.json" "$r-b30-32-shape-2026-11-26.json" >audit.json
python -W error warm_pair_exterior_support_gate_reference.py \
  results/warm-scalar-pair-late-support-bound-2026-11-26.json \
  --late-result "$r-b30-32-2026-11-26.json" \
  --end-shape "$r-b30-32-shape-2026-11-26.json" --end-count-audit audit.json
# Separate matched pair-only controls, n=20 with and without --pad 96:
python -W error scalar_pair_expanding_continuation_reference.py EOS.dat \
  "$r-b10-28-shape-2026-11-26.json" --stop 3 --steps 20 --quadrature 3 7 6 2
python -W error scalar_pair_expanding_continuation_reference.py EOS.dat \
  "$r-b10-28-shape-2026-11-26.json" --stop 3 --steps 20 --quadrature 3 7 6 2 --pad 96
```
