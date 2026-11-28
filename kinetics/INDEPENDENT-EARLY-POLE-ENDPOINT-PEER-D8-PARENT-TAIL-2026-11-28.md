# Parent-tail weight on the two D8 cold endpoints at fixed k

Hold the independently evolved fine and third D8 charge states at x=.05 and refit their transverse poles at k=.0015. At each of q=3.5e−6 and q=1.5e−5 split the actual radial source by retained parent energy, then compare the entire normalized cosine average with a single cos=.9 ray. The two rules (4,6,12) and (8,10,20) agree for every cumulative coefficient partition, using each full coefficient for normalization. The unsplit direct integral also closes. This tests where the computed signal lies within the archive, not whether the archive is physically complete.

| Chart | q | Parent E≥28, forward share over all cosines | Parent E≥28, forward share on cos=.9 |
| --- | ---: | ---: | ---: |
| Fine D8 | 3.5e−6 | 0.2980% | 1.5615% |
| Third D8 | 3.5e−6 | 0.5926% | 3.4552% |
| Fine D8 | 1.5e−5 | 0.01610% | 100% |
| Third D8 | 1.5e−5 | 0.02590% | 100% |

At q=1.5e−5 on that ray, E≥30 accounts for 81.698% of the fine D8 forward coefficient and all of the third D8 coefficient. The inverse is resolved as zero there on both retained supports. Yet E≥30 contributes only 0.00818% and 0.01327% of the two angle-averaged forward coefficients at that q. A tail-dominated ray therefore does not establish tail dominance of the integrated angle average. The complete bare forward angle coefficient on third relative to fine is about .445 at q=3.5e−6 and .552 at q=1.5e−5. These are two pointwise bare coefficients, before the scalar occupation factors, not a band-integrated or physical rate.

This does not contradict the [kinematic k scan](INDEPENDENT-EARLY-POLE-ENDPOINT-PAST-K-CAP-APERTURE-2026-11-28.md); that test measured what angles are allowed at each k, not how much source they carry. Both controls are needed before inferring sensitivity to an energy cutoff, a wider k band or an archive extension.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-cell-endpoint; d=2026-11-28
for g in fine third; do
  gzip -dc "$p-$g-d8-x05-state-$d.json.gz" > "/tmp/$g-d8-parent.json"
done
python -W error early_pole_endpoint_cold_parent_tail_partition.py \
  /tmp/fine-d8-parent.json /tmp/third-d8-parent.json
```
