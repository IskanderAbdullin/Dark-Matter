# Explicit momentum events witness the non-nested-grid shape difference

The previous complete fine/third D2 band comparison certified a sharp minimum compatible normalized TV from the two native panel ledgers. Here we construct equivalent evidence from the reported integral capacity cover, without relying on a solver objective in this step.

Select all fine panels not in the left half of the cover. Every third panel intersecting this physical set in a positive-length interval belongs to the right half of the cover, so its normalized mass gives an upper bound even when individual parts of crossing panels are unknown. The opposite construction selects all uncovered third panels. For all three response columns the same six disjoint fine intervals and seven disjoint third intervals result. The stated lower bounds therefore refer to explicit sets in momentum, not to an interpolation or half-interval crossing.

| Response | Fine mass on fine event | Third mass on same event, at most | Guaranteed fine minus third |
| --- | ---: | ---: | ---: |
| Forward number | 36.56149% | 10.32806% | 26.23343 pp |
| Inverse number | 36.96027% | 10.55284% | 26.40743 pp |
| Forward excess | 47.51821% | 16.02452% | 31.49369 pp |

The opposite event yields the same lower bound in the other direction. The exact edges, positive-length cover test, repeated mesh certificate and uncertain opposite masses are included in the receipt. These statements concern compatible finite held cold response shapes over the retained band, not the actual unbinned spectra or a whole-production prediction.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results; d=2026-11-28
python -W error early_pole_endpoint_cold_disjoint_witness_event.py \
 --fine "$p/fine-d2-endpoint-continuum-cold"-{1,2}.jsonl.gz \
 --third "$p/third-d2-endpoint-continuum-cold"-{1,2}.jsonl.gz \
 --states "$p/early-pole-cell-endpoint"-{fine,third}-d2-x05-state-$d.json.gz \
 --archive "$p/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json" \
 --cdf "$p/early-pole-cold-fine-third-d2-native-cdf-$d.json"
```
