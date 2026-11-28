# A small total source difference hides larger changes across individual scalar bins

Use the [matched fine screened comparison](INDEPENDENT-EARLY-POLE-SCALAR-SELF-MATCHED-SCREENED-2026-11-27.md) with the same initial prepared 97% parent, pristine screened starting weights in both arms, and sixteen scalar bins. This separately sampled short control ends at exposure 8e−9 with four macro steps. The only changed term is the held scalar shape-redistribution proxy. Save each scalar mode's actual cold number count and signed instantaneous decay/inverse source; compare the two arms bin by bin.

At the last sample, their scalar redistributions differ by 1.19386e−3 in the summed absolute bin counts, while the net source difference is +4.18162e−5. The summed absolute differences in individual bin source contributions are 3.76419e−3, roughly 90 times the net. Individual signed source differences span −.001269 to +.001027. At the earlier sample 2e−9 that absolute-to-net ratio is about 16, so the mode cancellation is not simply the result of one terminal rounding operation. The sum of the recorded mode differences agrees with the independent difference of total rates to 1e−15 absolute.

This calculation explains why the total source can move very little despite noticeable changes in the modeled scalar spectrum. The factors are specific to this initial prepared spectrum, time window and finite channels; they are not a result on cosmological production or a general relation between the shape and total abundance.

From `kinetics/`:

```sh
OPENBLAS_NUM_THREADS=1 python early_pole_joined_scalar_self_proxy_reference.py fine \
  --duration 8e-9 --steps 4 --hard-slow --screened --mixes .97 --mode-diagnostics
```

Raw: `results/early-pole-screened-scalar-source-mode-cancellation-2026-11-27.json.gz`.
