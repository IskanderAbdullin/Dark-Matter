# Matched scalar-self proxy: temporal refinement and a second prepared state

At fixed exposure 2e−8, repeat the [matched 97% fine screened control](INDEPENDENT-EARLY-POLE-SCALAR-SELF-MATCHED-SCREENED-2026-11-27.md) with sixteen rather than eight macro steps. In a separate test, use the third archived state with a refreshed scalar shape-proxy mass and vertex (both arms still reset their screened graph to the same pristine starting weights). All cases use sixteen scalar bins, angular order six, the same finite hard channels, and the two screened elastic terms.

| Archive and scalar treatment | Steps | Integrated gain without / with scalar self | Relative self effect |
| --- | ---: | ---: | ---: |
| Fine, held | 8 | 9.69149065e−10 / 9.69880488e−10 | +0.07547% |
| Fine, held | 16 | 1.00386127e−9 / 1.00441383e−9 | +0.05504% |
| Third, refreshed | 8 | 1.60861496e−9 / 1.60721737e−9 | −0.08688% |

On fine the no-self gain changes by +3.58% under time refinement, substantially more than the within-grid scalar term, so its small positive difference is not a continuum measurement. The fine decay source is initially negative and first sampled positive at 1.25e−9 for both sixteen-step arms (eight-step first positive sample: 2.5e−9). On third both are already positive at the initial sample. Its slightly negative scalar-self effect does not isolate either the archive or the response refresh because both differ from fine/held. Matched cross-controls are needed for that attribution.

In each pair the scalar proxy is optional and phenomenological. Maximum energy and entropy closure residuals over the compared samples are below 7e−15 and 5e−13. No claim is made about the continuum collision operator, a self-consistent in-medium vertex, or a final abundance.

From `kinetics/`, with one BLAS thread:
```sh
OPENBLAS_NUM_THREADS=1 python early_pole_joined_scalar_self_proxy_reference.py fine \
  --duration 2e-8 --steps 16 --hard-slow --screened --mixes .97
OPENBLAS_NUM_THREADS=1 python early_pole_joined_scalar_self_proxy_reference.py third \
  --duration 2e-8 --steps 8 --hard-slow --screened --mixes .97 --refresh
```
Full samples: `results/early-pole-scalar-screened-time-refresh-2026-11-27.json.gz`.
