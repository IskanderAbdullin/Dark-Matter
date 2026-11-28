# Forward and inverse have different retained angular windows

Extend the [parent-only cap scan](INDEPENDENT-EARLY-POLE-ENDPOINT-PAST-K-CAP-APERTURE-2026-11-28.md) to the inverse reaction on the two evolved D8 grids. In each charge cell reconstruct the same two numerical moments. The highest positive occupied parent shell ends at energy 32.00003906; the highest daughter shell ends at 26.00000017 on both grids. Since (E_mathrm{daughter}=E_mathrm{parent}+omega-E_s), the inverse channel has a parent-equivalent cutoff (26.00000017-(omega-E_s)), distinct from the forward cutoff. No statistical weight is assigned to any allowed angle here.

At q=3.5e−6 the sampled first openings in k are:

| D8 chart | Forward first opening | Inverse first opening | Inverse full cosine measure |
| --- | ---: | ---: | ---: |
| Fine | .00142967983 | .00145924756 | .00148273517 |
| Third | .00143557421 | .00146524552 | .00148881591 |

At k=.00145 and q=3.5e−6, the forward cap allows 88.532% of cosine measure on fine and 62.808% on third, while **both inverse caps give zero**. At q=1.5e−5 and k=.00150 the forward cap is fully open, but the inverse cap still allows only .79144 and .73209 of cosine measure. These statements are about the retained archived charge support on this past-spacelike branch, not counts, rates, unique roots over all k, or the full physical inverse reaction. A 4096-cosine projection and an independent radial energy inversion check these finite-support geometric statements.

Reproduction from `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-cell-endpoint; d=2026-11-28
python -W error early_pole_endpoint_d8_forward_inverse_cap_window.py \
  "$p-fine-d8-x05-state-$d.json.gz" "$p-third-d8-x05-state-$d.json.gz"
```
