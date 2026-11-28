# Actual G8 continuation through x=.02

From the exact accepted x=.01 state, G8 completes the paired 1024/2048
continuation. The scalar yield is 1.02181406787e−5, the relative paired
yield drift is 2.1847424e−6, and the maximum absolute continuity defect
is 1.69e−13. The new prefix independently verifies both earlier states.

At the same native S4 measurement mass, G8 total pressure is 0.0478002%
lower than G7. Relative total count differs by −2.0236e−9. The
between-cell count and per-particle terms are −0.0454695 and
−0.00233072 percentage points, respectively. Its negative pressure
sign is unchanged from the two earlier matched clocks. This is a
numerical result on the reported continuation, not a continuum
extrapolation or a full G8 history.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; d=2026-11-28
python -W error early_pole_endpoint_warm_g7_prefix_check.py --levels 7 8 \
 --native "$p"-{start,x01,x02}-state-$d.json.gz --native-chain "$p-chain-$d.json" \
 --reports6 "$p"-gauss7warm-{start,x01,x02}-$d.json.gz --states6 "$p"-gauss7warm-{start,x01,x02}-state-$d.json.gz \
 --reports7 "$p"-gauss8warm-{start,x01,x02}-$d.json.gz --states7 "$p"-gauss8warm-{start,x01,x02}-state-$d.json.gz
```
