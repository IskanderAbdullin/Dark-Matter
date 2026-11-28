# An imposed expanded clock isolates the effect of velocity rescaling

This extends the frozen ballistic stress test to an **imposed** conformal interval with a linear scale-factor clock. At each momentum, average the free speed $q/\sqrt{q^2+(M a)^2}$ uniformly over $a\in[1,A]$. The exact average is

$$
\bar v(q)=\frac{q[\operatorname{asinh}(MA/q)-\operatorname{asinh}(M/q)]}{M(A-1)}.
$$

An independent numerical clock integral checks this on each chart to 5.9e-15 relative. No background parameters are fitted. Use the same scalar distributions and common initial formal measurement mass as the original ballistic test.

The counterfactual equal-rms distance uses a chart-specific phase $b=u/\sqrt{\langle \bar v^2\rangle}$. For G3 and G4 evaluate each distribution at **both** chosen phases. Averaging the two routes through this 2×2 table exactly decomposes the difference into a fixed-phase shape term and a phase-rescaling term. The numbers at $u=1$ are:

| A | G4−G3 equal-rms loss | Fixed-phase part | Phase-rescaling part |
| ---: | ---: | ---: | ---: |
| 1 | +1.09983e-5 | −2.54975e-9 | +1.10008e-5 |
| 2 | +1.25895e-5 | −2.60738e-9 | +1.25921e-5 |
| 10 | +5.39229e-6 | −2.77944e-9 | +5.39507e-6 |
| 100 | +2.50829e-7 | −2.84805e-9 | +2.53677e-7 |
| 1000 | +8.93856e-9 | −2.84741e-9 | +1.17860e-8 |

The fixed-phase shape difference remains negative and small throughout this toy check, whereas the apparently much larger positive equal-rms result comes predominantly from the separate chart phase choices. Its size is strongly sensitive to the chosen expansion interval. Holding the phase common to both charts avoids that interpretation error. With no metric forcing, collision term or physical scale matching, this is not a cosmological transfer calculation.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; d=2026-11-28
python -W error -m unittest test_early_pole_endpoint_radiation_clock_streaming_probe
python -W error early_pole_endpoint_radiation_clock_streaming_probe.py \
 --native "$p"-x05-state-$d.json.gz --gauss "$p"-gauss{3,4}warm-x05-state-$d.json.gz
```
