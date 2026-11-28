# Strict thermal vector correction on two fixed invariant-mass strips

Use the author's [spectral master-integral release v1.1-spf](https://doi.org/10.5281/zenodo.10371215)
to evaluate the vector trace at fixed invariant mass in
[Jackson and Laine's convention](https://arxiv.org/abs/1910.09567):
the output columns are $M/T$, $k/T$, the Born trace over $T^2$,
and the strict NLO coefficient over $g_s^2 T^2$.
At each fixed $x=M/T$, let $u=(\omega-M)/T$. Integrate over
$0.0125\le u\le32$ with the same boost weights as the
[thermal Born reference](INDEPENDENT-HARD-THERMAL-QUARK-BORN-2026-11-26.md).
There is no integration over invariant mass or cosmic temperature
in these strip results.

For fixed illustrative $\alpha_s=0.30$, use $g_s^2=4\pi\alpha_s$:

| Fixed M/T | Strict NLO term / thermal Born count integrand, after boost integration | Same ratio for heat |
| ---: | ---: | ---: |
| 4.0 | +0.122094 | +0.118825 |
| 4.5 | +0.118750 | +0.116060 |

For $\alpha_s=0.20$ and 0.40 the strict correction scales linearly
in these fixed-coupling comparisons; this is not a theory uncertainty.
Linear interpolation or dropping alternate boost probes changes
these ratios by less than 1.0e-4 in absolute terms at $\alpha_s=0.30$.
The native Born traces independently agree with the analytic KMS
trace formula. At separate probes $M/T=5.5,7.5,10.5$, the
unintegrated NLO coefficients are also positive in the sampled
points, but no continuum behavior is inferred from that sample.

This is strict perturbation theory at selected timelike points,
without a light-cone matching prescription or lattice calibration.
Do not apply a common +12% factor to the entire thermal source,
a reheating threshold, or the [measured vacuum bands](INDEPENDENT-PARTIAL-HADRON-REHEAT-SOURCE-2026-11-26.md).
The grid is incomplete in invariant mass. Even the full narrow
range $4\le M/T\le4.5$ accounts for only about 7–8% of the
specified hard Born count in the two example radiation windows.
No NLO interpolation or application across that interval is
assumed here. Requiring $(\omega-k)/T\ge1$ retains 96.2% and
99.1% of the Born count on the two strips and changes their
relative NLO corrections to 0.12337 and 0.11906. A tenfold
tighter native tolerance at two $M/T=4$ probes ($u=0.65,2.2$)
changes the strict correction coefficients by at most 3.79e-6
relatively.
These controls do not estimate higher perturbative orders. The interpolation
comparison quantifies only this boost integration at the
stated external native precision, not its truncation uncertainty
or the uncertainty of QCD near the crossover.

## Reproduce

Obtain the upstream archive separately (checked SHA-256
fa19dcc5db840756b2009e31f64089d611ba259ebfe8f428f5cb3787d2da494a)
and GSL development headers. Our adapter does not vendor the
external GPLv3 solver. Put this repository's path in repo,
and the extracted spectral-source path in spec:

~~~sh
g++ -O2 -std=c++11 -I"$spec" -I"$spec/inc" \
    "$repo/kinetics/hard_vector_nlo_probe_driver_reference.cpp" \
    "$spec"/src/*.cpp -lgsl -lgslcblas -lm -o vector_trace_probe
python - <<'PY' >probe-input.txt
import math
for M in (4.0, 4.5):
    for u in (.0125,.06,.15,.35,.65,1.,1.5,2.2,3.1,4.4,6.5,9.,13.,18.,24.,32.):
        print(M, math.sqrt(u*(2*M+u)))
PY
./vector_trace_probe 3e-4 <probe-input.txt >strips.txt
OPENBLAS_NUM_THREADS=1 python -W error hard_vector_nlo_strip_reference.py \
    strips.txt results/hard-vector-nlo-focused-probes-2026-11-26.txt \
    --strict results/hard-vector-nlo-tolerance-check-2026-11-26.txt --eos EOS.dat
~~~

The focused probes use each $M/T\in\{5.5,7.5,10.5\}$ and
$u\in\{0.7,2,5\}$ with the same native tolerance. The archived
numeric probe outputs can be used directly for the boost analysis.
