# Enclose the frozen transverse pole using per-cell charge moments alone

The [optical audit](INDEPENDENT-EARLY-POLE-CELL-ENDPOINT-TRANSVERSE-AUDIT-2026-11-28.md) compares concrete reconstructing shapes. Here ask what the charged cell number and first energy moment guarantee *regardless* of its intrabin shape, for the same collisionless isotropic transverse response at k=.0015. Assume nonnegative charge weight and support within each declared energy interval. The argument does not use the extra Pauli cap, so it also covers positive atomic test measures. This is not a bound on collision kernels that sample a different, higher-order functional of the populations.

Write u=m²/E² and β=k/ω. The by-parts kernel per particle has the convergent positive-coefficient series
```text
h(E,β) = (1/E) sum_{n>=0} β^(2n) (1-u)^n [2+(2n+1)u] / [(2n+1)(2n+3)] .
```
A cell's number and energy fix the average of the affine secant of h. Bound the residual by the usual second-derivative interpolation error. For a term written as F(u)/E, its second derivative is (2F+10uF'+4u²F'')/E³; take explicit interval bounds on u and 1−u in each charged bin. For a tighter signed enclosure, split the **error bound** interval into twelve energy pieces without adding information about the distribution. Bound the residual of each series term about the original cell secant on every piece, then combine coefficient signs for 0≤β≤βmax. A positive lower bound on the zero-β gap gives βmax. The n>32 polynomial/geometric tail is bounded separately.

Solve the monotone pole equation using the resulting uniform field enclosure. The same bound applied to β derivatives gives a conservative separate residue interval; no correlation between its extremes and frequency extremes is assumed. The intervals are evaluated in floating point with checks and a numerical margin, not an interval-arithmetic audit of every rounding operation.

| Case | Relative frequency enclosure width | Relative residue enclosure width |
| --- | ---: | ---: |
| Original fine, one charge interval | 4.36e−4 | 6.52e−4 |
| Original third, one charge interval | 2.51e−4 | 3.65e−4 |
| Original fine, eight charge intervals | 6.46e−6 | 9.30e−6 |
| Original third, eight charge intervals | 3.76e−6 | 5.37e−6 |
| Evolved fine, two charge intervals | 1.06e−4 | 1.54e−4 |

The original exact step response, a continuous shape fit to both moments, and the positive atomic stress test placing each cell's full weight at its mean energy all fall within their applicable intervals. The atomic test is not a physical Pauli-limited state. On the evolved data the number-only step fit does not generally match the energy moment and is reported only as a comparison. The width of these intervals decreases with additional **known subcell moments**, not by imposing a preferred shape. Close agreement of two candidate optical closures is not by itself a rigorous rate bound.

From `kinetics/`, decompress the fine one-interval, fine two-interval and third one-interval x=.05 states and run:
```sh
OPENBLAS_NUM_THREADS=1 python -W error early_pole_endpoint_transverse_moment_envelope.py \
 FINE_X05.json FINE_D2_X05.json THIRD_X05.json
```
See the [receipt](results/early-pole-endpoint-transverse-moment-envelope-2026-11-28.json).
