# Propagate shell counts on a recorded metric without differentiating a histogram

Use the [direct scalar comoving shells](INDEPENDENT-WARM-DIRECT-SOURCE-CLOCK-AND-PSD-2026-11-26.md)
and the stored q64 spectrum as separate collisionless
test tracers on two prior Newtonian-gauge metric fixtures.
Both initial metric times are after the archived source
interval. Give all the tracers the same approximate early
adiabatic monopole, dipole and quadrupole. This assumption
does not model production perturbations.

The density numerator in the relativistic line-of-sight
formula contains q³ epsilon F(q) f'(q), where F is the
angularly propagated initial/source kernel. Move the
derivative onto q³ epsilon F by parts so the calculation
uses the positive number measure q² f(q)dq. The derivative
of speed with respect to log q is v(1−v²), and its time
integral gives the derivative of the travelled path.
This avoids directly differentiating bin noise. Independently
recover the result from the interior derivative and *both
endpoint jumps* of a piecewise constant q-number density.
Agreement is better than 2.58e−10 in amplitude; finite
differences check the analytic log-q response derivative
to 2.34e−9 relative on the sampled nodes.

| k (Mpc⁻¹), z=20 | Direct lepton scalar | Measured pion scalar | Direct yield-weighted mix | Held q64 scalar |
| --- | ---: | ---: | ---: | ---: |
| 2.6944 | .72132 | .64566 | .71508 | .65673 |
| 6.7360 | .46932 | .39298 | .46303 | .38676 |

The table uses a native CDM amplitude from the *same fixed
fixture*, not a measured transfer function. Refining time
and momentum quadratures changes results by at most
1.35e−5. Holding the same shell counts but using a uniform
log-q instead of q number measure changes them by at most
4.30e−4 in the tested representations. The unresolved
overflow bin above q/T0=18 contains below 1.5e−9 of either
source; this diagnostic places that small bin in [18,36]
and does not assume that interval is its strict support.

On these two modes, the difference between the direct mix
and the held scalar response comes almost entirely from
forcing on the recorded metric: the initial-only amplitude
differences are below 2.21e−7, compared with total scalar
differences .05835 and .07627. Neither is a large
matter response by itself. With the archived stable count
held, the scalar rest fraction is 2.6324e−6 for the archive
and 2.0992e−6 for the direct vacuum source. Treating the
other tracer as CDM, a same-count scalar replacement
changes the normalized matter amplitude by only
1.54e−7 and 2.01e−7 on these fixtures. Replacing the
scalar count with the actual independent vacuum yield as
well gives 3.06e−7 and 4.87e−7 in the same proxy.

These are controlled tracer comparisons. The metric is not
regenerated; neither the dark backreaction nor scalar
production perturbations, structure formation or total
thermal production is solved here.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
D=results
python -W error warm_direct_vacuum_fixed_metric_response_reference.py EOS.dat \
 "$D/warm-multishell-corrected-q64-fine-stage-48-2026-11-26.json.gz" \
 "$D/warm-direct-vacuum-scalar-psd-2026-11-26.json" "$D"/los-metric-off/*.npz
~~~

[Fixed-metric results and independent derivative checks](results/warm-direct-vacuum-fixed-metric-response-2026-11-26.json).

The [follow-up response controls](INDEPENDENT-WARM-DIRECT-SOURCE-RESPONSE-CONTROLS-2026-11-26.md)
compare fixed-peak production clocks, an alternate metric fixture, and removal of the rare low-q population.

The [matched native setup and positive-shell homogeneous bounds](INDEPENDENT-WARM-DIRECT-NATIVE-BASELINE-2026-11-26.md)
separate any subsequent perturbation response from small expansion changes.
