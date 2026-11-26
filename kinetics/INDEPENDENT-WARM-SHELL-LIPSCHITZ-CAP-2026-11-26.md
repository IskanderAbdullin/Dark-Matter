# A conditional regularity test for the light occupation cap

The [sharp two-moment audit](INDEPENDENT-WARM-SHELL-OCCUPATION-FEASIBILITY-2026-11-26.md)
constructs narrow Pauli-saturated profiles with the same two
saved light moments as very dilute profiles. Ask what explicit
additional information would suffice to rule out peaks of 1/2
without assuming a particular reconstructed shape.

For a shell ([a,b]), assume the true isotropic radial
occupation is (L)-Lipschitz **in physical momentum**:
(|f(p)-f(q)|le L |p-q|) inside this shell. A peak (f(x)ge c)
implies (f(p)ge[c-L|p-x|]_+). For every level between 0 and
(c), this cone intersects the shell in an interval of length
at least as large as the cone centered at (a). Among intervals
of that length, weighting with (p^2dp) makes the innermost
one cheapest. Integrate the levels to obtain the sharp
number-only criterion
[
 N ge N_{min}(L)
   = int_a^b p^2[c-L(p-a)]_+,dp.
]
The one-sided ramp at the inner edge attains this lower bound.
A smaller measured shell number excludes such a peak if the
slope bound holds. The actual kinetic energy is not imposed
here, so it could strengthen the test. No continuity across
shell edges is required. The (p^2) weight is important at
the origin: for a cone of width (c/L) in a shell with
(a=0), the allowance is (c^4/(12L^3)), not proportional
to (1/L).

For (c=1/2), the least critical slope over all shells and
all six matched saved times is (1.35244	imes10^5) GeV⁻¹
(q16) or (1.36684	imes10^5) GeV⁻¹ (q32). Both occur at
the innermost represented shell of the first continued frame.
These thresholds are 6.387 and 6.455 in unit-length linear
shell coordinates. A common physical bound strictly below
the corresponding slope would therefore exclude peaks of
1/2 throughout the represented light shells **at these
saved times**, without assuming the reconstructing ansatz.
For reference, the within-shell slopes of the archive's
specific entropy reconstruction have conservative
upper/critical ratios below (3.9	imes10^{-7}) (q16)
and (1.1	imes10^{-7}) (q32). This is not a derivative
estimate for any unknown physical distribution.

The missing step is proving a light-sector bound on
momentum regularity, or directly bounding peaks, under
the actual inverse, elastic, gain and loss dynamics.
The calculation does not cover exterior tails, the heavy
state, hidden intermediate times or the full future history.
Numerical verification includes exact massless controls and
independent quadrature of cones centered elsewhere in a shell.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error warm_shell_lipschitz_cap_reference.py
~~~

[Diagnostics](results/warm-shell-lipschitz-cap-2026-11-26.json).
