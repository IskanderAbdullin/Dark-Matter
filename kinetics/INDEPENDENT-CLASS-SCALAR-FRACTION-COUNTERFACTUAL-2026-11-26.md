# Counterfactual scalar fractions for the same archived momentum shape

In the [low-reheat tracer calculation](INDEPENDENT-CLASS-SCALAR-TRACER-2026-11-26.md),
the retained 86 eV scalar contributes only $2.53\times10^{-9}$ of
dark density. To separate the effect of its momentum distribution
from its tiny abundance, hold the shape fixed and assign it
counterfactual fractions between 0.2% and nearly 100%. This is
an artificial population assignment at fixed total dark density,
photon temperature and primordial spectrum. None of these
fractions is an output of the production calculation.

Use synchronous unmodified CLASS, with the same compensated
positive scalar PSD and the same adiabatic initial prescription.
The last case leaves $\omega_{\rm CDM}=10^{-8}$ to retain
the synchronous CDM anchor.

| Fraction of dark density reassigned to scalar | $P/P_{\rm CDM}$ at k=10 h/Mpc, z=0 | First bracketed matter-power half, z=0 (h/Mpc) |
| ---: | ---: | ---: |
| 0.2% | 0.98236 | Not in range |
| 1% | 0.91466 | Not in range |
| 10% | 0.40025 | 3.120 |
| 50% | Not quoted in the strongly suppressed region | 0.2252 |
| Nearly all | Not quoted in the strongly suppressed region | 0.1576 |

The expression $\ln(P/P_{\rm CDM})/f_\phi$ at k=10 h/Mpc
is −8.900 and −8.920 for the 0.2% and 1% cases at z=0.
This is a finite numerical response, not a direct computation at
the true abundance. The extreme suppressed tails of the
most scalar-rich runs are not convergence-certified and are
not used to support a physical exclusion.

Write $B$ for baryonic density and $D$ for total dark density
in the same units, and use a superscript 0 for the all-cold
comparison. The two factors are
$$
G=\frac{B\delta_b+D\delta_c}{B\delta_b^0+D\delta_c^0},
\qquad S=1+\frac{f_\phi D(\delta_\phi-\delta_c)}
                         {B\delta_b+D\delta_c}.
$$
At 1%, k=10 h/Mpc, z=0 the cold+baryon growth amplitude
(normalized to the cold comparison with the dark density
provisionally assigned to its cold contrast) is 0.95993.
Replacing the reassigned density by its actual scalar contrast
adds a weighting factor of 0.99630. Squaring the product gives
the power ratio 0.91466. Thus it would miss most of the effect to
weight the dilute scalar tracer by its fraction **without**
letting the cold+baryon component respond. This decomposition
uses late synchronous density contrasts and is compared with
the native total linear power only where $P/P_{\rm CDM}>0.05$.

## Independent raw momentum checks

Using the original unsmoothed file, 160 nodes with qmax 63.8377,
and the later fluid switch, the 1% and 10% tests give z=0
power ratios 0.914686 and 0.400350 at 10 h/Mpc. Across seven
probes from 0.1 to 10 h/Mpc and redshifts 0, 5, 20, the
maximum relative difference from the smooth test is 0.0392%
and 0.434%, respectively. The 10% half-power scale at z=0
moves from 3.12014 to 3.12082 h/Mpc. These comparisons vary
the numerical settings and the smoothing together.

Near the first half-power crossing, the nearly-all-scalar
case shifts from 0.157634 to 0.156522 h/Mpc at z=0 using
the independent raw, short-k run. The relative shift is
about 0.7% at all three redshifts. The 50% case shifts from
0.225237 to 0.225555 h/Mpc at z=0 (0.14%). Shape- and
grid-dependent differences elsewhere around the cutoff are
larger than the change in its half-power location.
The short raw tests are not
used to validate the much smaller, uncontrolled high-k power.

These are illustrations of a fixed shape, not a statement about
why this many scalars would be produced, observed Lyman-alpha
limits, a galactic fit, or the true nonlinear matter spectrum.
The adiabatic backward extrapolation remains a modeling choice
rather than a solution for production perturbations.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error class_archived_scalar_fraction_counterfactual_reference.py native fine --prepare
for c in 000 0002 001 010 050 all; do
  f=scalar_fraction_"$c"_smooth
  /path/to/class "native/$f.ini" >"native/$f.log" 2>&1
done
python -W error class_archived_scalar_fraction_counterfactual_reference.py native fine \
    --codes 0002 001 010 050 all
# After corresponding raw native configs have also run:
python -W error class_archived_scalar_fraction_counterfactual_reference.py native fine \
    --codes 001 010 --raw-codes 001 010
python -W error class_archived_scalar_fraction_counterfactual_reference.py native fine \
    --large-fraction-controls 050 all
~~~
