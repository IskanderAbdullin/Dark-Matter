# Detailed balance is a gate, not a rate, for the common charged aperture

Use the [shared bare positive charged cut](INDEPENDENT-FINITE-MASS-QED-SHARED-POLAR-DELTA-2026-11-27.md)
as a common-support diagnostic. At fixed neutral-scalar energy
$E_s$, both legs can be indexed by parent energy $y$: the daughter
has energy $y-E_s$. For the Born process and its inverse, let
$G=f_P(1-f_b)$, $L=f_b(1-f_P)$. Then $G-L=f_P-f_b$, the same
occupation difference occurring in the paired bare polar delta.
This equality does **not** identify the two diagrams with the rates
$G$ and $L$.

As an equilibrium control, use both charged legs at temperature
$T$. Their FD odds give $G/L=e^{-x}$ pointwise, where
$x=(E_s-\mu_P+\mu_b)/T$. For $x>0$, the scalar factor
$n=1/\operatorname{expm1}(x)$ makes $(1+n)G-nL=0$
pointwise. The common-temperature relation tests the gain/loss
part of detailed balance; the nonzero chemical-potential
differences here are algebraic fixtures, not conserved charges
of this model. The integrated ratio is the same for uniform
Born measure and the signed charged kernel with any *common*
regulator tested. It is not evidence of a physical positive
spectral density: the raw signed integrals can both be negative.
Four FD choices on each aperture give relative pointwise balance
errors at most $6.2\times10^{-15}$ and log-ratio errors below
$2\times10^{-15}$.

On the two frozen nonthermal $q=4.1\times10^{-6}$ apertures,
all represented subcells instead have (G>L>0). Unweighted
Born log loss/gain is (-1.250693) (fine) and (-0.836526)
(third). No nonnegative scalar occupation can neutralize this
channel alone, either pointwise or with the positive Born
measure. In fact its signed net divided by forward-plus-reverse
turnover is bounded below by $(G_{\rm int}-L_{\rm int})/(G_{\rm int}+L_{\rm int})$ for all nonnegative $n$ at this
fixed scalar momentum. This statement concerns a frozen
single channel: it does not decide the sign of the full
simultaneous collision operator.

The full bare charged kernel is positive near coincidence but
negative for $|y-y_0|>4.129\times10^{-6}$, compared with an
aperture width of about 0.0135. If it is used as a weight for
these statistical products, **both raw integrals are negative**
at guards $10^{-5}$ through $10^{-12}$ on all nine events.
Reversing the overall sign would not give a nonnegative
differential kernel either. The two one-sided archive cells
control each formal singular log coefficient. As a check on
the wide sign separation, extrapolating the exact step-cell
primitives would put the raw zeros below $10^{-468}$ in
archive units. Those numbers are not physical width estimates
or reliable extrapolations of the reconstructed spectrum. Direct
log-radial quadrature agrees within $8.0\times10^{-13}$
absolute with the closed archive integrals at guard $10^{-9}$.

The appropriate next gate is positivity and detailed balance of
a **matched** scalar self-energy/collision kernel including the
corresponding real, virtual and residue contributions, not
assigning this bare absorptive correction a standalone rate.
No physical neutral-scalar distribution, thermal width, or
cosmology prediction follows.

Replay with the same pinned archive spectra from `kinetics`:
```sh
OPENBLAS_NUM_THREADS=1 python -W error finite_mass_qed_aperture_statistics_reference.py
```
[Receipt](results/finite-mass-qed-aperture-statistics-2026-11-27.json).
