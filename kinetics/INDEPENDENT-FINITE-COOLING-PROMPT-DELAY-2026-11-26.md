# Testing the prompt-decay approximation in the finite cooling reference

Same 100 MeV vacuum parent, light fermion and 86 eV scalar as in the
[finite cooling-clock calculation](INDEPENDENT-FINITE-COOLING-REHEAT-PORTAL-2026-11-26.md).
Only the lifetime and daughter momentum shift are bounded here.

## Vacuum clock and momentum bounds

For a parent born with energy $E_i$, freely redshifting momentum only
increases its decay hazard relative to $\Gamma M/E_i$. For an
interval $\delta$ in $\ln a$, the elapsed time is at least
$\delta/H_{\max}$ on the prescribed clock. Thus survival up to that
interval is at most $\exp[-A/(E_i/M)]$, with
$A=\Gamma\delta/H_{\max}$. Averaging the pointwise bound
$\exp(-z)\le(n/e)^n z^{-n}$ over the positive parent birth
measure gives bounds directly from its second and fourth energy
moments, with no arbitrary rare-parent momentum cut. Parents born
within the last $\delta$ of $\ln a$ before the 1 MeV snapshot
are conservatively counted separately.

At fixed angle in the rest frame, the daughter comoving momentum
transverse to the parent trajectory is $aP\sin\theta$; its
longitudinal momentum is
$P\cos\theta\sqrt{a^2+(k/M)^2}+(k/M)E_*$, where $k$ is the
parent's comoving momentum. Consequently a short-delay decay changes
the daughter $p/T$ at 1 MeV by no more than
$(P/T_{\rm snap})(e^\delta-1)$ in norm compared to prompt
production, for the same rest-frame angle.

For the ratio-5 scenario, the minimum width-to-Hubble ratio is
2.54×10¹¹. At a delay of 10⁻⁹ in log scale factor, the fourth-moment
bound plus the very-late birth contribution leaves at most
8.28×10⁻⁹ of parents undecayed by the snapshot or outside the
short-delay class. Within that class both daughters' absolute
momentum shifts are below 4.81×10⁻⁸ in units of the snapshot T.
For the instantaneous case the analogous parent fraction is below
2.71×10⁻¹⁹. These are bounds for the resolved birth measure;
they check a vacuum lifetime approximation only, not in-medium effects,
availability of inverse channels, or neglect of the omitted processes.

## Reproduction

Use the same EOS input and JSON as the previous finite-clock calculation.

~~~sh
OPENBLAS_NUM_THREADS=1 python -W error finite_cooling_prompt_delay_bound_reference.py \
    EOS.dat results/finite-cooling-reheat-portal-2026-11-26.json
~~~
