# Separate soft scalar occupation from stable-fermion depletion

The [late stable floor](INDEPENDENT-WARM-REHEATING-POST-2MEV-ELASTIC-2026-11-26.md#a-conditional-stable-remaining-count-floor-with-arbitrary-elastic-redistribution)
grants its single occupation cap to both outgoing scalars
at all energies. This can grossly overstate the effect of
a low-energy population. Keep the same modeled anchors and
vacuum loss channels, but assume two caps: $F_s$ for
$E_\phi<E_*$ and $F_h\le F_s$ for $E_\phi\ge E_*$.
Both are pointwise bounds in the plasma frame, assumed at
all later times. Neither is derived from the two moments
of the recorded scalar distribution.

## Kinematic gate and comparison equation

Write $E$ for the total pair's energy, $s$ for its
invariant mass squared, and $\beta=\sqrt{1-s/E^2}$.
In the lab, the two scalar energies are
$E[1\pm\beta\sqrt{1-4m_\phi^2/s}\cos\theta]/2$.
Their maximum is at least $E/2\ge m_\chi$, so for
$E_*<m_\chi$ they cannot both be in the soft band.
Their minimum is at least
$E(1-\beta)/2=s/[2E(1+\beta)]\ge m_\chi^2/E$.
Thus a soft daughter costs pair energy at least
$E_g=m_\chi^2/E_*$, independently of scattering angle.
For example, production of a scalar below 10 keV
requires at least 40 GeV of lab pair energy here,
in the frame where the occupations are defined.
This is a conservative necessary
gate, not an estimate of the actual population of such pairs.

The two Bose factors multiply to at most
$(1+F_h)^2+(1+F_h)(F_s-F_h)\,1_{E\ge E_g}$.
If $n=n_++n_-$, the integrated opposing-charge pair
density is at most $n^2/4$. Its energy-weighted version
is $n_+\rho_-+n_-\rho_+\le n\rho_{\rm dark}$, assuming
the standard two-body Boltzmann collision integral.
Consequently the soft excess contributes at most
$2k_0(1+F_h)(F_s-F_h)n\rho_{\rm dark}/E_g$
to the loss of $n$, using the existing all-invariant
vacuum coefficient $k_0=y^4/(8\pi M^2)$.
The hard baseline contributes at most
$k_0(1+F_h)^2 n^2/2$.
The factor of 2 counts the loss of both fermions per event.

For the three charged-lepton final states, instead of
allowing every incoming pair to sit at the vector pole
independently of the available energy, divide the prior
vacuum rate bound by $E\ge\sqrt{s}$. The maximum of
$\sqrt{s}/[(s-V^2)^2+(V\Gamma)^2]$ is at
$s=[p+\sqrt{4p^2+3w^2}]/3$, where $p=V^2$ and $w=V\Gamma$.
Use $2k_E n\rho_{\rm dark}$ for the vector loss, where
$k_E=7.77837\times10^{-11}\ {\rm GeV}^{-3}$ includes all
three leptons. Take the stronger of this comparison
and the independent one using the old uniform *vector*
rate with the new scalar band cap. Retain the previous
uniform *scalar* cap comparison as a separate check.

With $b=a/a_{\rm anchor}$, the
[allowed comoving energy](INDEPENDENT-WARM-REHEATING-POST-ANCHOR-SCALAR-CEILING-2026-11-26.md#allow-the-specified-vacuum-charged-lepton-energy-input)
bounds $\rho_{\rm dark}b^3\le\rho_{\rm cap}$ even with
the specified positive vacuum lepton input. The photon
floor $H\ge H_\gamma/b^2$ then gives
$dX/du\ge -A X^2-B X$ for
$X=b^3n/n_{\rm anchor}$ and $u=1-1/b$. Here
$A=k_0(1+F_h)^2n_{\rm anchor}/(2H_\gamma)$ and
$B=(2\rho_{\rm cap}/H_\gamma)[k_E+k_0(1+F_h)(F_s-F_h)/E_g]$.
Reverse pair reactions can only add stable count.
The final retained fraction is at least
$e^{-B}/[1+A(1-e^{-B})/B]$, with the continuous
limit $1/(1+A)$ at $B=0$.
No fixed fermion or scalar shape is used after the anchor.

## Conditional result at the 80 keV visible anchor

| Assumed hard cap $F_h=1$ | Pair-energy gate | Assumed soft cap $F_s$ | Stable $\Omega h^2$ lower |
|---|---:|---:|---:|
| Split at 1 keV | 400 GeV | 10,000 | 108.25 |
| Split at 10 keV | 40 GeV | 10,000 | 18.59 |
| Split at 100 keV | 4 GeV | 100 | 108.46 |

The 10 keV case needs a soft cap of at least about
$3.49\times 10^4$ before this comparison can reach .12,
instead of about 170 if the same cap is assigned to
all modes. These are necessary values for the bound
to become inconclusive, not sufficient values for
depletion. Under the 80 keV energy budget, the vector
term has linear exposure .0799; the prior all-pair
quadratic exposure was 9.26. The uniform vector
comparison is stronger for the warmer illustrated
2 MeV anchor, where the energy budget is much larger.

The result remains conditional on the recorded anchor,
the caps, ordinary photon/entropy evolution, the
specified charged leptons and the selected stable
number-changing network. It does not bound a singular
condensate occupation, an unconstrained hard band,
an extra source of energy, other loss channels, or
hadronic rates. Large occupation in a cold component
is not by itself evidence for rapid conversion of
the massive stable population.

~~~sh
python -W error warm_prompt_reheating_cold_band_loss_reference.py \
 results/warm-prompt-reheating-late-shape-free-floor-2026-11-26.json \
 results/warm-prompt-reheating-post-anchor-scalar-ceiling-2026-11-26.json \
 results/warm-prompt-reheating-cold-lepton-energy-allowance-2026-11-26.json
~~~

[Two-anchor bounds and kinematic controls](results/warm-prompt-reheating-cold-band-loss-2026-11-26.json.gz).
The driver tests 16,000 seeded massive boosts, a
near-forward configuration within .00044% of the
analytic gate, the vector maximum on an invariant
grid and the count comparison against a direct ODE.
