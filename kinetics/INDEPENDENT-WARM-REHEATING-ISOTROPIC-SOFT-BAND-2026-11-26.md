# Remove the pointwise soft cap for isotropic stable-pair loss

The [two-band bound](INDEPENDENT-WARM-REHEATING-COLD-BAND-LOSS-2026-11-26.md)
uses a finite occupation cap even in the soft scalar band.
Here a different argument uses *only the number of soft scalars*
and makes no assumption about their pointwise occupation.
It assumes homogeneous isotropic and factorized incoming
stable-fermion distributions, resolved scalar phase-space
occupation, and the usual local vacuum Boltzmann operator.
It still needs a uniform bound $F_h$ on occupation in the
hard band at all times. Neither the kinetic continuation nor
an actual coherent condensate is modeled by this argument.

## Angular suppression at small outgoing momentum

Use $e$ and $k=\sqrt{e^2-m_\phi^2}$ for the energy
and momentum of the soft daughter; let $e<E_*$.
The incoming energies and momenta are $E_i,p_i$.
For total energy $E=E_1+E_2$ and total momentum
magnitude $P$, the outgoing momenta obey
$$
 E-P \le 2e\le 2E_*.
$$
Indeed $P\ge P_h-k$ by the triangle inequality,
and $E_h-P_h\le e-k$ because the other daughter has
$E_h>e$ and both have the same mass. Consequently
$E-P\le e+k+(E_h-P_h)\le2e$.
On the incoming side
$$
 E-P\ge \sum_i(E_i-p_i)
       =\sum_i m_\chi^2/(E_i+p_i).
$$
Each incoming energy is therefore at least
$E_{\min}=m_\chi^2/(4E_*)$, not just one high-energy leg.
Use $E_*<m_\chi/4$ so that $p_i\ge r E_i$, where
$r^2=1-m_\chi^2/E_{\min}^2>0$.

The boost relation gives $(s-2Ee)^2\le 4k^2(E^2-s)$.
Since $P\le E$, the support of $s$ lies inside
$[2Ee-2Ek,\,2Ee+2Ek]$, with width at most $4Ek$.
Also $s\le4Ee$, and the previous gate $E\ge m_\chi^2/E_*$
gives $P\ge L E$ for $L=\sqrt{1-4E_*^2/m_\chi^2}$.
At fixed incoming magnitudes, isotropy gives the normalized
relative-angle measure $ds/(4p_1p_2)$.

In the pair CM the energy exchanged on either Yukawa
fermion line vanishes. In the standard Hermitian norm the
operator for one such line has norm
$1/\sqrt{M^2+|\mathbf q|^2}\le1/M$: the spatial slash is
anti-Hermitian and its square is $-|\mathbf q|^2 I$.
Summing the two lines with spinor norm product $\sqrt{s}$
bounds the spin-averaged squared amplitude by $K_a s$,
with $K_a=4y^4/M^2$.
With the identical-final factor already included, the
lab cross section times Møller speed per CM polar cosine
is bounded by $K_a s\,v_\phi/(128\pi E_1E_2)$,
after averaging over azimuth; $v_\phi=\sqrt{1-4m_\phi^2/s}$.
The derivative from lab daughter energy to that cosine is
$2/(P v_\phi)$. Count either soft daughter (multiply by
two), and integrate over the allowed relative-angle
interval. The upper differential coefficient in $de$
is
$$
 \frac{K_a E e k}{8\pi L E_1 E_2 p_1 p_2}.
$$
The critical factor is $e k\,de=k^2\,dk$: this is exactly
the radial number-density measure, not an integral of
occupation against $de$ alone. It suppresses the vanishing
angular support as $k\to0$, so the bound stays finite even
for arbitrarily narrow resolved soft distributions.

After integration against the continuum occupation and the
single hard daughter factor $1+F_h$, the extra
angle-averaged coefficient is at most
$$
 (1+F_h)\frac{\pi K_a\,n_{\rm soft}}{2r^2 L E_{\min}^3}.
$$
Here $E/(E_1E_2p_1p_2)\le 2/(r^2 E_{\min}^3)$.
The baseline remains $k_0(1+F_h)^2$;
the replacement uses $f_{\rm soft}$ itself as an upper on
the nonnegative *extra* product, without assuming its supremum.

## Numerical comparison

Grant the entire allowed comoving dark energy to the count
$n_{\rm soft}\le\rho_{\rm cap}/m_\phi$ for the whole future,
even though the physical density dilutes. Keep the hard and
three-lepton comparisons above and solve the same count
inequality. For the modeled visible 80 keV anchor, with
$F_h=1$, the extra dimensionless quadratic exposure is:

| Hard/soft split | Minimum energy of **each** incoming stable fermion | Extra exposure upper | Stable $\Omega h^2$ lower |
|---:|---:|---:|---:|
| 1 keV | 100 GeV | 1.10e-16 | 130.86 |
| 10 keV | 10 GeV | 1.10e-13 | 130.86 |
| 100 keV | 1 GeV | 1.10e-10 | 130.86 |
| 1 MeV | 0.1 GeV | 1.15e-7 | 130.86 |

The 2 MeV anchor is less constraining with the conservative
all-pair vector comparison, but still gives .284 at the
1 MeV split and $F_h=1$. Controls compare the analytic
support interval and angular integral against independent
quadrature, and compare the spin-norm ceiling against the
exact numerical trace. This does not prove that $F_h=1$
is dynamically maintained. Coherent fields, anisotropic or
correlated beams, other stable loss channels, and initial
states with different dark energy are outside the proof.

~~~sh
python -W error warm_prompt_reheating_isotropic_soft_band_reference.py \
 results/warm-prompt-reheating-late-shape-free-floor-2026-11-26.json \
 results/warm-prompt-reheating-post-anchor-scalar-ceiling-2026-11-26.json \
 results/warm-prompt-reheating-cold-lepton-energy-allowance-2026-11-26.json
~~~

[Full comparison and angular controls](results/warm-prompt-reheating-isotropic-soft-band-2026-11-26.json).
