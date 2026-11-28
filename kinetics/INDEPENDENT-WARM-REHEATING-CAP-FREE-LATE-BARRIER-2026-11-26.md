# A conditional late fermion survival bound without scalar occupation caps

This extends the [isotropic soft-band bound](INDEPENDENT-WARM-REHEATING-ISOTROPIC-SOFT-BAND-2026-11-26.md)
by dropping the hard scalar occupation cap as well. It
still assumes the same conditional anchor, homogeneous
isotropic factorized populations, the resolved vacuum
Boltzmann collision operator, the previously enclosed
charged-lepton energy input and no additional stable
loss channels. The proof uses the ordinary Pauli
phase-space cap for the stable fermion per spin.

Expand the final scalar factors in the gross annihilation
rate as $1+f_{\phi,1}+f_{\phi,2}+f_{\phi,1}f_{\phi,2}$.
The constant term is bounded by the prior vacuum
coefficient. The soft part of the linear terms uses the
[angular bound](INDEPENDENT-WARM-REHEATING-ISOTROPIC-SOFT-BAND-2026-11-26.md)
with $F_h=0$ in its extra coefficient: the product term
is handled independently. Grant the whole available
scalar energy to the soft count for this inequality.

## Hard linear term

Give the stimulated scalar momentum magnitude $k$
and energy $e$ and call the incoming fermion momenta
$\mathbf p_1,\mathbf p_2$. The two identical linear
terms cancel the final-state symmetry division of two.
With $d\Pi_k=d^3k/[(2\pi)^3 2e]$, integrating out
the unstimulated scalar leaves
$2\pi\delta_+[(p_1+p_2-k)^2-m_\phi^2]$.
The spin-averaged bound $|\mathcal M|^2\le K_a s$
and $s\le4E_1E_2$ cancel the denominator $4E_1E_2$.
At fixed first fermion and scalar momentum, angular
integration over the second fermion gives a normalized
delta integral bounded by $1/(4p_2|\mathbf p_1-\mathbf k|)$.
The remaining normalized angular average of
$1/|\mathbf p_1-\mathbf k|$ is $1/\max(p_1,k)\le1/p_1$.
No angular cutoff is imposed on the distribution.

Thus the light count loss from the hard linear term is
at most $(\pi K_a/2) I_h J_+J_-$, where
$I_h=\int_{\rm hard}dn_\phi/e\le\rho_\phi/E_*^2$
and $J_\pm=\int dn_\pm/p$.
For either charge, two spin states with occupation
at most one imply $J_\pm\le C_P n_\pm^{2/3}$,
$C_P=(3\pi^2)^{2/3}/(2\pi^2)$: filling the smallest
momenta first maximizes this decreasing weight
at fixed count. It follows that
$J_+J_-\le C_P^2(n^2/4)^{2/3}$ for the total count density $n$.

## Quadratic term and a barrier without ODE integration

By reversal symmetry, the part of the forward loss
weighted by $f_{\phi,1}f_{\phi,2}$ is bounded by the
vacuum inverse rate on that scalar pair distribution.
The forward light-occupation product becomes a
multiplicative factor no larger than one in the inverse
angular integral. The forward spin-average factor of
four and final identical-particle division of two
become a factor of four and no final symmetry division
in the inverse cross section. Initial identical scalar
pair counting then divides by two; the loss of two
light fermions cancels it. With the spin-summed
bound, the light loss is at most $K_B\rho_\phi^2$.
Indeed the inverse cross section times speed divided
by the two incoming energies is bounded by
$2K_a\sqrt{1-4m_\chi^2/s}/(\pi s)$;
this peaks at $s=6m_\chi^2$, giving
$K_B=K_a/(3\sqrt{3}\pi m_\chi^2)$.
We drop inverse gains from the comparison entirely.

The scalar energy at any instant is no larger than
the total available dark energy minus the stable rest
energy. Write $X=b^3n/n_0$ and
$R=\rho_{\rm cap}/(m_\chi n_0)-X$. Thus
$\rho_\phi\le b^{-3}m_\chi n_0 R$.
On any candidate first descent from $X=1$ to $.8$,
use the largest possible $R$ and $X$ separately.
For the hard term the extra dilution factor $1/b$
relative to the clock $u=1-1/b$ can be discarded.
At the 80 keV anchor with the indicated 4.9 MeV split,
all four loss contributions including the specified
leptons add to at most $.17358$ per unit $u$, for
$0\le u\le1$. Therefore $X$ cannot reach $.8$ and
in fact remains at least $.82642$. In the adopted
conditional normalization this gives
$\Omega_\chi h^2\ge138.48$, with **no pointwise bound
on any scalar occupation**. The same simple barrier
is inconclusive at the warmer 2 MeV anchor.

The result does not justify the 80 keV anchor as a
full kinetic solution, apply to a coherent field or
correlated/anisotropic states, or allow an extra
energy source or stable number-changing mechanism.
The angular identity and inverse energy maximum
are independently checked numerically.

~~~sh
python -W error warm_prompt_reheating_cap_free_late_barrier_reference.py \
 results/warm-prompt-reheating-late-shape-free-floor-2026-11-26.json \
 results/warm-prompt-reheating-post-anchor-scalar-ceiling-2026-11-26.json \
 results/warm-prompt-reheating-cold-lepton-energy-allowance-2026-11-26.json
~~~

[Barrier ledger and independent checks](results/warm-prompt-reheating-cap-free-late-barrier-2026-11-26.json).
