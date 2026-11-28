# A daughter yield does not by itself determine the sign of Bose amplification

The [paired early Hartree run](INDEPENDENT-EARLY-POLE-PAIRED-HARTREE-MATRIX-2026-11-26.md)
keeps nonthermal fermion occupation in its own
bins. To check how much of that population lies
in relevant states, use the unweighted overlap in
energy for a single scalar mode, not a number-
weighted average over all fermions.

For an on-shell parent window $[L,H]$, measured
in temperature units, let
$I=\int_L^H f_P(E_P)dE_P$ and
$B=\int_{L-E_\phi}^{H-E_\phi}f_b(E_b)dE_b$.
The Bose-growth coefficient is proportional to
$I-B$. The spontaneous seed is instead
$\int_L^H f_P(E_P)[1-f_b(E_P-E_\phi)]dE_P$.
Even at $B=I$, the initial-parent FD monotonicity
implies that the seed is at least
$1-f_P(L)$ times its empty-daughter value.
Thus there is no contradiction if the growth
switches sign while the seed persists.

For a check independent of the ODE integrator,
bound the minimum total fermion count required
to obtain $B\ge I$. The count per $dE_b$
has weight $g_b E_b\sqrt{E_b^2-b^2}/(2\pi^2)$,
which increases in $E_b$. Of all allowed
$0\le f_b\le1$ distributions, filling the
lowest $I$ of the window at $f_b=1$
therefore achieves the required unweighted overlap
at least count. Integrating this extremum gives
$$
 Y_{b,\min}=\frac{g_b}{6\pi^2 S}
 \left\{[(L-E_\phi+I)^2-b^2]^{3/2}
       -[(L-E_\phi)^2-b^2]^{3/2}\right\},
 \qquad S=s_\mathrm{vis}/T^3.
$$
This is necessary, not sufficient, and is
independent of any chosen scalar occupation.
For the previously held $x=.05,m_\phi/T=.0011285$
clock at $p/T=4.1\times10^{-6}$, the sign
threshold is average daughter occupation .24820
over the accessible window and minimum total
daughter yield $1.0533\times10^{-5}$ for
the four-state convention. Direct numerical
integration checks the rearranged count.

At the dominant count mode of the actual fine
coupled endpoint, $B/I=.529$ and the seed is
.892 of the empty-daughter seed. The paired
daughter count exceeds the necessary minimum
there, but is not sufficiently concentrated to
reverse the gain at that mode. The analogous
coarse values are .562 and .884, so the exact
local overlap is less converged than the total
yield. The weakly populated higher-momentum
vacuum pivot can have $B/I>1$ while the seed
fraction remains close to one; modes with
negligible or uncovered parent support are
not assigned a sign classification.
No assumed thermal shape for the daughter is
used to claim these overlap results.

~~~sh
python -W error early_pole_daughter_pauli_bound_reference.py "$RUN_DIRECTORY"
~~~

[Endpoint overlap and necessary-count ledger](results/early-pole-daughter-pauli-bound-2026-11-26.json).
