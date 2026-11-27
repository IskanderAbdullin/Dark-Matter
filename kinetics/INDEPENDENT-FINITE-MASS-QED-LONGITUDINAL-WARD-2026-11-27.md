# Bare longitudinal checks for the massive scalar-current vertex

The [inclusive vacuum reference](INDEPENDENT-FINITE-MASS-QED-INCLUSIVE-VACUUM-2026-11-27.md)
combines the virtual Feynman-gauge amplitude, its external pole
residues, and real radiation. Here we check what the covariant
gauge-fixing term, proportional to K_mu K_nu, does before
the common loop integration. This is the free-propagator application
of the [Ward–Takahashi identity](https://doi.org/10.1007/BF02832514).

Let D_a(P) = slash(P)−m_a, S_a(L) = D_a(L)^−1 and L=P−K.
Then slash(K)=D_a(P)−D_a(L) on either line, even when the
two lines carry different equal-charge masses. Between external
on-shell spinors the longitudinal vertex reduces to the Born
scalar insertion: each neighboring internal propagator cancels
against an inverse propagator and the two resulting minus
signs multiply to plus.

The longitudinal external self-energy on either leg is

    slash(K) S_a(P−K) slash(K)
      = D_a(P) S_a(P−K) D_a(P) − D_a(P) − slash(K).

For a common even loop weight, the odd last term integrates
to zero. The first term projected onto the on-shell spinor is
quadratic in the off-shell displacement, away from regulated
singularities. The linear term supplies the common residue
coefficient −1, and the two half-residues cancel the longitudinal
vertex coefficient +1. Fixed-spatial-momentum derivatives in both
the parent rest frame and a boosted frame recover this cancellation
to 1.3e−10 after a two-step extrapolation. The unprojected matrix
identity agrees to 2.9e−14 absolute across the tested loop samples.

For a bare charged cut put N(L)=slash(L)+m_a on either energy
branch. The cut inverse annihilates N(L); at the current
vertex the adjacent longitudinal contraction therefore vanishes
between the external spinors for either of the two single charged
cuts. On a cut self-energy line the remaining numerator is
D_a(P) N(L) D_a(P). Its algebraic projection and first derivative
at fixed regular loop denominator vanish. The tested positive
and negative branches have projected vertex and self-energy
numerator defects below 3e−15 absolute.

The statements concern bare algebra with a common regulator.
They do not prescribe a thermal double-pole distribution, justify
the exchange of a singular limit with a derivative, or include
physical plasma longitudinal modes, which are not the covariant
gauge-fixing term. Real-time weights, spectral residues, and
medium real/virtual matching still need to be fixed.

Replay from kinetics:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_longitudinal_ward_reference.py
```
[Raw numerical receipt](results/finite-mass-qed-longitudinal-ward-2026-11-27.json).
