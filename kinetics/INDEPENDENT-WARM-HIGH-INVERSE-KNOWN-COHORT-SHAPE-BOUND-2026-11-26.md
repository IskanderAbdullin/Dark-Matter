# Bound the initial high inverse source from known cohorts without an intra-cohort shape

One-momentum [cohort activation](INDEPENDENT-WARM-REACTIVATE-KNOWN-COHORTS-2026-11-26.md)
did not constrain other distributions with the same birth
moments. At the 10.391 MeV anchor, bound the bare (L+phi	o H)
inverse energy rate for (p_Hge1) GeV from **those known
stored cohorts only**, conditional on the stated admitted support.

With the outer edge doubled to .604500 GeV, the other
incident particle's maximum energy leaves a minimum of about
.400 GeV for **each** incoming leg to cross the heavy cut.
This exceeds the q32 core edge, so both particles have to
come from the old exterior under this support hypothesis.
Use each cohort's proper collisionless birth redshift to
express the new minimum energy at its own birth. The
number above it is at most the observed birth energy
excess divided by the gap from the birth floor.

The full isotropic parent collision measure has coefficient
(C=g_HGamma M^2/(4pi^2P_*)) on
(dE_H,dE_phi=dE_L,dE_phi). Drop the exact allowed
angular restriction and final heavy Pauli factors upwards.
For each leg above threshold, the energy integral of
(f,dE) is at most (2pi^2/g) times the eligible number
divided by the minimum (pE) there. Multiplying these
two upper integrals and the maximum parent energy gives
an initial high-inverse energy bound of (3.78	imes10^{-19}),
in units of current retained dark energy times Hubble.
Different admitted support gives the documented different
upper limits, not the same universal bound.

As controlled examples, put two moment packets in each
cohort, with the upper packet at .52 or .58 GeV after
redshift. The lower packet sits above the birth floor and
the weights exactly preserve each birth number and mean
energy. Represent their resulting shell moments positively.
These synthesized profiles have actual initial rates of
(1.49	imes10^{-20}) and (6.64	imes10^{-21}) on the
same normalized scale, both below the moment-only bound
but far above the earlier empty-padding prompt sample.
These are alternatives, not physical estimates of the old tail.

This bound is for activation at one anchor with the stated
finite support and selected inverse reaction. It gives no
bound on untracked or beyond-support tails, regeneration,
extra visible or dark heavy sources, or cumulative use.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error warm_high_inverse_known_cohort_shape_bound_reference.py EOS.dat
~~~

[Diagnostics](results/warm-high-inverse-known-cohort-shape-bound-2026-11-26.json).
