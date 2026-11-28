# Enclose and resolve the direct cold scalar yield at the warm starting temperatures

The simpler [off-pole bound](INDEPENDENT-WARM-LATE-COLD-VISIBLE-OFFPOLE-2026-11-26.md)
still overstates the contribution above 25 MeV because
it does not require a compatible invariant mass and outgoing
angle in each total-energy interval. Use two independent
approaches on the same direct thermal charged-lepton,
vacuum heavy-pair model and the same equality cut p/T<.1.

## A positive kinetic fold

At each total incoming energy Q and invariant mass R,
the outgoing parent energy runs over Q/2 ± D, with
D=sqrt(Q²−R²) sqrt(1−4M²/R²)/2. Restrict invariant mass
to values permitting E_parent above the cold recoil gate,
and split its boundaries at the exact incoming mass
thresholds. For each parent energy use the exact open
unpolarized two-body probability to emit inside the cut.
Fold both the FD angular moments of the incoming bath and
their Maxwell–Boltzmann comparator. Integrating in the
alternate invariant-mass-first source chart gives the
same two local FD rates within 0.0261% at 35 and 50 MeV.
The two temperature/kinetic quadratures agree within
0.036% on the integrated FD values in this comparison.
The remaining incoming high-energy tail is guarded by
the already conservative vacuum cross-section envelope.

## A separate one-sided enclosure

Partition total energy, invariant energy and the parent's
lab boost angle. For an upper, bound the positive vacuum
kinematic factors at their ordered extrema, the resonance
factor at the closest invariant-energy endpoint or the
pole, and the cold probability at the largest permitted
value in the angular box. This allows a rectangular
superset of the physical integrations. For a lower keep
only boxes wholly within the allowed invariant range.
The extrema bound parent energy on both ends; the open
cold probability is unimodal, so its lower on an energy
interval is the smaller of its endpoint values. Bound the
visible FD factor from below using the least possible
individual incoming energy and the fixed lower bound on
their sum. This inequality applies before averaging over
any spin-correlated angular factors, which are nonnegative.
Integrate the exponential Boltzmann factor over each
total-energy interval analytically.

Both the source and the allowed injection cut are monotone
with the hotter endpoint in the adiabatic history. On each
temperature bin, use the opposite endpoint for H and the
entropy denominator, with their signs checked on every EOS
cubic. Below 12 MeV, retain the prior one-sided upper down
to .1 MeV; a lower need not include that positive contribution.
The enclosure is based on analytic bounds evaluated in
ordinary binary64, not formally directed interval arithmetic.

| Starting T (MeV) | Numerical FD / held comparison | FD interval lower / upper |
| --- | ---: | ---: |
| 25 | 1.324e−11 | 7.54e−12 / 2.22e−11 |
| 35 | 4.484e−5 | 2.96e−5 / 6.49e−5 |
| 50 | 0.18097 | 0.1277 / 0.2450 |

At 50 MeV the numerical incoming electron and muon shares
are about 50.37% and 49.63%; tau is negligible in this fold.
About 55.25% of the weighted charged-lepton events lie
between .4 and .6 GeV in pair invariant mass. The shares
are useful for testing which timelike bands matter, not a
substitute for a thermal hadronic spectral function.

As a distinct robustness check, if a spin-half parent in
this same vacuum two-body model is polarized, a linear
spin analyzer has directional weight at most twice the
unpolarized value by positivity. Thus changing only that
decay angular factor would multiply the upper by at most
two, still below 0.4900 of this reference yield at 50 MeV.
The positive lower in the table requires the unpolarized
model and is not extended to unknown polarization.

At 50 MeV the source is no longer unconstrained for this
model: it supplies a non-negligible but subdominant share
of the held cold comparison yield. The denominator is
the same q64 initial-population diagnostic from a different
history, not an independently measured abundance. The
low-temperature tightening does not remove channels from
the full model. Parents already present, late scattering,
stimulated decays, and hadronic production are not included
here, and the unpolarized decay is a stated model condition.
No claim about the complete final cold abundance follows
without those ingredients.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
b=results/warm-multishell-corrected-q64-fine-stage-48-2026-11-26.json.gz
python -W error warm_late_cold_visible_exact_recoil_fold_reference.py EOS.dat "$b"
python -W error warm_late_cold_visible_recoil_box_upper_reference.py EOS.dat "$b"
~~~

[Direct FD and MB folds](results/warm-late-cold-visible-exact-recoil-2026-11-26.json);
[separate interval enclosures](results/warm-late-cold-visible-recoil-box-2026-11-26.json).
