# Stationary neutral-operator sum rule and chemical balance

This finite Hamiltonian witness is an independent check of the
[stationary spectral gate](INDEPENDENT-EARLY-POLE-WILSON-STATIONARY-LINE-GATE-2026-11-26.md).
Two species share one oscillator and a two-level orbital space. The total
fermion number is one, with the parent and daughter labels distinct and
conserved by the toy bath. Diagonalizing the two charge blocks gives
paired nonnegative Lehmann weights for the operator that converts the
parent to the daughter. There is no collider, expansion or continuum
linewidth in the model.

At each line with emitted energy omega, inverse/forward equals
exp[(omega−Δμ)/T] at common equilibrium temperature. The total transition
weights independently sum to the parent and daughter probabilities,
respectively, for the unitary inter-block map. Both checks close below
2e−13. When the two dynamical blocks differ only by a constant energy
and the map is the identity, all emission weight stays on the single
line at .73, even though both species are coupled to the oscillator.
A species-dependent orbital half-shift .035 puts .25563 of emission
weight outside a .02 window around .73; on the identical dynamical
blocks, rotating only the probe by .3 radian puts .08434 outside the
window. These are stationary discrete sidebands, not decay widths.

A positive smoothing kernel applied separately to both spectra is
not enough to preserve the thermal sign at the displayed frequency.
For half-shift .035, at omega=.665 below the chemical gap .68, common
Gaussian smoothing of width .035 produces stimulated difference
−.06361 instead of the required positive sign. It produces −.11594
at the gate where the difference should be zero. Smoothing the
symmetric spectral core first and then applying the thermal factors at
the target frequency gives +.03080 and exact zero, respectively. This
is a thermally constrained plotting or remapping convention, not a
physical continuation of the discrete oscillator.

The two largest oscillator truncations agree within 1e−7 in emission
spread and off-window weight. A continuum calculation would have
to supply the neutral operator and its matched interacting states or
a conserving correlation kernel rather than infer it from finite
contour coherence.

Reproduce from kinetics/ with early_pole_stationary_neutral_operator_reference.py.
Raw data: results/early-pole-stationary-neutral-operator-2026-11-26.json.gz.
