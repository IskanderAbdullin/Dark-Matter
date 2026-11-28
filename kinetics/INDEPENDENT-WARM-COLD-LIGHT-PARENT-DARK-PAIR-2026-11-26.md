# Frozen dark pair creation in the parent modes that feed very cold light

The [fixed-light inverse-decay check](INDEPENDENT-WARM-NEAR-REST-LIGHT-INVERSE-DECAY-2026-11-26.md)
isolates parent energies near 0.260 GeV, even though the light
particle is almost at rest. Quantify an independent ingredient
missing from that frozen comparison: the represented positive
sources (phiphi	o Har H) and (Lar L	o Har H) at
the required parent energies.

Tag one heavy charge and its momentum. Integrate over its
would-be annihilation partner and over the incoming CM angle
and the outgoing lab-angle panels split at the actual scalar
or light shell edges. Use the tree amplitudes for
(Har H	ophiphi,Lar L) with the two inverse incoming
occupation factors. For the light-fermion final state, the
spin average is (y^4[((M+m)^2-t)/(t-m_phi^2)]^2);
the scalar final state uses the summed (t/u) spin amplitude
and the identical-scalar symmetry factor. Omit both final
heavy Pauli factors to overcount the positive inverse gain
at fixed archived daughter distributions.

An independent check integrates the tagged gain over one
heavy charge for a smooth synthetic scalar distribution
with finite support. It gives an event rate within
(2.3	imes10^{-5}) relative of the original scalar-first
two-body fold. On the warm saved q32 states the two finest
orders agree to better than 0.3% at the three sampled
parent energies around the rest-light recoil.

At the first 10.391 MeV frame, the q32 scalar-pair and
light-pair sources divided by a lower frozen decay damping
are at most about (1.48	imes10^{-32}) and
(7.77	imes10^{-32}) **at the sampled energies**. These are
the possible additions to the frozen heavy occupation from
these two represented positive gains before other losses.
The heavy occupation used in the previous local rate was
around (3	imes10^{-19}). For q16 both tagged gains are
exactly absent on represented support: twice the largest
retained incoming energy is below the tagged energy plus
the mass of the required heavy partner. This does not mean
all heavy production on q16 is absent.

These are frozen shape-dependent statements at two frames and
three parent energies each. They are neither uniform bounds on
unknown intra-shell/exterior spectra nor a calculation of
elastic gain, heavy transients or coupled light evolution.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error warm_cold_light_parent_dark_pair_source_reference.py
~~~

[Diagnostics](results/warm-cold-light-parent-dark-pair-source-2026-11-26.json.gz).
