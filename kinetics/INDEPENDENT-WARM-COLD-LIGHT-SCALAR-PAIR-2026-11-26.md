# Scalar pair conversion at zero tagged light momentum

The [fixed-light decay calculation](INDEPENDENT-WARM-NEAR-REST-LIGHT-INVERSE-DECAY-2026-11-26.md)
shows a local rate sensitive to the high-energy scalar support.
Independently test the two-body channel (Lar Lleftrightarrow
phiphi) at exactly zero momentum of a tagged light spin and
charge, on the same archived distributions.

Integrate over the other fermion charge (two spin states) and
the scalar final polar angle in the center-of-mass frame. For the
tagged particle at rest (s=2m(m+E_{ar L})). Split at the
actual represented outgoing scalar energy edges. For the inverse
gain allow partners above the represented fermion edge whenever
the two represented incoming scalars can create such a partner;
that final state has zero archival occupation in this diagnostic.
Keep the exact two scalar Bose factors and both fermion Pauli
factors. The ratio of gain and loss for a common FD/BE equilibrium
with the same chemical potential is one to about (4	imes10^{-14})
at the event level.

At the first 10.391 MeV checkpoint the q16 and q32 net light
rates are (-1.293	imes10^{-5}) and (-1.205	imes10^{-5})
in occupation per Hubble time. Both stay negative across the
six matched frames; refining the last two quadratures moves
net values by less than (3	imes10^{-5}) relative. Cropping
**only** the q32 scalar support to the matching q16 edge changes
the q32 estimate by about (6.1	imes10^{-5}) relative, not a
sign flip. Thus this local pair contribution contrasts with the
support-sensitive one-to-two reaction, despite using the same
archived scalar populations.

This test is fixed at the precisely rest tag and does not determine
the behavior over the full light distribution, unresolved exterior
population, or late history. No final yield is inferred.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error warm_cold_light_scalar_pair_reference.py
~~~

[Diagnostics](results/warm-cold-light-scalar-pair-2026-11-26.json).
