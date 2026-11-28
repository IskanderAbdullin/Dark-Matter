# The closed prompt grid removes a real recoil path on q16

The independent [tagged light decay calculation](INDEPENDENT-WARM-NEAR-REST-LIGHT-INVERSE-DECAY-2026-11-26.md)
found net q16 cold gain and net q32 cold loss for frozen
parent distributions. Does the archived prompt operator
actually apply those local rates? Compare the **cell-number
moment** of that operator on the saved shapes to the independent
tagged rate integrated with the same (p^2dp) measure in
the first light shell. Do not compare a point evaluation
at (p=0) to a projected cell average.

The existing closed prompt geometry admits an event only
if both daughters are represented. On q16 the upper scalar
energy edge is below the minimum allowed recoil energy at
small light momentum. Solving that minimum-energy crossing
independently gives (p_Lsimeq9.584) MeV at the first
10.391 MeV frame and (10.137) MeV at the last 10.137 MeV
frame. Consequently the native prompt calculation has no events
at all in the first 12 light shells in the first five frames
and the first 13 in the last frame. Those fully disconnected
shells account for 1.31–1.33% of represented light number
in the first five frames and 4.50% in the last. They are
disconnected for this operator, not for every reaction.
The independent kernel permits the scalar to leave the
represented region when a frozen heavy parent decays, so gives
a first-shell gain (4.4193	imes10^{-4}) per Hubble time
at the first q16 frame, instead of the native structural zero.
It uses a stationary parent surrogate, not a measured source.

There is no such recoil exclusion at the first light shell on
q32 in these six frames. Its first native cell-average net
rate is (-2.02181	imes10^{-5}) per Hubble time and matches
the independently integrated open-decay kernel. Across the six
q32 snapshots the two first-cell values agree to within
(1.5	imes10^{-9}) relative at the tested quadratures. This
cross-check does not vindicate closed support at other momenta,
later times, or in the scalar exterior.

The comparison explains why the frozen q16 prompt gain
cannot be inserted into the existing trajectory as if it
had already been evolved. A physical finite-grid extension
would have to route the escaping scalar count and energy,
partner accounting and possible re-entry, not just add its
light counterpart.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error warm_cold_prompt_support_grid_reference.py
~~~

[Diagnostics](results/warm-cold-prompt-support-grid-2026-11-26.json.gz).
