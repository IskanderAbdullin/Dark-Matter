# Opening the prompt boundary gives sharply different instantaneous spill rates

The [prompt geometry test](INDEPENDENT-WARM-COLD-PROMPT-SUPPORT-GRID-2026-11-26.md)
identified cold light shells disconnected on q16. Measure
the independent **full-angular**, stationary-parent on-shell
operator at all six matched warm endpoints near 10.4 MeV.
The incident scalar and light profiles are exactly the saved
retained shapes, with zero occupation outside for this control.
Unlike the closed solver, the outgoing state is allowed beyond
either edge. This extends the earlier
[fixed-snapshot outward kernel](late_prompt_exterior_snapshot_reference.py)
to the current matched checkpoint and adds a same-shape crop.

The heavy stationary abundance at each momentum is fixed by
the exact inverse/decay angular averages over the chosen
incident shapes and full physically allowed final angular range.
At the outer boundary the inverse occupation product vanishes
by construction. Count the positive births outside separately
for each species and include rest masses in the energy ledger.
These are local-in-time rates with a prescribed state. A positive
first crossing is not permanent escape, and neither a one-way
reservoir nor continuing frozen rates describe its later return.

At the first endpoint, normalized to each species' current
retained count per Hubble time, the q16 scalar/light rates are
about (6.34	imes10^3) and (5.47	imes10^3), while q32
gives (1.22	imes10^{-3}) and (2.15	imes10^{-3}).
The giant q16 dimensionless rates cannot be sustained at the
same shapes; less than (3	imes10^{-4}) of its fast gross
decays leave the domain, but that gross turnover is very large.
The q32 sum of outer **full energy** rates is about .01587 of the
current retained dark energy per Hubble time at the first frame.
Full energy rates and number rates are separate; none are yields.

Crop both incident species on the exact existing q32 shell edge
matching q16, keeping the full nonthermal core unchanged. The
q32 cropped scalar/light rates rise to (6.38	imes10^3) and
(5.51	imes10^3), within 1% of actual q16. This identifies the
finite-support boundary as the driver of the mismatch, not
different fitted shapes inside the common core. The disparity
persists at all six endpoints and is stable at the two
event-quadrature orders.

The next physical check is to evolve newly populated exterior
modes with both directions of transfer and verify closure of
both conserved species numbers and of energy. This frozen
comparison cannot determine the actual residual exterior
distribution, the coupled dynamics, or reuse of those modes.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error warm_post_checkpoint_open_prompt_flux_reference.py
~~~

[Diagnostics](results/warm-post-checkpoint-open-prompt-flux-2026-11-26.json.gz).
