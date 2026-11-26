# Explicit prompt re-entry prevents gross first crossings from accumulating

The [instantaneous open boundary check](INDEPENDENT-WARM-POST-CHECKPOINT-OPEN-PROMPT-FLUX-2026-11-26.md)
showed a very large initial one-way q16 birth rate beyond
its retained edge if the unrepresented states are empty.
Test what happens when new exterior populations can react
back rather than collect in a one-way bookkeeping sink.

Use the actual q16 10.391 MeV shape inside the original
boundary. Add the matched wider q32 grid, with both species
initially empty in the added shells. Also add a q64 support
control and a time-halved q32 sequence. Evolve only prompt
on-shell transfer in this finite represented system with
the radiation clock and physical momenta held fixed. Both
conserved species counts and total dark energy are computed
directly from the jointly evolving moments. The narrow
control, by design, has zero outside its original grid.

At (HDelta t=5	imes10^{-6}), the explicit wider calculation
has roughly (2.05	imes10^{-5}) of admitted scalar number and
(3.70	imes10^{-5}) of light number outside the old q16 boundary.
In contrast, extending the **initial frozen one-way** rates
linearly would give .0317 and .0273. This comparison shows
why initial crossing rates are not accumulated populations.
The positive coupled exterior fills rapidly and its occupation
then enters both inverse and decay statistics. It does not
demonstrate a unique microscopic re-entry history or count
individual crossing events.

This is a very short, fixed-background, no-external-source
transient, not the expanding warm trajectory. It does not
erase inherited unknown tails, prove convergence of all spectrum
moments, or establish the late abundance. All raw support/time
comparisons, outermost-shell counts and conservation defects
are in the diagnostic file.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error warm_prompt_reacting_exterior_transient_reference.py
~~~

[Diagnostics](results/warm-prompt-reacting-exterior-transient-2026-11-26.json).
