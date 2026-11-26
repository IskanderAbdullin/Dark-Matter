# The re-populated prompt tail resembles the actual wider checkpoint

The [explicit short re-entry test](INDEPENDENT-WARM-PROMPT-REACTING-EXTERIOR-TRANSIENT-2026-11-26.md)
started on the q16 checkpoint and populated explicit initially
empty exterior shells. Compare its resulting shell number and
kinetic-energy profiles with the saved **actual q32** checkpoint
at the same visible temperature, 10.391 MeV. Match physical
momentum edges exactly on the common core. Use a common q64
grid so both the q16-origin transient and a separately advanced
q32-origin control have room for further transfer.

The archived q32 tail above the old q16 boundary holds
(2.08167	imes10^{-5}) and (3.74778	imes10^{-5}) of
the represented scalar and light numbers, respectively, and
(1.16662	imes10^{-4}) and (1.52809	imes10^{-4}) of
their energies including rest mass. At (HDelta t=5	imes10^{-6}),
the respective q16-origin fractions are lower by 1.56% and
1.21% in number and by 1.59% and 1.15% in full energy.
After independent normalization of each tail's shell weights,
the (L^1) distances in tail number are .000995 and .000189,
and in tail kinetic energy .001131 and .000212. This comparison
tests more than a match of the total count.

In the separate transient starting with the actual q32 spectrum,
the already populated q16-to-q32 band remains close to its
initial shape. The new q64 populations beyond the **q32**
edge are (1.15	imes10^{-11}) and (2.22	imes10^{-11})
of the initial species numbers at this very short endpoint.
Both conserved species numbers and total dark energy are
audited independently.

This is a matched-support check of the isolated prompt
operator on fixed physical momenta and a fixed radiation
clock. The actual archive reflects a different, longer
calculation. Similar tail shapes here do not establish a
common history, full kinetic equilibrium, actual inherited
populations beyond q32, nor late abundance.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error warm_prompt_reentry_shape_audit_reference.py
~~~

[Diagnostics](results/warm-prompt-reentry-shape-audit-2026-11-26.json.gz).
