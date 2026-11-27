# Short coupled 4-to-2 sensitivity, not a chemical freeze-out trajectory

Couple the [reversible cold number-changing stencil](INDEPENDENT-EARLY-POLE-SCALAR-NUMBER-CHANGE-STENCIL-2026-11-26.md)
to the live parent, daughter, evolving Hartree pole and Friedmann
background, starting from a stored ratio-0.06 state at x=0.05.
The threshold squared-amplitude multiplier is an explicit nuisance
parameter, not a known in-medium suppression. Threshold geometry
and hot-node deposition are recalculated at the current pole.
The integrator stops if the cold-to-pole hierarchy gets too small.

The distinct number channels are checked: fermion daughters minus
removed parents still balance, but scalar count no longer follows
fermion daughter count. Its difference is checked against a
separately integrated 4-to-2 source. Total energy is still checked
against the integrated stress trace; the chemical transfer must
not remove energy from the tracked sector.

The centered-time short integrations pass energy and number
ledgers. At multiplier 10⁻⁴, between x=0.05 and x=0.0500001 on
the third grid, integrated chemical loss is 6.44 × 10⁻⁸ of
initial scalar number, while the decay channel adds 4.38 × 10⁻⁶.
These are short, frozen-input-start diagnostics, not a bound on
subsequent nonlinear growth of the hard occupation. The multiplier
10⁻⁹ test reaches x=0.0501; it is near the no-chemical limit here.

A single exploratory third-grid noncentered run at multiplier
10⁻⁷ reaches x=0.052 and gives 0.000235717 versus 0.000414781
in the same-grid pure-decay comparator. Integrated chemical
loss is 50.94% of the initial scalar number, partly offset by
15.68% net decay production. The chemical ledger mismatch is
3.2 × 10⁻⁸ of the starting count. This trajectory suggests that
accepting a small initial sink does not bound the later
reaction in the illustrative model; the hot-state occupation
feeds back. The longer result is not confirmed across momentum
grids or solver starts.

The first attempts with stronger coupling over a longer interval
failed with a step-size error in the off-center time coordinate.
Re-centering makes the very short tests reproducible, but no
longer nonzero-chemistry convergence claim is made here. A
threshold angular closure and rapid motion across hot grid cells
can affect the evolution; a more resolved number-changing kernel
and physical medium amplitudes are still required. Nothing in
this diagnostic supplies a final spectrum for CLASS.

~~~sh
python -W error early_pole_scalar_chemical_sensitivity_reference.py results/early-pole-spectrum-third-m06-2026-11-26.json.gz --strength .0001 --xend .0500001
python -W error early_pole_scalar_chemical_sensitivity_reference.py results/early-pole-spectrum-third-m06-2026-11-26.json.gz --strength 1e-9 --xend .0501
~~~

[Raw successful ledgers](results/early-pole-scalar-chemical-short-ledgers-2026-11-26.json).
