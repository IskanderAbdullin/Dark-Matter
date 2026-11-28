# Direct quantum rate corrections do not remove the earlier applicability stop

Take the earliest 99.98 MeV seed in the
[initialization test](INDEPENDENT-WARM-REHEATING-EARLY-SEED-SENSITIVITY-2026-11-26.md).
Fold the sampled quantum collision factors directly
along this path, with the same method as the
[separate 60 MeV initialization](INDEPENDENT-WARM-REHEATING-TABULATED-QUANTUM-2026-11-26.md).
Repeat after advancing on the new rates and refining
the table. All channels use the exact common-
temperature signed affinity, but the heavy channels
still approximate final heavy Pauli by the midpoint
of its local factor interval. The maximum width of
that interval is about .00151.

| Early rate representation, 99.98 MeV seed | Visible T where the .28 phase test stops | S at that point, per visible entropy |
|---|---:|---:|
| MB event rates, FD/BE moments | 86.0943 MeV | 7.0606e−5 |
| Tabulated quantum factors, 9 samples | 86.0732 MeV | 7.0668e−5 |
| Tabulated quantum factors, 13 samples | 86.0748 MeV | 7.0663e−5 |

The refined thermal light↔scalar detailed-balance
ratio agrees with the direct reverse/forward fold
to better than 3e−6 at all held calibration points.
The fine 13-point quadrature leaves an independent
maximum occupation-bracket defect below .00144,
recorded for mesh control rather than hidden.
The signed count pathway is therefore not
substantially changed by these measured quantum
corrections within this particular reduced thermal
calculation. The fixed phase-factor stop remains
a declared test boundary, not a phase transition.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error warm_prompt_reheating_tabulated_quantum_reference.py EOS.late.dat --early
~~~

[Both calibration grids and expanding histories](results/warm-prompt-reheating-early-tabulated-quantum-2026-11-26.json).
