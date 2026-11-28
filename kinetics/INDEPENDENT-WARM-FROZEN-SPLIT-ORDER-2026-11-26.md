# Compare three collision schedules on one frozen warm step

The [expanding q16 timing test](INDEPENDENT-WARM-EARLY-TIME-REFINEMENT-2026-11-26.md)
cannot attribute its whole error to the implicit collision solver.
Here hold the momentum geometry, elastic quadrature panels,
visible temperature and Hubble factor fixed at the same
mid-geometry from the 12 MeV q16 shape. Use the same signed
pair projection and explicit midpoint pair kicks. Write P for a
pair map and C for a prompt-plus-elastic collision map. Compare:

- P half, C full, P half.
- C half, P full, C half.
- P half, two C half maps, P half.

The last schedule is a control with the same number and sizes of
C maps as the reversed schedule. These are symmetric *operator
schedules*, not entirely second-order methods when a C map is
backward Euler. The schedules do not all use the same number
of pair evaluations; this is not a cost comparison.

For normalized interval 0.006 with backward Euler, the reverse
versus standard schedule changes the inclusive scalar count by
+6.89×10⁻⁵ at one subinterval and +1.83×10⁻⁵ at four.
The stable changes are −2.77×10⁻⁵ and −7.38×10⁻⁶.
Splitting the central C map into two halves instead does not
remove the discrepancy: the same reverse-versus-control scalar
changes are +6.85×10⁻⁵ and +1.81×10⁻⁵. In the 2→4 change,
the standard scalar count shifts by +2.08×10⁻⁵ and the reverse
one by +3.11×10⁻⁶. A smaller adjacent change in one frozen
test does not establish smaller full-trajectory error.

Changing the C maps to the tested L-stable two-stage solver
does not remove the order difference: at four subintervals
the reversed scalar shift is +1.88×10⁻⁵ and stable shift
−7.56×10⁻⁶. Thus the isolated benefit of that collision solver
does not by itself settle the splitting error.

The signed pair number and total energy audit includes new
exterior particles. They do not feed back into this frozen
calculation. There is no FRW motion within the interval,
production from an external bath, or feedback onto its
geometry and radiation clock. The default warm integrator
and the actual resumed q32 run are unchanged.

~~~sh
export OPENBLAS_NUM_THREADS=1
S=results/50MeV-from-empty-pair-quantum-shape-2026-11-26.json
python -W error warm_frozen_split_order_reference.py EOS.dat "$S" --method backward \
  --schemes pair-collision-pair collision-pair-collision >be.json
python -W error warm_frozen_split_order_reference.py EOS.dat "$S" --method backward \
  --schemes pair-two-half-collisions-pair >half.json
python -W error warm_frozen_split_order_reference.py EOS.dat "$S" --method sdirk2 \
  --schemes pair-collision-pair collision-pair-collision >sdirk.json
python -W error warm_frozen_split_order_audit_reference.py EOS.dat "$S" be.json half.json sdirk.json
~~~
