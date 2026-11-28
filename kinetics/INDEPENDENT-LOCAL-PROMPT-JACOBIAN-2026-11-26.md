# Parent-local assembly of the physical prompt Jacobian

The prompt integral is block-local in each parent lab-energy panel.
The four nonzero features of an angular row are the scalar count
and kinetic-coordinate weight and the corresponding two fermion
features. For the exact existing centered event integral, assemble
its analytic derivative within only the union of active columns
of each parent panel. Apply the moment projection there, then
accumulate in the global Jacobian. The physical value path,
kinematic support, and rate normalizations are unchanged.

This check uses the archived wide physical calculation immediately
after two accepted collision intervals, reconstructing the third
midpoint by redshift. It has 94,420 angular nodes in 1670 parent
groups, with 30.27 active columns on average compared with 92 total.
The assembled derivative agrees with independent application of all
92 old analytic directional derivatives to 5.6e-15 relative at the
actual state and at a perturbed state near the common reference.
The JVP finite-difference check is within 3.3e-11 and invariant-null
defects within 4e-16. The value is identical to the earlier evaluation.

On one single-thread measurement the old full matrix reconstruction
took 1.18 s, versus 0.27 s for parent-local assembly at the same
actual state; one-time panel compilation was also 0.27 s. These are
hardware diagnostics, not a general solver runtime guarantee.
The benchmark measures Jacobian construction, not how many
Newton iterations a strongly stiff rare exterior shell needs.

Prepare the [benchmark EOS files](INDEPENDENT-BENCHMARK-SETUP.md).
From `kinetics`:

```sh
OPENBLAS_NUM_THREADS=1 python -W error compact_prompt_jacobian_audit_reference.py EOS.dat EOS.late.dat results/wide32-embedded-vacuum-2026-11-26.json results/wide32-pre-third-collision-checkpoint-2026-11-26.json
```
