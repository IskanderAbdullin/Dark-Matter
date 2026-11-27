# A manufactured positive scalar step closes the nonlinear Hartree energy

Use the [exact endpoint Hartree gradient](INDEPENDENT-EARLY-POLE-SCALAR-HARTREE-DISCRETE-GRADIENT-2026-11-27.md) to test a nonlinear projected update, not just a rate ledger. From a closed finite stencil, keep the same incoming event cells and outgoing adjacent-cell brackets. Reweight the two convex outgoing interpolations by the available energy slack until each event direction conserves both the constant count vector and the endpoint Hartree gradient. Remove an event if its energy target falls outside both brackets combined. The nonlinear endpoint coordinates enforce a positive fixed-number simplex.

The manufactured mobility is a nonnegative coefficient frozen from the initial screened-free Born stencil on this fixture, using the rate/affinity ratio where defined. The event flux is proportional to the projected entropy secant. Consequently the entropy increment is independently predicted by the sum of mobility times squared projected affinity. This is a numerical discrete-gradient control with a frozen mobility, not a live physical matrix element.

Four positive manufactured runs across narrow and broad 16-cell scalar supports complete the implicit step without a rejected bracket. Scalar count is unchanged at the reported precision; absolute finite Hartree energy discrepancy stays below 3.4e−21. Entropy increases and the independent squared-affinity ledger coincides within 9e−13 relative. The largest scaled nonlinear residual is 1.1e−14. The adjustment in the individual leg's represented momentum magnitude can reach 4.57e−5 of maximum momenta on the broader comparison support (about 1.2e−11 on the narrow support). This demonstrates why exact count and energy alone do not validate an angular or momentum remap. No event-by-event vector momentum, physical live rate, or cosmological trajectory is inferred. The preceding [live held scalar proxy](INDEPENDENT-EARLY-POLE-JOINED-SCALAR-SELF-DECAY-PROXY-2026-11-27.md) does not use this manufactured update.

From `kinetics/`:

```sh
OPENBLAS_NUM_THREADS=1 python early_pole_scalar_secant_projected_fixture_reference.py
```

Raw: `results/early-pole-scalar-hartree-secant-projected-2026-11-27.json`.
