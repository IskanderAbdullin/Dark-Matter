# Common-affinity-centered stiff prompt return and expansion controls

The prior [manufactured FRW split](INDEPENDENT-JOINT-KINETIC-MODES-AND-IMPLICIT-2026-11-26.md) intentionally weakened the collision rate. Repeating that proposed macroscopic split with the original difference-of-products residual and the physical amplitude failed to reduce the implicit nonlinear residual. This follow-up avoids subtracting nearly equal capture/release products. It is still **one resonant channel on explicit finite support, with injection disabled after initialization**, not a complete kinetic cosmological model.

## Centered stiff collision residual

At fixed lab parent energy, factor the common equilibrium affinity from both capture and release. Write the angular nonequilibrium exponent as a deviation from the common thermal reference with the same two species counts and total energy. Instead of (B_i A-A_i B) computed directly, evaluate differences of the small deviations using expm1. The thermal reference is an arithmetic origin for the residual, not an imposed initial state or final equilibrium constraint. Evaluate the Jacobian of the same centered formula independently by finite differences. At small steps, direct and centered fluxes agree at about 1.3e-14 and the resulting implicit states at about 3.3e-15 relative in the chosen per-cell normalization. The new flux vanishes exactly when all deviation coefficients are zero.

The strictly populated test reaches a single step of .02 Hubble times at the unscaled prompt amplitude, with relative count/energy invariants below 2.5e-14 and positive entropy generation. This is not a claim of temporal accuracy for arbitrary spectra at such a step.

## Unscaled prompt channel on two controlled FRW starts

Start from a strictly positive manufactured shape near equilibrium, or initialize the warm number shells using the direct-lepton source histogram with constant occupation in each shell (so the initial within-cell kinetic moment is projected, not measured). The second control retains comoving q/T_initial from .161 to 15.87, dropping about .034% of scalar number and .156% of stable fermion number in the source histogram. Because it **deliberately omits the coldest cells**, it cannot establish the cold tail. Both outgoing and incoming support cuts are applied at each projected event. Update the same number/energy and pressure-work ledgers as before, including opposite heat and independently integrated entropy.

At a=1.0725, the 3- and 6-split warm-start calculations transfer .08093 and .08090 of the initial dark energy to scalars in comoving (a^4ho) units. Both stable counts and the total dark energy ledger (with pressure work) agree to relative below 3.2e-14 in these controls. Thus the earlier reduction of collision amplitude is not needed for this restricted prompt-only calculation. The warm result is not the full universe: no reacting continuing source, no below-pole angular redistribution in this time advance, and no numerical treatment of omitted extreme momentum populations.

## Separate nonlinear below-pole elastic check

The previously linear below-pole Gram is also checked in a joint nonlinear FD/BE event flux, with exact binwise partner energies, empty-bath-independent detailed balance, and the real s pole excluded. The centered near-equilibrium reference is exactly stationary by construction, and the nonlinear directional derivative matches the separate projected linear Gram to 2e-12. Event entropy matches the entropy-variable projection to roundoff, and all three collision invariants close. The finite perturbation used in this test produces an 8.6–8.9% nonlinear correction to the linearized rate. This independent flux test is not yet coupled to the full nonlinear FRW update.

## Reproduce
Use the same central visible EOS and generated raw spectra described in the earlier references.
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error kinetics/joint_resonant_centered_stiff_reference.py EOS_TABLE RAW_SPECTRA.json
python -W error kinetics/joint_resonant_frw_split_reference.py EOS_TABLE RAW_SPECTRA.json --physical-amplitude
python -W error kinetics/joint_resonant_frw_split_reference.py EOS_TABLE RAW_SPECTRA.json --physical-amplitude --warm-source-init
python -W error kinetics/joint_below_pole_elastic_nonlinear_flux_reference.py EOS_TABLE RAW_SPECTRA.json
```
