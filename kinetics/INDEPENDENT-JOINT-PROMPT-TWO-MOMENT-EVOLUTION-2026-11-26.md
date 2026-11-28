# Short joint nonthermal resonant evolution with separate number and energy

This is a closed-support, fixed-background short collision test at the sampled 50 MeV reheat to 45 MeV visible injection snapshot. The scalar and stable fermion spectra both evolve. It is not a reheating history, expanding prediction, coherent zero-mode calculation, or late relic abundance. Across these short intervals expansion and continued production are explicitly disabled. The only channel here is the previously tested, stimulated prompt resonant redistribution, not the whole dark sector. The largest integration time is 5e-8 Hubble times for this particular snapshot; H is a normalization constant.

## Reconstruct and update the two independent moments

Advance each finite shell's number and kinetic-energy moment, with a normal scalar BE and a Pauli-limited FD finite-shell shape reconstructed after each stage. Slopes are parameters of a local basis, not temperatures of individual shells. The separate support-bound trial uses cubic phase-space measure to check the FD packing energy even at tiny fill; subtracting nearly identical occupied interval endpoints initially generated a false rejection, so normalized packing means are integrated directly. At 25 and 33 shells the recovered moments of a single small trial step agree to 1.43e-13 and true event-energy conservation to 9.43e-15 relative to the scalar exchange.

The longer solver explicitly evaluates both reconstructed occupations at the same on-shell lab energy pairs. At each parent energy, its absorbing and releasing distributions share a single conservative prompt fixed point. Both projected shell rates follow from the same signed event flux; total scalar+fermion energy, and each stable count, remain invariant. The shared physical support excludes events where either constituent lies outside the selected box in both directions, rather than deleting late escaped particles. Admissibility is checked again after each stage; the solver backtracks instead of clipping moments or occupations. An aggressive-step stress test rejects 11 attempts.

The microscopic entropy production from event ratios independently matches the projected entropy derivative to 2.91e-14 in the manufactured non-equilibrium test. Continuum BE/FD detailed balance represented by the local shapes leaves only roundoff. The accepted evolution preserves both count ledgers and total energy to about 5e-16 relative in the sampled runs, and raises reconstructed entropy.

## Short times and numerical comparisons

At 1.5e-9 Hubble times the evolving scalar receives .01254 of the original total dark energy in the 25-shell refined Heun calculation. Halving the step changes this by about .075%, while the 33-shell projection changes it by about 3.5%. Higher parent/angle degrees do not change the shown digits. Reducing the shared momentum cutoff from p/T=90 to 17.58 changes this short transfer by about .012%, while the initially discarded number fractions are recorded separately.

At 1.5e-8 Hubble times the 25- and 33-shell controls give .04173 and .04302 of the initial dark energy transferred to scalars. A separate common-temperature, maximum-entropy solution on the same fixed box, with both counts and total energy fixed, would eventually have .0810 in this normalization. The tested long runs have achieved about .67–.70 of the available entropy increase, with nonzero differences in the scalar and fermion shell populations relative to that solution. This tests finite-box relaxation, not a momentum-independent kinetic clock or the full nonlinear cosmological limit.

## Reproduce

Generate direct-lepton raw input as before from the [central visible EOS](https://member.ipmu.jp/satoshi.shirai/DM2020/eos2020.dat):
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error kinetics/empty_lepton_pair_spectrum_reference.py EOS_TABLE > RAW_SPECTRA.json
python -W error kinetics/joint_prompt_two_moment_step_gate_reference.py EOS_TABLE RAW_SPECTRA.json
python -W error kinetics/joint_prompt_two_moment_evolution_reference.py EOS_TABLE RAW_SPECTRA.json --verify
python -W error kinetics/joint_prompt_two_moment_evolution_reference.py EOS_TABLE RAW_SPECTRA.json --height 1.5e-8 --dt 2.5e-10 --step 3 --verify
python -W error kinetics/joint_prompt_two_moment_evolution_reference.py EOS_TABLE RAW_SPECTRA.json --height 5e-8 --dt 4e-9 --verify
```
