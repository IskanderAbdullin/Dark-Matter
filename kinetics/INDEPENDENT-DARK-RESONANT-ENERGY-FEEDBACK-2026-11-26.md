# Actual event heat and a conservative dark-bath feedback benchmark

The [prior exterior study](INDEPENDENT-PARENT-FACTORIZATION-EXTERIOR-2026-11-26.md) identified a non-negligible scalar-to-fermion energy ratio but did not measure the heat. Here the resonant capture/prompt-reemission event is resolved in laboratory energy. Each scalar energy change is paired with the opposite heavy-partner energy change at the same fixed on-shell parent energy. We do not count the full on-shell pole twice. These are independent snapshot tests and a restricted two-temperature moment benchmark, not a converged evolving nonthermal momentum spectrum.

## Event-level audit

A parent-energy factorization integrates the incoming and outgoing scalar laboratory energies separately and preserves both heat weights, not just the destination number. Independent thermal tilted-integral evaluation agrees within 7e-14 in checked cases, and the reciprocal-energy transpose relation holds to roundoff in the common-temperature thermal test.

For the collisionless nonthermal injection snapshots at TR=20, 50 and 100 MeV with visible T=.9 TR, the instantaneous scalar heating from this channel is positive and respectively 2.58e4, 1.63e7 and 1.30e8 times the prescribed fermion energy per Hubble time. It acts **against** the decrease in fermion energy. This is an instantaneous rate under an unrelaxed distribution; multiplying it by a Hubble time is not a physical integrated heat. Its magnitude shows that, for these source-only distributions, a fixed collisionless bath is not a self-consistent forward evolution. It does not imply that an initially cold tagged scalar follows the same clock as the warm injected bulk.

The inferred instantaneous intermediate on-shell parent energy fraction relative to stable fermions is about 6.1e-8, 1.1e-4 and 3.9e-3 respectively. The bulk feedback clock remains below .001 times the estimated lab parent decay clock in these controls. These are useful minimum checks for prompt elimination but cannot replace comparison to all nonthermal eigenmodes or other channels.

## Energy-conserving two-temperature time test

Start from the measured number and mean energy of each prescribed injection spectrum and close them as scalar MB and fermion FD at separate temperatures. Then **stop injection** and conserve each stable count, with visible expansion and only the resolved on-shell resonant heat coupling. The transfer is positive to the cooler bath, using the exact sign of the conditional tilted-energy covariance, and the opposite heat appears in the fermion equation. Two implicit solvers and a separately constructed adiabatic common-temperature slow manifold check the strong-coupling evolution. The initial fast equilibration uses total energy, not an imposed average temperature. A direct solver initialization trial was too large for the initial heat stiffness and was corrected by limiting the first step using the initial derivative.

At TR=20, 50 and 100 MeV in this controlled closure the cumulative comoving scalar heat is .00403, .07777 and .08951 of the initial dark energy. The end gaps are below 5e-8 relative, and an independent zero-heat control retains .077, .338 and .360 temperature gaps. Energy and entropy ledgers close below 5e-11 in the displayed normalization, and the common manifold agrees with the resolved endpoint to about 4e-9. The heat clock in the controlled model substantially exceeds H. These tests validate bookkeeping and the stiff limit **conditional on the thermal closures**, not that actual distributions remain FD/MB during their initial heat transfer.

## Joint nonthermal collision projection

An instantaneous joint scalar/fermion operator partitions each exact lab energy interval at both species' bin edges. The scalar and fermion sources retain exact per-event energy, include the stimulated prompt-parent fixed point, and have separate number and opposite energy ledgers. The empty scalar/fermion and full fermion gates point into the admissible domain. The true weighted event heat test is consistent in sign and magnitude with the separate prescribed-fermion projection; varying the meshes changes the displayed transfer by 0.2–4.2%.

Assigning only a representative energy to each net **number-bin** update instead gives total energy defects of order 14–158% of the actual scalar heat in the tested meshes, depending on the case. This is not a violation of the continuum reaction; it is a warning against evolving only cell counts with that reconstruction. A positive energy-aware time update and converged reacting fermion distributions are still needed.

## Reproduce

Use the authors' [central visible EOS](https://member.ipmu.jp/satoshi.shirai/DM2020/eos2020.dat) without redistribution, and the earlier direct-lepton generated spectra:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error kinetics/empty_lepton_pair_spectrum_reference.py EOS_TABLE > RAW_SPECTRA.json
python -W error kinetics/resonant_on_shell_energy_exchange_reference.py EOS_TABLE RAW_SPECTRA.json
python -W error kinetics/conservative_dark_resonant_feedback_reference.py EOS_TABLE RAW_SPECTRA.json
python -W error kinetics/dark_resonant_thermal_slow_manifold_reference.py EOS_TABLE RAW_SPECTRA.json
python -W error kinetics/joint_nonlinear_prompt_return_moments_reference.py EOS_TABLE RAW_SPECTRA.json
```
