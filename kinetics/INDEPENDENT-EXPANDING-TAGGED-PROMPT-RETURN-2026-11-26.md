# Tagged-scalar transport through an expanding prescribed injection bath

This is a fixed-background linear test extending the frozen-injection snapshots. The other scalar modes do not stimulate scattering or decay, the fermion bath does not respond to the test scalar, and the invisible energy does not backreact on H. All selected fermion histories are obtained by direct lepton production plus prompt heavy fermion decay, using [the central SM EOS](https://member.ipmu.jp/satoshi.shirai/DM2020/eos2020.dat) and free redshift. The matched-FD option is an alternative shape with the same instantaneous number and mean energy, not an asserted equilibrium or mathematical bound.

## Numerics and external support

The comoving momentum grid is fixed. At each snapshot its physical edges use a=(s_R/s_SM)^(1/3), and the nonthermal below-pole scattering generator retains positivity under linear time interpolation. Column sums close exactly for tag number plus external support; an outer-grid event is an explicit accounting category, not physical disappearance. The time solver is checked by no interactions, isolated sink, isolated elastic motion, tighter integration, and a separate exponential midpoint product. The latter gives the cold benchmark survival fraction to better than 2e-7.

For 20 MeV reheat through 11 MeV visible temperature (a=1.845), a rest-bin impulse keeps about 99.45–99.55% of its tagged count inside the included grid across the baseline and extended support tests. The q/T_R=1.3 prompt-return test retains only .308 in its original momentum cell. About .907 lies in bins wholly below the diagnostic final p/T_SM=.005 cut. These are different measures from overall in-mesh survival. A deliberately warmer starting cell loses more; its exact fraction is grid and boundary dependent. For 50 and 100 MeV reheat the probabilities to leave an included *narrow* grid or encounter at least one inverse decay become large. The split between these two outcomes is strongly support dependent, and does not give a physical relic count.

## Distinguish inverse absorption from prompt return

A scalar captured in an inverse decay normally has a scalar descendant when the heavy on-shell fermion decays back. The code therefore compares one-event tagging (stop at first capture) with an explicit conserving prompt-return map. The latter integrates the occupied partner energy and the outgoing scalar angle with Pauli blocking of the released fermion. Each on-shell capture immediately returns exactly one scalar on or outside the mesh; the absolute number cannot vanish in this channel. Sampled rest widths after time dilation exceed H by more than 7e9. Elastic scattering uses a below-pole cut, so no independent full on-shell elastic pole is added to this return map.

At 20 MeV reheat, with support up to q/T_R=1.3, the rest-bin tag survives inside with probability .99476 in the prompt-return control; the mean explicit return count before outer loss is only .00406 in this finite period. At 50 and 100 MeV, using q/T_R up to 3 in a comparison, essentially all initial rest tags leave that finite mesh after repeated on-shell returns and off-pole scattering, with mean return counts about 5.2 and 20.5 respectively before exit. Rays leaving the mesh can physically return; these results do not establish either too few or too many scalars at later times. A true nonlinear or coherent zero mode is not part of the sampled finite-volume space.

## Prompt scalar injection as a prescribed dilute source

Insert the direct scalar source using exact decay-angle probabilities and independent time nodes. A coarse source trapezoid failed the cold count check; the refined independent source quadrature agrees with the previous full integrated Y within .12% for 20 MeV and .021% for 50 MeV reheat. The original impulse and the actual prompt source are distinct: at 50 MeV essentially none of the source launches below the diagnostic p/T_SM=.005 cut. With elastic reprocessing, the model makes a nonzero soft tail; its eventual numerical fraction still depends on boundary treatment. For example, with q/T_R=4, the absorbing return model gives 2.7e-17 of all injected tags in resolved soft bins while the censoring outer boundary gives 1.8e-8. Both label out-of-mesh injection and later first exits; neither is a physical lower or upper bound on the final full-system abundance. The latter also leaves tags in the source cell if their sampled event would exit, so it is only a boundary stress test. At this point even the larger finite-cell average occupation among the 50 MeV reacting controls is below 3e-4, but this does not bound the pointwise occupancy of any unresolved spike or imply absence of a condensate.

## Reproduce
The external EOS is not republished. Use `python -W error` and one BLAS thread in each command.
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error kinetics/expanding_tagged_injection_transport_reference.py EOS_TABLE SNAPSHOTS.json --create
python -W error kinetics/expanding_tagged_injection_transport_reference.py EOS_TABLE SNAPSHOTS.json --TR .02 --step 2 --matrix-check
python -W error kinetics/expanding_tagged_injection_transport_reference.py EOS_TABLE SNAPSHOTS.json --TR .02 --step 2 --support 1.3 --prompt-return
python -W error kinetics/expanding_tagged_injection_transport_reference.py EOS_TABLE SNAPSHOTS.json --TR .05 --step 2 --support 3 --prompt-return
python -W error kinetics/expanding_tagged_injection_transport_reference.py EOS_TABLE SNAPSHOTS.json --TR .05 --step 2 --support 4 --prompt-return --inject --reflect
python -W error kinetics/empty_lepton_pair_spectrum_reference.py EOS_TABLE > RAW_SPECTRA.json
python -W error kinetics/prompt_resonant_scalar_redistribution_reference.py EOS_TABLE RAW_SPECTRA.json
```

## Wider-support negative convergence check

At 50 MeV, extend the comoving upper gate from q/T_R=4 to 24 (25 base bins), retaining exactly the same source and bath histories. The direct outside-launch fraction drops from .139 to 3.16e-9. Elastic-only later exit also drops below 2.1e-7 of the full source, but the absorbing prompt-return iteration still accumulates a first-exit probability .99988. Meanwhile, the censoring-control terminal soft fraction rises from 1.76e-8 to 1.59e-6, and its maximum bin-average occupation rises from 2.8e-4 to .096. A wide absorbing boundary does not supply a converged result for the hot examples: many sequential returns can make the probability of ever crossing an artificial surface large even if most particles would later reenter. The censoring control also changes, so it does not replace real return from an extended exterior. These negative results are a reason to build and audit an explicit exterior-return model before interpreting the distribution, not to extrapolate either outcome.

Reproduce the wider check with `--TR .05 --step 2 --bins 25 --support 24 --prompt-return --inject --reflect` on the same prepared snapshots.
