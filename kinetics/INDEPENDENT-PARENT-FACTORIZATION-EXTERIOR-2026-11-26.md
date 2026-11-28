# Parent-energy factorization and exterior return for the tagged scalar

The [previous wide-boundary check](INDEPENDENT-EXPANDING-TAGGED-PROMPT-RETURN-2026-11-26.md) explicitly showed why a first exit from the grid cannot be equated with destruction or permanent escape. This follow-up factors the sequential on-shell operator by the physical parent lab energy and keeps actual transport into and back from a wider exterior. It is still a prescribed bath, a low-occupancy linear tag calculation unless otherwise noted, and not a relic prediction.

## Orthogonal parent-energy and detailed-balance check

At fixed allowed parent energy, incoming scalar absorption in cell i is proportional to the integral of f_chi over that cell's physical energy slice. The prompt decay into j is proportional to the integral of 1-f_chi over its slice, normalized over all allowed outgoing energies, including the last exterior. Fubini turns many independent angular returns into a single positive product integral. Include the laboratory phase-space Jacobian and total scalar width; do not also add the full on-shell pole in an elastic 2→2 rate. On overlapping nonthermal controls, refined total rates differ from the earlier independent angle calculation by at most .26%, while the largest refined active-row L1 change is below .09%.

In a genuine FD bath, multiplying the incoming tag by its Maxwell–Boltzmann equilibrium factor gives
`exp(-E_phi/T) f_chi(E_parent-E_phi) = exp(eta_chi-E_parent/T) [1-f_chi(E_parent-E_phi)]`.
Thus the exact equilibrium-weighted pair kernel is a symmetric positive Gram integral, with the closed mesh's expected stationarity and number conservation. Direct FD and factorized forms agree within 4e-16 relative without symmetrizing measured entries. A constant-within-cell trial distribution is not an exact equilibrium projection; its measured differential stationarity defect is 1.1–3.8% in the sampled finite shells. An optional MB within-cell trial checks projection sensitivity, not thermalization of the actual nonthermal bath.

## Stimulated prompt-parent candidate

With scalar stimulation included, at each fixed parent energy let A_i integrate f_chi f_phi and B_i integrate (1-f_chi)(1+f_phi) over cell i. The isolated prompt parent fixed point is f_1=A/(A+B), where A,B are the sums over slices. Its scalar transfer is proportional to `B_i f_1 - A_i (1-f_1)`, which conserves scalar count exactly in the complete prompt return and is nonnegative in an empty cell. Its continuous thermal BE/FD equilibrium and parent FD occupation are tested separately. An analytic Jacobian passes independent finite differences within 1.4e-9 relative in the sampled cases.

For three synthetic scalar shapes peaking at f_phi=.003, the instantaneous projected transfer differs from the linear tag limit by at most .047%; for shapes peaking at f_phi=1, it differs by roughly 9–10%. These are instantaneous controlled comparisons, not a bound on cosmological history or on an unresolved sharp mode. The stationary-parent algebra does not include other sources or loss processes for the parent, and this test is not a time evolution with full stimulated kinetics.

## Resolve motion outside the old inner mesh

At 50 MeV reheat, retain actual exterior momentum cells out to q/T_R=64, while keeping the previous q/T_R=24 surface only as a counter. Solving the evolving graph with actual transitions in both directions gives **8.56122 outward and 8.56122 inward crossings per injected tag** across that surface, but only 4.31e-6 of the tags reside beyond it at the endpoint. The count that ever reaches the distant q/T_R=64 boundary is about 2.5e-8 of the injected count. In this prescribed calculation, recording the first exit at q/T_R=24 as a loss discarded particles that usually returned. Counting repeated crossing events does not measure a probability above one.

At this support the linear terminal soft fraction under the diagnostic p/T_SM=.005 cut is 4.69e-8. This is not a physical abundance or a claim that the whole time-dependent scalar spectrum has converged. Resolving the exterior fixes an identified error with loss accounting; momentum-bin shape, source history, missing coupled kinetics and unresolved soft structure remain separate checks. Large external source bins are integrated by positive overlap probabilities; subtracting cumulative sums led to cancellation in some wide-support tail rows, and was rejected.

## Reproduce
Use the same [central SM EOS](https://member.ipmu.jp/satoshi.shirai/DM2020/eos2020.dat) and the preceding direct-lepton spectra or snapshot recipe. Do not redistribute the EOS from this repository.
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error kinetics/empty_lepton_pair_spectrum_reference.py EOS_TABLE > RAW_SPECTRA.json
python -W error kinetics/factorized_prompt_return_graph_reference.py EOS_TABLE RAW_SPECTRA.json
python -W error kinetics/stimulated_prompt_parent_elimination_reference.py EOS_TABLE RAW_SPECTRA.json
python -W error kinetics/expanding_tagged_injection_transport_reference.py EOS_TABLE SNAPSHOTS.json --create
python -W error kinetics/expanding_tagged_injection_transport_reference.py EOS_TABLE SNAPSHOTS.json --TR .05 --step 2 --bins 22 --support 64 --orders 5 4 3 12 --factorized --new-mesh --outer-diagnostic 24 --prompt-return --inject --reflect --return-only
```

### Extended-support and bath-backreaction caution

Retaining additional explicit support (q/T_R=96 and 128) keeps the ultimate far-boundary probability below 2e-13 in the sampled controls, while actual crossings of the interior diagnostic surface continue in both directions. The mean final-snapshot return probability of the entering resolved exterior flux rounds to one in both runs. This is a fixed-bath model check, not a mathematically rigorous error bound. The *number* of interior crossings changes from 8.56 to 3.29–4.47 across the selected meshes/time selections, and the soft fraction changes to 2.31–2.69e-8. Thus repeated physical reentry, but not the precise projected momentum distribution, is supported by these refinements.

The final scalar energy inside these wider grids is about .39–.40 of the energy of the independently prescribed stable-fermion spectrum. This ratio is not itself a measured energy transfer integral (the scalar source already brings energy). It does show why fermion feedback cannot be declared small merely because the dark sector is dilute relative to the visible radiation. Both fermion and scalar equations must exchange energy consistently for an actual spectrum or thermalization conclusion. Reproduce these extended checks by selecting q/T_R=96 (23 base bins, snapshot step 2) and 128 (23 base bins, step 3) in the above command.

At q/T_R=96, the MB within-cell trial gives a soft fraction 1.54e-8 and scalar-to-fermion energy ratio .367 rather than .401, with about .380 mean crossings in each direction. This is a projection sensitivity in the prescribed nonthermal bath, not a thermal-bath alternative. It gives another reason not to quote an exact crossing count or a converged soft fraction.
