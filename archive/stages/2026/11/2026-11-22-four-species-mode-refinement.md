# Stage 2026-11-22 — second-order Fermi-logit energy-grid refinement in coupled χ1χ2+real φ+one lepton

- **ID:** DM-KINETICS-20261122-FOUR-SPECIES-LOGIT-MODE-REFINEMENT
- **Active run:** started UTC 2026-11-22T16:22:01Z; current step follows the [cutoff–mesh independence stage](2026-11-22-four-species-cutoff-mesh-separation.md). Exact stage wall time not independently recorded; no asynchronous execution claim.
- **Parent Git SHA:** [d85b7cfff3fd95ca31276c1be17aea1a733709de](https://github.com/IskanderAbdullin/Dark-Matter/commit/d85b7cfff3fd95ca31276c1be17aea1a733709de), draft PR #1 on public research branch. Its real CI run [#45044674246](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/45044674246) was still in progress at last check; neither old nor current CI provides a physical dark-matter certification.
- **Project route:** [GOAL](../../../../GOAL.md), [AGENTS](../../../../AGENTS.md), [PLAN](../../../../PLAN.md), [kinetics spec](../../../../kinetics/CONTINUOUS-RESEARCH-SPEC.md), [scientific guardrails](../../../../kinetics/RESEARCH_GUARDRAILS.md), [NEXT_ACTION](../../../../kinetics/research-control/NEXT_ACTION.json).

## Research question, hypothesis, independent control

The prior 256-energy-node χ1χ2 source still missed the direct analytic-f full Pauli 5D weak entropy by 5.73%. A generic positive H theorem is insufficient to make a particular Q or S quantitatively correct. Prove the interpolation remainder of a C² Fermi log-odds L(E) and test whether at **identical incoming 5D phase-space, cutoff and angular quadrature** the collision entropy error shrinks approximately as h². Independently measure signed χ1,χ2 heat errors and preserve sign changes under refinement.

For Fermi logit L(E), |L−IhL|≤h² sup|L''|/8 and |expit(L)−expit(IhL)|≤h² sup|L''|/32. Verify the bound directly on three oscillatory physically allowed χ1,χ2,lepton logit curves at 64,256,1024 nodes and ten thousand evaluation points each. This is pure numerical interpolation theory and not an error bound for the full nonlinear collision or cosmological rate.

## Inputs and exact reproducible outputs

Synthetic Mχ1=.1,Mχ2=.02,mℓ=0, G=.01,mV=.1, real φ mφ=.001 in the coupled model, one common lepton. Fixed mesh endpoint40.1, incoming kinetic cutoff12, segmented radial12×Gauss3, relative mu6, CM log-transfer12, azimuth12. Oscillatory logit deviations +.3sin4E, −.2sin3E, +.25cos4E; neither daughter spectrum comes from an early-Universe integration.

[Complete scientific report](../../../../kinetics/FOUR-SPECIES-MODE-REFINEMENT-AUDIT.md), [new runner](../../../../kinetics/four_species_mode_refinement_gate.py), [three tests](../../../../tests/test_four_species_mode_refinement_gate.py), [machine data](../../../../kinetics/results/four-species-mode-refinement-2026-11-22.json). Commands:
    OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_four_species_mode_refinement_gate.py
    OPENBLAS_NUM_THREADS=1 python -m kinetics.four_species_mode_refinement_gate --full
Three new focused local tests **passed in 13.84 s**. No RNG, no external datasets; future real GH CI separately verified.

Reference direct unprojected analytic-occupation same-5D Qχ1=−1.2672292278252052e−11, Qχ2=+8.75060244836362e−12, collision entropy +1.1236203110330876e−10. With energy-mode intervals 256/512/1024/2048/4096, relative errors:
- Qχ1: −0.5030%, −0.0192%, −0.0178%, **+0.0021%**, −0.00054%.
- Qχ2: −0.2297%, −0.0745%, **+0.0049%**, −0.0056%, −0.00151%.
- H: **−5.7321%, −1.4830%, −0.4003%, −0.0966%, −0.02447%**.

Observed entropy absolute-error ratios for h/2: **3.8652,3.7047,4.1459,3.9462**, consistent with O(h²) for this fixed test. Rate Qχ1 and Qχ2 errors are NOT monotone (including actual sign flips), so do not call individual heat a stable Richardson extrapolation. Eventwise H, total energy, number conservation remain exact to ~1e−23 without hiding interpolation bias.

Separate fixed-grid bin1024 incoming-cut8/12/16/20 check: χ1 Q −1.18677/−1.26700/−1.27211/−1.27227 ×10^-11. Cut16→20 has a small relative change (~0.0127%) on these fixed angles; cut8→20 differs ~7.2%. This is only a tail test at specific low incoming and outgoing quadrature, not physical 5D rate convergence.

## Negative results / declared next gate

**Positive:** mathematical and independent numeric second-order convergence of Fermi-logit mode reconstruction in discrete Pauli/H statistics with an analytic-Fermi reference. **Negative:** old 256 modes are not adequate for <1% entropy error in oscillatory shape, even when all charge/H laws pass; 4096-node 0.024% entropy error does not certify full angular/radial/UV kinetic rate and time evolution. The chosen historic spin factor gχgℓ=4×4, electron/neutrino counting, untested plasma dynamics and gauge sector remain synthetic. No 86-eV φ freeze-out f(q), Γelastic(T), absolute ΩDM, BEC, galaxy/cosmological CLASS readiness.

Next: real strong-to-weak **joint** 5D convergence with fixed sufficiently fine energy-mode grid while varying physical p/k cutoff, nlocal, nμ, CM and azimuth at 3+ values; include mℓ=.5 and mχ1/.02. Independently resolve moving Yukawa φ low-energy and strict b=0 initial boundary. Only then transition to full physical finite SM/Zprime+FRW and observational gates. Preserve evidence and nonmonotone Q negative in Git.
