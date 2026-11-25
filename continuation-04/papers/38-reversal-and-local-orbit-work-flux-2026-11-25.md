# 38. Radial reversal and local current-work flux
2026-11-25 UTC. For spherical collisionless force L anticommutes with phase reversal R. Replace high-order residual energy flux by (EF(z)-reverse EF(Rz))/2 and phase flux by (TF(z)+reverse TF(Rz))/2. Initial attempt had tiny negative net boundary outflow -6.74e-12; enforce low-order one-sided vacuum flux at open phase and energy-domain edges. Keep positivity limiter and exact odd equilibrium source.

For cellwise constant residual occupation z, exact center-orbit power under piecewise shell potential is -de_eff*w*z*[δΦ(r_hi)-δΦ(r_lo)], since dt/dθ=T and phase count de_eff*w*T*TIME*dθ. Partition each orbit's work among its neighboring internal energy interfaces by the available half-widths; divide by center-energy gap to obtain conservative face-count flux. One-sided partition at edge, isolated energy rows explicitly logged. Retain parity phase flux; not a proof of full Hamiltonian consistency. No global energy fitting. The energy moment of this local allocation matches the independent center-orbit power to absolute 1.08e-19 (isolated contribution negligible); this does not imply exact discrete field work.

Same coarse .1 Gyr setup:
| variant | Eref | field energy | residual | N ledger |
|---|---:|---:|---:|---:|
| limiter |1.9263672e-6|-1.1376708e-6|7.8869641e-7|-1.11e-16|
| parity + vacuum |8.9546539e-7|-1.1406115e-6|-2.4514608e-7|0|
| local center-orbit work |1.2724524e-6|-1.1397383e-6|1.3271419e-7|2.22e-16|
Parity fine nphase64/nrad12: residual -3.3941922e-7, N 2.22e-16. All 20/0 steps; positive outflow after boundary fix (~4.8e-12 coarse). Work limited fraction 3.87e-8. Work Gaussian probe changes [7.6110493e-5,.0002348776,.0003975665,.0001738734,1.2738091e-5], cutoff .000429672, true-unbound .0000801165. Remaining error ~12% of field term: local current control improves but not yet quantitative energy conservation.

Code also includes experimental adjoint variant: evaluate discrete shell energy gradient from paper 36, pair its phase differences with infinitesimal limited reference-streaming flux, distribute opposite work to internal energy faces. Results pending at checkpoint. Reproduce python lowered_bose_halo_action_variants.py --mode parity (or work/adjoint); fine --nphase 64 --nrad 12. Local prototypes; no remote CI or physical convergence claim.
