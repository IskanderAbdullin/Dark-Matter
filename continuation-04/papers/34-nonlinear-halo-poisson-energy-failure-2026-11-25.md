# 34. Full nonlinear halo residual Poisson step: energy failure
2026-11-25 UTC. Common-L grid of paper 33, finite δ=.1 linear model, M1 K4 nL1 nphase32 nrad10 nz64 kphase4, raw C22. Advance residual Hamiltonian δΦ on full f with fixed reference streaming and evolving Poisson. Corner differences produce divergence-free canonical phase fluxes, exact particle number and L² ledgers incl boundaries; physical total energy not automatically conserved.

Prescribed δQ=.001 Q0 applied to f0: relative weighted source difference from independently integrated linear equilibrium force .59995, odd error .49588, spurious even .40856, spurious shell density source .16905 of linear-source L1. Direct full-f upwinding causes numerical diffusion of equilibrium even for arbitrarily small perturbing force. Zero δQ gives exactly zero force. Number and L ledger at rounding. Energy source .04501136 (unnormalized per Gyr), whereas equilibrium odd source gives zero.

At t=.1 Gyr, dt=.005 nominal:
| method | accepted / rejected | N error | Eref change per initial mass | δ-field energy | total energy residual |
|---|---:|---:|---:|---:|---:|
| full f |20/0|2.22e-16|1.74791290e-4|-7.45972417e-7|1.74045318e-4|
| δf plus independent equilibrium source |56/107|5.40945e-12|1.94032196e-6|-1.13187704e-6|8.08444914e-7|

Full-f error is 233 times the perturbation-field energy magnitude; cannot use for physical inference. δf evolution reduces it 215x but error still ~71% of field term; repeated positivity rejections and tolerance clipping explain number drift. δf formula is FV(δf,δΦ)+G δQ, includes nonlinear residual on δf; not a validated energy-conserving solver.

Gaussian probe changes (.5,1,2,5,10)R full-f [6.22639e-5,.00018791486,.00031670213,.00012609432,-5.55266e-6]; δf [7.567408e-5,.00023382332,.00039632855,.00017382525,1.2658836e-5]; frozen [4.831364e-5,.00014539036,.00025033282,9.34145e-5,-1.33056e-5]. Cutoff-born fraction .000445215/.000429666; unbound .0000802479/.0000802086. These are numerical comparisons, no physical density/lifetime claim.

Reproduce: python lowered_bose_halo_action.py ; repeat --deltaf. Local NumPy/SciPy, remote CI unverified. Next: positivity-preserving conservative flux limiter for equilibrium-subtracted source, then independent energy-work ledger. No change to cosmology/BEC/CLASS status.
