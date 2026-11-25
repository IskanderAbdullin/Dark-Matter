# 36. Exact discrete shell energy and split work ledger
2026-11-25 UTC. Let m=W(f-f0) be signed shell mass perturbation, D(r) the cumulative shell fraction row (linear within shell in r³). A=Ic^-1 [∫ DᵀD/r² dr + 11ᵀ/R] (same 12-point shell quadrature as Feedback). U=-mᵀAm/2; particle-energy gradient per unit phase count H_i=e_i-[WᵀAm]_i/V_i. For arbitrary two states U_b-U_a=-(m_a+m_b)ᵀA(m_b-m_a)/2 exactly. Use this identity to allocate the field work of the SSPRK2 middle step to separately accumulated C22 and residual-force increments. Reference streaming is measured from the actual remap, not an independent velocity approximation. Boundary ledger includes reference energy of lost particles; exterior perturbing field work on lost particles omitted (outflow only ~5e-12 of mass at .1 Gyr).

Coarse setup paper 35:
| component | reference energy change | perturbation field work |
|---|---:|---:|
| reference stream | -3.08e-17 | -1.1587558261e-6 |
| C22 | -6.42e-18 | +1.9686626037e-8 |
| nonlinear residual force | +1.9263672212e-6 | +1.3983924544e-9 |
All per initial mass. Sum closure per step <=9.03e-19, final closure <=7.43e-19, field agrees with independent Feedback to 1e-21. Finite difference energy gradient .000182659479024896 vs exact .000182659479821593. Thus main residual is a mismatch of discrete radial current work and residual energy-face transport, not C22; collision overlap contributes ~2.5% of remaining error. Max individual collision shell-mass L1 increment 1.59e-6 per initial mass.

Spatial/phase refinement M1,nphase128,nrad16,nz128,dt=.0025: 262144 phase states, 450754 events, 40/0 steps, mass ledger 2.22e-16, limited flux fraction 4.24e-9; Eref 2.1830142940e-6, field -1.3868796844e-6, residual 7.9613460959e-7. Almost unchanged absolute residual despite substantial phase refinement. Gaussian probes [7.987107e-5,.0002466662,.0004147492,.0001773358,1.78742e-5]. Therefore phase refinement alone does not fix energy consistency; coarse binding quadrature and force parity are next controls.

Reproduce: python lowered_bose_halo_energy_ledger.py; python lowered_bose_halo_action_limiter.py --nphase 128 --nrad 16 --nz 128 --dt .0025. Local runs, remote CI unverified. No physical convergence claim.
