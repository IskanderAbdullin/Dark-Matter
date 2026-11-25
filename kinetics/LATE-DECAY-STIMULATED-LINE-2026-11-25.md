# Late decay alternative: stimulated emission and inverse decay

Date: 2026-11-25. Continuation of fa5bb52. This is a diagnostic of one proposed microphysical mechanism, not an established relic model or an observational exclusion.

## Setup and reproducibility

A cold scalar A is produced with the previously computed collisionless PSD, then converts to a self-interacting scalar B and massless eta by A -> B + eta. Use mA=86 eV, mB=sqrt(mA^2-2 mA pstar), vacuum lifetime 10 Gyr, and pstar from 1e-9 to .1 eV. The reference backgrounds use H0=.674*100 km/s/Mpc, Om=.315, Or=9.2e-5. Inputs: parent_angle240_k00685_history.json and parent_angle320_sparse_k00685.json, from the matched-width kinetic work. They retain their original occupations and absolute normalization. Local nonrelativistic pilots truncate parents at 20 times the peak momentum (97.72/97.62% of number); hot tail and earlier production/relativistic epochs not evolved. Source code and raw snapshots are committed. Reproduce from repository root with .venv/Scripts/python.exe (Windows) or equivalent Python environment; commands in each script's CLI.

## Two deliberately incomplete eta-empty limits

late_decay_bose_screen.py: conserve A+B number using radial recoil probabilities P(qB|qA)=qB/(2 qA b). Ignore inverse, eta stimulation, B scattering and media. Original outputs had invalid initial geometry trial and bad age integration at the wide integration interval; retained as failed numerical evidence. Corrected clipped-cosine bins have maximum row error ~4.4e-16 and log-age integration. At pstar=.001 eV, 240/320 parent inputs, decay fraction ~.811/.814 at r=1+z=10001; first 1% at r~193861/185624, 50% at 40051/41244. At .01 eV, ~.795 at r=10001. These are not actual cosmological fractions. No general bound at finite recoil; eta inverse term can have either sign.

late_decay_instant_equilibrium.py: opposite instantaneous ideal BE+condensate closure for B, still omitting eta. Fix B number and kinetic energy, include angular BE enhancement and approximate condensate source Nc*fA(b). At pstar=1e-8/1e-6, 240/320 input grids give conversion 1% near r~1860-1890, 50% near r~314, 90% near r~81; final condensate ~99.93% of retained parent number. This is a closure, not evidence for physical formation of a condensate or thermalization. Number residual at machine precision; thermal and chemical assumptions remain uncontrolled.

## Exact invariant line gain

For isotropic A/B gas, massless eta, cubic amplitude g, vacuum decay Gamma_A=g^2*pstar/(8*pi*mA^2) in natural units. Collision term for an eta mode k is

C_eta = g^2/(16*pi*k^2) integral_{Emin(k)}^infinity dE_A [f_A*(1+f_B) + (f_A-f_B)*f_eta],

where f_B is evaluated at E_B=E_A-k and Emin=mA*cosh(ln(k/pstar)). Integrating the coefficient of f_eta over log k / H at frozen background yields exactly

tau_gain = 2*pi^2*Gamma_A*(n_A-n_B)/(H*pstar^3).

Proof: for fixed parent momentum the k interval is pstar*gamma*(1+-beta). Integrating k^-3 dk gives 2*p*E/(m^2*pstar^2). The analogous inverse-decay interval uses pstar_B=(mA^2-mB^2)/(2*mB), so mB*pstar_B=mA*pstar; hence the same coefficient multiplying n_B. This includes exact time dilation. Actual redshifting mode uses k=pstar*r/rc and q_min=mA*abs(1-(rc/r)^2)/(2*p0*rc). Independent split-bin quadrature and 12/24-point refinement agree within roundoff for the tested truncated PSDs. Along a mode with B=0 frozen and eta initially empty, f_eta=expm1(tau); huge gain is a failure diagnostic for freezing, not a physical number prediction.

At r=1e5, pstar=.001 eV, both inputs give tau~3.25e4; pstar=.01 gives ~32.5. At r=1e6 .001 gives ~3.30e5. For a mode crossing at r=1e8 starting at r=1e8 (partially traversed line), .001 gives ~1.61e7. This severe issue arises even after production; the ten-Gyr lifetime alone cannot set conversion time. Initial eta occupation is not consistently negligible. Finite Doppler width, expansion, medium dispersion and inverse reactions matter.

## Nonlinear local burst with inverse decay

late_decay_local_burst.py explicitly evolves A,B,eta with quantum net flux f_A*(1+f_B)+(f_A-f_B)*f_eta. Exact relativistic lab daughter kinematics and time dilation; finite cell projection conserves A+B and eta-B, with measured (not exact) energy error. Original implicit integration followed an unphysical negative branch when the seed was below absolute tolerance; saved as failed. Corrected analytical frozen-reservoir positive seeding at 1e-12 or 1e-14 of initial number and a resolved initial growth step. At r=1e5, pstar=1e-5 eV, two independent grids (100/120 and 160/240 parent/eta bins, GL 4/6) yield B=.4904/.4905 after .002123/.002128 seconds. Time/peak Doppler crossing ~7.6e-9; two number invariants at roundoff, energy error relative to initial kinetic energy ~9e-8/1.1e-6. No significant negative occupations. At pstar=.001, n100 gives .4901 after 1589 s, time/crossing=.00569, kinetic-energy projection error .00289. Fine-grid comparison and elastic feedback are separate tests. At pstar=.01 the local approximation cannot follow the entire Doppler passage: only ~2.3e-10 conversion in one peak crossing in the frozen-volume model; this is not a bound on the expanding line.

## Competing self-interaction and medium effects

late_decay_medium_screen.py evaluates existing independent NR C22 net moment integrals on B snapshots, sigma/m=.3. Integrated absolute fractional-moment rate screens for p4, N(q>2peak), N(q>6peak) are (3.74e-6,3.49e-6,7.48e-6) at pstar=1e-5, but (3.08,.597,6.72) at .001. Thus elastic feedback is important for the latter and not justified to neglect. This does not establish a microscopic scattering time for every mode. Initial `lambda*f` / collective effects also need checks; weak vacuum coupling alone is insufficient.

Illustrative NR coherent chemical potential mu=g_NR*nB, g_NR=lambda/(8*m^2). A toy instantaneous density cap nB<=delta/g_NR (delta=mA-mB) could block conversion at high density but is not derived from a nonequilibrium dispersion calculation. At sigma=.3, r=1e8, all-B mu/pstar~3578,35.8,3.58 for pstar=1e-5,.001,.01. Vacuum coherent 4->2 loss per particle/H ~38.7 at full B there. If B instantaneously follows that putative cap, full conversion would become allowed at r~6.54e6,3.03e7,6.54e7 respectively; integrated vacuum coherent 4->2 loss fractions (no feedback) ~6.77e-8,.00314,.630. Last is not small and invalidates ignoring number loss. B-only contact-pressure sound horizon to equality on this toy track ~6.94e-5,1.39e-4,1.78e-4 Mpc; not total sound horizon, transfer function or constraint. Kinetic pressure, hot tail, two-fluid response, condensate formation and eta in-medium mass/width excluded. Small recoil branch merits real dispersion+inverse kinetic treatment, not a claim of viability.

## Literature and next steps

- Decaying dark matter with self-interactions: https://doi.org/10.1103/PhysRevD.102.043024
- Halo response to decays: https://doi.org/10.1103/PhysRevD.111.103522
- Structure constraints on DDM: https://doi.org/10.1088/1475-7516/2015/07/009 ; https://doi.org/10.1103/PhysRevD.93.023510
- Satellite constraint context: https://www.ornl.gov/publication/milky-way-satellite-census-iv-constraints-decaying-dark-matter-observations-milky-way (not transferable directly to stimulated Bose conversion).
- Vacuum scalar number-changing normalization: https://doi.org/10.1007/JHEP08(2019)050 ; Bose-stimulated production: https://doi.org/10.1007/JHEP03(2025)202.

Priority: evolve the coupled B elastic, eta redshift/inverse and in-medium mass system through postproduction; test cutoff, coherent limit, energy conservation and chemistry before computing transfer/halo observables. Compare with soft mixing; do not fit abundance using these incomplete closures.

Raw evidence: kinetics/results/late-decay-stimulated-line-2026-11-25.json (including failed preliminary variants). Six regression tests for this and halo library passed.
