# Early heavy-state bootstrap and Hartree coherence: checkpoint

2026-11-25. Independent continuation from fa5bb52. Separate worktree avoids interfering with concurrently running late-decay and SPARC calculations. All bounds are conditional on their kinetic approximations and initial conditions.

## Occupancy bootstrap

For isotropic constant-amplitude relativistic 2->2 scattering, the angular reduction yields an interval length L<=2 p_external. For a fixed external species and unordered incoming pair, the gain coefficient is <= |M|^2/(64 pi^3 m_min) times I_i I_j (1+F_s+F_k), with I=int f dE, F the occupation supremum and E>=m_min. The quartic Bose term cancels exactly between forward and reverse; never bound it as an independent gain. At fixed number n and sup F, filling a sphere p_*=(6 pi^2 n/F)^(1/3) maximizes I and the inverse-energy variance V. Use I<=F p_*^2/(sqrt(m^2+p_*^2)+m), V<=min(n/m,F p_*^2/(4 pi^2)). This removes the earlier dependence on a frozen A spectral shape for the occupancy bootstrap. Integrated number source uses the bound sigma*v_M<=|M|^2/(32 pi E1 E2). All bounds deliberately ignore thresholds and reverse loss.

Before portal-off r=1e14, B-only diagonal B^4 collisions conserve B number; the coherent positivity bound is used with number<=2 Nref and the stipulated external A occupation cap. After vacuum projection, all six diagonal quartic channels included. For an assumed heavy cap F_H, integrate global direct/catalytic source and rotation, then use that number envelope to bound fixed-mode gains including HH->HH. A comparison inequality F'<=2 a sqrt(F)+G implies F<=(sqrt(F0)+int a+sqrt(int G))^2. Check first-exit at BOTH occupation caps. At seed=1e-8, sigma=.1,.3,1.25, vacuum postcut heavy caps close at 1e-4,1e-4,1e-3; maximum F_H upper estimates 3.865e-5,4.177e-5,8.003e-4; final number upper estimates 2.49e-15,1.86e-15,3.43e-15 of Nref. At larger tested seeds the deliberately coarse bound does not close; this is NOT an exclusion. Crucial caveat: the bounds on the B variance imply possible Hartree shifts >> the small vacuum A/B splitting, so the vacuum postcut assumption is not justified by the bootstrap. Do not treat this as a proof of early safety. Preportal assumptions remain external and offdiagonal collision kernels and 2<->4 are absent.

## Self-consistent positive Hartree coherence

For postportal common E=sqrt(p^2+m^2), H=(D_eff P_B+delta x sigma_x)/(2E), D_eff=D+lambda_B V_B/2, solve momentum-dependent complex amplitudes with small B amplitude scaled out. For constant x=seed the conformal energy (per comoving volume) is

 F_c = (D/2) J_B + delta x K + c J_total + lambda_B r^2 J_B^2/8,

where J_s=int (r/E) f_s d^3p_com/(2pi)^3, K=int (r/E) Re rho_AB, c=(delta x)^2/D. The local shifted quadratic form >=D f_B/4; at fixed x the explicit expansion derivative is nonpositive because r/E decreases monotonically and r^2 decreases. Hamiltonian evolution conserves the instantaneous energy at fixed r. Thus J_B<=4 F_c(initial)/D independent of coherent phases within this collisionless common-energy model. This is not a bound with arbitrary dephasing/scattering.

At sigma=.3, seed 1e-8, r=1e14->1e11, pure-A start, 240/320 production inputs with different quadrature: max B occupation 3.84e-6/3.91e-6, max Hartree shift/splitting 1.920/1.910, final heavy eigenstate fraction 3.143e-17/3.159e-17. Max norm residual <=1.9e-14, independent integrated conformal work residual over initial shift scale 3.02e-9/3.41e-10. The common-energy expansion parameter <=1.7e-12. Final Hartree shift/splitting .0383/.0451: reduced, not identically zero. Collision exposure evaluated on the coherent trajectory is ~1.01e-4 per B particle (conservative maximal rate), fixed-mode gain diagnostic ~7e-11. These are NOT a coupled collision solution or an a priori bound on redistributed B.

Alternative initial lower eigenstate and phase+ uniform B amplitude of delta x/D (more than the early coherent peak bound) at r=1e14->2e11 give max f_B 3.27e-6/5.98e-6 and final heavy eigenfractions 2.86e-17/1.83e-16. The phase test initially has Hartree shift/splitting ~3.9e4; positivity, norm, work checks still pass. Trimmed low momentum fraction below 5e-10 of number and 3e-9 of UR variance, explicitly logged. These tests support small conversion in the specified positive Hartree model, but omit primordial B, offdiagonal collisions, medium effects of the derivative portal after cutoff and later coherent evolution.

## Sources and reproducibility

- Relativistic isotropic Bose collision angular reduction and condensate: https://doi.org/10.1103/PhysRevE.93.032131 and https://doi.org/10.1103/PhysRevD.100.091501.
- 2026 BEST implementation: https://github.com/best-hep/best and https://arxiv.org/abs/2603.28848 (prior project benchmark found high-occupation MC conservation problems, so use deterministic angular bounds here).
- General density-matrix transport and coherent collisions: https://doi.org/10.1007/JHEP02(2024)217 and https://doi.org/10.1103/PhysRevD.100.083536. Their neutrino collision terms are NOT scalar B^4 collision terms.
- Run early_bose_bootstrap.py and early_hartree_coherence.py with the prior 240/320 frozen PSD outputs; full raw arrays in results/early-bose-bootstrap240-2026-11-25.json and results/early-hartree*.json.

Next: assess true coupled diagonal and offdiagonal collisions, direct extra species production, and UV protection. A subthreshold virtual mediator was screened separately as a possible way to reduce the splitting cost while suppressing halo conversion; verify the full time-dependent trigger and extra number-changing processes before calling it a candidate.
