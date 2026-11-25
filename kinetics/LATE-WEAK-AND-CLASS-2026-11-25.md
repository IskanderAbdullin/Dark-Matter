# Late weak-parent candidate and local CLASS validation (2026-11-25)

## Late weak-parent kinetic fit
M1=86 TeV, m2/M1=0.9, m_phi=86 eV; x0=0.1, initially empty parent, kappa_SM1=47, kappa2=0.03, kappa=0.005278180162600994. Width was fit on 160 cells and held fixed on 240/320 cells. Omega h^2=0.1207153994 and 0.1209727907, hot fraction q/T_actual>0.1=0.003565679 and 0.003556194. Nominal uncancelled seesaw mass ~0.05024 eV. Equilibrium thermal scalar mass/peak ratio ~0.00359 (contrast early weak case ~0.26); this is only a simple-Yukawa thermal diagnostic. Independent bath-only shape screen: peak TV distance 0.01818, energy discrepancy 2.04%, number discrepancy <0.1%. It excludes dark depletion and cannot establish closure in the stimulated regime. Unprotected loop scale ~28 MeV versus 86 eV; protective symmetry remains essential.

## Local CLASS-community 26.0.0
For the earlier strong thermal-parent candidate (not the late weak candidate), actual synchronous linear transfer was calculated at z=0,5 and k=0.003--20/Mpc on both 240 and 320 cell PSDs. At z=5 k=10/Mpc the ratios P/P_CDM are 0.989236364 and 0.988734742. Fine-mesh minimum sampled ratio 0.983762068 at z=5; no half-power crossing in range. Independent three-tracer partition drift <=0.44%. These conditional collisionless linear results do not establish viability or observational likelihood. Small apparent large-scale enhancement is not physically validated.

Thermal Bose reference comparison with archived official CLASS 3.4.3: half-amplitude crossing z=5 at k=0.42961044 vs 0.43275181/Mpc (0.726% shift); abundance and early radiation agree. One deeply suppressed thermal mode at k=0.78978243 has invalid power ratio; recorded as null, not a physical value.

## Negative results and next action
Late weak candidate PSD export failed transformed quadrature tolerance: third tracer moment errors -0.00297, -0.02920, -0.000316. No CLASS spectrum for that candidate. Repartition/rescale tracer or refine quadrature and check moments before solver. Full hierarchy alternate run failed because CLASS-community 26 renamed tolerance input; adapter patched to source-confirmed tol_perturb_integration, not rerun. Prioritize late weak candidate with momentum-resolved parent including depletion and self-consistent finite-T self energies. No observational viability claim.

Raw numerical records: results/late-weak-and-class-2026-11-25.json. Local environment and public research branch preserve all completed calculations. Main not updated (server-side ref error); workflow draft not activated (missing workflow scope).

## Quadrature follow-up

Both late weak-parent meshes and both partitions fail 120/240-node moment checks. All four pass at 480 nodes (max relative error 2.28e-5). Added explicit --bins 480 option with 960-node alternate check; analytic export passed both meshes. Recorded full checks in kinetics/results/late-weak-quadrature-2026-11-25.json and kinetics/results/late-weak-export-2026-11-25.json. The earlier failed coarse check is preserved. Next run actual CLASS at these validated settings; transfer and full hierarchy remain untested for this candidate.
