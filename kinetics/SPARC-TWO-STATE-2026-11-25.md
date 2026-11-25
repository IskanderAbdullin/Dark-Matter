# SPARC one/two-state equilibrium profiles: held-out radial diagnostic

2026-11-25. This is a deliberately conditional equilibrium profile comparison, not a full galaxy dynamics, observational likelihood or universal particle-mass constraint.

## Data and methods

Downloaded official SPARC tables:
- https://astroweb.cwru.edu/SPARC/MassModels_Lelli2016c.mrt (SHA256 9108994b12cc401b94a1768beca61c53ec354779385c9c9cc571049f3043244c)
- https://astroweb.cwru.edu/SPARC/SPARC_Lelli2016c.mrt (SHA256 5aa0501f6b0d881fa579030e315e7b5b6ef561a5bd3a07472f9929c7e5728243)

The mass-model table matches the earlier mirror exactly after line-ending normalization. All 175 galaxies/3391 points, selected quality <=2 and inclination >=30 deg: 153 galaxies/3168 points. Train on three of four ordered radii, hold out index mod 4=0; 2322 train/846 test. Radial points within each galaxy are correlated, so this is not a formal independent-sample likelihood. Also reproduce the earlier selected 20 dark-dominated galaxies (353/131 points) and repeat that selection with baryon/distance/inclination nuisance priors.

Reference: Lelli, McGaugh & Schombert, AJ 152,157 (2016), https://doi.org/10.3847/0004-6256/152/6/157 ; original Bose proposal Abdullin & Popov https://doi.org/10.1088/1475-7516/2021/11/055 (https://arxiv.org/abs/2107.06237) . The two-state assumption uses common temperature, chemical potential and potential with internal gap d=Delta m/T from the soft-mixing benchmark sigma=.3; at other m we use the stipulated scaling Delta m proportional to m^2.5, not a re-derived viable EFT. Density proportional to Li_3/2(exp(-alpha-y)) + Li_3/2(exp(-alpha-d-y)). One state omits second term. No finite truncation, condensate, two-fluid transport or halo formation. R is determined by m,V0,central fugacity; solve log-radius Poisson with RK4 table, interpolate log V^2. 14 independently integrated checks max relative V^2 error .001351 coarse, .000318 fine. At d=0 exact equivalence at m_two=m_one/2^(1/4), i.e. 72.3 eV for m_one=86 eV. This is an identifiability issue, not new particle-mass inference.

Fit two halo parameters for each galaxy at each *global* fixed m. Compare PIS and Burkert; separate classical isothermal limit (Lane-Emden). Fixed disk/bulge M/L .5/.7, gas signed; prior pilot uses Gaussian .5+-.1 and .7+-.14 plus quoted distance and inclination errors, +/-3 sigma, physical bounds. Distance rescales radii and baryonic V^2, incl rescales projected predicted velocities. Same nuisance choices for all models. No covariance, selection likelihood, systematic error floor, or population priors. Multistarts and explicit critical boundary; initial optimizer sometimes failed (notably two200 train1476 vs robust868 on pilot). The old Node code rerun gives PIS 550.7877/252.9488, Burkert 961.8378/411.5207 and one-state 86 1749.9502/616.8263, all reproduced to small library differences. Preserve failed/local-min results.

## Results

Quality 153 galaxies, fixed nuisance; chi-square sums (train/test):

| Model | mass eV | train | test |
|---|---:|---:|---:|
| PIS | - | 10008.24 | 3990.42 |
| Burkert | - | 12712.99 | 5069.50 |
| One | 74 | 14072.04 | 5400.11 |
| One | 86 | 13147.7 | 5110.60 |
| One | 120 | 12207.04 | 4831.46 |
| One | 200 | 12021.59 | 4826.30 |
| Two | 74 | 13138.96 | 5108.16 |
| Two | 86 | 12584.21 | 4917.92 |
| Two | 120 | 12113.68 | 4816.63 |
| Two | 200 | 12028.41 | 4832.85 |

All fits return finite feasible candidates; optimizer global optimality not certified. Mass scan selected by TRAIN objective only, not holdout. 200 eV is the best scanned mass for either Bose model, not a measured optimum or lower limit. At 200 eV no critical-boundary fits. High-mass limit tends to classical profiles, making mass weakly identifiable. Two-state 74 nearly maps one-state 86 as expected. At fixed 86, improvement of two-state over one-state is aggregate delta test chi2=-192.68, 76/153 galaxies better, 77 worse; paired galaxy-resampling per-point delta descriptive interval [-.656,.043]. At train-selected 200, two vs one aggregate +6.55, 65 better/88 worse. PIS vs both is better in aggregate; one200 minus PIS +835.88, two200 minus PIS +842.43. Galaxy-paired descriptive 95% bootstrap intervals for these per-point differences [.154,1.971] and [.167,2.013]. Not calibrated confidence intervals or a proof of universal preference; totals can be driven by individual galaxies.

Robust 20-galaxy pilot: one86 1750.83/617.12; two74 1758.03/616.12; two86 1358.57/479.74; one200 853.72/354.57; two200 868.21/359.75; PIS 550.79/252.95. With fine library+nuisance priors: one86 1528.87/536.96; two74 1528.20/538.53; two86 1112.64/374.82; one200 633.74/274.75; two200 647.04/284.75; PIS 500.09/233.84; Burkert 651.10/261.44. With nuisance priors on this small selected sample, paired intervals for one200/two200 relative to PIS include zero; do not overstate the fixed-nuisance result.

## Conclusion and next priorities

The original 86 eV calibration cannot be imported into the two-state kinetic equilibrium. Under these assumptions, adding the second state changes central density almost as a mass rescaling; the aggregate fits prefer less quantum-degenerate profiles on the scanned range, without proving any universal mass exclusion. Need full quality-sample nuisance scan, realistic finite halo, nonthermal or chemically incomplete composition, halo formation/transport, selection and independent-galaxy validation. Keep PIS and classical isothermal controls separate.

Reproduce with kinetics/sparc_two_state_pilot.py and kinetics/sparc_paired_diagnostic.py. Archived full per-galaxy fits and source hashes: kinetics/results/sparc-two-state-2026-11-25.json.

## Classical isothermal control

Quality sample Lane-Emden/classical alpha=36 control: train12085.156/test4901.923. Prior pilot: train669.298/test284.803. Quality one200/two200 minus classical holdout delta -75.624/-69.075, descriptive per-point paired intervals [-.289,.019]/[-.281,.021]. Thus there is a modest aggregate preference over the classical limit but no robust distinction under this resampling diagnostic; do not claim exact degeneracy at finite m. PIS remains better in aggregate. Prior pilot one200 improves classical by 10.05, while two200 is almost identical to classical. All classical per-galaxy fits archived with paired comparisons.
