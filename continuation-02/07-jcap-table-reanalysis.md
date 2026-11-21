# Source-table audit — How informative is the 86 eV preference?
**Scope**: transcribed using a parser from the original uploaded LaTeX source, Table “The best-fit parameters for the 22 LSB galaxies.” This is a cross-model comparison of **published reduced chi-squares**, NOT re-fitting observation points.

## Exact-descriptive findings
- 22 rows in the source table. Two galaxies — **F571-8** and **UGC-11648** — were excluded from the common-mass minimum scan, which therefore used **20**. Their 86 eV reduced chi-squares in the pure model are 3.01 and 6.43 respectively; a phenomenological Bose-star population (with imposed alpha=gamma) lowers these to 0.92 and 0.46 in the original article.
- **13/22** source rows report \(\nu\ge2.6\) when rounded, i.e. at or very near the upper Bose-degeneracy boundary \(\zeta(3/2)=2.612375...\). In **3/22**, optimum is \(\nu=0\), so the particle mass is not determined by that galaxy's isolated optimum. Of the remaining rows, 19 list a numerical “best-fit mass”; these cover 29–152 eV in the published table.
- At FIXED \(m=86\) eV, the source reduced \(\chi^2\) is strictly smaller than that of the **isothermal sphere for 12/22** and equal for **6/22**. For the selected 20 only, it is strictly smaller in 10/20.
- At \(m=86\) eV, the source reduced \(\chi^2\) is strictly smaller than **Burkert for 6/22**, and smaller than **pseudo-isothermal for 5/22** (as printed).
- Median per-row difference \(\chi^2_{r,86}-\chi^2_{r,IS}\) across the 22 is **−0.0055**; the corresponding median difference versus Burkert is **+0.0200**. They are descriptive differences of dimensionless reduced statistics, not meaningful aggregate likelihood deltas.
- These counts do not account for the number of independent observations per galaxy, degrees of freedom, estimated distance/inclination uncertainty, baryonic templates or covariance. Models have differing parameter numbers. In particular **one MUST NOT sum reduced chi-squares and call it a total chi-square**.

## Interpretation
The article's good fits show an interesting phenomenological representation, but this table alone cannot establish a statistically sharp universal particle mass or superiority to two-parameter empirical cored profiles. A large fraction of optima lying on the physical \(\nu\) boundary calls for proper profile likelihoods with boundary-aware uncertainty treatment and sensitivity to priors. Post-hoc exceptions treated with an extra component must be compared on the same footing as the rest.

## Next decisive observational analysis
1. Obtain ORIGINAL 22 Hα velocity points including covariances and original sample metadata; independently verify table values, model \(N_{\rm dof}\), and the exact global objective minimized. Avoid confusing these data with the different SPARC sample.
2. Separately download SPARC/Lelli+2016c mass models. Include gas, disk, bulge, stellar mass-to-light prior, distance/inclination nuisance. Source: https://astroweb.cwru.edu/SPARC/MassModels_Lelli2016c.mrt (public text, headers inspected); paper https://arxiv.org/abs/1606.09251.
3. Fit all galaxies, compare Burkert, PIS, NFW and Bose on joint log-likelihood and predictive held-out data. Pre-register exclusion and condensate criteria, rather than fitting exceptions after inspecting residuals.
4. Global \(m\) uncertainty must include the effect of exact boundary solutions, outliers and model mismatch; present posterior predictive checks and leave-one-out diagnostics.
5. Independently test the nonlinear relation between per-galaxy parameters \(\nu,R,V_0\) and shared \(m\), rather than treating \(m\) as an independently free third parameter for each galaxy.

**No actual SPARC fit was executed in this work**; the official dataset was located and its column schema verified, but downloading the full file to this runtime was unsuccessful. No new confidence interval for 86 eV is claimed.