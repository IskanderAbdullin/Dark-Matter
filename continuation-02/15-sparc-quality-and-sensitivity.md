# Official SPARC quality metadata and first sensitivity exclusion

Date 2026-11-21. Accessed original SPARC Table 1 through public plain-text GitHub mirror (SHA 8cdfb503c58bf5f4fc6321e976f7caec6e353bf3); official table https://astroweb.cwru.edu/SPARC/SPARC_Lelli2016c.mrt, paper https://arxiv.org/abs/1606.09251. Recheck source hashes against original before a formal paper.

Official quality flag: Q=1 high, Q=2 medium, Q=3 low. In 175 galaxies we find **99 Q1**, **64 Q2**, **12 Q3**; **12 inclinations <30 degrees**. A pre-specified suggested quality sample Q<=2 and inclination>=30 degrees contains **153 galaxies** (without yet applying a minimum number of observed radii).

The exploratory 20-galaxy selected pilot (note 13) contains **zero Q3** galaxies but one galaxy with inclination <30 degrees: **F563-1 (25 deg, reported error 5 deg)**. Removing this system, without refitting any other galaxy, gives an exact sensitivity subtotal for 19 galaxies:

| Model | Train chi² | Holdout chi² |
|---|---:|---:|
| Bose m=86 eV | 1744.965 | 613.837 |
| Burkert | 956.502 | 408.153 |
| PIS | 546.462 | 248.235 |

The ranking of total pilot scores is unchanged, but this subtraction alone does **not** validate overall robustness: all other target galaxies and stellar nuisance parameters remain untested under alternative selection criteria. An inclination limit is pragmatic, not a theorem. Distance/inclination measurements and quality flags were never used as nuisance priors in pilot fit.

For complete subsequent reanalysis:
- Use quality Q<=2 and inclination>=30 degrees as one preregistered primary sample, also examine all 175 without dropping failures;
- Include galaxy-level luminosity, distance D and e_D, inclination and e_Inc, quality flag, referenced observation literature, and corresponding covariances;
- Report selection effects and host galaxy distribution; keep the 20 pilot visibly distinct from the JCAP original 22;
- Profile or marginalize nuisance parameters rather than fixing disk and bulge mass-to-light values;
- Re-run likelihood and holdouts after changed selection, not just rescale the old summed chi-square.

Data and metadata are summarized in sparc-quality-metadata.json, and can be reproduced with code/sparc_quality.py using the original Table 1 text.
