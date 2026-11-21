# Continued investigation — critical mean-field limit, 22-galaxy re-audit, publication genealogy
Date: 2026-11-21. Status: new analytical + independently executed numerical work; **not an observational refit and not a new cosmological constraint**.

## Outputs
- [06-critical-mean-field-regularization.md](06-critical-mean-field-regularization.md): self-consistent finite-\(g\) central limit; mathematical proof, crossover and physical estimate.
- [07-jcap-table-reanalysis.md](07-jcap-table-reanalysis.md): reconstruction of the published 22-row table and assessment of the meaning of \(m=86\) eV.
- [08-related-primary-publications.md](08-related-primary-publications.md): 2019 and 2025 papers by the original author group + closely related 2012/2022 work.
- [research-log.md](research-log.md): research chronology, validation, unsupported steps.
- [numerical-summary.json](numerical-summary.json): calculated values and documented assumptions.
- [code/mean_field_regularization.py](code/mean_field_regularization.py): reproducible potential-based finite-interaction ODE.
- [code/parse_jcap_table.py](code/parse_jcap_table.py): parse source 22-row LaTeX table; requires original uploaded .tex.

Run: pip install numpy scipy mpmath; python code/mean_field_regularization.py; python code/parse_jcap_table.py /path/to/rv2_JCAP_041P_0721.tex

## Important note
Full scientific novelty relative to all literature is **not established**. SPARC full mass-model file has an independently verified public source (https://astroweb.cwru.edu/SPARC/MassModels_Lelli2016c.mrt) but was not downloaded into the execution runtime, and NO SPARC observational fit was run in this continuation. Do not treat the published reduced chi-square counts as tests of the original raw observations.

Original source: I. G. Abdullin and V. A. Popov, JCAP 2021 055, DOI 10.1088/1475-7516/2021/11/055. See the project root for the original five manuscripts and 100-paper bibliography.

- [09-thermal-streaming-benchmark.md](09-thermal-streaming-benchmark.md): independently computed number-weighted ballistic free-streaming length to equality, \(z=5\), and now for a specified 86 eV real-scalar thermal relic; comparison at 29, 50, 86, 152 eV; explicitly **not** a transfer-function constraint.


- [10-collisionality-and-asymptotic-tail.md](10-collisionality-and-asymptotic-tail.md): cross-section/relaxation plausibility benchmark, low radial optical depth, and independently verified \(r^{-2}\) infinite-mass halo tail.

- [11-lowered-bose-finite-halo.md](11-lowered-bose-finite-halo.md): positive finite-mass Bose King-like profile, 4 solutions and numerical code.

- [12-sparc-data-audit.md](12-sparc-data-audit.md): real SPARC data QC (175 galaxies, 3391 radii), including signed-gas issue and per-galaxy JSON.
- [13-sparc-20-galaxy-pilot.md](13-sparc-20-galaxy-pilot.md): initial fixed-86-eV comparison against Burkert/PIS on 20 galaxies, train/holdout, detailed caveats and Node code.
- [14-lowered-bose-normalization.md](14-lowered-bose-normalization.md): new phase-space-density and mass normalization after escape-energy truncation.
- [15-sparc-quality-and-sensitivity.md](15-sparc-quality-and-sensitivity.md): official SPARC quality Q, inclination and basic sensitivity to exclusion.
- [CONSOLIDATED-REPORT-RU.md](CONSOLIDATED-REPORT-RU.md): overarching audit, level-of-evidence table, go/no-go priorities and status.
