# SPARC: official source and observation constraints

Official database: https://astroweb.cwru.edu/SPARC/
Public mass-model ASCII catalog: https://astroweb.cwru.edu/SPARC/MassModels_Lelli2016c.mrt

The official catalog describes 175 late-type disk galaxies. Column order:
\`ID D[Mpc] R[kpc] Vobs[km/s] e_Vobs[km/s] Vgas[km/s] Vdisk[km/s] Vbul[km/s] SBdisk SBbul\`.

Important notes verified in the catalog header:
- \`e_Vobs\` is a random kinematic uncertainty only. It excludes systematic inclination correction uncertainty.
- \`Vgas\` contains a factor of 1.33 for helium and can have a *negative signed contribution*.
- \`Vdisk\` and \`Vbul\` correspond to mass-to-light \(M/L=1\,M_\odot/L_\odot\) at 3.6 microns, so their **squared** values must be multiplied by mass-to-light ratios.
- For real inference, gas pressure, vertical disk systematics, inclination, distance, shared stellar \(M/L\) priors, and data selection must be dealt with.

How to use the starter code:
\`python code/sparc_likelihood.py path/to/MassModels_Lelli2016c.mrt\` (summarizes observations).
\`python code/sparc_likelihood.py path/to/MassModels_Lelli2016c.mrt --galaxy F583-4\` (fixed-parameter example **NOT optimized or calibrated**).
\`cd code && python -m unittest -v test_sparc\` (four synthetic unit checks).

The plain-text catalog has 3391 observation rows by source line count (3416 lines including 25 header lines); code must *verify its own parsed count* when executed. In the research environment external binary/text file download was blocked; it was accessed and verified from the web but was **not** saved or fitted. No claimed SPARC best-fit mass or \(\chi^2\) exists from this project yet.
