# Scalar transport through nested finite boundaries

Keep the same 16 logarithmic central cells in [1e−6,1e−5] and the same interaction parameter, adding either three or six cells at the same logarithmic spacing to **each** end. Start from the archived scalar shell counts, transferred by overlap in q³, and hold the hard charged spectra fixed. The collision shape proxy uses a fixed initial mass and vertex for each separate support. Each support is closed at its own edges, so count loss from the interior represents redistribution into exterior bins, not number-changing reactions.

| Archive | Interior count change, 3 cells per wing | Interior count change, 6 cells per wing | Final central 97% witness source, narrow / widest |
| --- | ---: | ---: | ---: |
| Fine | −1.9143% | −2.2337% | −.00460722 / −.00422721 |
| Third | −1.7667% | −2.0171% | +.03934262 / +.03934029 |

At exposure 2e−7, enlarging the momentum support alters the late fine witness magnitude by 8.25%, without changing its sign; the third changes by less than 0.006%. This is a central number and source comparison on identical central bins, not a comparison of total production over different momentum ranges. The native archived central sources change from +.07864578 to +.07687222 (fine) and +.07004236 to +.06864617 (third). Central initial source values agree to within 1.6e−10 across supports. The energy difference stays below 9e−21 absolute, with exact fixed-count coordinates and monotone recorded scalar entropy.

The widest support is approximately [4.21697e−7,2.37137e−5]; 3.67% of the candidate scattering geometry still lies outside this graph and is excluded (a geometric fraction, not a physical rate bound). This control quantifies finite-boundary sensitivity for a held scalar-only proxy. There is no active feedback of the changes on the charged distributions, no convergence claim in the exterior, and no cosmological freeze-out inference.

From `kinetics/`:
```sh
for g in fine third; do
  OPENBLAS_NUM_THREADS=1 python early_pole_scalar_nested_support_reference.py "$g"
done
```
Raw: `results/early-pole-scalar-nested-support-2026-11-27.json.gz`.
