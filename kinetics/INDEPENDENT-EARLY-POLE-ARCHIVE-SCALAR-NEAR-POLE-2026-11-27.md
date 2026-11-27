# Disjoint guarded Yukawa annuli approaching the parent pole

Extend the [below-pole tagged vacuum fold](INDEPENDENT-EARLY-POLE-ARCHIVE-SCALAR-BELOW-POLE-2026-11-27.md) to four disjoint invariant annuli on either side of the parent mass. The inner excluded region is `abs(sqrt(s)/M - 1) < 1e-5`; the outer radius is .1. Both fermion legs remain in the same [0.45,5] radial band. Neither the real parent pole nor a finite-width medium enters these archived, frozen-input results.

The inverse-pole coordinate 1/(s−M²) resolves the enhanced angular weights in each finite annulus (the radial energy at q=0). An independent order of integration evaluates the rescaled density first at fixed x=s/M²−1, then integrates it in 1/x. The two routes agree within 6e−6 relative for the dominant exterior net of the two nearest annuli and within 0.003 for all gross terms on those same annuli. The tiny interior net at finite q is not assigned this precision. A direct Dirac trace checks points inside the nearest annuli to 2.3e−15 relative, and the FD/BE packet balance has a defect below 1e−14.

| Archive | Tagged q | Sum over guarded annuli, exterior net / Href | Largest single gross annulus / Href | Matched direct source / Href |
| --- | ---: | ---: | ---: | ---: |
| Fine | 4.1e−6 | −1.237e−13 | 5.162e−6 | +4.422 |
| Third | 4.1e−6 | −2.044e−13 | 9.113e−6 | +3.794 |

The central-partner net is smaller at this finite tagged momentum. The exact-rest controls are retained separately: their archived first-cell scalar occupation is extrapolated to q=0 and they do not describe a population of exact-zero modes. For both finite-q archives, the closest annuli increase the gross rate by roughly an order of magnitude when moved one decade inward. The regularized invariant density tends to a nonzero value at the pole; that is **not** the physical on-shell scattering rate. In particular the small guarded net is not an estimate or bound for the excluded neighborhood. The choice of real-intermediate-state separation is tied to the kinetic on-shell approximation; see [Ala-Mattinen et al.](https://doi.org/10.1103/PhysRevD.108.096034). A matched resonance/decay calculation must still account for the excluded width without double counting; no such result is claimed.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for g in fine third; do
  for o in 8 14; do
    python early_pole_archive_scalar_below_pole_finite_q_reference.py "$g" --q 0 4.1e-6 \
      --order "$o" --near-pole-bands > "$g-annuli-$o.json"
    python early_pole_archive_scalar_below_pole_finite_q_reference.py "$g" --q 0 4.1e-6 \
      --order "$o" --pole-probe 0 -0.002 0.002 -0.0002 0.0002 -0.00002 0.00002 -0.000002 0.000002
  done
  python early_pole_archive_scalar_pole_cross_reference.py "$g" "$g-annuli-14.json"
done
```
Results and threshold guards: `results/early-pole-archive-scalar-near-pole-2026-11-27.json.gz`.
