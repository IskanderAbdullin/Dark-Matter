# The held scalar self-scattering proxy changes the source differently across prepared inputs

Add an isolated shape-evolution test to the [frozen screened scalar stencil](INDEPENDENT-EARLY-POLE-SCALAR-VERTEX-RESPONSE-2026-11-26.md). The live scalar cells cover [1e−6,1e−5]. The scalar-only two-to-two stencil has a reciprocal positive weight, with one fixed initial dispersion and a large-N-inspired vertex held at the initial occupations. Both charged populations are held. At each sampled shape, evaluate the dressed cold decay/inverse source separately for the native step input and the 97% parent witness. For this diagnostic only, the scalar Hartree mass is reevaluated from the changed occupations. The decay source is not put back into the shape ODE.

| Held charged input | Fine source, initial → final | Third source, initial → final |
| --- | --- | --- |
| Native steps | .079252 → .078661 | .072599 → .070037 |
| 97% prepared parent | −.006823 → −.004613 | +.065176 → +.039395 |

Values use 16 scalar bins, ten angular/azimuth nodes, three hard subdivisions and final reference exposure 2e−7. The shape L1 changes are 0.836 and 0.718 of the initial scalar number; number drift is below 1e−10 relative. These L1 changes quantify redistribution between bins, not lost scalar number. Fixed-dispersion kinetic energy closes below 1e−10 relative; reevaluating the Hartree energy changes it by at most 1.3e−15 absolute across the tested numerical fixtures. Entropy increases. The initial shape clock itself is sensitive to momentum sampling, while comparing order six to ten changes final q16 source values by less than 0.2% of the respective magnitudes. A q32 fine check still has the negative witness source at both ends. The coarser q8 fine witness starts slightly positive and changes sign, illustrating why one should not infer the sign from that projection alone.

The held vertex is a proxy, not a controlled N=1 in-medium amplitude. This is a shape sensitivity measurement in a closed IR support, with no elastic feedback from the hard sector and no simultaneous decay, photon processes, or expansion. It is not a relic-density result or a correction to add unchanged to a coupled trajectory.

From `kinetics/` with one BLAS thread:

```sh
for g in fine third; do
  OPENBLAS_NUM_THREADS=1 python early_pole_live_scalar_shape_feedback_reference.py "$g" --bins 8 --order 6
  for order in 6 10; do
    OPENBLAS_NUM_THREADS=1 python early_pole_live_scalar_shape_feedback_reference.py "$g" --bins 16 --order "$order"
  done
done
OPENBLAS_NUM_THREADS=1 python early_pole_live_scalar_shape_feedback_reference.py fine --bins 32 --order 8
```

Raw: `results/early-pole-live-scalar-shape-feedback-2026-11-27.json.gz`.
