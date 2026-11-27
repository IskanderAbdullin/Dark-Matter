# Finite-time relaxation of the explicit parent in frozen occupations

The [pole factorization check](INDEPENDENT-EARLY-POLE-SCALAR-POLE-FACTORIZATION-2026-11-27.md) compares the archived and stationary parents. Retain the actual initial parent this time, and solve the single-channel rate equation analytically with the archived daughter and scalar populations held fixed. Use the same restricted incoming/outgoing daughter support as before. For each parent energy the stationary fraction is f_*, the fixed lab relaxation rate is r, and the surviving initial offset is `(f_initial−f_*) exp(−r t)`. Integrate its effect on the tagged scalar's gross and net rates in time, using a series for the integrated exponential near zero. An independent Radau integration matches the endpoint and integrated response of twelve manufactured decay/return controls.

| Archive, q=4.1e−6 | Initial net rate / Href | Final net rate, exposure 7.2e−5 | Integrated net vs an entirely frozen parent |
| --- | ---: | ---: | ---: |
| Fine | +4.42151 | +4.37711 | −0.503% |
| Third | +3.79437 | +3.75789 | −0.482% |

Both rates retain their initial positive sign. At exposure 2e−8 the changes to integrated net in this model are about −0.00014%; using the stationary parent from the start would instead give an almost vanishing negative net for this mode. At 7.2e−5, changing radial quadrature order 14→32 changes the reported percentage reduction by less than .001 percentage points in both archives. Larger sample exposures are supplied only as tests of this frozen-bath formula. Depletion of the other occupations, charge scattering, a self-consistent medium width, and cosmological expansion are not modeled, so these are not corrections to impose on the five-channel trajectory or estimates of final abundance.

From `kinetics/`:
```sh
for g in fine third; do
  for n in 14 32; do
    OPENBLAS_NUM_THREADS=1 python early_pole_archive_scalar_pole_factorization_reference.py "$g" \
     --order "$n" --q 0 2e-6 4.1e-6 8e-6 --exposure 0 2e-8 2.4e-5 7.2e-5 .001 .01 .1
  done
done
```
Raw: `results/early-pole-archive-scalar-parent-transient-2026-11-27.json.gz`.
