# The pole-core sign sensitivity survives a genuinely continuous daughter ramp

Replace the staircase approximations of the [two-moment ramp](INDEPENDENT-EARLY-POLE-SCALAR-LOCAL-POLE-CORE-2026-11-27.md) by a single exactly linear occupation in p³ on the same central ramp support. Evaluate it at the actual incoming and outgoing momenta in every fixed-invariant packet. Split at the ramp endpoints rather than at internal artificial subcell edges. Two small shoulder changes preserve number and energy over the same two native cells, with Pauli-admissible occupations. The constant is subtracted **only over the same local window** as before; neither the pole-normalization convention nor the signed observable changes.

| Archive, exterior-partner net / Href | R=2e−4, native steps / continuous ramp | R=2e−3, native steps / continuous ramp |
| --- | ---: | ---: |
| Fine | +1.156e−19 / −1.767e−17 | −5.208e−17 / −1.198e−16 |
| Third | +2.332e−19 / +1.767e−19 | +4.141e−18 / +9.824e−18 |

R is the radius in x=s/M²−1, tagged q=4.1e−6, and the displayed width parameter is zero after local subtraction. The finite-width control gamma/R=.1 does not alter either continuous result's sign; it is an illustration, not the actual parent width. Refining the invariant/angular order 18→28 and the window order 8→14 changes the continuous exterior net by less than .2% in these panels.

For the squared-pole numerator F, the symmetric second difference is G(x)=F(x)+F(−x)−2F(0). A slope jump at the center would make G proportional to x on one side of the symmetric window and give a logarithmic width dependence. By contrast the resolved continuous examples give an approximately constant G/x² from x=1e−5 to 2e−7: near −8.95e−14 (fine) and +7.88e−16 (third). This tests local smoothness over those resolved offsets, not the physical width limit; the larger archive-dependent result can be caused by different resolved curvature, not necessarily a cusp at x=0. A manufactured segmented hinge and an independent smooth polynomial test pass. The original finite-window split and abrupt phase-space shoulders remain modeling choices.

There is no physical scattering correction to insert on the basis of this comparison. It only shows that the sign sensitivity on fine is not explained away by using a coarse step representation of the ramp.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1
for g in fine third; do
  python -W error early_pole_archive_scalar_pole_continuous_ramp_reference.py "$g"
done
```
Raw: `results/early-pole-archive-scalar-continuous-pole-ramp-2026-11-27.json.gz`.
