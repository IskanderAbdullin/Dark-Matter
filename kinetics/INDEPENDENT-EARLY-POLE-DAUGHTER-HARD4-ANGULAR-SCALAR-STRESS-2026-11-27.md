# Scalar stress under the combined hard and angular refinement

Evaluate the [actual hard4 angular78 pair](INDEPENDENT-EARLY-POLE-DAUGHTER-ANGULAR-HARD-CROSS-2026-11-27.md)
with the same held-scalar stress, homogeneous dilation, instantaneous
affine-shear and hypothetical bare-mass controls used for the earlier
three fixtures. The new terminal pair changes the retained kinetic
pressure increment by +.51458%, close to +.51547% at hard4 cold32.
The relative change of the full conformal scalar pressure is only
+8.428e−9. With zero bare mass, the full scalar sector still
satisfies P=ρ/3; slow retained modes alone do not classify the
stress of the whole scalar population.

At the last represented hard support in the earlier frozen geometry,
the depth-four angular fixture gives w=.0090980,
(dP/dρ)_uniform=.0169522 and affine shear slope S/P=.33623.
For the frozen initial source, the first sign-change root is now at
bare mass .0995114 times the initial massless effective gap,
versus about .0994995 on hard4 cold32. These are diagnostic sweeps
with fixed populations or geometry, not a physical bare mass or
a simulated trajectory with changing mass.

## Where the instantaneous shear comes from

Keep the two terminal retained populations fixed, set bare mass
to zero, and multiply every prescribed exterior spectator density
by a common factor ζ. Re-solve the Hartree gap in each arm.
At the original ζ=1, the retained aperture supplies 99.99988% of
the scalar-density moment that sets the gap, but only .211%
of the total affine shear slope. Removing all exterior spectators
changes the absolute local-minus-native shear by less than .06%,
but makes the same paired signal roughly 474 times larger relative
to the total shear. At a spectator factor near .00211, the
retained and exterior contributions to the shear slope are equal.
At every tested factor the full scalar stress stays conformal.
An independent finite affine deformation checks the shear calculation
at each sampled factor.

This reweighting identifies why the ratio is small in the
original model. It does not approximate dynamical removal of
spectators, infer their physical abundance, or predict a
propagating stress, Jeans scale or transfer function.

Replay from the kinetics directory:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error \
  early_pole_daughter_fourth_scalar_stress_reference.py \
  results/early-pole-first-window-hard-depth4-cold32-initializer-2026-11-27.json.gz \
  results/early-pole-unsplit-hard4-angular78-{native,local-daughter}-2026-11-27.json.gz \
  results/early-pole-daughter-shared-cold-angular-sampling-2026-11-27.json
```
[Full receipt](results/early-pole-daughter-hard4-angular-scalar-stress-2026-11-27.json.gz).
