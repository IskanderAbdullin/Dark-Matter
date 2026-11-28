# Kinematic access to the retained on-shell channel

For fixed scalar and fermion momentum magnitudes, the maximum pair
invariant is attained at opposite directions. Across finite shell
intervals, the largest invariant is at the two upper boundaries;
the smallest is obtained from their closest rapidities in the
same direction. These exact bounds give the list of kinematically
open pairs independently of the prompt angular quadrature. In the
tested q16 and populated embedded q32 meshes, the resulting open
pairs match the sampled support at the audited values of b.

On the 12 MeV archive, the scalar lower floor for any retained
on-shell pairing is q/T_R=1.255 with a maximum retained q/T_R of
16, and q/T_R=0.629 with a maximum of 32. About 12.6% and 2.45%
of the two archived scalar counts lie below their respective
floors. These are counts in fixed snapshot shapes, not integrated
transfer probabilities. Redshifted copies of those same archives
are used for the later rows, not actual evolving distributions.

Even at the **upper** boundary of the first scalar bin, q/T_R=0.005,
the light-fermion partner at the 12 MeV snapshot would need q/T_R
at least about 4028. As both retained physical momentum limits
redshift, the maximum pair invariant decreases, so the first bin
never directly participates in this on-shell operator on q16 or
q32, nor even on hypothetical q64/q128 supports. Indirect transport
to this cold bin through the elastic channel is separate. The
argument does not bound unknown populations beyond the mesh or
exclude other interactions.

```sh
cd kinetics
for STEM in 12-vacuum-24 wide32-embedded-vacuum; do
 OPENBLAS_NUM_THREADS=1 python -W error retained_onshell_support_audit_reference.py results/"$STEM"-2026-11-26.json results/"$STEM"-shape-2026-11-26.json
done
```
