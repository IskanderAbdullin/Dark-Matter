# An apparent kinetic singularity from sharp occupation bins

The [previous archived aperture](INDEPENDENT-FINITE-MASS-QED-ARCHIVED-APERTURE-2026-11-27.md)
uses flat shell occupations in the real energy shift. A continuous
one-loop energy function need not have a bounded momentum derivative
with that choice. This matters if using group velocities in a
phase-space Jacobian, which the aperture check deliberately omitted.

Write A(r)=(r/E_r)f(r) and K(p,r)=log[(p+r)/|p−r|].
Away from an occupation edge, differentiate the log integral by
parts, using ∂pK=2/(p+r)−∂rK. The result is the regular
2A/(p+r) term plus the integral of K against the distributional
derivative of A, with the two hard endpoint jumps included. For a
jump ΔA at r₀, the derivative of the energy includes

    −g² m² ΔA / (4π² p E_p) × log[(p+r₀)/|p−r₀|].

An upward jump therefore produces a negative logarithmic velocity
singularity even when every sampled velocity is positive. For the
two upward archived parent edges near the aperture, the coefficient
at g=0.01 is −1.18e−10 and −1.96e−10. Extrapolating the local
finite remainder to the sign reversal would require a distance
from the edge of order 10 raised to minus 3.68 billion and minus
2.22 billion, respectively. This is a mathematical artifact of
the step reconstruction, not a resolved branch, physical pole,
or a justification for exact continuum monotonicity of that
step model. The tracked daughter jumps are downward.

To remove this numerical artifact, join cells in u=p³/3 with
C1 smoothstep transitions to common edge values and adjust the
plateau so the prescribed shell average stays exact. The transition
supports occupy at most the chosen fraction of each side of a
cell (1.5%, 5%, 12% here); the plateau correction preserves the
Pauli bounds. There is still a cutoff at the artificial hard
domain edge, outside the tested aperture. The five represented
edge tests have finite smooth slopes greater than .99 with
independent finite differences. Changing only the dispersion
input from flat to smooth shifts the full aperture center by at
most 3.29e−7 at g=0.005. For both archived shapes and two cuts,
the relative shift caused by the nonfactorizable log itself
changes by ≤0.014% of that shift. The overall smooth occupation
integrals of the compared aperture remain unweighted and are
not matched physical rates. The previously reported numerical
apertures are stable under this test without making a global
monotonicity statement for discontinuous cells.

Replay from kinetics:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_shell_edge_slope_reference.py
```
[Raw receipt](results/finite-mass-qed-shell-edge-slope-2026-11-27.json).
