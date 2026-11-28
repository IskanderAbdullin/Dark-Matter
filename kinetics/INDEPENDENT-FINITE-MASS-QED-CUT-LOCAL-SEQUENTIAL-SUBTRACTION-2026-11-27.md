# Resolve the two bare sequential residues before taking a local finite part

Use the [occupied spacelike poles](INDEPENDENT-FINITE-MASS-QED-NONNULL-PHOTON-WARD-POLES-2026-11-27.md).
An independent spinor resolution of the internal on-shell
projector factors each diagonal pole numerator into a scalar
Yukawa subprocess and a same-charge transverse-current subprocess.
For this pure scalar coupling the summed scalar subprocess has
a spin density proportional to the identity (anisotropy
below 3.56e−15). The product of the separately summed traces
divided by the two intermediate spins reproduces the modulus
squared of the full numerator within 1.71e−14 for all 48 roots.
This checks the residue, not a product of quantum kinetic rates
with widths and intermediate occupations already supplied.

Next select four actual spacelike slices with no occupation edge
inside the two-pole interval. Put $h=2q\cdot K$, $t=D_1/h$ and
$D_2/h=1-t$. Normalize the smooth numerator weights by the
mean of the two diagonal values at their roots (including
$dE_P/dt$). For a declared diagnostic imaginary regulator, remove
the two analytical $\pi f(0)/\epsilon$ diagonal terms in the
two denominators separately. Keep the interference coherent.
The finite diagonal boundary terms and interference have exact
polynomial moment integrals; the smooth numerator compression is
checked independently. A finite-width direct numerical integral,
with the same subtractions, agrees within 1.65e−11.

At fixed nonzero pole separation the equal-width ladder and
second-to-first width ratios .1, 2, 10 converge to the same
local limit. In the first sample its two diagonal pieces are
about −3.546 and −2.964, with interference +5.040, for a
combined local finite part −1.470 in the stated units. Shrinking
the integration window within the occupied cells moves the total
by about −1.008; the smooth edge integrals exactly reconcile
the change. Neither sign is a physical rate. Omitting the
interference or treating the local finite part as window-independent
would both fail this check.

The numerical widths do not preserve a dressed Ward identity
by construction. This is not a fully matched real-intermediate-state
subtraction prescription for the plasma. See
[the analysis of RIS schemes](https://arxiv.org/abs/2309.16615)
for the general distinction between removing resonant overlap
and defining an independent finite Boltzmann term. The limit
here assumes separated roots and smooth weights on the selected
slice; it does not cover merging poles or a pole at a histogram edge.

Replay from kinetics with archived spectra:
~~~sh
OPENBLAS_NUM_THREADS=1 python -W error finite_mass_qed_sequential_spin_factorization_reference.py
OPENBLAS_NUM_THREADS=1 python -W error finite_mass_qed_cut_local_sequential_subtraction_reference.py
~~~
[Spin factorization](results/finite-mass-qed-sequential-spin-factorization-2026-11-27.json),
[finite-part diagnostics](results/finite-mass-qed-cut-local-sequential-subtraction-2026-11-27.json).
