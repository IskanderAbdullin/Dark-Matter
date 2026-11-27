# Match the separate charged cuts to the familiar one-loop dispersive log

The [bare external self-energy](INDEPENDENT-FINITE-MASS-QED-BARE-EXTERNAL-SELF-ENERGY-2026-11-27.md)
isolates two signs of a charged bare loop. The
[earlier equilibrium control](INDEPENDENT-FINITE-MASS-QED-LOG-THERMODYNAMIC-2026-11-27.md)
uses the standard charged-pole dispersion without reconstructing
the two signs separately. This check bridges both calculations
for their **on-shell energy projection**, not the off-shell derivative.

Write E=√(p²+m²), e=√(k²+m²), and R=p k. For each branch
the angular inverse photon is J1±=∫du/K². Combining the
two log primitives at equal particle and antiparticle occupations,
the ratio inside their summed logarithm is

    [(E e)²−(m²+R)²] / [(E e)²−(m²−R)²] = (p−k)²/(p+k)².

Thus J1+ + J1− = log(|p−k|/(p+k))/R. With each charged
branch using the r/2 d(e) measure, the combined projected
energy shift is proportional to

    (2/E) ∫ dk (k²/e) f(k)
      − (m²/(p E)) ∫ dk (k/e) f(k) log[(p+k)/|p−k|].

The p=0 check evaluates the two original branches using
rationalized denominators and compares with the integrated
logarithmic limit. The null-photon term enters with the
thermal sign of the [photon self-energy convention](https://arxiv.org/abs/1110.2171),
so it adds 2(π²/6)/E before g²/(4π²) in equilibrium
temperature units. With Fermi–Dirac occupations, 25 mass and
momentum fixtures agree with the independent
[Grohs–Fuller equilibrium pole result](https://arxiv.org/abs/1706.03391)
to 5.7e−11 absolute in those stripped units. The two angular
branches are independently refined before merging.

For the frozen archived occupations treat the charged particle
and antiparticle as copies of the same tabulated function, and
do not assume a photon bath. Direct separate-cut integration
and independent shell logarithmic-weight quadrature agree
to 3.1e−11 absolute. At the two central archived parent
events, the log lowers the **charged-only factorizable shift**
by about 0.089%; the light charged daughter reduction is about
0.00067%. At zero parent momentum the reduction is about 0.58%.
These numbers compare frozen contributions, not the correction to
a decay rate.

A derivative of the complete dressed pole differs from
a fixed-spin off-shell derivative of the on-shell projection.
This identity provides a check on the dispersion input, not a
physical residue, screened soft pole, rate matching, or a new
claim that the archived state was in equilibrium.

Replay from kinetics:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_dispersive_cut_bridge_reference.py
```
[Raw numerical receipt](results/finite-mass-qed-dispersive-cut-bridge-2026-11-27.json).
