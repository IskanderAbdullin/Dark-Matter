# Coupled shifted counts retain opposing near-balanced Bose/Pauli events

Extend the [single-event shifted limit](INDEPENDENT-TWO-MASS-SHIFTED-HIGH-BOSE-EXTENT-2026-11-27.md)
to three manufactured reversible events sharing four Fermi
and three Bose cells. Each event has independent, positive normalized
weights for the initial fermion, final fermion and final boson.
The capacities are deliberately unity. No physical rates are used.

Let the Fermi counts be f+εz and the Bose counts c/ε+εy.
For a vector of shifted extents x, write
z=z₀+V_F x and y=W_B x, where each Fermi event column sums
to zero and each Bose column to one. With u=f(1−f) and
θ_j=Σ_i (W_B)ᵢⱼ/c_i, the limit in time s=t/ε is

    dx_j/ds = u − (V_Fᵀ z)_j / θ_j.

Thus the genuinely coupled limit is a small affine matrix flow,
computed independently by a matrix exponential. Its entropy
increment, in units ε², is

    θᵀ x − (||z||²−||z₀||²)/(2u).

For finite ε, evaluate the rates from the shifted Bose and Fermi
exponents and integrate the coupled nonlinear flow. Compare the
integrated signed entropy production to direct 135-digit endpoint
entropies using exactly normalized weights. Exact normalization in
the high-precision reference matters: even roundoff in a column
sum would contaminate the nearly cancelling linear entropy term
before division by ε².

At ε=10⁻²⁴, all four absolute float64 Fermi counts round to the
same before and after values, while the three event extents remain
resolved and have opposing signs. The two numerical extent
vectors agree with their independent analytic limits to below
4e−15 in maximum absolute error. The finite-ε scaled entropy
integral matches the independent endpoint subtraction within
3.3e−14 relative across the tested starts and ε levels.
This neither demonstrates a full spectral calculation nor
fixes the physical occupancy scale.

Replay from the kinetics directory with NumPy, SciPy and mpmath:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error two_mass_coupled_shifted_bose_reference.py
```
[Full receipt](results/two-mass-coupled-shifted-bose-2026-11-27.json).
