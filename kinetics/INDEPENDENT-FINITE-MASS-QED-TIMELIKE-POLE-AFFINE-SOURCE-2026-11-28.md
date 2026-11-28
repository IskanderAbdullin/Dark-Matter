# Fixed-angle timelike pole absorption as a function of an unknown population

Use the two [frozen timelike modes](INDEPENDENT-FINITE-MASS-QED-POLE-PAIR-APERTURE-2026-11-28.md) at selected momenta above the [past-null support opening](INDEPENDENT-FINITE-MASS-QED-FROZEN-TRANSVERSE-PLASMON-GATE-2026-11-27.md). The plasma pole occupation is an independent input **n**, not a fitted Bose law. This test evaluates the bare on-shell charged matrix and the archived Pauli factors at fixed scalar momentum q=4.1e−6 and fixed cosine zero, not an integrated physical rate or a medium-matched correction.

Write T=q−K, with T⁰<0, T²<0, and V=|T spatial|. At each parent energy E the on-shell daughter shell fixes the parent momentum along T to [E T⁰−(M²+T²−m²)/2]/V. Solve the quadratic for the positive aperture edge and integrate over E and the residual azimuth, splitting both parent and daughter occupancy edges. The two-particle phase measure after angular reduction is proportional to dE dφ/V. The reported unnormalized coefficients include the positive-frequency retarded residue in ω and the corresponding squared bare vertex, with a common phase-space factor, coupling constants and the remaining photon measure omitted. The transverse part sums both spatial polarizations; the longitudinal part contracts the conserved charge current with the Coulomb 00 propagator, not with the covariant gauge-fixing term.

Let A and B be the resulting nonnegative coefficients weighted by f_P(1−f_b) and (1−f_P)f_b. For scalar occupation n_s, after division by n_s the net of absorption and its inverse is

```text
n [(1+1/n_s) A − B] − B.
```

The code evaluates the signed difference separately rather than subtracting two rounded numbers. The table gives the **conditional balancing value** n⋆=B/[(1+1/n_s) A−B] at each checked aperture.

| Polarization | k | fine n⋆ | third n⋆ |
| --- | ---: | ---: | ---: |
| transverse | .0015 | 2.83126e−9 | 2.26304e−9 |
| transverse | .00256 | 4.94870e−8 | 4.47793e−8 |
| longitudinal | .0012 | 1.42158e−8 | 1.38667e−8 |
| longitudinal | .00256 | 5.65411e−7 | 5.40705e−7 |

If the external bath factors and geometry are held fixed, this particular coefficient is negative below its n⋆ and positive above it. The small values do not imply that the mode is populated: at n=0 there is only inverse emission in this diagnostic. The separate transverse samples at .0015 have parent threshold energy about 18 in archive units, strongly suppressing both forward and reverse coefficients. Do not compare the small n⋆ there as an absolute event flux. In a common-temperature Fermi/Bose fixture with shared charged chemical potential, the independent directional ratio gives the ordinary thermal n⋆ to 9.19e−13, as required by detailed balance.

The bare denominators have [analytic sign protection](INDEPENDENT-FINITE-MASS-QED-TIMELIKE-NO-SEQUENTIAL-ROOTS-2026-11-28.md); a rationalized evaluation avoids collinear cancellation. The spinor Ward identities are independently checked for each diagram and for their sum. The largest squared single-diagram residual relative to the Ward-reduced Born matrix is 2.08e−25; the coherent weighted residual relative to the polarized norm is 1.60e−17. The longitudinal time component and the equivalent spatial current projection give integrated norms agreeing to 9.32e−10. Independent finite differences check the analytic pole residue to 1.11e−9 relative. Joint energy/azimuth refinement changes the raw coefficients by at most 2.54e−7 relatively.

Reproduce the [receipt](results/finite-mass-qed-timelike-pole-affine-source-2026-11-28.json) from `kinetics/`:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p timelike-work
for g in fine third; do
  python -W error finite_mass_qed_timelike_pole_affine_source_reference.py "$g" \
    > "timelike-work/timelike-pole-$g-mesh0.jsonl"
  python -W error finite_mass_qed_timelike_pole_affine_source_reference.py "$g" --refine \
    > "timelike-work/timelike-pole-$g-mesh1.jsonl"
done
python -W error finite_mass_qed_timelike_pole_affine_source_receipt.py timelike-work
```

No external charged dispersion or vertex resummation, spacelike Landau cut matching, photon kinetic population, virtual contribution, all-angle integral, or total scalar rate is fixed by these conditional tests.
