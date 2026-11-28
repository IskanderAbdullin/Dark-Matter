# Replace the assumed Pauli ceiling with a light-particle count bound

The [one-use cold-decay estimate](INDEPENDENT-WARM-LATE-COLD-SCALAR-SINGLE-USE-2026-11-26.md)
included a conditional factor for an assumed pointwise
light fermion occupation f≤.5. For open unpolarized decay
this assumption can be replaced by the complete retained
plus recorded light count, without specifying the shape
of the tracked exterior phase. The cold scalar Bose factor
is still held at the archived shape.

For the light species, the count measure in energy is
w(E)dE with w(E)=4 E sqrt(E²−m²)/(2π²), increasing from
E=m. Given 0≤f≤1 and the same weighted count, the filled
lowest-energy interval maximizes the unweighted integral
of f. To see this, compare to the filled step f_* ending
at E_F. The product (f−f_*)(w−w(E_F)) is nonnegative
on either side; equality of the weighted integrals implies
that the unweighted integral of f cannot exceed E_F−m.
Here p_F=(6π² N/4)^(1/3) and E_F=sqrt(p_F²+m²).
An independent numerical check uses exact cell volumes
for 120 arbitrary profiles and five saturated steps.

In each archive the maximum blocked energy length is
1.3264e−5 GeV. At p/T=.1 the minimum allowed open angular
energy width is 3.1153 GeV, so Pauli suppression of the
vacuum-width lower bound is no more than 4.26e−6. Multiplying
the vacuum angular bound by (1+max f_cold)/(1−4.26e−6)
gives the same upper decay probability to two significant
digits, about 5.94e−5. The one-use cold supply upper is
4.70e−14 (q32) and 1.42e−14 (q64) of the already
represented p/T<.1 scalar count. The improvement is the
removal of an assumed pointwise Pauli ceiling, not a
new estimate of the exterior spectrum.

This uses isotropic occupations and a full open decay
angular interval. An untracked population would require a
separate bound on its count; repeated re-entry, changed cold
Bose occupation and limits on actual cumulative production
remain unresolved. This argument does not bound a
closed-box normalization that drops some outgoing angles.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error warm_late_cold_scalar_pauli_count_bound_reference.py \
 results/warm-corrected-q32-time48-step48-2026-11-26.json.gz \
 results/warm-multishell-corrected-q64-fine-stage-48-2026-11-26.json.gz
~~~

[Count budgets, blocking bounds and saturation tests](results/warm-late-cold-scalar-pauli-count-bound-2026-11-26.json).
