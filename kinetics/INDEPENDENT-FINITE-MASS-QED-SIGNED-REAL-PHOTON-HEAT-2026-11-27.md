# In the signed photon bridge, the bath loses heat on both archives

Keep the exact moving apertures and the fixed backgrounds of the
[signed photon number bridge](INDEPENDENT-FINITE-MASS-QED-SIGNED-SUPPORT-BRIDGE-2026-11-27.md).
Here calculate the actual **real photon energy moment** in the finite
lab energy band [1e−5, 2e−4]. Identical scalar gain does not imply
identical photon heating in the two crossed directions.

After absorbing the scalar Bose weights into $F_pm,R_pm$, define
$N_pm=F_pm-R_pm$ and $P_pm=F_pm+R_pm$. The real forward emission
rate is $(1+n)F_+-nR_+$ and the forward absorption rate is
$nF_--(1+n)R_-$. Their scalar sources add, but their real photon
energies have opposite signs. Multiplying by the photon energy
cancels the common $1/k$ in the real photon number measure. The
real photon heat density per unit photon energy is proportional to

$$
(1/2+n)(N_+-N_-)+(P_++P_-)/2 .
$$

For smooth signed weights its limit is $2T N'(0)+P_0$. This is
not the limit for the scalar source, and it is not assumed at
archive steps. At common equilibrium $N_pm=mp	anh(k/(2T))P_pm$
cancels the photon heat separately. The manufactured thermal
test gives less than 4.1e−13 absolute residual. Independent
direct real gain/absorption histories reproduce the grouped
heat density to within 9e−17 absolute of the unscaled measure.

| Archive | Photon heat / positive Born scalar energy source | Leading-subtracted scalar number / Born |
|---|---:|---:|
| fine | −5.52538e−5 | −2.13606e−6 |
| third | −1.19093e−4 | +6.94187e−7 |

Both ratios use $g=.005$ at one cold scalar momentum. All five
finite subwindows cool the photon bath in both archive controls,
even though the scalar comparison switches sign. Two-mesh
changes of the photon ratio by subwindow are below 1.2e−14
on this fixture. This heat is the finite **real** photon
moment, not a change in the total bath energy or a matched
finite virtual or dressed plasma rate. The adjacent
scalar number column uses only its common leading subtraction.

Replay from kinetics with pinned archived spectra:
~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_signed_photon_heat_reference.py
~~~
[Fine receipt](results/finite-mass-qed-signed-real-photon-heat-fine-2026-11-27.json),
[third receipt](results/finite-mass-qed-signed-real-photon-heat-third-2026-11-27.json).
