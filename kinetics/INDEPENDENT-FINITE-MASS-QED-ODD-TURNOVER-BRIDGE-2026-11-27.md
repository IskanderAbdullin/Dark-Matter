# The odd part of the locally subtracted real source is turnover, not net

Use the [signed photon bridge](INDEPENDENT-FINITE-MASS-QED-SIGNED-SUPPORT-BRIDGE-2026-11-27.md)
with exactly common photon energies and one-charge normalization.
At a fixed scalar occupation, fold its Bose factor into the
non-photon forward and reverse weights $F_\pm,R_\pm$. Define
$N=F-R$, $P=F+R$. For the two signs, the photon factors give
$(1+n)F_+-nR_+ + nF_--(1+n)R_-$. Subtract $(1+2n)N_0$
with the same Born eikonal baseline. The result per unit photon
energy is

$$
\frac{(1/2+n)(N_++N_--2N_0)+(P_+-P_-)/2}{k}.
$$

For smooth signed weights its limit is $T N''(0)+P'(0)$;
this Taylor result is **not** imposed across histogram edges.
The equilibrium identity $N(k)=-\tanh(k/(2T))P(k)$ makes both
terms cancel exactly for a common bath. Manufactured thermal
controls give a maximum absolute residual of 2.13e−14;
a separate smooth nonthermal fixture tends to −0.10.

All 12 actual signed bridge controls reproduce the separate
forward-and-reverse calculation within 7.2e−10 relative when
using the same archived Born baseline. Replacing the odd
turnover by odd net changes the sampled results by up to
13.03% relative. This provides a sensitive regression check.

The finite integral over the same [1e−5, 2e−4] bridge is
+6.94187e−7 of Born in the third archive, compared with
−2.13606e−6 in the fine archive at g=.005. Each of the five
third-archive subwindows is positive. Its largest two-mesh
drift per subwindow is 3.28e−8 of Born before the common
gauge factor. Only the leading soft common term is subtracted;
these are not physical inclusive signs or a finite virtual
rate.

Replay (from kinetics):
~~~sh
OPENBLAS_NUM_THREADS=1 python -W error finite_mass_qed_odd_turnover_reference.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_integrated_signed_bridge_reference.py third
~~~
[Turnover audit](results/finite-mass-qed-odd-turnover-2026-11-27.json).
[Integrated third receipt](results/finite-mass-qed-integrated-signed-bridge-third-2026-11-27.json).
