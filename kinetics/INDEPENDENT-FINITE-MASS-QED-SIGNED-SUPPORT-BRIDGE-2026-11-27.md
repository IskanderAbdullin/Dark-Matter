# Pair both real photon signs and their leading subtraction before integration

The [original signed soft calculation](INDEPENDENT-FINITE-MASS-QED-PARENT-INTEGRATED-SIGNED-SOFT-2026-11-27.md)
examined lab photon energies up to 1e−5. Extend the exact
moving-aperture comparison through 2e−4 at the same cold scalar
momentum and on the same two archives. Include *all* translated
archive edge crossings inside the wider parent aperture.
At the 1e−5 interface the old limited-support evaluator
and the broadened one reproduce the signed kernels to
1.81e−18 absolute. The upper bridge bound is a checked range
of the numerical chart, not a screening scale.

The charged aperture equation is monotone here: above parent
energy .22 the phase derivative has a lower bound exceeding
.00043 from the scalar energy term and the maxima of the
longitudinal photon and scalar-momentum derivatives.
Check the outer sign brackets at .22 and 6, and constrain
Newton iterates inside them. Both signs retain their
moving boundaries and the exact FD factors. The control
would not apply unchanged once the absorption derivative
can change sign.

For each photon energy $k$, write $A^+$ and $A^-$ for the
respective exact emission and absorption weights before
putting in photon Bose factors. Keep separate $G$ and $L$,
and let $A^0$ be the independent zero-photon eikonal
baseline. The local leading-subtracted numerator over
the logarithmic photon measure is
$(G^+-G^0,L^--L^0)+n_\gamma(A^++A^--2A^0)$.
Thus the large leading terms are removed before integrating
the photon range. This is only their universal leading
subtraction, not the finite thermal virtual contribution.

The local signed net per unit photon energy, normalized
by the Born source times its eikonal angular factor
and *before* the common coupling factor, is:

| Photon energy | Fine archive | Third archive |
|---:|---:|---:|
| 1e−5 | −3880.59 | +1318.21 |
| 4e−5 | −4635.81 | +1500.21 |
| 2e−4 | −5118.27 | +1684.60 |

The opposite signs occur after the leading subtraction;
they are comparisons between archive occupations, not
between independent meshes of the same archive. The
largest absolute drift in either unnormalized $G,L$
residual component across the two tested meshes is 4.94e−6
over the six sampled photon energies. No finite virtual or
thermal external-leg term, self-consistent width, or
physical inclusive sign is inferred.

Integrating only this finite bridge on the fine archive gives
−2.13606e−6 of the Born net at g=.005. Each of the five logarithmic
subwindows contributes negatively. The largest two-mesh drift per
subwindow, before applying the gauge factor, is 1.21e−6 relative to
the Born reference (about 7.7e−13 after that factor). The analytic
zero-photon eikonal baseline agrees with its numerical replay within
2.42e−17 absolute. The integrated result is for this finite band
and leading subtraction only; it does not include what lies above
2e−4 or the finite virtual parts.

Replay with the pinned archive spectra in kinetics:
~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_signed_support_bridge_reference.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_integrated_signed_bridge_reference.py fine
~~~
[Integrated fine receipt](results/finite-mass-qed-integrated-signed-bridge-fine-2026-11-27.json).
[Pointwise receipt](results/finite-mass-qed-signed-support-bridge-2026-11-27.json).
