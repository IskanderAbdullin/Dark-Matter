# Broadening only the shell constraint does not preserve detailed balance

The [common-aperture statistics gate](INDEPENDENT-FINITE-MASS-QED-APERTURE-STATISTICS-2026-11-27.md)
tests exact Born kinematics. To test a **possible implementation error**
at finite width, fix the scalar energy $E_s$ and give the on-shell
mismatch $\delta=E_{P,\mathrm{pole}}-E_{b,\mathrm{pole}}-E_s$
a manufactured Lorentzian weight on a finite positive-energy strip.
Use equal bath temperature, vanishing chemical potentials, and
full FD and BE factors.

If the occupations still use the two pole energies, the ratio of
forward to reverse statistical weights is $e^{-\delta/T}$.
Equivalently the net divided by forward-plus-reverse turnover is
$-\tanh(\delta/(2T))$ at each mismatch. Sharing the same
positive Lorentzian between forward and reverse does not fix it.
If instead the FD factors use conserving frequencies
$E_P'=E_{P,\mathrm{pole}}-\lambda\delta$ and
$E_b'=E_{b,\mathrm{pole}}+(1-\lambda)\delta$, balance is exact for
any $\lambda$. The specific frequency routing in this control
does not derive the interacting spectral functions.

The pole-frozen error integrated over the strip is generally
linear in the width parameter $\gamma$, not quadratic. Write
$F$ for the smooth net including the phase weight, so $F(0)=0$,
and $H_0$ for the zero-mismatch turnover. For endpoints
$[-a,b]$ with $a,b>0$, the limiting ratio of integrated
net/turnover divided by $\gamma$ is

$$
 \frac{F'(0)\log(b/a)
 +\int_{-a}^{b} [F(\delta)-F'(0)\delta]\,d\delta/\delta^2}
 {\pi H_0}.
$$

A symmetric window cancels the term odd in mismatch but not the
even term from the changing FD weights and phase factor.
At (T=0.1) in the fine-mass geometry, the width-linear
coefficients are (+0.31866) (symmetric uniform phase),
(-0.79790) (same symmetric strip with sloped phase)
and (-1.97449) (asymmetric strip, uniform phase). The sign
change does not indicate any physical heating. It is a diagnostic
of the broadened-shell shortcut under manufactured weights.

Twelve controls use the archived scalar and charged energy
geometry but equilibrium test occupations. Independent
Lorentzian angle and subtracted real-energy quadratures agree
below 1.4 × 10⁻¹² in net/turnover across the width ladder.
All conserving routes are consistent with pointwise zero net.
This is compatible with the equilibrium thermal frequency
relation emphasized in the [Kadanoff–Baym framework](https://pure.mpg.de/rest/items/item_1798467_1/component/file_1798466/content).
It determines neither the plasma spectral shape nor the
normalized radiative scalar rate.

Replay from kinetics with the two pinned spectra:
```sh
OPENBLAS_NUM_THREADS=1 python -W error finite_mass_qed_dressed_frequency_balance_reference.py
```
[Receipt](results/finite-mass-qed-dressed-frequency-balance-2026-11-27.json).
