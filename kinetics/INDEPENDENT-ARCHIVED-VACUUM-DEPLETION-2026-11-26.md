# Late vacuum scalar-pair annihilation cannot close the archived count budget by itself

This is a conditional comparison after the 50 MeV reheating benchmark has
reached either its 12 or 9.3 MeV archived endpoint. Keep the measured
stable fermion-plus-antifermion yield at that endpoint. Afterward, assume
charge symmetry, adiabatic visible entropy, standard radiation as a floor
to the Hubble rate, and only vacuum χ₂χ̄₂ → φφ depletion. Drop **all**
further production and inverse reactions to give this vacuum channel
the most optimistic opportunity to remove the surplus.

Use the [upstream arbitrary-momentum envelope](https://github.com/IskanderAbdullin/Dark-Matter/blob/90c141722376dab64e4cd10422d545474817787c/kinetics/PHI-BOSE-PORTAL-VACUUM-DEPLETION-ENVELOPE-2026-11-24.md)
for the unpolarized scalar Yukawa benchmark. Both exchanged momentum
squares obey (t,u\leq0), so each spin-averaged diagram is bounded by
(y^4 s/[2(M^2-t)]) or its u-channel counterpart. Bounding their sum
by Cauchy and integrating the identical-particle angular measure gives
(K=y^4/(8\pi M^2)=3.97887\times10^{-12}\,\mathrm{GeV}^{-2}).
The bath-frame Møller factor relative to the center-of-mass rate is
(s/(4EE')\leq1), so this upper bound needs no thermal assumption.

The total stable yield (Y=(n_++n_-)/s_{\rm vis}) obeys
(dY/dt\geq -K s_{\rm vis}Y^2/2) because each event removes two
particles and (n_+n_-=(n_++n_-)^2/4). Thus its surviving fraction is
at least (1/(1+K Y_{\rm start}\int s_{\rm vis}dt/2)). Use the
tabulated radiation for one integration and independently bound
the integral with (H\geq H_\gamma=A_\gamma T^2). In the latter case,
entropy conservation yields an integral proportional to
(int_0^{T_f}(3g_s+Tg_s')dT=2\int_0^{T_f}g_s\,dT+T_fg_s(T_f)).
Bound this by (3T_fg_{s,\max}); no late-time matter is needed for
this conservative photon-only bound.

| Archive endpoint | Fixed-count Ωh² | Lower Ωh², photon bound | Lower Ωh², tabulated radiation |
| --- | ---: | ---: | ---: |
| 12 MeV | 88.18 | 8.497 | 17.54 |
| 9.3 MeV | 88.44 | 10.70 | 21.44 |

Even the loose photon-only lower bound exceeds 0.120 by over a factor of
70. This is **not** a lower bound if stimulated emission, other
annihilation/number-changing channels, extra entropy, or changes to the
reheating benchmark are allowed. If both accessible final scalar
occupations are uniformly bounded by F throughout this interval, the
maximum enhancement is at most (1+F)². Reaching 0.120 requires at
least F≈7.85/9.07 (optimum using the photon bound), or F≈12.50/14.35
(on the radiation clock). These are necessary thresholds for such a
uniform cap, not measured occupations or sufficient conditions.

At the two archived endpoints, **within the retained incoming χ₂
support only**, each final scalar has (E_\phi\geq
s/[2(E_{\rm tot}+|\mathbf P|)]\geq m^2/E_{\rm tot,max}).
This gives approximately 1.14/1.46 MeV, independent of angle and
well above the very cold scalar bins. At least one of the two scalar
daughters has lab energy ≥20 MeV by energy conservation. This does not
limit stimulation from exterior incoming momenta or from a later
state with changed supports. The next relevant calculation is the
signed stimulated rate on actual distributions, with exterior control.

## Reproduce

Prepare the [positive-temperature EOS and zero endpoint](INDEPENDENT-BENCHMARK-SETUP.md).
From `kinetics`, run:

```sh
OPENBLAS_NUM_THREADS=1 python -W error archived_vacuum_depletion_gate_reference.py EOS.late.dat results/conditional-relic-normalization-2026-11-26.json results/12-quantum-48-shape-2026-11-26.json results/early-9p3-quantum-shape-2026-11-26.json
```
