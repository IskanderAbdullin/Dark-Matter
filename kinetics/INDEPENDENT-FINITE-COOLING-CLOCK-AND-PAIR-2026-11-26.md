# Finite-cooling background clock and the post-1 MeV vacuum pair diagnostic

The [positive prompt distributions](INDEPENDENT-FINITE-COOLING-PHASE-SPACE-2026-11-26.md)
are inputs. All bounds below refer to these represented counts and
moments, not omitted interactions or a full reheating solution.

## Expansion at fixed visible temperature

If the represented distributions redshift without collisions from
the 1 MeV anchor, \(n=Y s_{\rm vis}\) and
\(p^2=p_{\rm anchor}^2(s_{\rm vis}/s_{\rm anchor})^{2/3}\).
Both \(E-m\) and \(p^2/(3E)\) are nonnegative, monotone concave functions of
\(p^2\). Bin chords and Jensen upper bounds therefore enclose the
homogeneous energy and pressure, including a conservative lower
value for any final overflow bin. Adding only the represented dark
energy to the same visible-clock comparator gives
\(H/H_{\rm vis}=\sqrt{1+\rho_{\rm dark}/\rho_{\rm vis}}\).
The upper-envelope mean-momentum representatives obey
\(d\rho/d\ln a=-3(\rho+P)\) in the independent derivative check.

| Clock quantity | Instantaneous | Finite, \(T_{\max}/T_R=5\) |
| --- | ---: | ---: |
| \(\Delta H/H_{\rm vis}\) at 1 MeV | \(2.9702\times10^{-7}\) | \(2.9388\times10^{-7}\) |
| \(\Delta H/H_{\rm vis}\) at 30 keV | \(1.12294\times10^{-5}\) | \(1.12293\times10^{-5}\) |
| Fractional shortening of 1 MeV → 30 keV proper time | \(7.48722\times10^{-6}\) | \(7.48712\times10^{-6}\) |

Visible-only time between those endpoints is 1465.80 s in the chosen
table. The correction compares the same post-anchor entropy clock,
with and without represented dark energy. It is dominated by the
common stable rest mass at lower temperature. It is not a nuclear
network, an observational constraint, or a derived \(\Delta N_{\rm eff}\).

## Pair conversion of the frozen stable measure

For symmetric Dirac particles and antiparticles, the vacuum forward
event rate is \(n_\chi^2\langle\sigma v_{\rm M}\rangle/4\), where
\(n_\chi\) counts **both** charges and spin states. Two scalars per
event give the undepleted first-order optical depth
\(\int (n_\chi\langle\sigma v_{\rm M}\rangle/2H)d\ln a\).
The collision calculation independently refines incoming angular,
outgoing scattering angle, clock and within-bin moments. Using the
bin RMS or two bin endpoints with the same moments changes the
result by ≤0.20%. The threshold \(p\)-wave coefficient is derived
from the expanded \(t/u\) spin amplitude; the full calculation uses
the exact amplitude. This diagnostic does not evolve or deplete the
two distributions.

From 1 MeV to 30 keV, the forward-only depth is \(4.66\times10^{-6}\)
for instantaneous reheating and \(2.58\times10^{-6}\) for the finite
cooling history. The corresponding first-order scalar rest-count
proxies are \(\Omega_s h^2=2.41\times10^{-12}\) and \(1.33\times10^{-12}\).
The former is comparable to its **very small** prompt-only count;
the latter is around 0.060% of its larger (still small) prompt-only
count. None of these values measures the reaction during
production or supplies a final scalar abundance.

At the 1 MeV anchor, a separate reconstruction uses both smoothing
widths and direct final-state lab-angle integration. It recovers
the vacuum event rate to better than 10⁻¹⁴. On these trial shapes
the forward Bose correction is below numerical resolution and the
inverse-to-forward ratio is below 10⁻¹⁸. Absolute inverse numbers
are tail-sensitive and are not a distribution-free limit. This
probe does not include the spectra produced later by the same
collision term.

## Conformal interval in the modeled history

Writing \(a_R\) for the radiation matching scale, the exact modeled
pre-radiation conformal increment is
\(2[1-\sqrt{a(T_{\max})/a_R}]/(a_R H_R)\).
The post-matching increment is integrated independently in log T
and log a. For \(k=10h/{\rm Mpc}\), \(k\Delta\tau\) between the
prescribed start and the 1 MeV snapshot is \(5.14\times10^{-4}\)
(instantaneous) and \(6.49\times10^{-4}\) (finite). These are
increments on the included intervals, not absolute conformal times
or a perturbative matching through the omitted heating branch.

## Reproduction

Use the full external EOS file and previously published artifacts.

~~~sh
OPENBLAS_NUM_THREADS=1 python -W error finite_cooling_collisionless_expansion_reference.py EOS.dat \
    results/finite-cooling-reheat-distribution-2026-11-26.json
OPENBLAS_NUM_THREADS=1 python -W error finite_cooling_post_anchor_pair_diagnostic_reference.py EOS.dat \
    results/finite-cooling-reheat-distribution-2026-11-26.json
OPENBLAS_NUM_THREADS=1 python -W error finite_cooling_signed_pair_anchor_reference.py EOS.dat psd \
    results/finite-cooling-class-psd-export-2026-11-26.json
OPENBLAS_NUM_THREADS=1 python -W error finite_cooling_horizon_clock_reference.py EOS.dat \
    results/finite-cooling-reheat-portal-2026-11-26.json
~~~
