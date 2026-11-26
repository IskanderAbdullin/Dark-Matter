# Passive pair emission during finite cooling, with its scalar recoil

The [prompt history](INDEPENDENT-FINITE-COOLING-PHASE-SPACE-2026-11-26.md)
is the unmodified source. This check accumulates its stable particles
before each temperature-panel midpoint and then uses that instantaneous
positive spectrum in a **first-order, spontaneous forward** pair tally.
No pairs deplete or heat the input distribution; no inverse, stimulation,
additional charged bath, transfer through heating, or dark rescattering
is included. This is not a lower or upper bound on the complete coupling.

For equal particle and antiparticle densities, the event count is
\(n_\chi^2\langle\sigma v_\mathrm{M}\rangle/4\); the two emitted particles
give \(2R\). To compare their counts at the 1 MeV anchor, integrate
\(2R(a)(a/a_f)^3dt/s_f\). In the cooling branch the clock and volume
are \(H\propto T^4\), \(a\propto T^{-8/3}\); in the subsequent branch
the supplied visible entropy is conserved. The independent count
calibration against the prompt-only result is checked in both branches.

The outgoing direction is integrated with the exact positive spin
amplitude in the center-of-mass frame, including its azimuth around
the incident axis. At fixed scattering cosine \(x\), the projection
on the boost direction has mean \(c_0 x\) and variance
\((1-c_0^2)(1-x^2)/2\), where \(c_0\) is the cosine between incident
and boost axes. Writing \(C+D\cos\phi\) for the corresponding lab
energy, the kernel integrates both the radial square and number
over the exact clipped azimuth. The analytic primitive and the
independent rotated lab-angle test agree below \(7\times10^{-16}\)
in their documented absolute/relative residuals. The two temperature,
source grid, and angular controls are retained separately.

| Prescribed history | New pre / post scalar Ωh² | New pre / post scalar RMS \(p/T\) at 1 MeV | Prompt-only → mixed scalar RMS |
| --- | --- | --- | --- |
| Instantaneous | 0 / \(3.02\times10^{-10}\) | not applicable / 7.38 | 6.70 → 7.37 |
| Finite \(T_{\max}/T_R=5\) | \(3.57\times10^{-10}\) / \(1.53\times10^{-10}\) | 2.87 / 7.12 | 0.604 → 2.047 |

For the finite example the recoil adds about 22.8% to the tiny
prompt scalar count but carries about 92.9% of the combined
squared-momentum sum. Consequently the earlier
[prompt-only tracer](INDEPENDENT-FINITE-COOLING-SCALAR-TRACER-2026-11-26.md)
should not be interpreted as the scalar spectrum with pair production.
The first-order maximum stable loss rate is below 0.0016 per Hubble
time on the sampled trajectory, and the new scalar rest count is
tiny relative to the nominal stable Ωh²=0.120. Neither check proves
that a full reactive spectrum has this shape.

## Resolved positive recoil and the conditional mixed shape

The clipped azimuth version retains positive counts and the squared
moment in each recoil bin. The tested CDF probes move by less than
4.0×10⁻⁴ over the two complete numerical controls. For the finite
example, the recoil-only fraction below p/T=1 is 0.1204 and the
prompt-only fraction is 0.9198. The number-weighted mixture gives
0.7713. The new particles account for 18.58% of the two-component
number (different from the 22.8% addition relative to prompt alone).

In the same prescribed background, the weighted chord/Jensen
brackets on the combined collisionless mean distance from 1 MeV
to equality are 20.496–20.646 Mpc (instantaneous) and 4.308–4.337
Mpc (finite). The finite prompt-only interval was 2.415–2.430 Mpc.
Total scalar rest Ωh² is still just 3.05×10⁻¹⁰ or 2.742×10⁻⁹
in this first-order comparison. These are geodesic shape
statistics, not a total-matter perturbation result. The optional
positive exporter reconstructs each component independently and
preserves the total count and second moment; a transfer test must
check the resulting synthetic shape separately. For the finite mixed
shape the external dilation/energy/pressure screen selects 320 and
224 transformed nodes, with worst residual 9.94×10⁻⁹ and 5.59×10⁻⁸
under the two smoothing choices. This screen configures the solver;
it does not predict a mixed transfer result.

## Reproduction

~~~sh
OPENBLAS_NUM_THREADS=1 python -W error finite_cooling_during_production_pair_reference.py \
    EOS.dat results/finite-cooling-reheat-portal-2026-11-26.json
OPENBLAS_NUM_THREADS=1 python -W error finite_cooling_during_production_pair_reference.py \
    EOS.dat results/finite-cooling-reheat-portal-2026-11-26.json --shape
OPENBLAS_NUM_THREADS=1 python -W error finite_cooling_pair_mixture_audit_reference.py \
    EOS.dat results/finite-cooling-reheat-distribution-2026-11-26.json \
    results/finite-cooling-during-production-pair-shape-2026-11-26.json
OPENBLAS_NUM_THREADS=1 python -W error export_finite_cooling_pair_mixture_reference.py \
    EOS.dat results/finite-cooling-reheat-distribution-2026-11-26.json \
    results/finite-cooling-during-production-pair-shape-2026-11-26.json psd_mix
OPENBLAS_NUM_THREADS=1 python -W error prepare_finite_cooling_scalar_tracer_reference.py \
    mix_native --psd psd_mix --manifest results/finite-cooling-pair-mixture-export-2026-11-26.json
~~~
