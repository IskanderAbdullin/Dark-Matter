# Crossed scalar absorption is a separate, infrared-finite channel

In addition to photon-assisted decay and charged-pair annihilation,
cross the photon and the scalar in the vacuum amplitude to obtain
$P+s\leftrightarrow b+\gamma$. This is **not** the previously checked
$P+\gamma\leftrightarrow b+s$ photoproduction process. Use the same
frozen massive parent, daughter and scalar spectra, with a Bose
photon bath at temperature 1 in archive units. No antiparticle
proxy is required for this charge assignment. The earlier
[restricted real-photon comparator](INDEPENDENT-EARLY-POLE-YUKAWA-REAL-RADIATIVE-WITNESS-2026-11-27.md)
already includes the same reaction under the daughter-plus-photon
label. Here its integration support and normalization are checked
independently and extended, not counted as a new reaction.

At given incoming parent momentum, set $B=P+q$, $W^2=B^2$.
The photon energy allowed in the bath frame starts at
$(W^2-m_b^2)/(2(B^0+|\mathbf B|))$. Since
$W\ge M_P+m_s$ and $B^0-|\mathbf B|\ge E_s-|\mathbf q|$,
it is always at least
$(1-m_b^2/(M_P+m_s)^2)(E_s-|\mathbf q|)/2$.
At scalar momentum 4.1e−6 the respective lower bounds are
4.32474e−4 and 4.30074e−4 on the two archives. This excludes
the channel from a lab photon window below 1e−5 without a
regulator subtraction.

Crossing the photon reverses both charged-photon invariants,
but no charged spin crossing sign is needed. Verify the analytic
azimuthal average against direct angle integration and polarized
spin traces (angle residual below 5.7e−12, direct spin residual
below 2.9e−9). The Ward check is below 4.1e−15. The radial
quadrature splits transformed daughter occupation edges and
uses the positive photon-energy variable to resolve Bose
enhancement. In the fixed scalar clock it uses the same
unaveraged spin and phase conventions as the other single-charge
real comparators. At common temperature, the charged plus photon
statistical weights obey gain/loss ratio $\exp(-E_s/T)$, giving
zero scalar net to better than 2e−16 of forward plus reverse
turnover after applying the Bose scalar factors.

For photon temperature 1 and $g=0.005$ the isolated contribution
relative to the positive Born source is −4.65037e−7 (fine)
and −5.45103e−7 (third). The last two meshes differ by under
1.8e−8 relative in gain, loss and net. The frozen parent
occupation has no support beyond the integration cap 56;
the positive gain tail is bounded by 5.74e−10 before coupling
factors using the daughter support and the exponential photon
tail. With the previous photon and charged acceptance windows restored,
independently computed gain and loss match that comparator
after its factor of two for conjugate charges is removed
(maximum relative difference 1.64e−5). The earlier window
contains just 39.54–39.76% of the full signed net at this mode.
Use the full value in place of the restricted one, not their sum.
This is a fixed-mode collision comparison, not a
dressed thermal width, a matched radiative correction, a
self-consistent cosmological evolution, or a second charge
assignment.

Replay from kinetics with the two archived spectra:
~~~sh
OPENBLAS_NUM_THREADS=1 python -W error finite_mass_qed_scalar_absorption_reference.py
~~~
[Receipt](results/finite-mass-qed-scalar-absorption-2026-11-27.json).
