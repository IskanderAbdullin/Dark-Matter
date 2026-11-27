# Crossed charged-pair annihilation has a hard photon threshold at cold scalar q

Complement the vacuum photoproduction reference with
(P+ar b\leftrightarrow s+\gamma) at the same unequal masses.
Cross the daughter and photon in the photoproduction invariant
formula, retaining the sign from the incoming antifermion
spin sum. Independently trace the two diagrams with denominators
(-2P\cdot k) and (-2\bar b\cdot k), where (k) is the outgoing
photon. Both polarization sums are unaveraged and have no
coupling factors. The single charge assignment does not include
its charge conjugate.

On 60 finite-mass events the crossed invariant result is
positive and agrees with the direct trace within
(3.1\times10^{-13}) relative. Ward contraction and polarization
gauge shifts are consistent with zero change. Full FD/BE
forward and reverse weights agree at common temperature on
the same events, including in moving bath frames. This verifies
the equilibrium weights and crossing geometry for the channel,
not its share of a matched finite-temperature rate.

At fixed bath-frame scalar momentum (q), write
(E_s=\sqrt{m_s^2+q^2}) and let (c) be the cosine between scalar
and photon. Since the pair invariant is
(m_s^2+2k(E_s-qc)), the photon threshold is
[
 k_{\min}(c)=\frac{(M+m)^2-m_s^2}{2(E_s-qc)}.
]
At archived (q=4.1\times10^{-6}) the best-case thresholds are
1.60257 and 1.61143 in the fine and third bath energy units, at
(c=-1). At the same (q), the head-on to codirectional
increase is less than 1% in either archive. This channel cannot
contribute to the previously tested (10^{-5}) lab photon
soft window at that scalar momentum. It remains kinematically
accessible at hard photon energy; the threshold by itself
does not measure the incoming pair occupation or compare its
rate with real emission, absorption, or the Born source.

Replay from kinetics with the two pinned spectra:
```sh
OPENBLAS_NUM_THREADS=1 python -W error finite_mass_qed_crossed_annihilation_reference.py
```
[Receipt](results/finite-mass-qed-crossed-annihilation-2026-11-27.json).
