# Restricted hard timelike thermal-quark diagnostic

This comparison isolates **statistics**, not the physical QCD
corrections. Choose three free massless quarks with
$N_c\sum Q_f^2=2$, $4\le x=\sqrt{s}/T\le12$,
$0.0125\le u=(\omega-\sqrt{s})/T\le32$.
Temperatures are 0.25–0.42 GeV, so the invariant masses in
this test are below 5.04 GeV and below the 10 GeV mediator.
No extrapolation through the pole, low-invariant-mass region,
or lower-temperature hadronic phase is made.

Use the free thermal vector trace in [Jackson and Laine, Eq. (2.4)](https://arxiv.org/html/1910.09567#S2.E4).
At fixed four-momentum the vacuum-normalized timelike trace is
$1-(2T/k)[\log(1+e^{-k_-/T})-\log(1+e^{-k_+/T})]$,
$k_\pm=(\omega\pm k)/2$. The KMS factor $n_B(\omega)$
times this trace equals the angular average of the two incoming
Fermi occupation factors. The independent identities check
agrees to 3.2e-15; the uncut MB boost moments
also agree with the analytic Bessel moments independently.
Compare to the same restricted charge-weighted vacuum mass
ratio with MB thermal weights, all other conventions unchanged.

| Radiation temperature slice (GeV) | FD Born optimistic scalar count expressed as Ωh² | Matched MB |
| --- | ---: | ---: |
| 0.25–0.340 | 0.0240771 | 0.0260588 |
| 0.25–0.34862 | 0.0269128 | 0.0291253 |
| 0.25–0.420 | 0.0552675 | 0.0597567 |
| 0.30–0.34862 | 0.0146635 | 0.0158627 |

Thus Boltzmann weights overestimate this particular free
count contribution by about 8% relative to FD. Quadrature
refinement moves the responses by at most 6.5e-8
fractionally. This is an optimistic count in which every direct
stable daughter and each prompt scalar is treated as if it could
end up as a final scalar. The heating proxy assumes new energy
stays in the dark sector and redshifts no faster than radiation.

**Do not add this test to the [measured sub-bands](INDEPENDENT-PARTIAL-HADRON-REHEAT-SOURCE-2026-11-26.md):**
the spectral support overlaps. It neither measures the QCD
in-medium spectrum nor supplies its error bars. It neglects the
interacting trace, nonzero quark masses, inverse evolution,
and dark/visible thermal feedback. This is not a new reheating
threshold or a prediction of the scalar abundance.

```sh
OPENBLAS_NUM_THREADS=1 python -W error born_thermal_quark_hard_portal_reference.py EOS.dat
```
