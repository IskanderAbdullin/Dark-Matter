# Instantaneous expansion contribution of the warm archive

Reconstruct the retained massive stable and scalar distributions
of the updated [warm reference](INDEPENDENT-WARM-COUPLED-B10-2026-11-26.md).
Use the original q≤16 trajectory, not the later q≤32 comparison.
Independently restart the 1.206 MeV checkpoint to exactly 1.000 MeV
at 12 and 20 steps. Integrate the represented energy and pressure
directly from the cell occupancies. For the exterior cohorts we
have their number and an energy interval, not an inferred shape.

Suppose a cohort with rest density $R_0=\sum_i m_i n_i$ and
kinetic-energy bounds $K_{0,-}, K_{0,+}$ is carried collisionlessly
from its ledger anchor to expansion factor $b\ge1$. For any mass,
$$
 \frac{K(p)}{b^2}\le K(p/b)\le\frac{K(p)}b.
$$
Thus its subsequent physical energy lies between
$R_0/b^3+K_{0,-}/b^5$ and $R_0/b^3+K_{0,+}/b^4$.
Its pressure is between one third and two thirds of its kinetic
energy, with the additional upper limit $\rho/3$. These inequalities
do not assume a thermal or monoenergetic exterior. A numerical
two-mass redshift control tests the interval against distributions
with both relativistic and nonrelativistic momenta.

Add the retained and exterior energy interval to the same
prescribed visible radiation density at fixed temperature. The
instantaneous fractional expansion change is
$\sqrt{1+\rho_{\rm dark}/\rho_{\rm vis}}-1$.

| Photon temperature | Fractional H increase, conditional on this archive | Dark pressure / dark energy |
| ---: | ---: | ---: |
| 12.00 MeV | about 5.047e-5 | 0.27283–0.27285 |
| 3.434 MeV | 7.5366e-5–7.5387e-5 | 0.17627–0.17645 |
| 1.206 MeV | 1.43590e-4–1.43635e-4 | 0.08108–0.08128 |
| 1.000 MeV, 20-step restart | 1.65901e-4–1.66116e-4 | 0.06706–0.06765 |
| 0.4138 MeV, independent later branch | 3.57597e-4–3.58283e-4 | 0.02456–0.02532 |

At 1 MeV, the 12- and 20-step endpoints change the expansion
interval by less than 3e-10 in absolute fractional H. The
exterior contribution to $\rho_{\rm dark}/\rho_{\rm vis}$ is
between 2.45e-6 and 2.88e-6 there. At the same state, omitting
those tails changes a rate per logarithmic expansion interval
through $1/H$ by at most about 1.44e-6 relatively. This is an
instantaneous clock check, not a bound on a re-evolved trajectory. A stable fixed-count proxy
far above the present reference density can therefore coexist
with a much smaller instantaneous expansion change in this
early, radiation-dominated archive. No observational limit
is applied here. The changing dark pressure-to-energy ratio
also prevents treating this entire component as extra radiation
with a single constant $\Delta N_{\rm eff}$.

As an explicit **collisionless extrapolation only**, redshift the retained
0.414 MeV distributions without further reactions to 100, 50 and
30 keV. This gives instantaneous fractional H shifts about
0.001548, 0.003210 and 0.005340, respectively. At 30 keV the
dark pressure-to-energy interval is about 0.000986–0.00118.
The extrapolation deliberately freezes the comoving particle
counts; it is not another coupled calculation nor a prediction
of primordial light-element abundances.

The comparison carries older exterior cohorts as collisionless
objects and does not feed their missing pressure/energy back
into the earlier collision solves. It is not a full evolution
with neutrino decoupling, nuclear abundances or an expansion
likelihood.

~~~sh
export OPENBLAS_NUM_THREADS=1
r=results/warm-coupled-pair
python -W error coupled_warm_pair_continuation_reference.py EOS.dat \
  "$r-b10-28-shape-2026-11-26.json" --stop 1.208342787706954 \
  --steps 12 20 --save-snapshots 1MeV-shapes.json >1MeV-rows.json
python -W error warm_archived_expansion_exterior_envelope_reference.py EOS.dat \
  results/50MeV-from-empty-pair-quantum-2026-11-26.json \
  results/50MeV-from-empty-pair-budget-2026-11-26.json \
  results/50MeV-from-empty-pair-quantum-shape-2026-11-26.json \
  --stage "$r-b35-32-2026-11-26.json" "$r-b35-32-shape-2026-11-26.json" \
  --stage "$r-b10-28-2026-11-26.json" "$r-b10-28-shape-2026-11-26.json" \
  --branch 1MeV-rows.json 1MeV-shapes.json \
  --branch "$r-b30-32-2026-11-26.json" "$r-b30-32-shape-2026-11-26.json" \
  --free-future .0001 .00005 .00003
~~~
