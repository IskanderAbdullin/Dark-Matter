# Which direct terms still prevent a light-occupation proof?

The [heavy-mode barrier](INDEPENDENT-WARM-HEAVY-OCCUPATION-BARRIER-2026-11-26.md)
does not establish the light fermion cap. As a separate
diagnostic, bound the cumulative positive creation into a
tagged light mode from three direct pair channels after
the warm 10.442 MeV anchor. The bound discards all losses.
It does not include the heavy decay/inverse redistribution
or light elastic gain, so it cannot prove a light cap.

- The scalar inverse vertex uses the same all-momentum
  angular delta and energy split as the heavy barrier,
  with the light threshold and heavy exchanged fermion.
  With one tagged outgoing spin, the sum over the other
  spins gives at most twice the unpolarized initial spin
  average, and the incoming scalar symmetry only reduces
  the estimate. Finite scalar energy alone bounds this
  part; no pointwise scalar cap is imposed.
- For heavy pairs producing light pairs, use the crossed
  invariant bound and a Fermi-ball radial moment with
  heavy per-spin cap 1/2. Scale the available heavy
  energy from zero up to the full dark-energy upper.
  All included charge channels are overcounted, as in
  the heavy creation check.
- Vacuum charged-lepton and selected charged-pion
  pairs use Boltzmann energy factors and the subpole
  upper for a tagged light mode. For a lepton invariant
  above the split, reserve half of the Boltzmann exponent
  to pay its threshold and use the remaining half for
  the outgoing particles. The selected massive pion
  band has no supported high invariant above the split.

For q32, the integrated positive scalar and visible pair
allowances are at most .20001 and .03203 in occupancy
units. The coarse heavy-to-light allowance is 15.39 if
the heavy species could use the entire upper
dark-energy budget. That is a loose bound, **not** an
occupation prediction; it is inconclusive even before
adding the missing gain terms.

| Uniform upper on heavy energy / full available dark-energy envelope | Direct-pair-only dose upper |
| --- | ---: |
| 0 | .23204 |
| .01 | .26474 |
| .02 | .31452 |
| .04 | .44015 |
| .05 | .51239 |

The direct subtotal reaches .50 around a heavy energy
allowance of 4.83%. That fraction is a target for a
tighter comparison, not a computed heavy abundance
or sufficient validation of the actual light cap.
One must also leave room for initial light occupations,
one↔two redistribution, elastic gains and any added
source outside this model.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_light_direct_pair_dose_reference.py EOS.dat \
 "$D"/warm-isolated-fermion-survival-bound-2026-11-26.json \
 "$D"/warm-signed-fermion-with-charged-pions-2026-11-26.json \
 "$D"/warm-isolated-recycling-energy-envelope-2026-11-26.json
~~~
