# Joint live decay and pair moments share one conserved energy

Combine the [decay-cell reference](INDEPENDENT-EARLY-POLE-LIVE-PAULI-GAUGE-MOMENT-DECAY-2026-11-26.md)
and the [hard photon-pair reference](INDEPENDENT-EARLY-POLE-LIVE-RECIPROCAL-PAIR-MOMENTS-2026-11-26.md)
in one evolving state. The same photon and two fermion fields,
local number/energy distributions and reciprocal gauge functional
serve both channels. A selected scalar mode has a variational Hartree
energy, not a separately prescribed energy gain. Each channel
contributes to the same count, energy and entropy ledgers.

From one starting state, vary the imposed strength of pair
production and inverse annihilation together. The comparison has
two incident photon samples and the relative weights below are
test settings, not couplings or cosmic time.

| Pair multiplier | Decay enabled | Scalar count gain after exposure .14 |
| ---: | :---: | ---: |
| 0 | yes | .0006234747 |
| 1 | no | 0 |
| .5 | yes | .0005716821 |
| 1 | yes | .0005381282 |
| 2 | yes | .0004899186 |

In this example pair dynamics reduces the scalar increment by 13.7%
at multiplier one relative to decay alone. An additive calculation
would miss the shared occupation factors. This is not a bound:
it need not have the same sign at a different bath state or across
the hard/cold integration domain.

The adjusted events close the total energy to 6e−17 or better, and
the step-by-step event entropy production agrees with the net entropy
change. Reverting just the events to bare energies while retaining
the interaction energy gives a +3.8893e−9 drift, independently matched
by the collisional work ledger. Total photon plus charged-particle
count, integrated photon exchange and integrated scalar count are
separately audited. Refining from 10 to 20 outgoing quadrature
nodes changes the scalar gains by less than 5e−14.

Next unresolved input: replace the two incident photon samples with
converged hard integration and a matrix element, and verify the
combined model under expansion. There is still no physical rate,
gauge threshold or CLASS-ready spectrum.

Reproduce from kinetics/ with
early_pole_live_joint_pair_decay_moments_reference.py --order 10 or --order 20.
Raw: results/early-pole-live-joint-pair-decay-moments-2026-11-26.json.gz.
