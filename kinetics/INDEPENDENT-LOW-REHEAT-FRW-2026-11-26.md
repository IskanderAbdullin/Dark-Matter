# Low-reheat FRW control near the prescribed stable-count root

An initially empty finite-shell run at T_R=9.0190699 MeV evolves the
physical vacuum visible lepton source, prompt daughter injection,
collisionless redshifting of the retained dark species, and the
visible energy withdrawn by both retained and exterior daughters.
The "vacuum" rows do not evolve any dark collisions or inverse visible
reaction. The grid retains q/T_R in [0.8,16], on 18 cells. Use the
[official low-temperature EOS extract](INDEPENDENT-BENCHMARK-SETUP.md)
with the positive knots at and above 1 MeV. Integrating the interpolated
entropy to construct thermodynamically compatible pressure changes
the radiation density by up to 0.137% relative to separately
interpolated tabulated radiation in these tests.

| Endpoint | Steps | Retained stable Ωh² proxy | Adding exterior source | Independent full vacuum source |
| --- | ---: | ---: | ---: | ---: |
| 3.157 MeV | 24 | 0.118238 | 0.119771 | 0.119965 |
| 1.127 MeV | 24 | 0.117739 | 0.119258 | 0.120000 |
| 1.127 MeV | 48 | 0.118276 | 0.119809 | 0.120000 |

The finest all-support ledger is 0.160% below the independent
stable-source integration; until the 48-step mesh the early source
pulse is less accurately integrated. About 1.28% of the inclusive
stable yield lands outside the retained grid. The corresponding scalar
proxy is 2.28e-12, still 3.47% below the independent 2.36e-12; the
endpoints do not resolve the rapid early heavy-pair source equally well.
The final visible temperature is lower by 3.56e-8 at the same scale
factor. The vacuum source still below 1.127 MeV adds only 8.63e-15
to the prescribed stable density proxy through 1 MeV. This is not a
late-loss bound, scalar dominance, complete coupled dark evolution,
or a reheating/BBN viability result.

Refining the vacuum runs to 96/192 intervals gives all-support stable
proxies 0.1199465/0.1199810 at 1.127 MeV, the finest within 0.016% of the
independent integral. The corresponding scalar proxies are
2.3426e-12/2.3582e-12, the finest within 0.23%. Stability against mesh
refinement does not validate other omitted processes.

A matched **state-dependent** calculation now reaches 3.157 MeV on both
12 and 24 intervals. It includes signed visible FD exchange, prompt
return, below-pole elastic heat, and visible energy feedback. On the
finer mesh the inclusive fixed-count stable proxy is 0.1197706215
rather than 0.1197706226 in the matched vacuum control (relative
change about -9.6e-9). The agreement compares only these implemented
operators. The new stable-pair Yukawa annihilation to scalars, other
number-changing reactions, and exterior re-entry are not included
in this evolution. No quantum run to 1.127 MeV is claimed.

The cold signed-shell portal now incorporates the visible Boltzmann
factor on the **forward** term before subtracting the independently
bounded inverse term. It does not form a standalone exp(q/T) and does
not drop inverse events when exp(-q/T) has underflowed. At 1.127 MeV,
quadrature nodes extend above q/T=8900; vacuum direct and heavy shell
moments agree with the existing positive integrator, the direct
thermal cancellation is verified, and an artificial inverse-only
tail with both incoming energies above 1 GeV remains nonzero. The
artificial tail tests arithmetic, not a physical distribution.
The finite-box equilibrium diagnostic also brackets its fermion
chemical potential above the rest energy at cold temperatures.

## Reproduce

From `kinetics`, set `EOS.late.dat` to the official extract. Run:

```sh
export OPENBLAS_NUM_THREADS=1
cmd=(python -W error joint_resonant_frw_split_reference.py EOS.late.dat --eos-lower-knot-limit .001 --physical-amplitude --empty-pulse .00002 --empty-reheat .0090190699 --visible-source --visible-backreaction --support-low .8 --empty-cells 18 --source-only --matched-vacuum-order)
"${cmd[@]}" --empty-stop .35 --steps 12 24 > matched.json
"${cmd[@]}" --empty-stop .125 --steps 12 24 48 > deep.json
"${cmd[@]}" --empty-stop .125 --steps 96 192 > refined.json
python -W error joint_resonant_frw_split_reference.py EOS.late.dat --eos-lower-knot-limit .001 --physical-amplitude --empty-pulse .00002 --empty-reheat .0090190699 --visible-source --visible-backreaction --support-low .8 --empty-cells 18 --empty-stop .35 --quantum-visible-source --with-elastic --vacuum-scaled-residual --allow-empty-shells --steps 12 24 > quantum.json
python -W error low_reheat_frw_count_check_reference.py EOS.late.dat matched.json quantum.json deep.json refined.json
```
