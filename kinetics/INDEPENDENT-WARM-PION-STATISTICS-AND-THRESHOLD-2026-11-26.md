# Condition the two-pion result on Bose statistics and explicit missing-strip caps

The [measured-band model](INDEPENDENT-WARM-PION-BAND-AND-SOURCE-HEADROOM-2026-11-26.md)
uses Maxwell–Boltzmann weights as a defined vacuum comparison.
Those weights do **not** upper-bound Bose populations.
If the corresponding vacuum spectral kernel is used as an
effective reaction of two charged massive pions, the true
Bose population at fixed T is at most its MB weight times
`B(T)=1/(1-exp(-mπ/T))` per incoming pion. Use B(T)² for the
two-body input, giving a finite factor at 50 MeV and a
rapidly vanishing excess during cooling. This does not
turn the measured inclusive FSR convention into a complete
finite-temperature radiative inverse channel.

The unmeasured interval between `4mπ²` and .10 GeV² is a
separate hypothesis. Use the vacuum two-body
`Rπ(s)=(1/4)(1−4mπ²/s)^(3/2)|Fπ(s)|²`
and supply a constant cap on the latter factor, first 1 then
4. The cap is **not a fitted or measured band**. Integrate its
positive upper cells without overlap with KLOE. A numerical
cap equal to 4 means |Fπ|≤2, not a fourfold increase in |Fπ|.

For a separate conditional bound on charged-pion elastic
heating, use the spin-zero conserved current
`J=pπ+p′π`, with spacelike `J²=4mπ²−t` and the dark
fermion spin average
`4(k·J)²+t J² ≤4(s−mχ²−mπ²)²`.
The equality before dropping the last term uses `q·J=0`
and equal incoming and outgoing dark masses. Multiply the
same positive vector bound by the imposed **spacelike**
form-factor cap, which is independent of the timelike one.
Two charged spin-zero degrees of freedom replace 12 charged
lepton degrees of freedom. Bound the massive second-energy
moment by `24 T⁵ Γ(5,mπ/T)/Γ(5)`, and include both incoming
and final pion Bose factors. These two ratios decrease as
the clock cools, so their 50 MeV values give an all-time
upper for the additionally allowed pion heat.

At the 50 MeV start the pair Bose multiplier is at most
1.13496. Per unit spacelike |Fπ|² cap, the added pion
heat exponent is at most 16.060% of the **already small**
charged-lepton heat exponent, not a 16% energy increase.

| Matched comparison | Scalar Ωh² upper |
| --- | --- |
| MB effective measured band | 0.02761518 |
| Same band with two-input Bose ceiling | 0.02823471 |
| Bose plus threshold and elastic caps =1 | 0.02835030 |
| Bose plus threshold and elastic caps =4 | 0.02869973 |

This covers only the specified effective charged two-pion
channel and capped elastic kernel. It does not put a
model-independent bound on pion form factors, thermal
hadronic rates, other hadrons or photons, or new numbers
of light scalar particles. Neither these moment inequalities
nor the measured window determine a relic abundance.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_pion_bose_and_threshold_reference.py EOS.dat KLOExSec.data \
  "$D"/warm-energy-weighted-recycling-kernels-2026-11-26.json
~~~
