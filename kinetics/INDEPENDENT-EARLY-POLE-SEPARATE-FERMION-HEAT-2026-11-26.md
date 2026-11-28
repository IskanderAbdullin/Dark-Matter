# Separate fermion temperatures alter the stimulated window

Relax the instantaneous common-temperature assumption of the
[shared photon decay comparator](INDEPENDENT-EARLY-POLE-PHOTON-DECAY-COMPARATOR-2026-11-26.md).
Each Dirac species has its own conserved moments and instantaneous
FD fit within that species. Exact projected decay events remove
parent energy E_scalar+E_daughter and add the two daughter energies.
Both hard gauge pair sources transfer heat to or from a finite
number-budget photon reservoir. The auxiliary interspecies heat
source moves the daughter energy toward its common-temperature
value at the given two counts and combined energy, without
altering either count. Its rate is a freely chosen comparator.

When the temperatures differ, the sign of the stimulated
per-event factor depends on the parent energy:
the Fermi gain condition is
(β_d−β_P)E_parent+β_P μ_P−β_d(μ_d+E_scalar)>0.
For a given scalar mode we therefore integrate the signed
occupation difference over the exact fermion aperture rather
than using a chemical-potential difference alone. Reported
crossings are the first resolved negative **window averages**,
not the time when a total number source vanishes.

For the third-grid benchmark, the initial daughter and parent
temperature fits are .37183 and 1.03735 in reference units.
At g=.005, the cold-mode average crosses near x=.050191 with
no imposed cross-species heat exchange, near .050317 at rate 100
per reference expansion time, and near .050962 at rate 1500.
The last case has near equal final temperatures; refinement
of the hard pair quadrature gives .050970. The rate-100 other
momentum grid gives .050308. At g=.004 and heat rate 100 the
average has not crossed by .052. These are relative sensitivity
checks inside the chosen hard-cut, fast-within-species model.
They neither fix a physical heat rate nor assert net decay
stops.

~~~sh
python -W error early_pole_separate_fermion_heat_reference.py results/early-pole-spectrum-third-m06-2026-11-26.json.gz --exchange 100
~~~

[Raw stress–energy and occupation-window checks](results/early-pole-separate-fermion-heat-2026-11-26.json.gz).
