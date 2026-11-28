# Bound late thermal charged-lepton production into the equality cold band

A production test should target the later physical momentum
cut, rather than treating p/T as fixed across visible entropy
release. In the same [prescribed ballistic background](INDEPENDENT-WARM-LATE-COLD-SCALAR-BALLISTIC-2026-11-26.md),
the maximum injection momentum that redshifts into p/T<.1 at
equality is .1 T_photon,today / a_injection. Use this as the
cut for the existing [thermal charged-lepton vacuum bound](INDEPENDENT-EMPTY-LEPTON-INJECTION-2026-11-26.md).
This only includes creation of a heavy pair by thermal
charged leptons through the prescribed portal and immediate
decay. It does not count an already present heavy reservoir,
indirect or hadronic production, rescattering, feedback or
the growth of Bose occupation.

The required initial pair energy is at least the recoil
gate for one parent plus the rest energy of its partner.
The bound replaces FD by its exponential upper and uses a
global vacuum cross-section maximum. Within each bin the
source upper is largest at the hotter end: the gate rises
as the allowed injection momentum falls during expansion,
and the exponential tail shrinks. Bound H and visible entropy
below using photons only at the colder endpoint. We verify
that the adiabatic temperature clock is monotone on every
EOS interpolation segment and refine the enclosing interval
sum, rather than extrapolating a late numerical rate.

| Starting photon T (MeV) | Minimum initial pair energy/T at start | log10 upper / archived held cold yield |
| --- | ---: | ---: |
| 7.39405 | 327.75 | −118.89 |
| 10.442 | 167.14 | −50.38 |
| 12 | 127.59 | −33.69 |
| 25 | 30.80 | +6.04 |
| 50 | 8.40 | +14.09 |

All rows end at .1 MeV. The archive denominator is the
same corrected q64 equality-selected scalar yield,
4.18022e−9, only for comparison; it does not come from
these five hypothetical starting temperatures. The bound
is strongly restrictive for this late direct source, but is
too loose to answer the same question at 25 or 50 MeV.
A positive logarithm in those rows is not a production
measurement or evidence for an early source. This test
does not bound sources below the tested terminal
temperature or any route that creates scalars at larger
comoving momenta and later scatters them into the cold band.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error warm_late_cold_visible_vacuum_dose_reference.py EOS.dat \
 results/warm-multishell-corrected-q64-fine-stage-48-2026-11-26.json.gz
~~~

[Yield envelopes and partition refinement](results/warm-late-cold-visible-vacuum-dose-2026-11-26.json).
