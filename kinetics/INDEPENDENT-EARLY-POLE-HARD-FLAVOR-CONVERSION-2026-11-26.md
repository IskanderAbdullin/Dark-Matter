# Hard pair conversion changes the charged chemical budget

A finite additional Abelian channel is
`P + Pbar ↔ b + bbar`. It does not overlap the previously
checked [elastic and Compton subsets](INDEPENDENT-EARLY-POLE-SAME-FLAVOR-REDISTRIBUTION-2026-11-26.md)
or the photon pair channel. Use the single vacuum vector
s-channel, with both parent charge states and both daughter
charge states labelled and the charge distributions assumed
equal within each flavor. Select `s≥1` and a finite window
on all four momenta. The [massless Born normalization](https://edu.itp.phys.ethz.ch/hs10/ppp1/Chapter5.pdf)
gives `g⁴(1+cos²θ)/(32πs)` for the CM cosine differential.
Finite mass changes both recoil and the spin polar polynomial.

The initial phase space includes one labelled P/Pbar
orientation with two spin states per charge. The reverse
term is evaluated on exactly the same kinematic states.
A second half for net reciprocity would be incorrect here:
the outgoing flavor is not integrated as a separate initial
reaction in this forward-oriented computation. Each event
changes parent/daughter counts by −2/+2 and transfers the
full incoming pair energy from parent to daughter. A fixture
at common temperature and common chemical potential has
vanishing net transfer; different chemical potentials at
that same temperature activate conversion.

## Frozen and live initial controls

On [.45, 8], after splitting every incident occupation
edge, vacuum massless integrations give daughter number/energy
rates +28.14/+70.58 (fine), +28.36/+71.14 (third), in
`T_ref³` and `T_ref⁴` per reference log-scale interval.
Independently integrating the inverse orientation gives
−28.10/−70.50 and −28.33/−71.06 for parent. This all-edge
reversal test agrees at around 0.12%, whereas coarser
unsplit incident grids can give a few percent mismatch.
The finite bare masses change the corresponding coarse
forward rates by about 0.1%, with no sign reversal. In the split
archived calculation the gross inverse event flux is less than
0.7% of the gross forward one; this is not a near-equilibrium
chemical state.
Expanding the common box to [.30, 10] changes the coarse
rates to about 31 and 77–80.

The live moment projection uses massless weights on [.45, 5],
then remaps one final radial momentum to exact common
variational energy and checks the outgoing momentum triangle.
No photons are touched by this operator. The combined
charged count, energy and nonnegative entropy production
are audited. The finer live initial daughter number/energy
rates at the reference coupling are +23.57/+52.61 and
+23.80/+53.20, close to the independent finite-window freeze
values +23.63/+52.76 and +24.39/+55.00 (fine and third).

In the same live box, the existing photon-pair packets supply
+29.37/+28.75 parent counts (fine/third), while the conversion
removes −23.34/−23.58. It offsets about 79%/82% of that
specific modeled parent supply. Including the much smaller
scalar decay count source leaves an initially **positive**
parent balance of +5.95/+5.10, not a net parent depletion.
The daughter gain increases from about 127 to about 150.
All numbers here are in the same reference clock and live
cell support, but the photon-pair and conversion angular gates
are their own finite model domains.

## Evolving finite comparisons

At the standalone fixed-scale exposure `.0003`, with photon
pairs on but without screened elastic scattering or live
Compton, adding conversion changes the cold scalar gain by
−5.42% (fine) and −6.68% (third); without photon pairs the
changes are −5.20% and −6.28%. In the shorter joint test
(`2e-6`, twelve split steps), with all four previously tested
channels and the active unlike-polarization refresh enabled,
adding conversion changes the gain by about −.037% and −.035%
on the two native step states. The same-flavor soft pole
is still only a finite subset with frozen response. All the
joined gains remain positive, with variational energy and
entropy ledgers closed. The effects at different intervals and
channel settings are not interchangeable.

These are finite hard subsets with frozen exterior and
an early reference clock. They do not establish complete
chemical equilibration, a late scalar spectrum, or the
cosmological transfer function.

Reproduce from `kinetics/` with one BLAS thread:

```sh
python early_pole_hard_flavor_conversion_reference.py --level 3 --incident-edges --without-fixtures
python early_pole_hard_flavor_conversion_reference.py --level 3 --incident-edges --without-fixtures --reverse
python early_pole_archived_live_flavor_conversion_reference.py third --level 1
python early_pole_archived_live_flavor_conversion_reference.py fine --evolve
```

Raw: `results/early-pole-hard-flavor-conversion-2026-11-26.json.gz`.
Joined-channel data: `results/early-pole-live-joined-charged-channels-2026-11-26.json.gz`.
