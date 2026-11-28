# Screened interspecies heat in the plasma-frame transfer coordinate

The preceding [strict vacuum-hard test](INDEPENDENT-EARLY-POLE-HARD-DISTINCT-FERMION-ELASTIC-2026-11-26.md)
retained only (s≥1) and (|t|/s≥.15). That cut is not a
plasma-frame momentum cutoff: the rest-frame exchange satisfies
(Q²=ω²-t). This check keeps the same incoming (s≥1), the reversible
radial box [.45, 8] on both legs, and the distinguishable two-flavor
elastic channel. It tests, but does not complete, the soft part.

## Amplitude and independent phase-space coordinates

With the momentum transfer axis along z, spin-summed conserved
currents give a time-time product (L), a mixed time/transverse
product (I), and a transverse trace (T). The squared Born
kernel, divided by (g⁴), has the terms
(L|D_0|² + T|D_T|² - 2 I Re(D_0 D_T*)).
The vacuum denominators (D_0=1/Q²) and (D_T=-1/t) reproduce
the exact invariant expression, including finite masses, to better
than (8×10^{-13}) in the randomized kinematic audit. The screened
kernel is reversible within (5×10^{-13}) on the same controls.
For the collisionless dressed propagator both complex polarization
components are kept, with the relative sign of the dissipative parts
and their interference. For the general scope of such approximations
see [Arnold, Moore and Yaffe](https://arxiv.org/abs/hep-ph/0209353)
and the [mass-dependent QED example](https://doi.org/10.1103/PhysRevD.107.116016).
Both masses enter the frozen polarization model. The direct
collision quadrature uses massless external kinematics.

The CM formulation divides each azimuth into constant support
intervals using both the outgoing occupation edges and the exact
(Q) boundaries. It can combine vacuum above a given (Q) and
screened below it without overlap. But insensitivity to the matching
boundary alone does not demonstrate convergence: the archive
reconstruction and the resolution where the outgoing energies cross
its narrow features also matter.

The independent direct formulation starts with the plasma-frame
three-vector exchange measure and integrates out the two incident
polar angles with the on-shell energy delta functions. This leaves
`dp dk Q dQ d(ω/Q)/(32π⁵)` before the relative-azimuth average,
with Dirac population factors and the factor one half for the full
reciprocal domain. At fixed (Q,ω,p,k) the squared kernel is at most
quadratic in the cosine of incident relative azimuth. Its three
angular moments are integrated analytically under (s≥1), or both
(s≥1) and the optional vacuum hard cut. This provides an independent
normalization and phase-space test and needs no infrared floor in
(|t|/s). It still has a stated lower bound in (Q).

Continuous archive reconstructions are either monotone PCHIP through
the native shell centers, or the derivative of a PCHIP cumulative
number in (p³/3), which preserves the number of every native shell.
The latter is checked to stay within physical Pauli bounds. For
piecewise constants, quadrature panels on both sides of each sharp
incident edge must be resolved on the (Q) scale; panels with only
Gauss points in the cell interiors can completely miss the transfer
at very small (Q). A uniform-temperature Fermi fixture with
different chemical potentials gives zero net screened heat;
heating only the parent gives positive daughter heat on the same domain.

## Results

In the direct vacuum check, imposing `|t|/s≥.15` again gives
30.61 for fine, agreeing with 30.34–30.84 from the earlier CM runs
despite the new coordinates and independent normalization. Numbers
below use `T_ref⁴` per reference log-scale interval. The frozen
Debye moment gives `m_D/T_ref ≈ .00288` for either archive.

| Direct screened subset | Fine | Third |
| --- | ---: | ---: |
| Cell-number-preserving reconstruction, `Q ≤ .32` | 217.65 | 214.33 |
| Sharp occupations with edge layers, `Q ≤ .04` | 118.73 | 115.91 |

The first row starts at `Q/T_ref=5e-5`. Its contribution below
`.001` is around .14 for fine and .14 for third; a further fine
extension from `5e-5` down to `1e-5` changes that by about one
millionth. Omitting the incident edge layers gives virtually zero
sharp-shell heat below `.005` on both grids. Resolving the layers
instead gives a nonzero result close to the continuous
cell-number-preserving reconstruction: on fine, for `.001 ≤ Q ≤ .02`,
77.51 with refined edge layers versus 77.57 with the refined
continuous representation. The temperature-balance residual in
the manufactured equilibrium fixture is below `5e-9` in these units.

In contrast, for the raw piecewise-constant CM estimator the
`Q=.04` to `.16` matching variation is below 0.1%, yet changing
the CM quadrature shifts the fine result from roughly 120 to 282.
Matching-boundary stability alone does not bound the integration
error. The extra controls isolate a substantial part of that error
as resolution of occupation migration across incoming cell edges.

For the cell-number-preserving case across the full allowed Q
range, splicing screened below `.08` and vacuum above gives
333.37 and 328.71. Moving the splice over
[.04, .16] changes either by less than 0.1%. The result remains
conditional on the archived representation, the radial box, and
the incoming invariant cut; this number is not a total gauge
transport rate.

This is an instantaneous frozen comparison on one channel with
hard radial and invariant cuts, not a fully matched kinetic theory.
Other gauge processes, relaxation of the spectra and expansion are
not included. A soft-only one-particle diffusion moment cannot be
substituted for this species-to-species energy exchange. No scalar
freeze-out or transfer-function input follows.

From `kinetics/`, with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1`:

```sh
python early_pole_frozen_screened_elastic_amplitude_reference.py
python early_pole_frozen_elastic_cm_transfer_reference.py --level 3 --without-fixtures
python early_pole_frozen_elastic_direct_q_reference.py --form cell-conservative --level 1 --window 5e-5 16
python early_pole_frozen_elastic_direct_q_reference.py --grid third --form steps --layer --window 5e-5 .04
python early_pole_frozen_elastic_direct_q_reference.py --form steps --window .32 16 --lower-tau .15
```

Raw: `results/early-pole-screened-elastic-transfer-2026-11-26.json.gz`.
