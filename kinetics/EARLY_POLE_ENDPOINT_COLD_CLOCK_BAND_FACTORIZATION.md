# The explicit scalar population accounts for about 1% of the final-clock band change

At each of S2 and S4, hold the independently evaluated coarse and fine cold-band kernels. The raw q panels never straddle a scalar cell and the native scalar geometry is identical on the two clocks for a given level. Divide each original forward coefficient by its own `1+f` and each inverse coefficient by its own positive `f`; replace that explicit factor with either clock's occupation without interpolating or redoing any integral. The diagonal reconstructs the actual measured bands. Compare the four positive entries in log units by averaging both replacement orders.

| Cold band | Explicit population share S2 | Explicit population share S4 |
| --- | ---: | ---: |
| Forward number | 1.1161% | 1.1146% |
| Inverse number | 1.0264% | 1.0246% |
| Forward kinetic excess | 1.1021% | 1.0943% |

The remaining approximately 99% is carried by the combined bare response kernel. That kernel includes the scalar dispersion and phase-space measure, fitted photon energy and residue, and charged-occupation closure. This calculation **does not** further decompose those terms or establish what causes their differences during evolution. It only separates the multiplicative scalar population from the already computed whole-band response on two final-step trajectories. Both accepted q meshes are checked on the crossed as well as native entries.

Use the same axis-specific q-panel arguments and state order as `early_pole_endpoint_cold_scalar_band_clock_contrast.py`, with `--matched results/early-pole-cold-scalar-band-clock-contrast-2026-11-28.json` for `early_pole_endpoint_cold_clock_band_factorization.py`. The argument binds the four original exact state hashes and native gross integrals.

## Separately tested optical and charged factors at two cold points

At q=3.5e-6 and 4.7e-6, independently re-integrate all eight positive corners of the scalar-gap, photon-energy/residue, and charged-occupation replacements for both S2 and S4. Restore the explicit scalar occupation factor analytically. Average the log changes over all six orders of the three bare terms; the final path exactly matches the physical corners of the computational fixture. Separate direct coarse and fine point evaluations recover these native endpoint ratios within 4.4e-13 relative; paired angular quadratures for every corner agree to about 6.1e-13 or better.

The scalar-gap term accounts for 98.71–99.03% of the local forward log change and 98.81–99.11% of the local inverse log change. The explicitly factored photon energy and residue together have magnitude close to 0.001% of that log change; charged closure is also negligible at these sampled q. This is direct support for the scalar-gap explanation at those two points, but the full-band test above distinguishes only **explicit occupation** from the combined kernel. A full-band subdivision of the bare terms has not been computed.

With the same four ordered states and `--matched` receipt, use `early_pole_endpoint_cold_clock_point_factors.py` (no panel or optical-reference arguments). Crossed intermediate corners are simply held numerical response evaluations, not self-consistent trajectories.
