# Spacelike collisionless noise on closed patches away from charged sequential roots

A spacelike photon cut has an occupied on-shell intermediate-state aperture in the [earlier local test](INDEPENDENT-FINITE-MASS-QED-CUT-LOCAL-SEQUENTIAL-SUBTRACTION-2026-11-27.md). Here deliberately restrict the charged parent energy to two closed intervals `[1.18,1.26]` and `[1.85,2.25]` at the same three representative spacelike external geometries. Both chosen patches lie away from the charged roots. This is a positive finite support/control calculation, not a cut subtraction, matching recipe, or total scalar rate.

For fixed K and q the angular parent shell defines p_parallel linear in E. Its lab z component is a linear center plus a term proportional to p_perp cosφ. Thus the possible D1=(P+K)²−M² values at each E are bounded by two functions linear in E plus or minus a constant times p_perp(E). The squared p_perp is an upward quadratic with positive discriminant, so its square root is concave on the positive interval. Each envelope derivative is monotone, and the global extrema come from an endpoint or one derivative zero. Using D1+D2=2 q·K provides the second denominator interval. This bounded test proves a minimum absolute denominator of 2.82766e−4 across every angular and energy point in the selected supports, rather than pitching an arbitrary small-width cut on a resonant integral.

Use the leading finite-mass frozen retarded [transverse cut](INDEPENDENT-FINITE-MASS-QED-FROZEN-TRANSVERSE-CUT-2026-11-27.md) or [Coulomb density cut](INDEPENDENT-FINITE-MASS-QED-LONGITUDINAL-SUM-2026-11-28.md). For each channel, its same leading collisionless noise-to-response ratio T_eff(u) determines a symmetric source. In positive-frequency real spectral normalization the two continuum coefficients are

```text
D^> = ρ [T_eff(u)/ω + 1/2]
D^< = ρ [T_eff(u)/ω − 1/2].
```

Both are positive at the selected small frequencies, and they reconstruct the same retarded commutator ρ. The symbol T_eff isn't an imposed background temperature. This calculation uses only the charged-sourced continuum; it adds neither a homogeneous incoming population nor any timelike pole. With separately evaluated bare directional currents A,B and archived scalar occupancy n_s, the unnormalized signed expression divided by n_s is
```text
(D + A/n_s) D^< − B ρ,  where D = A−B.
```
It deliberately omits the common contact factors and k,ω and angle integration. In the thermal control the exact pair ratio D^</D^>=exp(−ω/T) balances the independently integrated Fermi pair and Bose scalar exactly; the nonthermal archive is not assigned that law.

The 24 measured patch and polarization samples agree upon joint angular/energy mesh refinement to 9.64e−12 relative. The coherent Ward check is below 1.22e−19 in squared norms, and the two equivalent longitudinal charge/current projections agree to 1.72e−10. The independent common-temperature detailed-balance factor closes to 2.23e−16. Although the sampled conditional patch sources are positive, none of these figures constrains the sign after adding other geometries, resonant regions, virtual terms or a consistently matched medium calculation.

Run `python -W error finite_mass_qed_spacelike_off_root_patch_reference.py` from `kinetics/` for the [receipt](results/finite-mass-qed-spacelike-off-root-patch-2026-11-28.json).
