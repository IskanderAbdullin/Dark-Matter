# Removing the adjacent donor cells has only a small initial elastic effect

The [finite held chart](INDEPENDENT-EARLY-POLE-DAUGHTER-DONOR-SUPPORT-TIME-2026-11-27.md)
found that removing the two adjacent allowed donor cells, while
restoring the same full retained daughter number and bare energy,
raised the mean cold contrast by 2.12016e-4. Isolate the initial
elastic direction using the same fixed hard geometry and cold seed.
For each donor support, also reset the two active cells to native
while keeping the fitted outside compensation fixed. This
compensation-only control does not restore the full two moments
after the reset; it tests the initial slope at the reset input.

| Initial directional comparison, no-neighbors minus all-outside | Across three probes |
|---|---:|
| Retuned local, both elastic channels | +6.838e3 to +7.080e3 |
| Compensation-only, both channels | +5.920e3 to +6.178e3 |
| Difference of the two rows | +8.61e2 to +9.18e2 |
| Retuned local, unlike channel alone | +6.798e3 to +6.840e3 |

The donor difference is 0.221–0.229% of the original
uncompensated local-minus-native initial elastic slope difference;
the compensation-only control already gives most of it. The two
retuned initial rates differ by only +1.288e-8. Slope and finite
mean have different units, so their relative sizes do not quantify
a causal contribution to the finite contrast. The probe spread
and directional additivity defects are retained in the receipt.

From `kinetics/`, run
`python early_pole_daughter_initial_aperture_flux_reference.py results/early-pole-joined-fine-216-memory-2026-11-27.json.gz --include-retuned --include-donor-variant`.
The earlier `--include-retuned` run replays byte-for-byte, with
an optional held local snapshot checking its hard and cold input.
This is a finite directional diagnostic in the archived model,
not a relic result.
