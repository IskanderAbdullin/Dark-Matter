# Most of the warm pressure-source contrast appears before changing the occupations

The [independent warm W2 evolution](INDEPENDENT-EARLY-POLE-ENDPOINT-SCALAR4-WARM2-EVOLUTION-2026-11-28.md) lowers reconstructed terminal pressure by 6.23%, while total yield barely changes. Separate three frozen comparisons at each matched checkpoint: the actual S4 state on its original scalar grid, exactly the same constant-in-cell S4 spectrum copied to the W2 subcells with identical charges, and the independently evolved W2 state on the W2 grid. Hold the mass and expansion normalization at the S4 sample's values. None is an additional evolving solution.

Project the signed collision count source of a scalar node onto its constant-in-cell reconstructed moment. The pressure weight in a cell is the volume average of `q²/sqrt(m²+q²)` under q² dq; for the tail-only pressure weight, the integral is cut at q=.01 while the cell volume denominator is unchanged. Compute forward, inverse, and signed moments separately. The energy-projection, charged/scalar pairing, 4→8 collision quadrature, and 32→64 within-cell weight comparisons all pass. There are no expansion or varying-mass terms in this diagnostic derivative.

| x | Frozen warm-grid change, total pressure source | Additional change from the independent evolved state on W2 | Combined change |
| --- | ---: | ---: | ---: |
| .01 | −6.25996% | +1.72e−8 relative | −6.25996% |
| .02 | −6.29100% | +1.68e−7 relative | −6.29099% |
| .04 | −6.20798% | +1.19e−6 relative | −6.20787% |
| .05 | −6.15754% | +2.00e−6 relative | −6.15736% |

Thus in this controlled instantaneous comparison the geometry already produces nearly the entire observed source contrast at fixed occupations. The total scalar number sources at the same points differ only at the level of a few parts in a million. This does not establish an exact decomposition over the historical trajectory, control higher warm refinements, or show which numerical chart gives the physical pressure. The measured pressure source is a functional of the reconstructed finite-grid numerical state.

Run `early_pole_endpoint_warm_pressure_source_control.py` with the matched archived S4 and W2 exact states, and `--x .01`, `.02`, `.04`, or `.05`. The [final receipt](results/early-pole-endpoint-warm-pressure-source-x05-2026-11-28.json) gives all the absolute and signed coefficients, including the warm share.
