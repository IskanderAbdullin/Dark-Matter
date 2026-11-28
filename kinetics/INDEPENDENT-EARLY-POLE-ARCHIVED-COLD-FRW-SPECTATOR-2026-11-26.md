# Conserving redshift of the archived live cold window and its exterior

The [fixed-scale cold window](INDEPENDENT-EARLY-POLE-ARCHIVED-COLD-LONG-REFERENCE-2026-11-26.md)
does not move the resonance relative to the archived hard momentum
cells. This comparison extends the
[earlier toy expansion ledger](INDEPENDENT-EARLY-POLE-LIVE-JOINT-HARD-PACKETS-FRW-2026-11-26.md)
to the actual initialized charged archives and the same eight common
scalar intervals. Near-aperture cells use the actual archive splits.
A conformal scalar gap is still assumed.

Exterior charged cells do not evolve their comoving occupations,
but their massive comoving free energy, reciprocal field moments
and trace are recomputed from their exact archived support at each
scale. The photons outside the live hard cells retain a frozen
initial BE spectrum and enter the total energy, including below the
hard field cutoff. Only the prescribed hard portion enters the
gauge moment. The scalar spectators contribute to the same live
gap as the evolved scalar bins. No exterior contribution is silently
counted as an active collision source.

Two clocks use the same adopted visible radiation reference:
one ignores the dark contribution to H, the other rescales all
collision rates by 1/sqrt(1+R_dark/R_visible) with the explicitly
represented dark energy. The pair vacuum weights remain the initial
ones, and the hard pair prefactor and vacuum scalar normalization
follow the short reference expansion. This is a test of the
represented energy, not a complete Friedmann cosmological history.
In this initialization R_dark/R_visible≈.0514, which gives H
about 2.5% larger than the visible-only reference.

| Archived state | Clock | Decay-only scalar gain | Combined scalar gain |
| --- | --- | ---: | ---: |
| Fine | Visible reference | 6.24464e−6 | 6.02544e−6 |
| Fine | Conditional full H | 6.09157e−6 | 5.88287e−6 |
| Third | Visible reference | 5.77917e−6 | 5.53968e−6 |
| Third | Conditional full H | 5.63864e−6 | 5.41056e−6 |

The reference scale interval is ln(a/a₀)=8e−5. The separately
integrated redshift trace matches the comoving energy change to
within 6e−16 in these short controls, with the collisional work
and common number independently tested. The clock correction
slows both rates and has little effect on relative pair suppression
in this short window, about 3.4–4.1%. It does not validate the
physical scalar abundance or a late spectrum.

Reproduce from kinetics/ using
early_pole_archived_cold_frw_spectator_reference.py --grid fine/third
--expanding-only. Raw:
results/early-pole-archived-cold-frw-spectator-2026-11-26.json.gz.
