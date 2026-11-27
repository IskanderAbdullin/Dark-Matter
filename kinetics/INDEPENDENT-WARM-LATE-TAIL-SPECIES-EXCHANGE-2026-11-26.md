# Separate the two daughter currents and scan the hard-packet sign response

In the late [cut-aligned .68 history](INDEPENDENT-WARM-LATE-TAIL-TIME-REFINEMENT-2026-11-26.md)
and [.80 control](INDEPENDENT-WARM-LATE-ALIGNED-HARD080-RELAXATION-2026-11-26.md),
measure the scalar and light-fermion contributions to the
same signed convex excess above C=.714821 GeV. At each
parent energy the total formation and decay rates cancel
under the stationary prompt closure. Thus the excess can
be assigned directly to daughter emission minus daughter
consumption, separately for each species. The angular
sum reproduces the independently conserved total excess
current in three quadrature controls.

| Case and Hdt | Scalar current | Light current | Total current |
| --- | ---: | ---: | ---: |
| .68 start | +1.547e-12 | +1.462e-10 | +1.478e-10 |
| .68 at 6.25e-8 (locked 24-step run) | +1.010e-12 | −3.272e-13 | +6.832e-13 |
| .68 at 6.625e-8 | +9.810e-13 | −1.779e-12 | −7.983e-13 |
| .68 at 1e-6 (full finer path) | −1.242e-13 | −5.102e-13 | −6.344e-13 |
| .80 start | −9.184e-8 | −8.571e-7 | −9.489e-7 |
| .80 at 1e-6 (full finer path) | −2.869e-8 | −1.559e-10 | −2.884e-8 |

The rates are fractions of retained dark energy per Hubble
time. At the .68 total reversal the scalar contribution is
still positive; the light loss has grown enough to exceed it.
The .80 initial loss is mostly light, whereas the terminal
negative current is mostly scalar. These are changes of
the energy **in the tail**, not a net loss of particle
number from either full species.

## Sensitivity of the frozen sign to the hard-packet center

Keep the old receipt counts and exact mean birth energies,
hold the pad at 1.05 GeV, and vary the hard packet center
across 26 frozen choices. The sign near the threshold
depends on the chosen physical shell width, even though
those original moments are held fixed. In the 5 MeV
half-width family the sign bracket is [.710,.712] GeV;
in the 2.5 MeV family it is [.712,.714] GeV.

Four local crossing estimates and checks on either side,
using the same projected packet construction:

| Half-width (MeV) | Local sign change near center (GeV) |
| ---: | ---: |
| 5 | .711600 |
| 2.5 | .713603 |
| 1 | .714576 |
| .25 | .714831 |

The two daughter momentum cuts are .714821 and .714542 GeV.
The narrow-shell result is near the scalar cut, but this
finite family does not prove that the transition equals
either cut in a zero-width limit. These roots concern a
frozen onset and do not assign a terminal state, sign or
lifetime to the unknown exterior distribution.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
a=results/warm-corrected-anchor-signed-floor-2026-11-26.json
b=results/warm-multishell-corrected-q64-fine-stage-48-2026-11-26.json.gz
python -W error warm_late_hard_packet_cut_sign_scan_reference.py EOS.dat "$a" "$b"
python -W error warm_late_hard_packet_cut_center_refinement_reference.py EOS.dat "$a" "$b"
python -W error warm_late_species_cut_exchange_reference.py EOS.dat "$a" \
  results/warm-late-cut-aligned-hard068-test-step0-2026-11-26.json.gz \
  results/warm-late-aligned-turning-sample14-2026-11-26.json.gz \
  results/warm-late-aligned-turning-sample15-2026-11-26.json.gz \
  results/warm-late-aligned-hard068-terminal-6-2026-11-26.json.gz \
  results/warm-late-aligned-hard080-test-step0-2026-11-26.json.gz \
  results/warm-late-aligned-hard080-terminal-6-2026-11-26.json.gz
~~~

[Species-resolved currents](results/warm-late-species-excess-transfer-2026-11-26.json);
[26 centers on each width](results/warm-late-frozen-hard-center-scan-2026-11-26.json.gz);
[local crossing controls](results/warm-late-frozen-hard-center-roots-2026-11-26.json).
