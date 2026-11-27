# Direct archive occupations remove most apparent hard-aperture dependence

The [curvature diagnostics](INDEPENDENT-EARLY-POLE-DAUGHTER-LOCAL-FRACTION-2026-11-27.md)
implicate the two-moment occupation reconstruction in the initialized
local rate response. For a separate check, hold each mesh's native
packet geometry, native parent occupation, and scalar seed fixed.
At the same sampled daughter momenta, evaluate the raw archive
step occupation and the derivative of the monotone cumulative
number interpolation, and compare their effect with what the
corresponding fitted daughter states give. Only packets whose
daughter is in the fixed original interval are changed.

With fixed parent and scalar occupations, the paired packet
difference is proportional to minus their sum times the daughter
occupation difference. This form avoids subtracting two potentially
large gain/loss totals; it also independently reproduces the
previous fixed-geometry curvature within 1e-16.

| Hard depth × cold modes | Direct archived occupations | Two-moment fitted occupations | Fitted minus direct |
|---|---:|---:|---:|
| 1 × 8 | -0.0116528081 | -0.0155111946 | -0.0038583866 |
| 4 × 8 | -0.0116526849 | -0.0116709653 | -0.0000182804 |
| 1 × 32 | -0.0122152928 | -0.0162510964 | -0.0040358036 |
| 4 × 32 | -0.0122151925 | -0.0122359324 | -0.0000207398 |

The discrepancy between fitted and direct evaluations falls by
about 211× at eight cold modes and 195× at 32 cold modes. The two
native packet meshes give direct rate differences within 1.24e-7
at a fixed cold discretization. Direct means evaluating the two
archived occupation prescriptions, not knowing the true within-shell
distribution. These calculations retain the finite set of packet
samples and do not estimate an error bound for the physical source.
They neither remove the cold-mesh sampling effect nor permit use
of these frozen occupations as an independently evolving conservative
closure.

Run `early_pole_daughter_archive_pointwise_frozen_packets.py` on
the same four fixtures named in the local-fraction report. In each
fixture the raw file `results/early-pole-spectrum-fine-m06-2026-11-26.json.gz`
supplies the original daughter steps and their cumulative number
interpolant. The JSON receipt includes every mode and daughter cell
sampled in this comparison.
