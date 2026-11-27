# Fewer above-cut daughters can coexist with growing excess energy

Use the [three archived near-cut histories](INDEPENDENT-WARM-LATE-NEAR-CUT-SHAPES-2026-11-26.md)
at their common fixed energy cutoff C=0.7148213914435242 GeV.
Count particles with E>C and separately accumulate energy
weighted by (E−C) for the same daughters. The resident
moments independently satisfy E_above = C N_above + W.
All held species number is used as the normalization for
the count; W is normalized by the full held dark energy.

The two within-species angular count currents are negative at
all nine tested states (Hdt=0, 3e−8 and 1e−7 in each case).
In contrast, .710 / .005 and .713 / .0025 have positive
total W current at the first two samples. For the narrower
.713 test, the initial light count current is −1.199e−7 of
its species count per Hubble time, while the light W current
is +2.531e−9 of the total dark energy per Hubble time.
At the endpoint only 8.13% of the initial light above-cut
count remains in that case, but its light W is 1.878 times
the initial light W. The broad .713 comparison also loses
above-cut counts, but its light W falls rather than rising.

Two independent angular quadrature orders resolve the
signed rates relative to the gross ledgers. The sum of the
species' W current agrees with the native projected rate.
The occupation integration across the cut is split and
checked at 16 and 24 nodes per cell. The rates are derivatives
at selected states, not interpolated event totals between them;
the reported endpoint count change comes from the saved spectra.
These are closed, fixed-background prompt-only controls,
not estimates of an absolute relativistic population or
indicators of cosmological clustering.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MALLOC_ARENA_MAX=2
# From kinetics/, extract the exact archived phases.
python - <<'PY'
import gzip, json
from pathlib import Path
for case in ("710-3", "713-3", "713-narrow-3"):
    source = Path(f"results/warm-late-near-cut-{case}-2026-11-26.json.gz")
    out = Path(case)
    out.mkdir(exist_ok=True)
    for name, body in json.loads(gzip.decompress(source.read_bytes())).items():
        (out/name).write_text(json.dumps(body, allow_nan=False))
PY
python -W error warm_late_threshold_count_current_reference.py 710-3 713-3 713-narrow-3
~~~

[Count, full energy, excess, and signed/gross angular controls](results/warm-late-threshold-count-current-2026-11-26.json).
