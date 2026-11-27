# The locally subtracted real extension crosses null transfer smoothly

Combine the independently normalized [positive emission](INDEPENDENT-FINITE-MASS-QED-BARE-EMISSION-MONOTONE-BRIDGE-2026-11-27.md)
and [negative absorption](INDEPENDENT-FINITE-MASS-QED-ABSORPTION-CROSSES-NULL-TRANSFER-2026-11-27.md)
charts over the photon band [.0002,.00256]. Split the integration at
both angular endpoints where the absorption transfer first and last
crosses its null surface. The Born baseline is the exact zero-photon
eikonal result. Positive gain, inverse, and photon Bose factors are
kept separate before applying the [local leading subtraction](INDEPENDENT-FINITE-MASS-QED-ODD-TURNOVER-BRIDGE-2026-11-27.md).
The bare real photon energy moment has opposite relative signs
for emission versus absorption and needs a different combination.

| Archive | Leading-subtracted scalar source / Born | Real photon heat / Born scalar energy |
|---|---:|---:|
| fine | −6.10562e−5 | −5.73444e−4 |
| third | +1.86128e−5 | −1.39799e−3 |

These are single-charge frozen comparisons at one cold scalar mode,
with temperature 1 and g=.005 in archive units. The sign of the
subtracted scalar number differs between archives, while the real photon heat is negative
in every tested subinterval of both archives. The largest absolute
difference of a band ratio between the two refinement rules is
5.81e−8. Shared absorption checks at .00016, .00019, .0002 agree
within 1.08e−5 with the earlier signed parent chart. This is not
a physical inclusive correction: finite virtual terms and dressed
dispersion/residue contributions are not included.

Replay from kinetics:
~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_signed_transfer_match_reference.py
~~~
