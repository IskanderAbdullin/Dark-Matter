# Resolve the changing contributions across travelled-distance bands

The [distance ceilings](INDEPENDENT-WARM-LATE-COLD-TRAVEL-SOURCE-2026-11-26.md)
settle the smallest path bands but were too loose at 1 Mpc.
For the numerical estimate, find the momentum crossing
at the actual injection time for each quadrature node.
Fold the known charged-lepton FD kernel and the measured
charged-pion KLOE BE window into that moving cutoff.
All sources are vacuum direct creation only, with
the common 12–50 MeV expansion and no later collisions.

| Maximum path to reference equality | Measured-window pion + lepton numerical yield / held count |
| --- | ---: |
| .3 Mpc | 3.413e−7 |
| .5 Mpc | .002613 |
| 1 Mpc | .41667 |
| 2 Mpc | 1.8000 |

The numerical .3 Mpc value is well below the earlier
conditional upper for the full considered window. These
are different cumulative populations; the rapid relative
increase at longer paths does not remove the deficit at
shorter paths. The 2 Mpc numerator exceeds the *held
archive comparison* for this cut, not an observational
abundance; scatter-out, scatter-in and the origin of the
held particles are not inferred from this vacuum fold.
Two quadrature controls and a separate cold-cut root
bracket document the numerical result.


## Within-bin KLOE shape check

Holding the KLOE spectrum at the official edge-center
values between the outer centers and table boundaries,
repeat the representative p/T<.1 yield with linear
interpolation and shape-preserving splines in cross
section or log cross section. Relative to fixed bins
they shift the yield by +1.90e−5, +1.44e−5 and
−7.55e−5, respectively. This is only a numerical
sensitivity to three specified curves, not a statistical
envelope over possible true within-bin shapes.
The stored constant-bin value is recovered exactly.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
D=results
python -W error warm_late_cold_exact_travel_bands_reference.py EOS.dat \
  "$D/warm-multishell-corrected-q64-fine-stage-48-2026-11-26.json.gz" KLOExSec.data \
  "$D/warm-late-cold-distance-source-2026-11-26.json" \
  "$D/warm-late-cold-pion-distance-2026-11-26.json"
python -W error warm_late_cold_kloe_bin_shape_sensitivity_reference.py EOS.dat KLOExSec.data \
  "$D/warm-late-cold-visible-kloe-recoil-2026-11-26.json"
~~~

[Exact folded distance bands](results/warm-late-cold-exact-travel-bands-2026-11-26.json);
[KLOE interpolation sensitivity](results/warm-late-cold-KLOE-bin-shapes-2026-11-26.json).
