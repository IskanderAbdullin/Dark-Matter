# Frozen streaming at the mature warm endpoints: phase and reconstruction controls

Use the exact common-mass x=.05 G8, G9 and mixed300 occupations. Leave all unrefined cells in their identical within-cell convention. In refined warm cells compare the stored positive nodal populations with a different positive shape having the same archived number and mean energy. These two choices preserve the same archived number and mean energy but are not a bound on all possible unresolved shapes.

The test imposes a uniform scale-factor interval [1,A] in conformal time with no forcing or further collisions. At an assigned phase b integrate the angular loss 1−j0(b vbar), where vbar is the mean speed over that interval. Compare at the phase chosen once from G8, or rescale the phase independently in each chart to give the same RMS value u. A symmetric difference of crossed cases splits the latter into its fixed-phase shape and phase-rescaling parts.

In the positive smooth representative, at A=1 and u=1, the G9−G8 contrast is +4.6243e-10 at common phase, versus -1.99537e-6 when individually RMS matched. The larger negative effect is from rescaling, not from a larger fixed-phase shape difference. At A=100 and u=1, the common-phase G9−G8 contrast has opposite signs in the two allowed representations (-9.6037e-11 for the stored atoms, +5.1673e-10 for the smooth fit). At A=1000 the equal-RMS contrast also changes sign (+1.2781e-10 versus -1.6900e-9). The archived pressure ordering alone therefore does not determine these frozen phase-dependent effects.

All cases are evaluated on two integration meshes, with an independent numerical speed average for the analytic prescribed clock. These very small effects are for a dimensionless frozen test, not a transfer function, observational inference, or cosmological free-streaming distance.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
d=2026-11-28; p=results/early-pole-fine-charge2-scalar4; w=results/early-pole-warm
python -W error early_pole_endpoint_warm_three_chart_stream.py "$p-x05-state-$d.json.gz" \
 "$w-mixed300-actual-prefix-x05-$d.json" "$w-single-g9-upgrade-$d.json" \
 "$p"-{gauss8warm,gauss9warm,mixed300warm}-x05-state-$d.json.gz
```
