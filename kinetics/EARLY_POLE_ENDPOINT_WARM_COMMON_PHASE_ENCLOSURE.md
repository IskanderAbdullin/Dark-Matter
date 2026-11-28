# Conditional common-phase streaming signs at x=.05

Take the exact mature G8, G9 and mixed300 archived occupations, at the same native reference mass as before. Use the prior frozen test's G8 phase for u=1, once per imposed interval A, for **all** shapes and charts. Each refined warm cell can have any nonnegative population preserving its own number and mean excess energy. Cold cells remain in their archived constant convention.

At A=1 the outer conditional intervals resolve the same ordering for every admissible warm population (G9 above mixed300 above G8). The G9−G8 difference is in [+4.5105e-10,+4.7380e-10]; mixed300−G9 is in [−5.5204e-11,−3.2454e-11]. The comparison at A=2 still resolves both refinements above G8, but the mixed300−G9 sign can reverse. At A=10, 100 and 1000, constructive positive populations give **both signs for every comparison**. For A=10 the explicit G9−G8 witnesses are approximately −4.3694e-10 and +1.4454e-9. These are smaller than the absolute loss and depend on the fixed phase; the sign result does not depend on choosing the preceding smooth representative.

To construct the enclosure, parameterize each cell by normalized excitation energy. On a finite mesh the extreme mixtures consist of at most two grid energies; each one is feasible for the unmodified kernel. Bound the interpolation remainder by the chain rule: the first and second derivatives of 1−j0(beta v) are bounded by beta/2 and beta-squared/3 (in absolute value), via the angular cosine integral. The speed-derivative bounds use both ends of the prescribed mass interval [M,MA] and the momentum shell. In the energy coordinate the q derivatives are H E/q and -H² M²/q³, where H is that shell's energy width. The normal linear interpolation remainder yields the outer intervals. Independent cold quadrature and the full old nodal result are checked; the enclosures are ordinary floating-point, not directed-rounding certificates.

The interval A and phase carry no cosmological calibration. This is conditional frozen streaming, not a physical transfer function or a convergence result for the underlying interacting dynamics.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
d=2026-11-28; p=results/early-pole-fine-charge2-scalar4; w=results/early-pole-warm
python -W error -m unittest -v test_early_pole_endpoint_warm_three_chart_streaming_hull
python -W error early_pole_endpoint_warm_three_chart_streaming_hull.py "$p-x05-state-$d.json.gz" \
 "$w-mixed300-actual-prefix-x05-$d.json" "$w-single-g9-upgrade-$d.json" \
 "$w-three-chart-radiation-streaming-$d.json" "$p"-{gauss8warm,gauss9warm,mixed300warm}-x05-state-$d.json.gz
```
