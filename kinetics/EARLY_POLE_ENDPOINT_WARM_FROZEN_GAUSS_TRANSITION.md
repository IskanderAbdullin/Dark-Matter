# At fixed occupations, a fifth Gauss point reverses the warm source correction again

Lift the exact same piecewise-constant S4 occupations onto each positive volume-Gauss order without changing any cell's count or any charged node. Use the same scalar mass, collision clock, and response, evaluated on paired internal radial meshes. Compare with the saved warm Volume64 reference; its unsigned Volume32/64 discrepancy relative to total net source is 1.20e-5 or less at both endpoints.

| Frozen endpoint | G2→G3 net-source change | G3→G4 net-source change | G4→G5 net-source change |
| ---: | ---: | ---: | ---: |
| .01 | +0.09238% | −0.12395% | +0.11220% |
| .05 | −0.03244% | −0.14335% | +0.14804% |

The percentages divide each difference by the same Volume64 net warm pressure source at that endpoint. The signed G5-minus-reference discrepancy is +0.04575% and +0.05452%; corresponding sums of absolute cell deviations are 0.04767% and 0.05502%. The G5 rule reproduces pressure-volume polynomial moments through degree nine to 3.78e-14 on the actual grid and keeps all nodes and weights positive. The largest G3→G4 frozen contribution comes from q=[3.96368,5.72884], also the highest-impact parent in the independent evolved pressure attribution. At x=.05 this parent alone contributes −0.11511% of the total reference instantaneous source. The G4→G5 effect moves its largest contribution to q=[5.72884,8.28009]. In the frozen tests nearly all three adjacent pressure-source differences are forward, not inverse.

This is an instantaneous S4-occupation control. It supports proceeding to an independent G5 trajectory without assuming G5's evolved correction has the same sign or that the pressure has reached a continuum limit. In particular, even the G2→G3 signed frozen response at .05 disagrees in sign with the accumulated evolved pressure difference at that epoch.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; d=2026-11-28
for i in 01 05; do
 python -W error early_pole_endpoint_warm_frozen_gauss_transition.py \
   "$p-x$i-state-$d.json.gz" "results/early-pole-endpoint-warm-source-cells-x$i-$d.json"
done
```
