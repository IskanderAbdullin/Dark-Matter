# The fifth warm Gauss order raises pressure after the fourth lowered it

An independent G5 history evolves the same cold cells and D2 charged geometry, and five positive volume-Gauss nodes in each of the 22 complete warm cells. The state has 233 scalar, 448 daughter and 488 parent nodes. All five reports and endpoint states are archived. Their paired time-steps, continuity and positive weights pass the restart chain receipt, with SHA-256 links binding the actual previous fine states.

At common held reference mass, G5 minus G4 relative pressure changes are +0.112041%, +0.113468%, +0.110519%, and +0.112598% at x=.01/.02/.04/.05. At .05 this nearly offsets the preceding G4/G3 decrease: G5/G3 remains −0.012933%, and G5/geometric-W4 −0.503181%. The G5/G4 number change is only +3.64e-9 in relative fraction. The final pressure difference remains [+0.112590%,+0.112606%] under the archived warm count-and-energy conditional reconstruction. This is evidence against inferring monotone warm spatial convergence from G2 through G5, not a measurement of an asymptotic order.

The independently checked cell decomposition at .05 gives +0.118903 percentage points from changed cell populations, counteracted by −0.006305 points from pressure per particle. The two cells selected previously by the frozen split diagnostic contribute about 97.8% of the signed net G5/G4 increase. This comparison describes the two completed trajectories; it neither evolves the 219-node targeted alternative nor identifies a dynamical cause.

For the chosen matched-energy positive reconstruction, G5/G4 changes at formal probe mass 100 are +0.29010% in pressure and +0.81735% in pseudo-pressure. These are frozen stress probes only: a sign at a heavy formal mass must not be inferred from number and mean energy alone. The physical history, actual expansion and full perturbation transfer are not represented by these probes.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; d=2026-11-28
python -W error early_pole_endpoint_gauss_warm_chain_receipt.py --order 5 \
 --reports "$p"-gauss5warm-{start,x01,x02,x04,x05}-$d.json.gz \
 --states "$p"-gauss5warm-{start,x01,x02,x04,x05}-state-$d.json.gz
python -W error early_pole_endpoint_warm_gauss_order_contrast.py --levels 4 5 \
 --lower "$p"-gauss4warm-comparison-x{01,02,04,05}-$d.json \
 --higher "$p"-gauss5warm-comparison-x{01,02,04,05}-$d.json \
 --stress results/early-pole-frozen-stress-sweep-w124-gauss{4,5}-$d.json
for x in 01 02 04 05; do
 python -W error early_pole_endpoint_warm_gauss_transition_attribution.py --levels 4 5 \
  "$p"-x$x-state-$d.json.gz "$p"-gauss{4,5}warm-x$x-state-$d.json.gz
done
```
