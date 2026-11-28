# The coarser S4 panels cannot realize thirteen S8 reference means without splitting

Project the archived S4 and S8 scalar populations onto the same 31 S4 conditioning cells in the retained cold band. In each cell combine the two candidate number and excess totals before forming the mean, always using the held S4 reference mass to define excess. Both chart kernels remain independently fixed. At the pure S4 reference both chart panel families can match all 31 number/mean constraints. At each nonzero tested S8 reference weight (0.25, 0.5, 0.75, 1), the original S4 panel basis fails 13 of the requested means, while S8 passes all 31. The largest pure-S8 miss on S4 is 1.96% of the requested cell mean. This is an actual obstruction of the old finite positive panel basis, not a statement about arbitrary distributions in those physical cells.

The highest-number shared contiguous coarsening that supports both endpoint vectors on both bases retains 22 parents. Without forcing those merged boundaries on the response display, the three-field joint loss in fixed original S4 channel scales ranges from 0.000111924 to 0.000112803 through the reference blends. On all 31 original S4 parents with the pure S4 reference the comparable loss is 0.000133596. Both q meshes agree to 2.03e-13 on the reported optima. For this newer pair, the LP proposes candidate simplex weights and signed support weights; the published bracket comes from projecting them and computing independent lower/upper values directly. The largest such bracket width is 1.12e-9 even if nominal solver dual diagnostics are less tight. Ordinary floating-point, not directed rounding.

A direct geometry check produces a plan adding one edge in each of 13 old S4 panels, preserving exactly the retained 31 display bins common to S4/S8. The proposed 84-panel support can represent both reference endpoints and hence their blends on all 31 conditioning cells with chart-specific positive weights. It calls for 26 fresh panel evaluations per mesh (replace the original panel by its two children). This note includes **no** new reaction coefficients for those children; none may be filled in from the parent integral.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/fine-charge2-scalar; s=results/early-pole-fine-charge2-scalar
base=results/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json
python -W error early_pole_endpoint_cold_shared_reference_scan.py --levels 4 8 \
 --left "$p"4-endpoint-continuum-cold-{1,2}.jsonl.gz --right "$p"8-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --states "$s"{4,8}-x05-state-$d.json.gz --archive "$base" >/tmp/s4-s8-reference.json
python -W error early_pole_endpoint_cold_shared_reference_reintegration_plan.py \
 --left "$p"4-endpoint-continuum-cold-2.jsonl.gz --right "$p"8-endpoint-continuum-cold-2.jsonl.gz \
 --scan /tmp/s4-s8-reference.json
```
