# First-window finite hard-aperture refinement reverses direction at depth four

Keep the joined duration (3.3333333333333335e-7), 16 macro steps,
cold seed, packet quadrature, screening update and retained hard
endpoints fixed. Subdivide only the two hard cells per charged species
in the initial cold-source aperture. The photon support is unchanged;
the two charged supports each increase from 20 cells to 24, then 32.
Initial species numbers match within 2e-14, and each initial cold
source matches the separately evaluated static packet receipt.

| Aperture depth | Native mean | Continuous mean | Continuous below native |
|---|---:|---:|---:|
| 1 | 0.06694243427 | 0.06528878623 | 2.47025% |
| 2 | 0.06698094021 | 0.06542048071 | 2.32971% |
| 4 | 0.06680149350 | 0.06507936701 | 2.57798% |

From depth two to four the native mean decreases by about 0.27% and
the continuous mean by about 0.52%. This reverses both depth-one to
depth-two movements. The initial continuous source was already much
more stable between depths two, four and eight, so the static source
alone does not certify convergence of the joined trajectory. The
maximum recorded elastic step energy defect at depth four was
1.6e-15, below the existing 5e-15 budget; all six trajectories pass
the original closure gates.

Use `early_pole_cold_screened_split_reference.py` with
`--aperture-subdivision 4` and the same 16-step fine first-window
settings, then compare the six raw outputs using
`early_pole_first_window_aperture_depth_compare.py`. The receipt
does not establish an aperture limit, an error bar, or a physical
production rate. In particular the changed hard support also affects
the elastic updates even at fixed quadrature order. The useful next
control is a matched time-step refinement at fixed depth four.
