# Accuracy gate for finite native momentum sampling of the archived PSD

The accurate [native homogeneous CLASS calculation](INDEPENDENT-CLASS-PSD-BACKGROUND-AUDIT-2026-11-26.md)
used 4160 momentum nodes. A faster perturbation pilot need
not inherit that accuracy. Evaluate the actual exported stable
f₀(q) files on the fixed finite-q scheme-3 grids used in the
pilot (qmax=20). Check physical normalization, normalized
q-moments, and the dilation identity
∫q^(3+j) f₀'(q)dq = −(3+j)∫q^(2+j) f₀(q)dq
when the boundary term is negligible. The derivative is taken
from a separate cubic spline, not from CLASS's finite difference;
this diagnoses quadrature aliasing even with a smooth derivative.

| Representative grid and export | Count error | Largest normalized moment error | Largest dilation identity residual |
| --- | ---: | ---: | ---: |
| 64, standard log-smoothed export | 0.678% | 0.685% | 0.774% |
| 128, fine log-smoothed export | 0.209% | 0.210% | 0.541% |
| 1040, fine | 0.0211% | 0.0210% | 0.103% |
| 4160, fine | <5e-10 relative | <5e-10 relative | 7.3e-7 relative |

Changing the node count can produce non-monotonic errors because
the samples move relative to shell edges. Rescaling Ω to a common
value in CLASS may remove a count mismatch, not automatically the
dilation/perturbation mismatch. A completed coarse pilot is
therefore not a precision transfer result unless these numerical
effects and the integration convergence have been checked.
Tested [official manual quadrature convention](https://github.com/lesgourg/class_public/blob/v3.3.4/tools/quadrature.c).

```sh
OPENBLAS_NUM_THREADS=1 python -W error archived_class_manual_momentum_gate_reference.py class_psd_base/stable.dat class_psd_fine/stable.dat
```
