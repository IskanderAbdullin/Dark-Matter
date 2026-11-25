# Direct daughter-energy finite-volume kernel — 2026-11-25

## Scope and derivation

The previous [fixed-bath benchmark](VERIFIED-QUANTUM-BURST.md) had a 21.8% spread in mean momentum; the parent/angular quadrature and barycentric remapping may smear the stimulated peak. Same idealised reaction and masses as before (massless daughters, g1=g2=4, gphi=1, M1/T=.1, external thermal parent, no FRW). See Yin, JHEP 05 (2023) 180, https://doi.org/10.1007/JHEP05(2023)180.

With p the mother momentum, E=sqrt(p²+M²), q=(E+p mu)/2, r=E-q, the rate element g1 M/(2 pi²) (p²/E) dp dmu/2 becomes g1 M/(2 pi²) dE dq = g1 M/(2 pi²) dq dr, restricted to q r >= M²/4. This exactly removes angle and parent quadrature. Energies are in T and time in inverse rest width.

On the pair of daughter momentum cells define R_ij as the allowed area times g1 M/(2 pi²), and A_ij as the integral of the same measure weighted by the parent Fermi factor 1/(exp(q+r)+1). The cell rate is A(1-f2)(1+fphi)-(R-A)f2 fphi. Both number equations use the same rate; exact number difference and invariant boundaries follow. Split the integration at the threshold intersections; integrate the Fermi factor analytically in r, with Gauss-Legendre in q. The independent vacuum gain integral checks normalization.

Two closures are compared: integrated A (exact for constant cell occupations) and balanced A=R f_FD(qbar_i+qbar_j), where qbar is the exact shell-volume mean. The latter exactly preserves discrete thermal detailed balance; its vacuum gain bias vanishes under refinement. Neither claims exact per-event physical energy projection; the bath is external and the cell energies approximate true energies. This is a separate solver, not a modification of the earlier benchmark.

## Executed results

| scheme | cells | n/T³ | mean q | rms q | hot fraction q>.1 |
|---|---:|---:|---:|---:|---:|
| integrated | 100 | .367817 | .00354628 | .0506400 | .00100655 |
| integrated | 200 | .365019 | .00299373 | .0504627 | .00101386 |
| integrated | 400 | .364475 | .00275551 | .0504066 | .00101545 |
| integrated | 800 | .364371 | .00269054 | .0503904 | .00101576 |
| balanced | 200 | .363321 | .00266660 | .0503977 | .00101579 |
| balanced | 400 | .364086 | .00266784 | .0503880 | .00101584 |
| balanced | 800 | .364277 | .00266806 | .0503855 | .00101585 |
| balanced, extended cutoff + tighter time | 400 | .363992 | .00266766 | .0503892 | .00101585 |

Mean q changes by ~2.41% for integrated 400→800, ~0.0085% for balanced 400→800, and ~0.84% between the two 800-cell variants. Compared to earlier mean q=.00387 on its finest grid, the new mean is substantially smaller. The balanced closure has a 0.0177% vacuum rate bias at 800; its normalized equilibrium residual is <8e-16. Integrated residual reduces about fourfold on each doubling. Maximum daughter number mismatch <7e-14. All final states nonnegative and fermionic occupations <1. The exact shell average of q² is used, rather than squaring the average q. RMS remains nearly 19 times the mean: a cold number distribution is not equivalent to negligible hot tail.

The [code](daughter_energy_fv.py) and [selected actual local outputs/commands](results/daughter-energy-fv-2026-11-25.json) are reproducible with Python 3.13.5 / NumPy 2.3.5 / SciPy 1.17.0 (single BLAS thread). Eight runs succeeded. Convergence of number, mean, rms and threshold fraction at this single time is substantially improved; pointwise PSD, epoch dependence, cutoff/step independence as separate effects, physical IR regulation, and cosmological applicability are **not** established. No CLASS run.

## Checkpoint

Use the balanced closure with independently checked direct phase-space area and exact discrete detailed balance for a first radiation-dominated FRW fixed-bath model: q=p/T with T proportional to 1/a, x=M/T, H=H(M)/x²; d/dx = (Gamma/H(M)) x times the collision kernel (which itself depends on x). Quantify freeze-out, abundance and residual hot tail, including a reheating-temperature sensitivity; do not treat the external bath as a closed energy budget.
