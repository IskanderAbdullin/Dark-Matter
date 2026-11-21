# Two-component numerical experiment (not observed galaxy fits)

Recomputed the spherical Poisson–Bose ODE numerically with **unchanged** \(V_0\), \(R\), central noncondensed density normalization \(x_n(0)=1\), and central Bose-star ratio \(\alpha=12.8\), dispersion ratio \(\gamma=12.8\). The comparison deliberately only changes \(\alpha,\gamma\) and solves for the new total potential. These high values are a **phenomenological example**, not a preferred particle-mass/coupling prediction.

## Results at \(\nu=0.19\)

| \(\xi=r/R\) | Pure \(u=V/V_0\) | With Bose stars | Ratio |
|---:|---:|---:|---:|
| 0.1 | 0.05770 | 0.19821 | 3.4350 |
| 0.5 | 0.28487 | 0.48406 | 1.6992 |
| 1.0 | 0.54865 | 0.60688 | 1.1061 |
| 3.0 | 1.23104 | 1.11613 | 0.9067 |
| 6.0 | 1.52907 | 1.45653 | 0.9526 |
| 10.0 | 1.57373 | 1.55713 | 0.9895 |

The model increases central rotation speed strongly, but produces slightly **lower** speeds at intermediate dimensionless radius. This is a self-consistent redistribution under the combined potential, not a negative Bose-star mass.

## Results at exact critical degeneracy, \(\nu=\zeta(3/2)\)

| \(\xi\) | Pure \(u\) | With Bose stars | Ratio |
|---:|---:|---:|---:|
| 0.1 | 0.05654 | 0.19720 | 3.4880 |
| 0.5 | 0.26070 | 0.46533 | 1.7849 |
| 1.0 | 0.47290 | 0.53774 | 1.1371 |
| 3.0 | 1.00142 | 0.88332 | 0.8821 |
| 6.0 | 1.32467 | 1.23775 | 0.9344 |
| 10.0 | 1.46422 | 1.43389 | 0.9793 |

At the critical central limit with \(\alpha=12.8\), the independent local density-slope calculation gives analytical \(-d x_n/d\xi=2.0579426748\) and numerical \(2.0579299027\) at \(\xi=10^{-5}\). The difference is compatible with finite-radius truncation at this order.

**This does not validate formation of Bose stars, a physical value of \(\alpha,\gamma\), or the observed galaxy-rotation data**. Run the source-provided solver and this comparison locally; cosmological and baryonic tests remain outstanding.
