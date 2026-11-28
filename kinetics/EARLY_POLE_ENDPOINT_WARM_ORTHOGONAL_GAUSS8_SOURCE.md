# A better-conditioned positive volume rule for G7 and G8

On each finite scalar shell use the affine coordinate u in [0,1], with volume density proportional to (a+(b−a)u)². Form the Gram matrix and the matrix for multiplication by u in an orthonormal shifted Legendre basis. Ordinary Legendre integration of order p+2 is exact for these polynomial integrands at node order p. The generalized symmetric eigenvalues are the Gaussian nodes; squared projections of the constant vector are positive weights. Check every positive interior node and all normalized monomial moments through degree 2p−1 independently. The original G2 through G6 execution paths remain unchanged.

On the actual fine warm shell chart, p=7 and p=8 give 277 and 299 scalar nodes, with maximum moment residuals about 5.89e−15 and 1.44e−14. On warm parent shells the worst 2-norm condition numbers of the monomial Gram matrices at p=7,8 are about 5.93e8 and 1.92e10; in the Legendre basis the corresponding maxima are about 2.01 and 2.03. At zero lower edge, separate Gauss–Jacobi roots agree to 1.9e−14 in normalized node location and 2.0e−15 in normalized weight. The wide-to-narrow positive shell tests and independent lower-order reproduction are in the unit suite.

Both new held warm source rules were directly evaluated at the stored x=.01 and .05 S4 snapshots. Differences are divided by that snapshot's admitted Volume64 net source. The known G2–G6 subresults reproduce exactly.

| x | G6 native warm-cell L1 | G7 L1 | G8 L1 | G7−G6 net source difference | G8−G7 net source difference |
| ---: | ---: | ---: | ---: | ---: | ---: |
| .01 | 0.0422235% | 0.0263033% | 0.0244703% | +0.0592024% | −0.0478693% |
| .05 | 0.0529980% | 0.0342180% | 0.0286017% | +0.0527999% | −0.0489210% |

The evaluated L1 relative to Volume64 decreases through G8 at both stored snapshots, but the alternating net change persists. The .05 signed residual relative to Volume64 is +0.0213334% for G7 and −0.0275876% for G8. These are instantaneous frozen-source tests. They do not guarantee monotone evolved pressure, spatial or temporal convergence, or physical transfer.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results
for n in 01 05; do
 python -W error early_pole_endpoint_warm_frozen_gauss_transition.py \
  "$p/early-pole-fine-charge2-scalar4-x$n-state-$d.json.gz" \
  "$p/early-pole-endpoint-warm-source-cells-x$n-$d.json" --through 8
done
```
