# Local equilibrium conditioning depends on lab-energy spread even within the invariant pole window

At common temperature and zero chemical force, construct the positive
two-bath moment matrix for event number and lab heat. Its normalized
conditioning depends on the variance of lab event energy. Partition
the independent positive quadratures in invariant mass, using ±30
and ±100 widths around the vector pole. Within and between region
contributions exactly reconstruct the original variance.

At T=.6882487043 GeV, the spread of lab energies within the ±100-width
pole window contributes about 79.0% of the total variance. The actual
normalized condition number is 531.5. Keeping that region's number
weight but replacing all of its lab event energies by the region's
own mean raises it to about 2523. Replacing them by the invariant
pole energy instead gives about 3182. The outside regions remain
unchanged in both controls. At .30 GeV the in-window lab variance
contributes less than 0.007% and the same comparison changes the
condition number only slightly.

These are deliberate conditioning counterfactuals. They show why
an invariant-energy pole location cannot be substituted for the
distribution of transferred lab energies. They do not establish
a physical abundance, freeze-out time, or observational impact.

```sh
cd kinetics
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error fd_portal_equilibrium_variance_conditioning_reference.py
```
