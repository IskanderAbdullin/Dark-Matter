# Finite-displacement time refinement with the exponential FD filter

Compare the filtered and unfiltered discrete-gradient updates against a tight DOP853 integration of the **unmodified** fixed-reservoir moment ODE. Use the same three original positive/off-equilibrium thermal starts, the same horizon and 2/4/8/16 steps. The filtered observed error-halving ratios range from 4.0002 to 4.1295. Both collision phase-space orders reproduce the eight-step forces to within 2.50e-10.

At T=0.6882487043 and eight steps the filtered maximum moment error is 1.4382e-5 of its initial component scale, versus 7.0076e-5 unfiltered. At T=0.05 the effect goes the other way (6.4160e-6 versus 5.4773e-6). This illustrates that the filter is not uniformly more accurate away from equilibrium. At the eight-step resolution the ratio of modified to physical entropy rates was 0.99888 to 0.99989 across the tested steps. The discrete balance is separately checked; it does not make those rates equal.

```sh
cd kinetics
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error fd_portal_reservoir_filter_finite_order.py
```

This comparison concerns a finite test interval at fixed visible temperature, not a full kinetic or cosmological history.
