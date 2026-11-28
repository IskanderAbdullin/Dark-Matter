# Time-dependent chemical lag is much smaller than its instantaneous algebraic value

Start both baths at exactly the same temperature and zero chemical force. Evolve the nonlinear vacuum timelike FD model over log(a/a_start)=0.05, and separately propagate the linear forcing system along the common-entropy equilibrium-temperature chart. The immediately stationary algebraic solve is a third comparison, not a dynamical trajectory.

At initial visible T=0.55 GeV, terminal eta is -0.008319 in the nonlinear history and -0.008324 in the time-dependent linear calculation, versus -0.3010 from imposing instantaneous balance. At 0.68825 GeV the three values are -0.006349, -0.006400 and -0.03571. The linear history retains the slow mode's finite response time. These are restricted short runs from selected initial states, not a statement about prior thermal history or decoupling.

The nonlinear phase-space resolutions agree in both force coordinates within 2.43e-12. The two interpolation resolutions for the linear trajectory and the independent comoving entropy ledger pass the executable gates.

```sh
cd kinetics
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error fd_portal_expansion_forced_transient.py eos-range.dat
```

Use the same local EOS normalization as in the closed tabulated-bath benchmark. The linear result does not extend the vacuum timelike model to a complete reaction network or a relic-density estimate.
