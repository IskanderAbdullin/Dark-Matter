# Symmetric discrete-gradient FD collision step

For two FD baths with fixed visible temperature, this local update evolves dark number and energy density. Along the chord between the old and new moment vectors, recover the two thermal force coordinates and average them with a five-point Gauss rule. The collision source is evaluated once at the averaged forces in the discrete equation. In moment coordinates the gradient of `L = rho/T_SM - S` is `(alpha/2, beta/T_SM)`; in the exact chord integral this yields `L_new - L_old = -dt * sigma`. The executable verifies the chord pairing using a separate 12-point rule, the rate equation, positive entropy production, and the telescoped balance. The chord and averaged force are unchanged by exchanging the endpoints; this gives the time symmetry of the implicit equation, not a claim of uniqueness for arbitrary time steps.

Against a DOP853 integration of the differential moment equations, 2, 4 and 8 steps gave observed error-halving ratios of 4.001 to 4.114 across the three sampled starts. The maximum reported scaled entropy-balance error was 6.0e-15, and the maximum force difference for the independent phase-space mesh at four steps was 2.53e-10. The largest eight-step relative moment error was 7.01e-5, illustrating that the entropy check does not measure integration accuracy.

Reproduce from `kinetics`:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error fd_portal_discrete_gradient_reservoir_reference.py > /tmp/fd-discrete-gradient.json
diff -u results/fd-portal-reservoir-discrete-gradient-2026-11-28.json /tmp/fd-discrete-gradient.json
```

This checks a sampled local numerical method, not a full expanding-universe trajectory, particle abundance, or unconditional global convergence or positivity.
