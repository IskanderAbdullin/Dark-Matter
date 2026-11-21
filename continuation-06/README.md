# Continuation VI — energy-state geometry near Bose-halo escape

This sixth continuation examines the **single-particle density of states at the finite tidal-energy boundary** of a self-consistent lowered Bose halo. The key question from continuation V is whether an energy-space escape-flux model can support a physically meaningful finite halo, rather than arbitrarily choosing a truncated f(E).

- [Theoretical derivation, comparison, universal shallow-limit](01-energy-state-measure-and-universal-scaling.md).
- [Numerical energy-shell quadrature and comparison](code/phase_density_edge.py).
- [Independent universal shallow-cutoff scaling test](code/dos_scaling.py).
- [Static numerical results](results/phase_density_edge.json) and [small-cutoff results](results/dos_scaling.json).

**Key result:** at the same mathematical \(\eta=4\) normalization, the quadratic-lowered halo has ~19 times larger \(\mathcal G(E_t)\) than the linearly lowered halo, despite nearly equal total dimensionless mass. This is not an actual escape-rate prediction. A separate shallow-limit result is \(\mathcal G_\eta(E_t)\sim157.3291029\eta^2\).

New tests are included in the public GitHub Actions smoke workflow and the weekly extended numerical checks. All statements distinguish static orbit-energy geometry from dynamical Bose collisions, tides, stellar baryons, cosmological production, observations and nonlinear stability.
