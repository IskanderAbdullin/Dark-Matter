# Local equation of state of fragmented classical scalar (2026-11-25)

Extend vacuum-seeded lattice runs (PHI-BOSE-QUANTUM-SEED-2026-11-25.md) by 16 periods, sampling every 8 substeps. Define K=⟨u_z²⟩/2, G=⟨|∇u|²⟩/2, M=a⟨u²⟩/2, V=b⟨u⁴⟩/4. Isotropic mean pressure P=K−G/3−M−V and E=K+G+M+V. For periodic box exact continuum identity d⟨u u_z⟩/dz=2K−2G−2M−4V; thus E−3P=2M−d⟨u u_z⟩/dz. When the time-average boundary term vanishes, w=(2G/3+V)/E∈[0,1/3], independent of whether fragmented. In pure quartic massless theory the virial-averaged w=1/3 regardless of nonlinear redistribution. This is not an exact relation at an arbitrary instant or a cosmological evolution solution.

Actual averaged results (K,G,M,V normalized to time-averaged E):

| A | N | z interval | w | w_hom exact | boundary-corrected w | K | G | M | V | relative virial residual |
|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
|2.72085|48|767.829–884.831|.302287|.296782|.301898|.614197|.110840|.0465199|.228443|−.0012669|
|120|64|519.137–637.797|.334009|.333312|.332939|.604353|.187955|.0000237364|.207668|−.0011362|

The raw quartic w slightly exceeds 1/3 owing to finite-time boundary and numerical residual; not a super-radiation equation of state. Finite-time corrections + ~10⁻³ virial residual limit interpretation. Mixed-field w increases by ~0.005 in this window relative to the exact homogeneous action solution. Max sampled relative energy error .000891 and .000920. Average k>2 gradient share .3171 and .5764; Hartree nonzero-mode action-weighted rms k=1.5068 and 1.7943, proxy only because interaction significant. Peak |u|=1.6009 and 1.6054 (sampled), far below coherent fermion zero-crossing at T20 for these examples; this does not bound all times, rare configurations, or earlier epochs. Pressure local at fixed a and finite cutoff; no final density or f(q), no proof of BEC. Source 189c513; reproduce with python phi_bose_lattice_eos.py CHECKPOINT --periods 16 --output NEW_CHECKPOINT, using previous quantum-seeded runs. No remote CI checked.

Next physics: expansion with a controlled closure or classical field at accelerated expansion as a separate convergence/adiabaticity diagnostic, then quantum scattering/matching when occupations become small. Never infer count conservation from virial alone.
