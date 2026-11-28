# What common pair equilibrium would mean at the colder archive densities

The [same-moment separate fits](INDEPENDENT-WARM-PAIR-PHASE-UNCERTAINTY-2026-11-26.md)
are not at common thermal/chemical equilibrium. For a different
controlled comparison, keep only the two aggregate invariants
of the scalar-pair channel: total represented scalar-plus-light
number and total represented energy. At each actual archive solve
for a common FD/BE temperature and chemical potential. External
old cohorts enter the real Hubble normalization but are not
thermalized into this reservoir. This is a closed equilibrium
comparator, not a claim that the archived phase shapes reach it.

At the q32 7.394049 MeV checkpoint, the common fit has dark
temperature 7.08331 MeV and common chemical potential −53.5372 MeV.
The equilibrium would keep 89.264% of its archived light number
and 91.580% of its light total energy. At q64 7.754805 MeV the
analogous values are 7.37588 MeV, −55.8840 MeV, 90.599% of number
and 92.608% of light energy. Taking the original cut or a 1.07 GeV
thermal cut changes these fitted ratios by less than 4e-10 absolute.
Both count and energy are matched independently using two total-moment
quadratures, with separate finite-difference affinity checks.

At the broad hypothetical common equilibrium, use the
[number/heat tangent construction](INDEPENDENT-WARM-COMMON-PAIR-NUMBER-HEAT-2026-11-26.md)
with small affinities appropriate to the more negative chemical
potential. Normalize the rates with the radiation and total dark
energy at the real archived anchor. The slow and fast pair modes
are 0.812 H and 5.03 H for q32, and 0.923 H and 5.54 H for q64.
The very fast rates in the previous denser 17 MeV control therefore
cannot be transferred to these colder dilute reservoirs. In
particular the slow mode is not fast compared with expansion in
this comparator. It does not follow that a real kinetic solution
already near these snapshots relaxes at these tangent rates:
the actual baths need not remain separately thermal and the
prompt, elastic and external interactions are not this matrix.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error warm_pair_common_moment_equilibrium_reference.py EOS.dat \
  results/warm-multishell-trajectory-lifted24-q32-2026-11-26.json.gz \
  results/warm-pair-common-chart-q64-terminal-2026-11-26.json.gz --relax
~~~

[Matched total moments, equilibrium targets and rates](results/warm-pair-common-moment-equilibrium-2026-11-26.json).
