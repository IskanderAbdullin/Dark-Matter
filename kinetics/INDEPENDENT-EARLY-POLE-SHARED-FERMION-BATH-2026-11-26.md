# Equal charges do not selectively thermalize the decay parent

Consider a neutral real scalar with off-diagonal Yukawa coupling
to two Dirac species carrying the same Abelian charge. Equal
charges are required for this neutral-scalar Yukawa term to
be gauge invariant. If the Abelian bath is populated, gauge
momentum exchange can act on **both** hard fermions. It does
not automatically justify holding only the parent FD.

For a static geometry test, add the common hard asymptotic
squared mass g²T²/4 to both fermions. The aperture center is
proportional to M_parent²−M_daughter²+m_scalar² and is unchanged
by the common shift. Its width changes through the parent
mass in the discriminant. An independent vacuum-massless HTL
normal-branch solve checks the hard asymptotic shift. If a
hypothetical species-dependent split is added, the center moves.
At x=.0075, a 1% difference of the common squared shift and
g=.30 closes this particular aperture in the static toy,
whereas the common shift alone does not. The split is not an
assumption of the equal-charge model.

The soft angular energy separation can be much smaller than a
single charged-line damping scale even when this common-mass
aperture still survives. As an explicitly parametric check,
g=.001, .003 and .01 give a width-proxy to angular-splitting ratio
near 0.067, 0.51 and 4.47 in the stored .06 cases. This is not
a decay width, a gauge-invariant two-line kernel, or a
Lorentzian ansatz for the real spectral function. For the
hard thermal and transport approximations see
[Arnold, Moore and Yaffe](https://arxiv.org/abs/hep-ph/0209353);
for the non-exponential single-fermion damping caveat see
[Blaizot and Iancu](https://arxiv.org/abs/hep-ph/9607303).

A separate live comparator imposes kinetic relaxation of
both fermions while preserving each particle number and their
combined energy. The thermal target has a common temperature
and two chemical potentials. At x=.052 and rate 100 per
reference expansion time, third-grid scalar yield is .000431450
under joint kinetic relaxation, vs .000420711 for parent-only
and .000414781 without it. The other grid gives .000439328,
.000425758, and .000419545. These surrogate rates are not
derived gauge rates; this check does not add the common gauge
pole to the dynamic calculation.

For two FD spectra at equal temperature, the stimulated part
has positive sign precisely if μ_parent−μ_daughter > E_scalar,
independent of parent energy in the allowed window. Equal
chemical potentials give a negative stimulated contribution,
not zero spontaneous production. Holding the shared target
temperature and parent chemical potential fixed, filling the
thermal daughter to the zero-stimulation threshold would
require about 20 times its snapshot number at x=.05, or
about 2.5–2.6 times by x=.10. That is a population test,
not a time to equilibrate. Gauge pair sources, radiation
backreaction and the consistent in-medium decay kernel are
not included.

~~~sh
python -W error early_pole_shared_gauge_bath_kinematic_reference.py results/early-pole-spectrum-third-m06-2026-11-26.json.gz
python -W error early_pole_joint_fermion_thermal_reference.py results/early-pole-spectrum-third-m06-2026-11-26.json.gz --kind joint
python -W error early_pole_fermion_population_gate_reference.py results/early-pole-spectrum-third-m06-2026-11-26.json.gz
~~~

[Dispersion and width-proxy checks](results/early-pole-shared-gauge-bath-kinematics-2026-11-26.json.gz);
[joint live and chemical population checks](results/early-pole-shared-fermion-bath-controls-2026-11-26.json).
