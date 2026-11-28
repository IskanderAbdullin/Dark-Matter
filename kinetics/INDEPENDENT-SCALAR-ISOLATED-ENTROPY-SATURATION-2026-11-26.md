# Isolated scalar elastic equilibration does not force condensation

A divergent ideal Bose occupation at zero momentum
can occur at saturation without most particles
being cold, as in the [excited-state bound](INDEPENDENT-SCALAR-SATURATED-THERMAL-COLD-ORDER-2026-11-26.md).
This check instead asks whether an underoccupied
modeled scalar gas could **reach** saturation by
isolated, number-conserving kinetic equilibration
and expansion. It does not assume that scalar self
scattering is actually fast, or that the rest of
the dark sector stops exchanging energy or particles.

## Critical entropy decreases across the relativistic transition

Use $x=(E-m)/T$, $r=m/T$, and the saturated
occupation $f=(e^x-1)^{-1}$. The excited-state
radial number measure in $x$, up to temperature
factors, is
$$
 (x+r)\sqrt{x(x+2r)}\,f(x)\,dx.
$$
For $r_2>r_1$ the ratio of these density-of-states
factors is decreasing in $x$; the normalized
number measure shifts toward smaller $x$. The
entropy per occupied particle at this $x$ is
$$
 h(x)=x+(e^x-1)[-\ln(1-e^{-x})].
$$
It is increasing because $h'(x)=e^x[-\ln(1-e^{-x})]>0$.
Therefore the average entropy per excitation along
the saturation curve decreases monotonically from
$4\zeta(4)/\zeta(3)=3.60157$ to
$(5/2)\zeta(5/2)/\zeta(3/2)=1.28378$.
An ideal zero-entropy condensate only decreases the
entropy per **total** particle below the saturation
ratio. This is independent of how fast the elastic
collision process tries to establish thermal equilibrium.

## Apply to the modeled states, not to a coupled history

The two-temperature diagnostics give scalar entropy
per particle 5.039 at visible 2 MeV and 4.985 at
80 keV. Both exceed the maximum critical value.
If scalar number is conserved and the gas is
isolated afterwards, reversible expansion keeps
the entropy per particle fixed and ordinary elastic
relaxation cannot decrease it. This scalar alone
cannot equilibrate into a saturated or condensed
ideal gas under those assumptions. Slow or absent
kinetic relaxation supplies no equilibrium route
either; a driven nonthermal field is not tested.

At fixed scalar count the 80 keV starting specific
entropy would need to fall by at least 27.8% to be
compatible with any saturated equilibrium. At fixed
entropy the scalar count would need to rise by at
least $1.3841\times$; neither change is sufficient
to demonstrate condensation. If a **hypothetical**
half-abundance state inherited that same entropy per
particle, the illustrative 34.7% zero-momentum
fraction in the excited-state comparison would
require at least 52.8% entropy export at fixed count,
even granting the largest critical entropy for its
remaining excited gas. This last number is not a claim
that the modeled starting point has the half count.

Heat or number exchange with stable fermions,
new scalar production, nonequilibrium condensation
or dense interactions that change the dispersion
are outside the isolated comparison. The code checks
the microscopic entropy integral against
$(\rho+P-\mu n)/(T n)$ and the analytic relativistic
and nonrelativistic saturated limits.

~~~sh
python -W error scalar_isolated_entropy_saturation_reference.py \
 results/warm-prompt-reheating-post-2mev-elastic-2026-11-26.json.gz \
 results/scalar-saturated-thermal-cold-order-2026-11-26.json
~~~

[Entropy and saturation ledger](results/scalar-isolated-entropy-saturation-2026-11-26.json).

## How much was actually exchanged on the stored coupled path?

Reconstruct the ideal FD entropy of the stable species
alongside the scalar BE entropy, multiplying each
entropy per particle by its recorded yield. With
constant visible comoving entropy these two ratios
can be compared as comoving entropies. Between the
modeled 2 MeV and 80 keV anchors on the fine feedback
arm, the stable component gains 23.88% of its own
starting entropy, the scalar loses .816% of its own,
and the total rises 1.919%. The scalar entropy per
particle falls 1.071%, well short of the isolated
saturation gap. The two independently calibrated arms
give the same qualitative changes, and the total
rises across every stored adjacent pair.
This checks the sign and size of the sampled
approximate thermal history, not an H-theorem
for its MB collision closure. It does not constrain
an alternative cooling or number source.

~~~sh
python -W error scalar_coupled_entropy_archive_audit_reference.py \
 results/warm-prompt-reheating-post-2mev-elastic-2026-11-26.json.gz
~~~

[Separate species entropy ledger and both arms](results/scalar-coupled-entropy-audit-2026-11-26.json.gz).
