# A heavy-fermion occupation barrier from light occupations and source budgets

The [signed warm lower comparison](INDEPENDENT-WARM-SIGNED-FERMION-TIME-COMPARISON-2026-11-26.md)
requires a per-spin cap for both the light and heavy
fermions. This step asks whether the cap for the heavy
species can follow from the one for light. It does not
establish the light cap, or bound elastic scattering
that redistributes the heavy fermions.

## Pointwise decay and inverse decay

For a heavy mode with occupation (f_H), the signed
vacuum decay/inverse term divided by (Gamma M/E_H)
has local angular integrand
`(f_L−f_H) f_φ − f_H(1−f_L)`.
If (f_L≤δ) in all daughter directions, then exactly
at (f_H=δ) this is at most `−δ(1−δ)`, regardless of
scalar occupancy. Free redshift does not break a
uniform cap. Therefore a heavy occupation initially
below δ stays below if additional positive gain is
less than `δ(1−δ) Γ M/E_H` at every mode. This is a
pointwise barrier, not a comparison of total rates or
an assumption that the heavy species is thermal.

For vacuum visible pair injection, bound the incoming
thermal weights by the Boltzmann weight of the two
outgoing heavy particles. Apply detailed balance to
bound injection into the tagged heavy mode by
`exp(−E_H/T) n_H,MB K_vis`, with the heavy partner's
single-charge two-spin Boltzmann density. The global
vector coefficient uses the exact maximum of the
Breit–Wigner shape (s/[(s−m_V²)²+(m_VΓ_V)²])
and the supported charged-pion allowance as a separate
coefficient. At (T<M), (E_H exp(−E_H/T)) is
decreasing from the rest point. The ratio to the
decay margin is therefore maximized at the rest
point and highest included temperature.

## Bound internal pair production into the same mode

For a tagged heavy momentum (p), eliminate its
unobserved outgoing partner and the incoming angle
using the energy delta. The remaining angular factor
averages as (1/max(k,p)) for an isotropic incident
mode with momentum k. Use 32 as a deliberate
spin/label and light-charge overcount, so the
estimate does not rely on a favorable symmetry
factor.

- **Two incident scalars.** When both incoming energies
  exceed M/2, their required radial integrals are
  controlled by scalar energy. If one is softer, its
  partner needs energy at least (M²/E_soft).
  The hard radial tail brings a factor (E_soft²/M⁴).
  If (p≥M/4), (E_soft/max(k,p)≤2); otherwise
  each incoming energy obeys (E≥(E_H−p)/2>3M/8)
  by four-momentum conservation, so the same ratio is
  also at most 2. Thus the soft expression is controlled
  by scalar energy squared even for arbitrary scalar
  occupation. Use the earlier full dark-energy upper
  instead of requiring a known scalar share.
- **Light fermion pairs.** For the crossed massive
  transition, (t,u≤−(M²−m²)²/s); it follows from the
  positive product of the two CM propagator endpoints.
  Bound the scalar-exchange matrix element by a
  constant plus a term quadratic in s, and then
  (s≤4E_1E_2). The two remaining radial integrals
  can be upper-bounded by (J/m) and
  (mJ+π²n_charge), with (J=∫p f_L dp).
  At fixed charge density the maximum J under the
  occupancy cap is the filled low-momentum ball.
  Enlarge its density to total dark energy divided by m.
  This no longer uses the recorded light momentum shape.

For q32 at δ=1/2, the known positive visible and internal
production allowances together are at most 2.483×10⁻¹²
of the modewise decay margin. The result also includes a
separate diagnostic on five frozen retained daughter
archives. At 10.442 MeV the larger sampled stationary
parent occupation for decay/inverse alone is below
(5.93×10^{-12}); the smallest response relative to
nominal visible H on the probed grid is above (10^{12}).
Neither sampled value is used as a cap or a bound for an
unrepresented daughter distribution.

Thus the specified **number-changing heavy production**
cannot break a heavy cap if the light cap and the energy
allowance hold and there is no unbounded additional
driving. Heavy elastic momentum redistribution, other
vertices and the time-dependent light occupations must
still be checked. A fast decay is not used as a substitute
for a complete signed kinetic evolution.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_heavy_decay_occupation_barrier_reference.py EOS.dat \
 "$D"/warm-signed-fermion-with-charged-pions-2026-11-26.json \
 "$D"/warm-early-q{16,32}-time32-checkpoint-2026-11-26.json.gz \
 "$D"/warm-coupled-pair-b{35-48,10-42,30-32}-shape-2026-11-26.json > /tmp/heavy-visible.json
python -W error warm_heavy_internal_source_barrier_reference.py \
 "$D"/warm-signed-fermion-with-charged-pions-2026-11-26.json /tmp/heavy-visible.json
~~~
