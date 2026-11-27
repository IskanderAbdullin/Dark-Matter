# Use a heavy-energy dose to bound prompt Bose growth and a light count to bound inverse loss

The [boosted-parent limits](INDEPENDENT-WARM-LATE-COLD-PARENT-RATE-2026-11-26.md)
use an archived maximum cold Bose occupation for the
instantaneous band-integrated rate. To avoid silently
assuming that ceiling persists, fix the background and
track one scalar momentum. In the isotropic two-body
operator, the pointwise gain divided by (1+f_s) cannot
exceed 4 Gamma M²/(2 p* p E_s) times the unweighted heavy
occupation integral between its kinematic endpoints. For
E_H≥L(p) the weight in the heavy energy integral is at
least (2/pi²)L² sqrt(L²−M²). Divide any allowed eligible
heavy energy density by this lower weight to bound the
unweighted occupation integral. Integrating the resulting
differential inequality controls ln(1+f_s) without assuming
that f_s remains at a fixed value.

A declared finite parent energy cap creates an exact
inaccessible interval at the bottom of the band. Above
that floor, bound each momentum panel using the minimum
p E_s from its left end and the minimum heavy weight
from its right end. The minimum allowed parent energy
decreases with p below the rest-frame daughter momentum,
as is clear in the rapidity form, so the combined
endpoint inequality holds throughout the panel. Fourfold
refinement checks the bound without treating a sampled
maximum as an envelope.

At .1T with equal incoming daughter caps of 7.5 GeV
and a parent known to stay below their summed energy
(15.000027 GeV), the inaccessible scalar interval extends
to p=1.60001e−4 GeV. An upper over the accessible band
is 1.1084e5 GeV⁻⁴ times Gamma times the allowed heavy
energy density. Written in units of the 7.394 MeV archive,
this is 1.4920e7 times the integrated eligible energy
fraction times Hdt. Increasing the Bose factor by 1%
through this channel therefore needs a dose above
6.66e−10 in those normalized units. The frozen 3.3 and
6.5 GeV reactivations have eligible energy fractions of
only 1.63e−19 and 8.41e−20. If those were **sustained
upper bounds** during Hdt=1e−7, the corresponding maximum
Bose exponents would be 2.44e−19 and 1.26e−19. They are
benchmarks, not verified through-time caps.

## A separate loss envelope and exact pointwise samples

The inverse loss divided by f_s depends on the unweighted
light occupation integral above the required recoil.
Monotonicity of the light number weight converts the
eligible light count into an upper bound on that integral.
Use the recorded birth floors and means without specifying
their shape; no untracked light or newly populated hard
phase is allowed. At the 7.5 GeV daughter cap, if those
archived light-count bounds persist over Hdt=1e−7, the
fractional inverse loss at any held p/T≤.1 is below
8.02e−13 (q32) and 3.21e−13 (q64). This conditional bound
does not assume a value of the scalar occupation.

Independent exact pointwise folds of the hard 3.3 and
6.5 GeV reactivations show why a negative band-integrated
rate does not determine every sub-band's sign. At 3.3,
the frozen p/T=.095 sample has gain and zero loss, but
p/T=.1 has a negative current about −2.41e−6 times its
local occupation per Hubble time. At 6.5, p/T=.045 has
a very small positive source, while .05 is negative
about −8.63e−7 per Hubble time in the same local units.
The majority of very soft samples are zero or have
negligible positive gains; no interpolation or persistence
of these signs is assumed. The signed local result is from
both occupation products, not the dose bound.

The bounds exclude redshift, energy injected by other
channels, or any evolution that violates the stipulated
heavy or light budgets. They are not a full late spectrum
or proof about the growth of a condensate.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
a=results/warm-corrected-q32-time48-step48-2026-11-26.json.gz
b=results/warm-multishell-corrected-q64-fine-stage-48-2026-11-26.json.gz
python -W error warm_late_cold_pointwise_bose_dose_reference.py EOS.dat "$b" \
 results/warm-late-cold-scalar-reactivated-prompt-2026-11-26.json
python -W error warm_late_cold_pointwise_prompt_spectrum_reference.py EOS.dat "$b"
python -W error warm_late_cold_light_loss_bound_reference.py EOS.dat "$a" "$b"
~~~

[Heavy dose envelope](results/warm-late-cold-pointwise-bose-dose-2026-11-26.json);
[inverse light-count envelope](results/warm-late-cold-light-loss-bound-2026-11-26.json);
[pointwise signed spectra](results/warm-late-cold-pointwise-prompt-spectrum-2026-11-26.json).
