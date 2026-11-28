# Homogeneous pressure susceptibility is not the pressure ratio

Take the terminal scalar distributions and prescribed spectators from the
[frozen bare-mass test](INDEPENDENT-EARLY-POLE-DAUGHTER-BARE-SCALAR-STRESS-2026-11-27.md).
Dilate momenta by s and number densities by s³ without changing comoving
occupations, bare mass b or coupling λ. Primes below denote log-s
derivatives. Write I = Σ n/E, A = Σ n/E³, E² = p² + m² and
m² = b² + λI/2. Implicit differentiation gives

    m′ = (λ/2)(2I+m²A)/(2m+(λ/2)mA)
       = m − 2b²/(2m+(λ/2)mA).

Thus 0 ≤ m′/m ≤ 1, whence each E′/E is in [0,1] and I′/I is in [2,3].
The logarithmic slope of the kinetic pressure lies in [4,5]
(when nonzero), and that of V = λI²/8 in [4,6]. With the stationary
energy identity ρ′ = 3(ρ+P), the homogeneous derivative cₐ² = P′/ρ′
for w=P/ρ therefore satisfies

    4w/[3(1+w)] ≤ cₐ² ≤ min(2w/(1+w), (1+9w)/[9(1+w)], 1/3).

The second and third upper bounds follow by differentiating
ρ−3P=b²I. In particular cₐ² exceeds w at every 0<w<1/3 here.
The trace form, direct derivative, fourth-order fixed-bare dilation
and nonlinear scaling ρ(s,b)=s⁴ρ(1,b/s), P(s,b)=s⁴P(1,b/s)
agree on the accepted paired fixtures (hard depth 1 × 32 cold bins,
hard 1 × 78 angular-refined bins, hard 4 × 32 cold bins).
A separate exact p=0 fixture gives cₐ² = 2w/(1+3w) even though
its single-particle kinetic pressure vanishes.

| Frozen control | w | cₐ² |
|---|---:|---:|
| Zero bare mass | 1/3 | 1/3 |
| Last local daughter aperture, earlier frozen geometry | .30640–.30684 | .31930–.31953 |
| Last whole represented hard support, earlier frozen geometry | .009098–.009101 | .016952–.016958 |
| Selected bare-mass threshold | .001 | .0016603 |

These are responses to uniform dilation of a **held scalar-only**
distribution, not scalar propagation speeds, effective Jeans scales,
cosmological transfer functions or homogeneous evolution including
further production. They neither close higher kinetic moments nor
justify substituting this derivative into a perturbation equation.

## Replay from this PR

Run from the kinetics directory. The depth-one input is only the
first native snapshot metadata used by the scalar reconstruction,
not a rerun of the original preparation. The script validates exact
initial counts and momenta against the paired held histories.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python - <<'PY'
import gzip, json
from pathlib import Path

p = Path("results")
cases = json.loads(gzip.decompress(
    (p / "early-pole-daughter-bare-scalar-stress-2026-11-27.json.gz").read_bytes()))
for name, item in zip(("h1q32", "h1angular78", "h4q32"), cases, strict=True):
    Path(f"/tmp/scalar-stress-{name}.json").write_text(json.dumps(item))
PY
a=results/early-pole-first-window-hard-depth1-cold32-scalar-stress-input-2026-11-27.json.gz
b=results/early-pole-first-window-hard-depth4-cold32-initializer-2026-11-27.json.gz
s=early_pole_daughter_scalar_adiabatic_derivative.py
python -W error "$s" "$a" results/early-pole-unsplit-cold32-{native,local-daughter}-2026-11-27.json /tmp/scalar-stress-h1q32.json
python -W error "$s" "$a" results/early-pole-unsplit-hard1-angular78-{native,local-daughter}-2026-11-27.json.gz /tmp/scalar-stress-h1angular78.json \
    --angular-mesh results/early-pole-daughter-shared-cold-angular-sampling-2026-11-27.json
python -W error "$s" "$b" results/early-pole-unsplit-hard4-cold32-{native,local-daughter}-2026-11-27.json /tmp/scalar-stress-h4q32.json
```

The [combined receipt](results/early-pole-daughter-scalar-adiabatic-2026-11-27.json.gz)
contains both arms at each threshold and aperture ceiling, the exact p=0
unit fixture and the independent numerical checks.
