# A scalar source must supply both count and cold phase-space density

The [80 keV ledger](INDEPENDENT-WARM-REHEATING-ANCHOR-VIABILITY-LEDGER-2026-11-26.md)
separates the stable excess and scalar count deficit.
This complementary check asks how cold a proposed scalar
replacement could be at its **required** number density,
given a uniform upper bound on resolved occupation.
It neither derives such a bound nor proposes a source.

For one real species with $a_0=1$ and comoving momentum
$q$, require $0\le f(q)\le F$ and a fixed comoving count
$n_c=\int q^2dq\,f/(2\pi^2)$.
Among all such spectra, putting $f=F$ on
$0\le q\le q_*$ and zero elsewhere minimizes every
increasing average cost of $q$, with
$q_*=(6\pi^2n_c/F)^{1/3}$. To see this, subtract the
filled spectrum from any other spectrum. The difference
is nonpositive below $q_*$ and nonnegative above it
with zero total count. Multiplying by the difference
between the increasing cost and its value at $q_*$
makes both contributions nonnegative. The same count
argument gives the maximal fraction below a proposed
momentum cutoff, $\min[1,(q_{\rm cut}/q_*)^3]$.
Angular anisotropy cannot improve these geometric
number bounds if the cap applies in every direction.

In a fixed background, the collisionless comoving
path from the anchor to equality is strictly increasing
with $q$: the single-particle speed
$q/\sqrt{q^2+a^2m_\phi^2}$ grows with $q$ at every time.
Hence the filled sphere also minimizes the mean
comoving path. Here it is evaluated in the same
fixed background used for the independent ballistic
controls, keeping identical total background matter
when the scalar fraction changes.

| Scalar target \(\Omega_\phi h^2\) | Minimum mean path for \(F=1\) | Necessary cap to allow mean path \(\le .1\) Mpc |
|---:|---:|---:|
| .060 | 3.44 Mpc | \(1.97\times10^5\) |
| .120 | 4.16 Mpc | \(3.95\times10^5\) |

For the half target and $F=1$, at most .003873 of
the particles lie below $q/T_{\gamma,0}=.1$.
Putting at least half the count below that cut
requires $F\ge129.09$. Putting half below .01
requires $F\ge1.29\times10^5$. These are illustrative
coldness targets, **not** limits inferred from data.

The stored modeled 80 keV scalar shape has fitted
rest occupation .512 and mean reference ballistic
path 5.68 Mpc, but its conditional rest-density proxy
is only .0149. The target calculation uses the required
count, not a relabeling of that thermal fit as a new
source. A large allowed occupation alone does not
populate the required momenta, and the ideal filled
sphere need not be reachable by the actual dynamics.
The comparison supplies a useful target for candidate
production spectra, not a Boltzmann transfer function,
coherent-condensate result or observational exclusion.

~~~sh
python -W error scalar_cold_occupation_phase_capacity_reference.py \
 EOS.late.dat \
 results/warm-prompt-reheating-post-2mev-elastic-2026-11-26.json.gz
~~~

[Geometric and ballistic controls](results/scalar-cold-occupation-phase-capacity-2026-11-26.json).
