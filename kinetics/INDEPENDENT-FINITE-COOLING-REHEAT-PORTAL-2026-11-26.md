# Vacuum source with a prescribed cooling phase before radiation domination

The [instantaneous vacuum reference](INDEPENDENT-LOW-REHEAT-VACUUM-NORMALIZATION-SENSITIVITY-2026-11-26.md)
matches the stable fixed-count target near 9.019 MeV. Retain
the same FD visible source, physical masses, and vacuum
final-state portal. Ask how much this number changes if
production also occurs during a finite cooling phase.

Use a prescribed $T\propto a^{-\alpha}$, $H\propto a^{-\beta}$
before the radiation segment, joining continuously in $H$
to the radiation comparison at $T_R$. The standard illustrative
matter-dominated cooling choice is $\alpha=3/8,\beta=3/2$.
The [Giudice, Kolb and Riotto reheating analysis](https://arxiv.org/abs/hep-ph/0005123)
motivates distinguishing a maximum bath temperature and the
end of reheating; [Bernal, Deka and Losada](https://doi.org/10.1103/PhysRevD.111.055034)
discuss dependence on a parametrized cooling rate. We do not
solve either model. In particular our idealized join matches
the radiation $H$ at $T_R$ and is not an inflaton/radiation
equality calculation.

For physical number source $R(T)$ with $a_R=1$, the dilution
and expansion factors give a pre-radiation contribution
$$
 Y_{\rm before}=
 \frac{1}{\alpha H_R s_R}\int_{T_R}^{T_{\max}}d\ln T\,
      R(T)\left(\frac{T_R}{T}\right)^{(3+\beta)/\alpha}.
$$
Redshift momenta from each production time and propagate the
two prompt decay daughters with their exact second energy
moments. The code checks the source moment independently by
explicit production and decay angle integration, and checks
the cooling Jacobian both in $\ln T$ and in $\ln a$.

| $T_{\max}/T_R$ | Matched $T_R$ (MeV) | Fraction of stable number made before radiation | Prompt-only scalar Ωh² | Scalar RMS momentum / photon T at 1 MeV |
| ---: | ---: | ---: | ---: | ---: |
| 1 (instantaneous) | 9.01907 | 0 | 2.36e-12 | 6.696 |
| 1.5 | 6.61791 | 0.812 | 4.03e-12 | 3.709 |
| 2 | 6.36378 | 0.852 | 4.78e-11 | 2.080 |
| 3 | 6.26421 | 0.866 | 5.65e-10 | 1.032 |
| 5 | 6.23618 | 0.870 | 2.23e-9 | 0.604 |

At ratio 5, $T_{\max}=31.18$ MeV and the prompt scalar number
is larger than in the instantaneous case but still tiny compared
with 0.120. The stable RMS momentum/T changes from 4.457 to
3.292. These comparisons keep the conditional stable mass
density target fixed, not the end temperature. At a fixed
maximum temperature of 10, 20 or 40 MeV, the matched end
temperatures are 6.60417, 6.25805 and 6.23163 MeV. At fixed
ratio 3, changing only $\alpha$ from 0.25 to 0.375 and 0.50
moves the root from 7.01710 to 6.26421 and 5.13565 MeV.
The differences between assumptions are far greater than
the numerical root refinement in each assumption. The largest
relative root change in the ratio grid is 6.9e-9; independent
angular and power-law-clock checks agree to about 2e-15.

Rest-energy and RMS/Jensen bounds on the combined dark density
at 1 MeV give about 5.80e-7–5.94e-7 of the visible energy for
these examples. The pre-produced density at the join is below
1.12e-7 of the visible energy where there is a pre-stage.
These checks are local to the stated surfaces; they do not
establish the approximation throughout reheating.

For this collisionless continuation from 1 MeV, the concavity
of speed as a function of $p^2$ bounds mean speed by
$p_{\rm rms}/\sqrt{p_{\rm rms}^2+m^2}$. In the reference
radiation-plus-matter geometry this gives upper bounds on
the remaining mean comoving displacement to equality. For
the instantaneous and ratio-5 examples, respectively, they
are 20.98 and 3.21 Mpc for the scalar, 0.247 and 0.182 kpc
for the stable species. These are **kinematic upper bounds**,
not actual mean displacements, transfer functions or structure
constraints. They do not include motion before 1 MeV.

The reference neglects the rising-temperature branch,
in-progress bath thermalization, dark collisions and inverse
rates, particle/visible feedback, changes to the portal and
post-freeze entropy injection. It is a conditional source and
spectrum-moment comparison, not a BBN or observational bound.

~~~sh
OPENBLAS_NUM_THREADS=1 python -W error finite_cooling_reheat_portal_reference.py EOS.dat
~~~
