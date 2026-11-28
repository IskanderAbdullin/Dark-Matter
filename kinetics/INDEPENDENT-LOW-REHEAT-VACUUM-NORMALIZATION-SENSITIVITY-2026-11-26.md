# Sensitivity of the low-reheat vacuum abundance normalization

In the prescribed-radiation vacuum reference, the stable count
matches $\Omega_\chi h^2=0.120$ near $T_R=9.019$ MeV. Compute the
local logarithmic derivative directly from the positive integrand
at the upper reheating endpoint, and cross-check against finite
differences of the integrated yield. In the fixed spectral
normalization, changing the portal mixing $\varepsilon$ multiplies
all counted sources by $\varepsilon^2$, so the matched
reheating root has
$$
 \frac{d\ln T_R}{d\ln\varepsilon}
   =-2\left(\frac{d\ln\Omega_\chi h^2}{d\ln T_R}\right)^{-1}.
$$
This comparison keeps the width and other couplings fixed.

| Stable mass | Nominal reheating root | Local d ln (stable Ωh²) / d ln T_R | T_R root at 0.9 times mixing | At 1.1 times mixing |
| ---: | ---: | ---: | ---: | ---: |
| 18 MeV | 8.88810 MeV | 4.5066 | 9.32014 MeV | 8.52487 MeV |
| 20 MeV | 9.01907 MeV | 4.7902 | 9.43108 MeV | 8.67195 MeV |
| 22 MeV | 9.18187 MeV | 5.0650 | 9.57810 MeV | 8.84731 MeV |

At the nominal 20 MeV mass and mixing, changing $T_R$ by −1% or
+1% changes the fixed-count stable proxy to 0.11434 or 0.12584.
The prompt-only scalar proxy at the nominal matched root is
$2.36\times10^{-12}$. It does not include secondary scalar pair
production from the [coupled low-reheat result](INDEPENDENT-LOW-REHEAT-QUANTUM-BRIDGE-2026-11-26.md).
The quadrature comparison changes the fitted roots by at most
1.1e-8 fractionally, while the model parameter choices change
them by much more. This matters when comparing two low-temperature
cutoffs or portal normalizations: the vacuum root is a tunable
reference point, not a measured reheating temperature.

The calculation uses FD initial leptons, vacuum final states,
prompt heavy-pair decay counts and prescribed radiation from
instantaneous reheating to 1 MeV. Perturbing the stable mass
changes its production threshold and mass-to-abundance factor;
it does not change the fixed mediator width or the mass of the
heavy parent. Visible feedback, inverse depletion, dark reactions
and other nonthermal initial conditions are not included.

~~~sh
OPENBLAS_NUM_THREADS=1 python -W error low_reheat_vacuum_root_sensitivity_reference.py EOS.dat
~~~
