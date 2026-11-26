# Short-time relaxation of the resolved exterior prompt current

The earlier q32 archive showed a large gross outgoing and returning
rate but did not evolve it. Here evolve the actual embedded q32
FD/BE entropy shells around that same initial snapshot, **without
expansion or visible source**, with the retained stationary-parent
prompt integral only. The q16 boundary is internal, not absorbing.
The fixed H normalization includes tabulated radiation and retained
dark energy at the start of the short control.

The ensemble ODE uses each shell's positive number/kinetic
susceptibility to convert the projected moment flow into derivatives
of its two entropy coefficients. A separate analytic ODE Jacobian
includes the derivative of this susceptibility, not just the collision
Jacobian. Its finite-difference directional check is about 1.13e-10.
With Radau, both represented species counts and energy are conserved
to better than 1.4e-13 relative across the archived checks, and entropy
increases.

The scalar signed exterior current starts at -178.49 N_phi H (net
inward), falls to -3.22 N_phi H after 3e-8 Hubble times, and to -0.143
N_phi H after 3e-7 Hubble times. The exterior scalar fraction changes
from 2.2135e-5 to 2.1024e-5. The coldest integrated scalar shell count
stays at about 1.91e-8 of scalar number. The transient therefore cannot
be estimated by holding even the initially small gross-current
difference fixed for a long interval. This is a short fixed-background
prompt-only response, not a cosmological evolution or an abundance result.

## Numerical controls

The prompt evaluation centers each parent at a local affine pivot,
using direct positive sums when angular occupations are widely separated.
The projected flux, Jacobian-vector products and conservation checks
agree with the previous implementation. An independent elastic-cache
check also avoids computing the unused global Gram in a cache-only
build, and represents each elastic panel in its active local columns.
The pairwise version computes a stable antisymmetric exponential
difference. The 5635-panel test recovers the original flux and analytic
JVP at floating-point precision. Single measured times are hardware
diagnostics, not general performance guarantees.

## Reproduce

Prepare the [external EOS data](INDEPENDENT-BENCHMARK-SETUP.md).
Use the earlier raw empty-bath spectrum file for the isolated numerical
controls. From `kinetics`:

```sh
export OPENBLAS_NUM_THREADS=1
python -W error compact_elastic_panel_reference.py EOS.dat empty.json
python -W error pivoted_prompt_reference.py EOS.dat empty.json
python -W error prompt_exterior_fixed_time_reference.py EOS.dat EOS.late.dat results/wide32-embedded-vacuum-2026-11-26.json results/wide32-embedded-vacuum-shape-2026-11-26.json --extent 3e-7
```
