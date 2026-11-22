# Exact CM daughter-cell angular overlap in the conservative Z′→ν+antiν event graph (2026-11-22)

## Specific numerical hypothesis

The first conservative **time-evolving** energy-lattice quantum event
graph [audit](VECTOR-NEUTRINO-REVERSIBLE-EVENT-GRAPH-AUDIT.md)
preserves charge/energy/Pauli/Bose/entropy but has large finite-grid
bias relative to an independent continuum weak 1↔2 integral.
Part of the error could be a **nonconservative angular endpoint
weight**, rather than defective quantum gain/loss. Test this with a
different daughter energy lattice and exact partial-cell intersection,
holding the same vector energy midpoints, dispersions, physical
Pauli/Bose event flux and independent continuum source.

Pre-declared pass criterion: a flat event quantum factor must integrate
the full physical moving-vector polar angle with probability exactly
one on *every* allowed vector shell, not merely approximately on
averaged rates. Compare with a naive center-hit rule on exactly
the same new lattice and quantify separate N,E,S errors against the
independent continuum Gauss 144×144. Fail closed on rate error,
no posthoc conservation rescaling and no claim of full convergence.

## Geometry and positive microreversibility

Take E_V,j=j h, p_j=sqrt((jh)^2−mV²). Neutrino/antineutrino
cell centers are k_i=(i+1/2)h, k_l=(l+1/2)h with
j=i+l+1, exactly conserving E_V,j=k_i+k_l.
Each neutrino energy cell is [ih,(i+1)h].

The exact on-shell CM polar range for fixed vector p_j is

    kmin=(E_j−p_j)/2,  kmax=(E_j+p_j)/2,

and the **analytic overlap length** for daughter cell i is

    Delta k_{j,i} =
       max[0, min((i+1)h,kmax) − max(ih,kmin)].

The nonnegative reaction weight is

    W_{j,i,l}=wV,j Γ0 (mV/E_j) Delta k_{j,i}/p_j,

with wV,j=3 p_j E_j h/(2π²), neutrino g=2,
antineutrino g=2 and the same F−R as the prior
conservative solver. Because the cells partition [0,E_j],

    sum_i Delta k_{j,i}=kmax−kmin=p_j

and therefore **for every vector mode**

    sum_i W_{j,i,l}=wV,j Γ0(mV/E_j).

This is an exact algebraic normalization and a mathematical
difference from the center-hit approximation Delta k=h for any
cell whose midpoint lies within the kinematic range.
The exact overlap does not rescale the total rate after summing
events and does not modify occupations to repair conservation.
Each on-shell elementary reaction still preserves C,L,E and
increases quantum entropy under the previously verified
bounded implicit extent.

## Measured independent weak numerical comparison

Same synthetic unequal nonthermal ν/antiν and 0.75 BE vector
spectra as the initial stage, mV=.4, Emax=4, Γ0=T=1.
Vector E midpoint weight remains wV,j~E_j p_j h.
Reference independent continuum Gauss 144×144:
dNV=+0.0058088005769424,
dEV=+0.0115779873168571,
entropy=+0.00144991988902934.

Both variants below use **identical half-offset daughter
energy bins**, so their comparison isolates the
partial angular-cell overlap.

| h | Center-hit N error | Exact overlap N error | Center-hit E error | Exact overlap E error |
|---:|---:|---:|---:|---:|
| .40 | +3.3531% | −0.7119% | +6.7802% | +4.2202% |
| .20 | +5.2037% | **+0.0824%** | +5.0406% | **+2.2292%** |
| .10 | +1.7525% | +0.2056% | +2.7410% | +1.1547% |
| .05 | +0.4030% | +0.1640% | +1.0860% | +0.5922% |

For the flat angular factor at h=.2,
the center-hit shell normalization error reaches
**+34.164%**, whereas analytic overlap gives
max absolute shell error <1e−15. At h=.4/.1/.05
the overlap shell residual stays below 4e−15.
This is a **verified mathematical fix** to
angular normalization of the piecewise constant
cell representation, not yet a continuum-ready
thermal production source.

Entropy-source relative errors at h=.2 are
+5.30465% (center) versus +0.80588% (overlap).
At h=.05 they are +0.63633% versus +0.28158%.
The overlap construction is therefore beneficial
for the tested N/E/S weak moments and exactly
removes a substantial angular alias.

## Critical negative finding and caveats

Even with **exact** angular-shell probability,
the finite-vector energy/radial measure, moment
collocation and UV cutoff still bias the weak energy
moment by +2.23% at h=.2 and +0.59% at h=.05,
while number h=.05 remains +0.16%. The sign
and magnitude do not establish an asymptotic
h² regime from the chosen limited grids.

The first lattice and this second half-offset lattice
have different daughter grid origins. Do NOT compare
their errors as though only an angular weight changed:
the fair center-hit-versus-overlap comparison is
within this stage. Likewise an exact discrete energy
ledger describes the **discretized** mode energies,
not automatically the continuum integral of E f(E).
The mass threshold vector p=0 cell, cellwise
neutrino energy weights, cutoff and physically
nonzero ν masses were not continuously resolved.
Time stepping remains first order in the separate
implicit events, with the original dt-bias.
No active weak ν↔e±γ collisions, vector/dark pair
channels, χ1/χ2→real φ Bose, H/FRW or physical
mediator parameter constraints have been solved.

## Reproduction and next decision

- [Executable overlap vs center reference](vector_neutrino_shell_cell_overlap.py)
- [Four focused tests](../tests/test_vector_neutrino_shell_cell_overlap.py)
- [Compact grid scan](results/vector-neutrino-shell-overlap-2026-11-22.json)
- [Permanent stage](../archive/stages/2026/11/2026-11-22-vector-neutrino-shell-overlap.md)

    python -m pytest -q tests/test_vector_neutrino_shell_cell_overlap.py
    python -m kinetics.vector_neutrino_shell_cell_overlap

NumPy, SciPy, deterministic Gauss and exact interval geometry;
mV=.4, Emax=4, h=.4,.2,.1,.05, Γ0=1; all values
arbitrary scaled scientific fixtures. New actual
GitHub Scientific CI must be checked, not inferred
from local prototype or documentation. Primary
quantum weak context: [Froustey/Pitrou/Volpe 2020]
(https://arxiv.org/abs/2008.01074).

Next priority: conservative vector-energy and
neutrino cell-integrated 0th/1st moment measures,
nonnegative microreversible remapping and a
coupled higher-order time solver with independent
N/E/S spectra and 3+ h/dt/cutoff controls. Then
physical ν–e±γ weak plasma and χ1/χ2/φ + Z′
gauge-consistent Bose kinetics. **No physical χ
Γ_el, 86-eV ΩDM, frozen f(q), BEC, or CLASS-ready
claim follows.**
