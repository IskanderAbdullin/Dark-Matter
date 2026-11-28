# Positive energy-assigned nodal stencil for early paired decay

This independent check modifies the numerical stencil, not the early
particle model. It uses the same held FD parent, massless fermion
daughter, positive small quartic and prescribed background as the
[paired Hartree control](INDEPENDENT-EARLY-POLE-PAIRED-HARTREE-MATRIX-2026-11-26.md).
The scalar momentum nodes and conservative Hartree variation are unchanged.
Both runs use a scalar momentum floor of $10^{-9}T$.

For scalar energy $E_i$, project the exact geometric daughter aperture
onto nonnegative piecewise-linear hat weights `w_ij` at daughter nodes
`e_j`. Their sum and first moment reproduce the clipped aperture's uniform
measure and its energy moment. Then evaluate the held parent FD occupation
at each **discrete** channel's energy `E_i+e_j`. Each channel loses that
energy from the reservoir and supplies the same energy to the two
daughter distributions. This is an algebraic collision energy identity;
the bath itself is not evolved, and the FD-weighted physical continuum
moment is still approximated. A receiving node can lie outside the
exact aperture that contributed to its hat. At nonzero grid spacing,
this is a positive quadrature, not exact physical kinematics.

The forward and reverse terms at each node use the same massless
fermion occupation, so both faces $f_b=0,1$ are inward. With thermal
occupations at the discrete energies and zero chemical potentials,
the channel-level detailed-balance identity holds. Scalar and fermion
number increments use exactly the same signed weights. At an
independent fixed aperture test, geometric width and energy moment
errors are below $4\times10^{-15}$ and $2\times10^{-13}$;
thermal detailed balance agrees to $2\times10^{-15}$ relative.
Compared to exact FD aperture integration, the largest empty-daughter
seed error over the selected active rows is 0.159% coarse and 0.0413%
fine. The row selection omits rows below $10^{-7}$ of the maximum.

At $x=0.05$, integrated scalar yields are $4.43529\times10^{-4}$
and $4.38025\times10^{-4}$ (coarse and fine), about 1.26% apart.
The fine value is 0.073% above the matched piecewise-constant
overlap result $4.37706\times10^{-4}$; momentum refinement matters
more for this yield than the change of stencil at the tested fine
resolution. Pair-count discrepancies remain below $2\times10^{-15}$
in both nodal trajectories. None of this fixes the held-reservoir,
nonequilibrium-width or expanding-bath assumptions.

~~~sh
python -W error early_pole_paired_nodal_energy_reference.py --static-only
python -W error early_pole_paired_nodal_energy_reference.py
python -W error early_pole_paired_nodal_energy_reference.py --scalar 130 --fermion 190
~~~

Raw trajectories: [coarse](results/early-pole-paired-nodal-coarse-2026-11-26.json.gz),
[fine](results/early-pole-paired-nodal-fine-2026-11-26.json.gz).
