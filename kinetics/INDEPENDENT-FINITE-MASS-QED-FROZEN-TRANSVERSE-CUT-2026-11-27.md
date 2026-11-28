# Below the transverse pole gap the frozen bath still has a spacelike cut

Extend the [frozen transverse pole gate](INDEPENDENT-FINITE-MASS-QED-FROZEN-TRANSVERSE-PLASMON-GATE-2026-11-27.md)
to the retarded transverse spectral density for spacelike photon
four-momentum. The same isotropic collisionless charged bath,
coupling and antiparticle proxies are held fixed. This is the
retarded spectral weight, not an assumed thermal photon population
or a scalar collision integral.

For positive frequency with $u=\omega/k<1$, the Landau imaginary
part is nonzero only for charges with speed $v>u$. The retarded
logarithm fixes its sign. On a histogram, integrating the charged
derivative by parts gives a positive radial primitive proportional
to $u(1-u^2)\int_{p_{\min}}^\infty p f(p)dp$, where
$p_{\min}=m u/\sqrt{1-u^2}$. The lower boundary term vanishes for
the transverse projection. The signed jump formula and independent
positive cell form agree to better than 9e−17 of the Debye mass
squared. This prevents numerical branch choices from generating
negative spectral weight.

Integrating $\omega\rho_T/\pi$ over positive frequencies gives
the transverse pole weight plus the cut weight. The sum agrees
with one within 7.1e−14 on the sampled momenta. In this particular
collisionless response the cut carries about 3.27% of that
canonical weight at momentum .0008, and 7.93% at .00256.
The part below the held scalar energy is respectively about
3.27% and .70–.71% of the canonical weight. These are spectral
weights of a propagator, not fractions of a reaction rate.

The receipt also reconstructs spacelike absorption points below
the scalar energy with positive spectral density and nonzero
charged archive occupations inside the kinematic aperture.
Removing the lightlike vacuum pole therefore does not remove all
low-energy electromagnetic support. Turning this observation into
a rate requires the non-null collision contraction and an
appropriate statistical photon correlation; in particular the
nonthermal archives do not establish an equilibrium fluctuation–
dissipation relation. Longitudinal modes and collisional widths
are outside this check.

Replay (kinetics):
~~~sh
OPENBLAS_NUM_THREADS=1 python -W error finite_mass_qed_frozen_transverse_cut_reference.py
~~~
[Raw receipt](results/finite-mass-qed-frozen-transverse-cut-2026-11-27.json).
