# Conservative FRW production with two finite daughter masses

The [new solver](two_mass_bath.py) uses rationalized energy differences, analytic full/empty cell classification and split quadrature on boundary-crossed cells. The invariant two-body measure contains eta=2D/x²; the width is the actual channel rest width, not a fixed microscopic coupling. Both daughter energies and shell pressure traces enter total energy and Friedmann expansion. The inherited massless radiation expression is explicitly renamed as a momentum proxy for massive chi2.

## Executed results

Windows Python 3.9.7 / NumPy 1.26.4 / SciPy 1.13.1, one BLAS thread. Geometry: 1820 independently sampled cells (1192 nonzero), maximum reference error normalized by rectangle area 2.67e-8, orientation error 1.73e-11, previous one-mass kernel error 1.27e-14.

| m2/M | omega_phi h² | fraction p/T>.1 | mean p/T |
|---:|---:|---:|---:|
| 0 | .118234 | .439894 | 1.03844 |
| .5 | .160733 | .136321 | .254650 |
| .707107 | .169487 | .021382 | .024806 |
| .9 | .174701 | .000473217 | .000455870 |

At m2/M=.99 BDF failed with minimum-step error. This is retained as a negative numerical result, not physical exclusion. The completed rows passed balance, directional Jacobian, daughter-number, Pauli and continuity checks; maximum absolute continuity residual about 5.5e-6. Mesh convergence and matched-abundance width tuning remain open.

[Actual selected outputs](results/two-mass-bath-2026-11-25.json). The local project has been moved into the explicitly authorized R:/Workspace/projects/Dark-Matter-Research; its isolated venv was recreated and verified.

Yin JHEP05(2023)180, sections IV.1/IV.4 (https://doi.org/10.1007/JHEP05(2023)180), motivates mass degeneracy and eventual return of RHN energy to the SM bath. At m2/M=.9, the equal-number stable fermion rest-density estimate is omega_chi2 h²~1.57e8; the scalar abundance is also above target at this fixed width. Therefore no massive stable daughter scenario has been established. Next: implement detailed-balanced effective removal/inverse creation of chi2 with independent number-source and energy closure; retune physical width only after its effect is understood.
