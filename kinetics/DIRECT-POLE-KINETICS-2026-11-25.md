# Direct cubic mediator: velocity-dependent elastic evolution

The earlier 1% tail-change diagnostic survives the velocity-dependent calculation for the fixed reference candidate. The 220-bin runs change N(q>2 q_peak) by **−0.81022%** and **−0.81060%** on the independent 240- and 320-cell production spectra. This is a homogeneous elastic result at fixed originally simulated abundance and a chosen illustrative cut.

| Production cells | Elastic bins | Time slices | Relative N(q>2 q_peak) change |
| ---: | ---: | ---: | ---: |
| 240 | 160 | 1 | −0.786987% |
| 240 | 160 | 8 | −0.786987% |
| 320 | 160 | 8 | −0.812457% |
| 240 | 220 | 8 | −0.810217% |
| 320 | 220 | 8 | −0.810603% |

The worst observed 160→220-bin shift is 0.0233 percentage points, below the [0.05-point comparison target](DIRECT-POLE-KINETIC-PLAN-2026-11-25.md). The fine-mesh shifts at q>6 q_peak are +0.2200%, +0.1853%, and at q>20 q_peak are +0.1383%, +0.1491%. Differences between time averaging and dividing into eight ordered exposure intervals are below 7e-8 percentage points for N(q>2 q_peak) on the tested 160-bin input. The archived fixed-contact calculations on 160 bins changed N(q>2 q_peak) by −0.81890%, −0.80120%. An isolated difference between the old and new cross sections should be interpreted on equal grids. Relative discrete number and kinetic-energy drift are at most 9e-16 and 2.3e-16. The weighted Bose detailed balance defect is below 8e-16 relative to the dilute forward scale. The greatest negative number fraction from accepted integration states is below 2e-23.

## Kinetic calculation and its scope

For comoving particle momenta in peak units, angular reduction gives a uniform measure dP over the intersection of the two pair triangles. The square of the relative speed in these units is 2(p1²+p2²)−P². Integrate the pole factor b²/[b+2(p1²+p2²)−P²]² analytically over the allowed P, where b=a[Vref/(r vpeak,0)]². The coefficients are positive and multiply both directions of each conservative graph event. The integration uses Friedmann exposure proportional to r⁴ x⁴ d ln a/H. The pre-existing event coefficients in d(p²) carry the required factor of eight relative to using d(p²/2) for the three energy measures. The independent Maxwell relative-speed test checks absolute rate normalization against direct integration; the other test checks the angle primitive against direct quadrature, including the contact limit.

All represented grid velocities at z=100 are below 0.00363c. In exact kinematics s−4M²=|p1−p2|²−(E1−E2)², and |E1−E2|≤vmax|p1−p2|. Hence the error in this numerator from the NR limit is at most vmax² relative to |p1−p2|². For the original stabilizing square, the neglected contact and crossed amplitude is bounded by g²[1/mS²+(s−4M²)/mS⁴]. With s−4M²≤4M²vmax²/(1−vmax²) as a conservative bound, the omitted amplitude is less than 1.4e-5 of the s-channel amplitude over the represented z≤100 range. These are kernel-level estimates. They do not bound a full interacting history before z=100, a primordial mediator component, forward medium dispersion, or cosmological perturbations.

Reproduce from the repository root with `python kinetics/direct_pole_kinetics.py --mesh 240 --n 220 --slices 8 --nodes 16 --output /tmp/direct-pole-240.json` and repeat with `--mesh 320`. Setting `DM_ELASTIC_COMPILED=1` uses the existing optional Numba backend. Reference data and outputs are [summarized here](results/direct-pole-evolution-summary-2026-11-25.json) with full [240-cell](results/direct-pole-evolution-240.json) and [320-cell](results/direct-pole-evolution-320.json) evolutions. Local suite at this checkpoint: 75 passed, 1 skipped.

Independent methodological references: [Chu, Garcia-Cely & Murayama, velocity-dependent resonant scattering](https://doi.org/10.1103/PhysRevLett.122.071103) and [Yang & Yu, angular-dependent gravothermal transport](https://doi.org/10.1088/1475-7516/2022/09/077). Their observational halo analyses are outside the homogeneous quantum evolution here.
