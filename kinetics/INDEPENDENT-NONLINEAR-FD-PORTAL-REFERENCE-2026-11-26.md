# Conditional nonlinear lepton-portal FD trajectories

This extends the preceding local finite-width FD gate. The restricted model contains one symmetric 20 MeV Dirac pair species plus an internally thermal relativistic scalar and a consistently reconstructed tabulated visible EOS. It includes direct electron, muon and tau channels only. It omits chi1, scalar particle-number evolution and nonthermal states, complete hadrons, in-medium propagation and separate real mediator dynamics. These are numerical controls, not a completed physical relic history.

## Thermal state and conservation

For each log scale factor integrate actual pair number, the dark temperature, the SM temperature, and independent conservation ledgers. Forward timelike events remove two fermions; the inverse creates two. The visible bath and dark bath receive exactly opposite heat. Charged-lepton elastic scattering adds heat but no particle source and produces nonnegative total entropy; all modeled flows use the same sign convention. Verify the timelike calculation with an independent set of entropy coordinates and finer invariant and angle quadrature as well as no-coupling and initial-displacement controls.

Starting at common T=0.6882487043 GeV and zero FD chemical potential, the finite-width timelike test reaches T_D/T_SM=0.5116243 at visible T=10 MeV versus 0.4754995 without reactions. The corresponding a^3 n/n_init values are 1.1044313 vs 1, and n/s_SM is 0.01442046 vs 0.01294367. Here n counts both charge sectors. Initial dark temperature ratios .8 and 1.2, with mu_D/T_D=-.4 and +.35, lead to .53913 and .48105 at the endpoint. The actual number can increase or decrease; equilibrium at one instant does not erase dependence on the initial conditions. The common start repeated with evolving entropy coordinates and raised quadrature orders changes the ratio by less than 9e-9 and normalized count by less than 1.8e-8. The separate classical control matches the actual initial number and temperatures, not the FD chemical potential; its terminal ratio .5058071 is not a Pauli-only comparison.

## Independent elastic geometry and combined trajectory

Transform the elastic phase space to pair invariant energy u, total laboratory energy L and incoming polar angle along the pair boost. With exact Møller relative velocity the transformed measure is P times p_star squared divided by u, with P the lab total three-momentum magnitude and p_star the CM magnitude. The outgoing polar angle uses the same boost axis, with relative azimuth integrated independently. This implementation does not depend on the upstream incoming-momentum code. On the selected controls, a separate incoming-coordinate calculation agreed within 2.2e-5. The included public test also checks angular/invariant refinement, thermal cutoff and an independent finite difference of the equilibrium elastic conductance. The conductance interpolator is positive and exact in heat direction at equal temperatures. An early log-linear temperature interpolation failed the 11% accuracy gate at the coldest off-grid point; the default is now shape-preserving interpolation in log temperature, with the old interpolation retained only as a trajectory sensitivity.

Adding finite-vector charged-lepton elastic heat shifts the common-start terminal T ratio by about +0.00280 and count by about +0.000116 in a^3n/n_init units. These are restricted sensitivities, not a bracket for neglected channels or a determination of physical kinetic decoupling. Both pair creation and changes in H or bath temperatures can change the final counts; the added elastic term itself does not create pairs.

## Reproduce

Use the central three-column [Saikawa-Shirai EOS data](https://member.ipmu.jp/satoshi.shirai/DM2020/eos2020.dat) (not redistributed). All test data were generated in the cloud with no manual Actions.
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error kinetics/nonlinear_timelike_mb_trajectory_reference.py EOS_TABLE MB_GRID.npy --build-grid --nodes 321
python -W error kinetics/nonlinear_timelike_fd_trajectory_reference.py EOS_TABLE MB_GRID.npy
python -W error kinetics/boost_axis_elastic_heat_reference.py
python -W error kinetics/nonlinear_combined_charged_lepton_fd_reference.py EOS_TABLE ELASTIC_GRID.npz --build
python -W error kinetics/nonlinear_combined_charged_lepton_fd_reference.py EOS_TABLE ELASTIC_GRID.npz
```
