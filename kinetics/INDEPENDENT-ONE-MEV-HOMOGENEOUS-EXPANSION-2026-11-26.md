# Homogeneous expansion increment from the low-reheat test at 1 MeV

Carry the reproduced early low-reheat quantum state to photon
temperature 1 MeV with the signed pair-only continuation and q64
support. Evaluate the energy density from each retained momentum
cell as E=√(m²+p²); do not treat either species as strictly massless.
Use the reference visible-radiation clock, not a neutrino transport
or nucleosynthesis network.

Retained scalar and stable energies are about 1.09e-10 and 5.85e-7
of visible radiation. Add a conservative upper for the earlier
positive exterior cohorts by allowing them to dilute in number
but lose no particle energy; use the independent redshifted upper
bound for the later exterior tails. Including these uncertainties,
the fractional Hubble increase from the modeled dark energy is
2.926e-7 to 2.964e-7 at the same photon temperature. This is
equivalent in energy to (3.53–3.57)e-6 of one ideal massless neutrino
flavor at Tν=Tγ, **not** a prediction of observable ΔN_eff.

Visible thermalization after low reheating, neutrino phase space,
light-element abundances, and omitted interactions or extra
particles are not evaluated. The small direct dark density cannot
establish BBN viability by itself.

Use the same published early quantum snapshot:
```sh
export OPENBLAS_NUM_THREADS=1
python -W error scalar_pair_expanding_continuation_reference.py EOS.late.dat results/low-reheat-early-pair-quantum-shape-2026-11-26.json --stop 3.1734356669329107 --steps 12 24 --pad 64 --save-snapshots one_mev_shape.json > one_mev.json
python -W error low_reheat_mev_expansion_gate_reference.py EOS.late.dat results/low-reheat-early-pair-quantum-2026-11-26.json one_mev.json one_mev_shape.json
```
