# Long wavelength and vanishing pair chemistry do not commute

Interpolate the [chemically equilibrated](INDEPENDENT-CHARGED-LEPTON-LOCKED-SOUND-DAMPING-2026-11-27.md)
and [number-conserved](INDEPENDENT-CHARGED-LEPTON-FROZEN-CHEMICAL-SOUND-2026-11-27.md)
static controls. Multiply the complete timelike kernel (forward and
inverse number, energy, and momentum) by ε, leaving the elastic
kernel and background state unchanged. This maintains detailed
balance and total conservation; it is not a change to just one
row of the relaxation matrix. Normalize both time and wave number
to the active elastic heat coefficient divided by dark enthalpy,
the same clock at every ε.

In entropy coordinates write C(ε)=C_el+ε C_pair, with both matrices
positive semidefinite. At ε=0 there are three slow conserved
directions. For ε>0 only two remain. Let v_c be the unit direction
in the first null space orthogonal to the second. Degenerate
perturbation theory gives the small chemical decay rate
γ_c=ε (v_cᵀ C_pair v_c)+O(ε²). If v_s is the equilibrium acoustic
direction, the same almost-conserved chemical mode yields a
leading contribution to the small-k attenuation:

    D(k→0,ε) = lim_{k→0}(−Re λ/k²)
    ε D(k→0,ε) → |v_cᵀ B v_s|² / (v_cᵀ C_pair v_c).

Thus a coefficient diverging as 1/ε is not an infinite damping
rate at fixed wavelength. Its expansion shrinks to wave numbers
small compared with γ_c. The order of limits also changes the
sound speed. For k≪γ_c the sound is chemically equilibrated.
For γ_c≪k, while elastic exchange is still fast, the sound
is effectively number-conserved.

| T (GeV), h_V/h_D=10 | Equilibrium speed | Number-conserved speed | Limit of ε D at small k |
|---:|---:|---:|---:|
| .015 | .57502401 | .57543066 | 9.39306e−5 |
| .050 | .57711659 | .57712240 | 1.57685e−6 |

The numerical scan varies ε from 10⁻³ to 10⁻⁶ and compares the
fixed-k generator with the independent projected singular
coefficient. At the smallest ε it reproduces both acoustic speeds
on opposite sides of the chemical scale and agrees with the
predicted εD and γ_c/ε to relative error below 2e−6.
This is a static local dispersion check with a controlled process
multiplier, not the location or shape of any cosmological transition.

Replay from the kinetics directory:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error \
  charged_lepton_finite_chemical_sound_crossover.py \
  results/charged-lepton-portal-combined-momentum-2026-11-27.json
```
[Receipt](results/charged-lepton-finite-chemical-sound-crossover-2026-11-27.json).
