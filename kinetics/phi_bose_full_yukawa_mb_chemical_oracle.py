"""Independent classical Maxwell-Boltzmann 5-field chemical-relaxation oracle.

This is NOT the FD/BE physical quantum spectrum: here all five
species including real 86eV φ use ideal dilute MB populations,
n_i=g_i m_i² T K2(m_i/T)/(2π²), chi_i=n_i/T.
Forward equilibrium rates:
- separate χ1 and anti χ1 1↔2 decay, g=2 each:
  gamma_dec = nχ1 * Γrest * K1(m1/T)/K2(m1/T);
- both χ+antiχ→2φ positive spin/identical-scalar vacuum rates:
  gamma_pair_i = Gondolo-Gelmini MB integral from INDEPENDENT
  vectorized physical crossed-fermion amplitudes, with correct m_in/mv.
The shared 4x5 reaction matrix yields one exact net Dirac U1 zero
chemical mode, four positive finite relaxation eigenmodes for an
equilibrium-shaped fixed-T FIVE-field MB thermal ansatz.
No physical SM portal/thermal mass/FRW/86eV relic or Bose condensate.
"""
import json
import math
import numpy as np
from scipy.special import kv
from kinetics.phi_bose_yukawa_hubble_gate import rest_width
from kinetics.phi_bose_pair_independent_vectorized import reaction_density_MB
from kinetics.phi_bose_yukawa_charge_audit import stoichiometry
from kinetics.phi_bose_physical_cutcell_weak import M1,M2,MPHI

MPLANCK=1.2209e19
U1=np.array([1.,-1.,1.,-1.,0.])

def ideal_MB_number(T,mass,g):
    if T<=0 or mass<0 or g<=0:
        raise ValueError("Physical positive MB T/mass and spin count")
    if mass==0:
        return g*T**3/math.pi**2
    return g*mass**2*T*kv(2,mass/T)/(2*math.pi**2)

def equilibrium_MB_four_channel_rates(T,angle_order=128):
    if T<=0:raise ValueError("Positive temperature")
    gamma_decay=ideal_MB_number(T,M1,2)*rest_width()*(
        kv(1,M1/T)/kv(2,M1/T))
    heavy=reaction_density_MB(T,angle_order,
                              m_in=M1,m_virtual=M2)["gamma_forward_GeV4"]
    light=reaction_density_MB(T,angle_order,
                              m_in=M2,m_virtual=M1)["gamma_forward_GeV4"]
    return np.array([gamma_decay,gamma_decay,heavy,light])

def chemical_MB_eigenmodes(T,gstar,angle_order=128):
    if T<=0 or gstar<=0:raise ValueError("Positive input")
    gamma=equilibrium_MB_four_channel_rates(T,angle_order)
    masses=(M1,M1,M2,M2,MPHI)
    spin=(2,2,2,2,1)
    n=np.array([ideal_MB_number(T,m,g) for m,g in zip(masses,spin)])
    chi=n/T
    S=stoichiometry(True).astype(float)
    L=(S.T@np.diag(gamma/T)@S)/np.sqrt(np.outer(chi,chi))
    L=(L+L.T)/2
    eigen=np.linalg.eigvalsh(L)
    H=1.66*math.sqrt(gstar)*T*T/MPLANCK
    null=np.sqrt(chi)*U1
    defect=float(np.linalg.norm(L@null)/
                 (max(eigen)*np.linalg.norm(null)))
    return dict(T_GeV=T,gstar=gstar,H_RD_GeV=H,
                gamma_four_channels_GeV4=gamma.tolist(),
                particle_number_GeV3=n.tolist(),
                susceptibility_GeV2=chi.tolist(),
                chemical_eigenvalues_GeV=eigen.tolist(),
                chemical_eigenvalues_over_H=(eigen/H).tolist(),
                true_Dirac_U1_relative_matrix_zero_defect=defect,
                scope="conditional five-species ideal MB fixed-T fugacity linearized reaction network, not quantum full momentum solver/SM/FRW",
                class_ready=False)


def audit():
    cases=[chemical_MB_eigenmodes(T,g)
           for T,g in ((.05,10.75),(.1,17.25),(.3,61.75),
                       (1.,61.75),(10.,106.75))]
    return dict(schema_version=1,model="Yukawa y=.001 masses .1/.02/86e-9 GeV, two Dirac+anti plus real phi, four reactions, ideal MB around equilibrium",
                results=cases,class_ready=False)

if __name__=="__main__":
    print(json.dumps(audit(),indent=2,allow_nan=False))
