"""Chemical-affinity consistency of the FULL real-scalar Dirac Yukawa network.

Uses only reaction stoichiometry and exact FD/BE odds; it DOES NOT
integrate a physical 2->2 collision rate or assert thermalization.
For a REAL scalar, particle/antiparticle chemical potentials are opposite
when only vector Dirac number remains. Pair reactions force mu_phi=0 at
common-temperature full chemical equilibrium; then mu_chi1=mu_chi2=mu_U1.

Earlier mu1=-.15,mu2=.05,mu_phi=-.2 is a valid fixed point ONLY
of the restricted particle decay chi1<->chi2+phi, not the full network.
"""
import json
import math
import numpy as np
from scipy.special import expit
from kinetics.phi_bose_yukawa_charge_audit import (
    SPECIES,CHANNELS,stoichiometry,
)


def chemical_affinities(mu):
    """A_r=sum(mu_initial)-sum(mu_final)=-S@mu, in GeV."""
    v=np.asarray(mu,float)
    if v.shape!=(5,):
        raise ValueError("Expected chi1, anti-chi1, chi2, anti-chi2, phi")
    return (-stoichiometry(True)@v).tolist()


def charge_symmetric_potentials(mu1,mu2,muphi):
    """Unbroken net U(1) permits mu(anti)= -mu(particle)."""
    return np.array([mu1,-mu1,mu2,-mu2,muphi],float)


def paired_fermion_bose_weights(e2,eant,k1,k2,
                               mu2,muanti,muphi,T=1.):
    """Exact microreversible quantum forward/inverse occupation factors.

    Valid only for on-shell energy equality e2+eant=k1+k2, and
    strictly subcritical BE mu_phi<min(k1,k2). No matrix element included.
    """
    if T<=0 or abs((e2+eant)-(k1+k2))>2e-12:
        raise ValueError("Invalid temperature or off-shell event")
    if min(k1,k2)<=muphi:
        raise ValueError("Bose chemical potential >= scalar event energy")
    f=expit(-(np.array([e2-mu2,eant-muanti]))/T)
    b=1/np.expm1((np.array([k1,k2])-muphi)/T)
    F=float(np.prod(f)*np.prod(1+b))
    R=float(np.prod(1-f)*np.prod(b))
    return F,R


def event_entropy(F,R):
    """Quantum pair contribution (F-R) ln(F/R)>=0 for positive microevents."""
    if F<=0 or R<=0:
        raise ValueError("Only interior positive occupations allowed")
    return (F-R)*math.log(F/R)


def audit():
    restricted=charge_symmetric_potentials(-.15,.05,-.2)
    exact=charge_symmetric_potentials(.05,.05,0.)
    f,r=paired_fermion_bose_weights(.1,.1,.09,.11,.05,-.05,-.2)
    return dict(schema_version=1,species=SPECIES,channels=CHANNELS,
                restricted_particle_decay_fixed_point_mu_GeV=restricted.tolist(),
                full_network_affinities_at_restricted_point_GeV=chemical_affinities(restricted),
                full_chemical_eq_one_example_mu_GeV=exact.tolist(),
                full_chemical_eq_affinities_GeV=chemical_affinities(exact),
                pair_forward_inverse_example={"F":f,"R":r,"F_over_R":f/r,
                                               "predicted_exp_affinity":math.exp(.4),
                                               "positive_entropy":event_entropy(f,r)},
                only_continuous_full_U1=[1,-1,1,-1,0],
                limitations="Four topology network and equilibrium odds only; not actual 2to2 Yukawa amplitude/thermal rates, hidden SM portal, nonthermal time evolution, BEC or freezeout",
                class_ready=False)


if __name__=="__main__":
    print(json.dumps(audit(),indent=2,allow_nan=False))
