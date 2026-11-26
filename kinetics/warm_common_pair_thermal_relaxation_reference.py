"""Convert the two-bath response into physical linear relaxation eigenvalues."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import expit,roots_legendre

from late_radiation_background_reference import MPL,TabulatedRadiation
from three_species_decay_onsager_reference import m,mp
from warm_light_occupation_stage_audit_reference import read


def susceptibility(mass,T,mu,g,BE,n,cap=1.07):
    x,w=roots_legendre(n)
    edges=np.r_[0.,np.geomspace(1e-6,cap,28)]
    h=np.diff(edges)/2
    p=((edges[:-1]+edges[1:])[:,None]/2+h[:,None]*x)
    E=np.hypot(mass,p)
    u=(E-mu)/T
    f=1/np.expm1(u) if BE else expit(-u)
    vol=h[:,None]*w*p*p*g/(2*np.pi**2)
    v=np.stack((np.ones_like(u),u),axis=-1)
    a=(vol*f*(1+f if BE else 1-f)).ravel()
    vv=v.reshape(-1,2)
    C=(vv.T*a)@vv
    assert np.all(np.linalg.eigvalsh(C)>0.)
    return C,float(np.sum(vol*f*E))


def modes(req,eos,n):
    T,mu=req["T_GeV"],req["mu_GeV"]
    Cb,Eb=susceptibility(mp,T,mu,1,True,n)
    Cf,Ef=susceptibility(m,T,mu,4,False,n)
    D=np.linalg.inv(Cb)+np.linalg.inv(Cf)
    d,v=np.linalg.eigh(D)
    R=(v*np.sqrt(d))@v.T
    J=np.array(req["extrapolated"]["richardson_response"])*req["equilibrium_forward_event_number_energy"][0]
    np.testing.assert_allclose(J,J.T,rtol=2e-10,atol=0.)
    J=(J+J.T)/2
    K=R@J@R
    l,vec=np.linalg.eigh(K)
    target=np.linalg.eigvals(D@J)
    np.testing.assert_allclose(np.sort(target),l,rtol=1e-13,atol=0.)
    assert np.all(l>0.)
    visible=TabulatedRadiation(eos,T).state(1.)[1]
    H=np.sqrt((visible+Eb+Ef)/3)/MPL
    Cjoint=np.block([[np.linalg.inv(Cb),np.zeros((2,2))],[np.zeros((2,2)),np.linalg.inv(Cf)]])
    flux=np.block([[J,-J],[-J,J]])
    z=np.linalg.eigvals(flux@Cjoint)
    assert max(abs(np.sort(z)[0:2]))<1e-13*max(l)
    np.testing.assert_allclose(np.sort(z)[2:],l,rtol=3e-13,atol=0.)
    return dict(quadrature=req["quadrature"], susceptibility_order=n,scalar_susceptibility=Cb.tolist(),fermion_susceptibility=Cf.tolist(),
                dark_energy_density_GeV4=[Eb,Ef],radiation_density_GeV4=visible,H_GeV=float(H),decay_rates_GeV=l.tolist(),
                decay_rates_over_H=(l/H).tolist(),decay_times_in_Hubble_time=(H/l).tolist(),
                slow_fast_normalized_coordinates_columns=vec.tolist(),
                naive_diagonal_flux_rates_over_H=(np.linalg.eigvalsh(R@np.diag(np.diag(J))@R)/H).tolist())


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("response",type=Path)
    a=p.parse_args()
    r=read(a.response)
    output=[modes(z,a.eos,n) for z in r for n in (10,18)]
    for x,y in zip(output[::2],output[1::2],strict=True):
        np.testing.assert_allclose(x["decay_rates_GeV"],y["decay_rates_GeV"],rtol=1e-11,atol=0.)
    print(json.dumps(dict(scope="Thermal-reservoir tangent modes from the scalar-pair response. Visible temperature equals synthetic bath T only for H normalization; not a true evolving kinetic spectrum.",cases=output),indent=2,allow_nan=False))
