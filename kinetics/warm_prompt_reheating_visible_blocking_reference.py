"""Signed FD visible-to-dark rates on held two-temperature states."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import expit,roots_legendre

from empty_lepton_pair_spectrum_reference import Source
from three_species_decay_onsager_reference import m,M,mp


def channel(tvis,td,eta,mass,n):
    inv,eng,ini,fin=((10,28,20,26),(16,46,32,42))[n]
    src=Source((inv,eng,ini),tail=72.,dark_mass=mass)
    Q,P,V,B,A=src.event_moments(tvis,include_boltzmann=False)
    x,w=roots_legendre(fin)
    E0=Q[:,:,None]/2+P[:,:,None]*np.sqrt(V)[:,None,None]*x/2
    E1=Q[:,:,None]-E0
    # All visible FD angle factors are in B and the rank-two moment A.
    mu=w*(B[:,:,None]+A[:,:,None]*(3*x*x-1)/2)/2
    assert np.min(mu)>=0.0
    Bf=np.sum(mu*expit(E0/td-eta)*expit(E1/td-eta),axis=2)
    f=Bf*np.exp(-Q/tvis)
    r=Bf*np.exp(2*eta-Q/td)
    u=B*np.exp(-Q/tvis)
    return np.array([[v.sum(),(v*Q).sum()] for v in (u,f,r,f-r)])


def run(path):
    g=json.loads(path.read_text())
    if "arms" in g:
        samples=[dict(visible_T=v["T_visible"],dark_T=v["T_dark"],phase_peak=v["peak"])
                 for v in g["arms"][0]["rows"]]
    else:
        samples=g["history"][-1]["samples"]
    out=[]
    for i,x in enumerate(samples):
        if i not in (0,3,7,11,14):
            continue
        tv,td=x["visible_T"],x["dark_T"]
        eta=np.log(x["phase_peak"])+np.array((m,mp,M))/td
        levels=[]
        for n in (0,1):
            L=channel(tv,td,eta[0],m,n)
            H=channel(tv,td,eta[2],M,n)
            assert np.all(L>=0.) and np.all(H>=0.)
            levels.append(dict(light=L.tolist(),heavy=H.tolist(),
                               light_visible_net_over_empty=(L[-1]/L[0]).tolist(),
                               heavy_visible_net_over_empty=(H[-1]/H[0]).tolist(),
                               total_net_energy_over_empty=float((L[-1,1]+H[-1,1])/(L[0,1]+H[0,1]))))
        for k in ("light_visible_net_over_empty","heavy_visible_net_over_empty"):
            np.testing.assert_allclose(levels[0][k],levels[-1][k],rtol=8e-5,atol=0)
        out.append(dict(index=i,T_visible_GeV=tv,T_dark_GeV=td,controls=levels))
    # Manufactured exactly common neutral chemical state cancels, with
    # finite quantum occupations on both sides of the lab event.
    eq=channel(.06,.06,0.,m,1)
    assert np.max(abs(eq[3]/eq[0]))<2e-15
    return dict(scope="Signed vacuum e/mu/tau vector source on held thermal states. Full FD initial and "
                      "final occupancy and exact fixed-Q reverse relation. No hidden channels or update "
                      "to the evolving dark distribution in this control.",
                rows=out,common_neutral_equilibrium_net_over_empty=(eq[3]/eq[0]).tolist())


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("continuation",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.continuation),indent=2,allow_nan=False))
