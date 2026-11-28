"""Time-resolved lower comparison under a fermion cap, not a scalar cap."""
import argparse
import json
from pathlib import Path

import numpy as np

from three_species_decay_onsager_reference import m,M,y
from warm_fermion_survival_with_leptons_reference import subpole


def read(p):
    return json.loads(Path(p).read_bytes())


def setup(N,E,K,clock,steps):
    b=clock["b_stop"]; h=clock["floor_min"]; Hph=clock["photon_anchor_GeV"]
    t=np.linspace(0,np.log(b),steps+1)
    # Under the known radiation lower floor, exact upper integrals in
    # each cell. After b_stop use only freely redshifting photons.
    integ=np.stack([np.r_[np.diff(-np.exp((2-p)*t))/h,np.exp((2-p)*t[-1])]
                    /((p-2)*Hph) for p in (3,4,5)],axis=1)
    eta=np.geomspace(1.008,4.,96)[:,None]
    cuts=np.linspace(6.,9.9,40)[None,:]
    m_i=np.array([m,M]); m_x=m_i[::-1]
    Cs=np.max(2*y**4/(np.pi*m_x**2*m_i**4))/np.sqrt(1.-eta**-2)
    Ch=np.max(2*y**4/(np.pi**3*m_x**2*m_i**2))*eta**2
    coup=.5*N*(K+subpole(cuts))
    J=.5*(3*np.pi**2*N)**(2/3)
    ws=2*Cs*J*E*E/(4*N)
    wh=2*Ch*(J*2**(-2/3))**2*E/N
    deficit=2*E/(cuts*N)
    return integ,coup,ws,wh,deficit,eta,cuts


def comparison(model,cap):
    weights,coup,ws,wh,deficit,eta,cuts=model
    Y=np.ones((len(eta),cuts.size))
    x=cap**(1/3)
    ws,wh=ws*x,wh*x*x
    for A,B,C in weights:
        r=np.cbrt(Y)
        # The exact decreasing comparison solution loses no more over
        # this cell than the rates at its start, evaluated at the previous
        # lower starting value. Scalar comparison preserves ordering.
        Y=np.maximum(0.,Y-A*coup*Y*Y-B*wh*r**4-C*ws*r*r-2e-14)
    R=Y-deficit
    i,j=np.unravel_index(np.argmax(R),R.shape)
    return max(0.,float(R[i,j])),float(eta[i,0]),float(cuts[0,j])


def run(fermion,leptons,isolated,coarse):
    F,L,I,C=map(read,(fermion,leptons,isolated,coarse))
    output=[]
    for a,b,c,d in zip(F["cases"],L["cases"],I["cases"],C["cases"],strict=True):
        assert a["q"]==b["q"]==c["q"]==d["q"]
        n=a["core_start_fermion_number_GeV3"]; o=a["conservative_initial_fermion_Omega"]
        E=b["available_energy_upper_GeV4"]; K=L["coefficient"]
        clock=c["bounds"][1]["clock"]
        data=[]
        for steps in (320,1000,2800):
            s=setup(n,E,K,clock,steps)
            at_half=o*comparison(s,.5)[0]
            assert at_half>F["target"]
            rows=[]
            for x in (0.,.0001,.001,.003,.005,.01,.05,.1,.3,.5):
                z,eta,cut=comparison(s,x)
                rows.append(dict(cap=x,lower_Omega=float(o*z),eta=eta,cut_GeV=cut))
            data.append(dict(time_lower_cells=steps,fermion_Omega_lower_at_half_per_spin=at_half,
                             entire_signed_cap_range_excludes_target=True,values=rows))
        old=d["finite_fermion_cap_needed_to_evade"]
        val="fermion_Omega_lower_at_half_per_spin"
        assert all(d1[val]<=d2[val] for d1,d2 in zip(data,data[1:]))
        output.append(dict(q=a["q"],old_additive_threshold=old,controls=data))
    return dict(scope="Homogeneous isotropic warm signed network, both dark fermions "
                      "capped per spin, scalar occupancy unbounded. Time-dependent "
                      "loss comparison with finite Euler lower cells; not an actual evolution.",
                target=F["target"],cases=output)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("fermion","leptons","isolated","coarse"):
        p.add_argument(k)
    a=p.parse_args()
    print(json.dumps(run(a.fermion,a.leptons,a.isolated,a.coarse),indent=2,allow_nan=False))
