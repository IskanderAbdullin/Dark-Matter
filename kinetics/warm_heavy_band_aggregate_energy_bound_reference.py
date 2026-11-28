"""Aggregate finite-band comoving heavy energy from a positive inverse-source bound."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre
from joint_visible_heavy_adiabatic_reference import parent_branch

from three_species_decay_onsager_reference import M,m,mp,P,Es


def read(p):
    return json.loads(Path(p).read_bytes())


def inverse_bound(E,delta,cut,n_e,n_split):
    assert 0<cut<M*P/mp
    H=np.hypot(cut,M)
    emax=(H*Es+P*cut)/M
    emin=(Es**2+(mp*cut/M)**2)/emax
    walls=np.linspace(emin,emax,n_e+1)
    l,u=walls[:-1],walls[1:]
    kl=np.sqrt((l-mp)*(l+mp)); ku=np.sqrt((u-mp)*(u+mp))
    D=M*M-m*m-mp*mp
    A=np.sqrt(D*D-4*m*m*mp*mp)
    elo=np.maximum(m,(D*D+4*m*m*kl**2)/(2*(D*u+ku*A)))
    ehi=H-l
    # Optimize the light radial integral for fixed scalar energy by
    # filling the lowest allowed E first. Bound light energy allocation
    # above by the larger endpoint of each density-split cell.
    p0=np.sqrt(np.maximum(0.,(elo-m)*(elo+m)))
    frac=np.linspace(0.,1.,n_split+1)
    load=(3*np.pi**2*E/delta)*frac[1:,None]/elo[None,:]
    pf=np.cbrt(p0[None,:]**3+load)
    top=np.minimum(ehi[None,:],np.hypot(m,pf))
    measure=delta*np.maximum(top-elo[None,:],0.) * ((top+elo[None,:])/2+u[None,:])
    # The scalar may use all of the remaining density in an arbitrary
    # narrow shell. Use the left split endpoint and scalar denominator.
    B=measure*(1-frac[:-1,None])*(2*np.pi**2*E)/(kl*l*l)[None,:]
    j,k=np.unravel_index(np.argmax(B),B.shape)
    return float(np.max(B)),[float(l[k]),float(u[k])],[float(frac[j]),float(frac[j+1])]


def normalization_check():
    d=M*M-m*m-mp*mp
    root=np.sqrt(d*d-4*m*m*mp*mp)
    def f(x): return .04*np.exp(-x/.028)
    def g(x): return .10*np.exp(-x/.033)
    errors=[]
    for cut in (.05,.16,.38):
        HH=np.hypot(M,cut)
        b=(HH*Es+P*cut)/M
        a=(Es**2+(mp*cut/M)**2)/b
        x,w=roots_legendre(24)
        pieces=np.linspace(0.,cut,70)
        r,s=(pieces[:-1,None]+pieces[1:,None])/2,np.diff(pieces)[:,None]/2
        PP=(r+s*x).ravel()
        E=np.hypot(PP,M)
        H=parent_branch(E,f,g,angle=24)
        left=(2/np.pi**2)*np.dot((s*w).ravel(),PP**2*H[:,0]*H[:,-1])
        cuts=np.unique(np.clip(np.r_[np.linspace(a,b,500),HH-m],a,b))
        de=np.diff(cuts)[:,None]/2
        eps=((cuts[:-1,None]+cuts[1:,None])/2+de*x).ravel()
        k=np.sqrt((eps-mp)*(eps+mp))
        Emin=np.maximum(m,(d*d+4*m*m*k*k)/(2*(d*eps+k*root)))
        Emax=np.minimum((d*eps+k*root)/(2*mp**2),HH-eps)
        span=np.maximum(Emax-Emin,0.)[:,None]/2
        sample=Emin[:,None]+span*(1+x)
        right=M/(np.pi**2*P)*np.sum((de*w).ravel()*f(eps)*
                    np.sum(g(sample)*(sample+eps[:,None])*w*span,axis=1))
        errors.append(float(abs(left/right-1.)))
    assert max(errors)<5e-5
    return errors

def run(pions,barrier,low):
    X,Y,Z=map(read,(pions,barrier,low))
    out=[]
    for a,b,c in zip(X["cases"],Y["results"],Z["cases"],strict=True):
        assert a["q"]==b["q"]==c["q"]
        E=a["controls"][-1]["new_available_energy_upper_GeV4"]
        cases=[]
        for r,s in zip(b["cases"],c["cases"],strict=True):
            delta=r["cap"]
            assert delta==s["light_cap"]
            ext=delta*(1-delta)*r["total_bounded_input_over_outgoing_decay_margin"]
            rows=[]
            for old in s["checked_subbands"]:
                pc=old["comoving_cut_GeV"]
                h=np.hypot(M,pc)
                Q=[]
                for n1,n2 in ((400,40),(1200,120),(2400,240)):
                    u,e,share=inverse_bound(E,delta,pc,n1,n2)
                    # Q_inv/(Gamma M)=g M/(4 pi^2 P)*u; Q_other/(Gamma M)
                    # integrates a uniform per-mode bound over this band.
                    val=h/(1-delta)/E*(M*u/(np.pi**2*P)+2*pc**3*ext/(3*np.pi**2))
                    Q.append(dict(scalar_energy_cells=n1,density_split_cells=n2,
                                  below_cut_heavy_energy_upper_over_E=float(val),
                                  worst_incoming_scalar_energy_GeV=e,worst_light_energy_fraction=share))
                rows.append(dict(cut_GeV=pc,pointwise_bound=old["below_cut_energy_fraction_upper"],
                                 aggregate=Q,
                                 strongest_bound_for_below_cut=min(1.,old["below_cut_energy_fraction_upper"],
                                              Q[-1]["below_cut_heavy_energy_upper_over_E"])))
            cases.append(dict(light_cap=delta,controls=rows))
        out.append(dict(q=a["q"],rows=cases))
    return dict(scope="Initially empty finite comoving heavy band and no unbounded heavy "
                      "elastic exchange. Exact inverse kinematic thresholds and light Pauli "
                      "cap plus a total energy allowance, scalar per-mode occupancy unbounded. "
                      "No bound on the complementary above-cut band.",
                exact_phase_normalization_relative_defects=normalization_check(),cases=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("pions","barrier","low"):
        p.add_argument(k)
    a=p.parse_args()
    print(json.dumps(run(a.pions,a.barrier,a.low),indent=2,allow_nan=False))
