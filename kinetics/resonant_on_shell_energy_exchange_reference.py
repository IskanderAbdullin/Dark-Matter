"""On-shell scalar heat ledger and independent exact thermal-tilt quadrature."""
import argparse,json
from pathlib import Path
from functools import lru_cache
import numpy as np
from scipy.special import roots_legendre, expit
from three_species_decay_onsager_reference import M,m,mp,P,Es,width
from factorized_prompt_return_graph_reference import graph,C
from tagged_scalar_nonthermal_bath_reference import Bath
from sm_dark_entropy_exchange_reference import StandardModel,MPL

@lru_cache(None)
def quadr(n):
    return roots_legendre(n)

def thermal_heat(Tphi,Tchi,aphi,achi,order=(70,42),tail=55.):
    """Resonant elastic heat Q_phi for MB scalar and FD fermion. Opposite heat on chi."""
    assert Tphi>0 and Tchi>0
    if Tphi==Tchi:return 0.
    z,w=quadr(order[0])
    x,v=quadr(order[1])
    th=max(Tphi,Tchi)
    a=(z+1)*np.sqrt(tail)/2
    E1=M+th*a*a
    W=th*a*np.sqrt(tail)*w
    P1=a*np.sqrt(th*(E1+M))
    top=(E1*Es+P1*P)/M
    bot=(Es*Es+(P1*mp/M)**2)/top
    # Integration in log scalar energy handles ever wider moving boosts.
    ll=np.log(bot)
    k=np.log(top/bot)
    E=np.exp(ll[:,None]+k[:,None]*(1+x)/2)
    wt=k[:,None]*E*(v/2)
    b=wt*expit((E1[:,None]-E)/Tchi-achi)
    total=b.sum(axis=1)
    mean=(b*E).sum(axis=1)/total
    d=(Tphi-Tchi)/(Tchi*Tphi)
    y=E-mean[:,None]
    u=d*y
    num=np.empty_like(mean)
    near=np.max(abs(u),axis=1)<.9
    # Subtract the zero covariance at equal temperature analytically.
    num[near]=np.sum(b[near]*y[near]*np.expm1(u[near]),axis=1)
    # Distant temperature pairs retain a bounded and stable direct average.
    mx=np.max(u[~near],axis=1,initial=-np.inf)
    num[~near]=np.sum(b[~near]*y[~near]*np.exp(u[~near]-mx[:,None]),axis=1)
    s=np.zeros(len(mean));s[~near]=mx
    raw=-C*np.sum(W*np.exp(aphi+achi-E1/Tchi+d*mean+s)*num)
    assert raw*(Tchi-Tphi)>=0.,(Tphi,Tchi,raw)
    return raw

def test(eos,spectra):
    output=[]
    for T in (.018,.045,.09):
        for ratio in (.65,.93,1.,1.0000001,1.12,1.5):
            q=thermal_heat(T*ratio,T,-8.,-9.)
            ref=thermal_heat(T*ratio,T,-8.,-9.,order=(110,70),tail=67.)
            err=abs(q/ref-1) if ref else abs(q)
            assert err<2e-7,err
            if ratio==1.: assert q==0.
            output.append(dict(Tchi=T,Tphi=T*ratio,Q_phi=q,quad_refinement=err,
                               entropy_density_per_time=q*(1/(T*ratio)-1/T)))
    # Independent two-angle energy matrix using MB scalar shell trial.
    E=[]
    for T in (.018,.045,.09):
        ed=np.r_[0.,np.geomspace(.5*mp,30*T,26)]
        R,Ein,Eout=graph(ed,T=T,eta=-9.,trial_T=.8*T,energy=True)
        re=thermal_heat(.8*T,T,0.,-9.)
        q=(Eout-Ein).sum()
        error=abs(q/re-1)
        assert error<5e-6,error
        R,Ein,Eout=graph(ed,T=T,eta=-9.,trial_T=T,energy=True)
        rec=np.linalg.norm(Ein[:,:-1]-Eout[:,:-1].T)/np.linalg.norm(Ein[:,:-1])
        assert rec<1e-13
        E.append(dict(T=T,event_energy_reverse_relative=rec,tilt_vs_event_relative_error=error))
    # Mixed nonthermal snapshots, using the actual injected scalar masses per cell.
    mix=[]
    for r in spectra["rows"]:
        if r["TR"] not in (.02,.05,.1) or abs(r["TF"]/r["TR"]-.9)>1e-6:continue
        T=r["TF"]
        b=Bath(r,eos)
        edge=np.array(r["scalar_differential_bins"])*T
        g,u,v=graph(edge,T=T,bath=b,energy=True)
        V=np.diff(edge**3)/(6*np.pi**2)
        phi=r["stable_yields"][1]*eos.state(T)[2]*np.array(r["scalar_normalized_counts"])
        f=phi/V
        q=f@np.sum(v-u,axis=1)
        q_sep=f@v.sum(axis=1)-f@u.sum(axis=1)
        assert abs(q/q_sep-1)<1e-12
        physical_qf=f@((u-v).sum(axis=1))
        # Verify the heat with higher parent and outgoing-angle degrees and
        # a deliberately different projection of the measured scalar spectrum.
        _,a2,b2=graph(edge,T=T,bath=b,energy=True,order=21,theta=24)
        refine=f@np.sum(b2-a2,axis=1)
        mask=np.unique(np.r_[np.arange(0,len(edge),2),len(edge)-1])
        e2=edge[mask]
        z=np.add.reduceat(phi,mask[:-1])
        v2=np.diff(e2**3)/(6*np.pi**2)
        _,a3,b3=graph(e2,T=T,bath=b,energy=True)
        coarse=(z/v2)@np.sum(b3-a3,axis=1)
        assert abs(refine/q-1)<.02 and abs(coarse/q-1)<.2,(T,q,refine,coarse)
        W,A,B,father=graph(edge,T=T,bath=b,components=True,parents=True)
        ff=np.r_[f,0.]
        capture=A@ff
        released=B@(1+ff)
        steady=capture/(capture+released)
        dE=W/C
        par=(2/np.pi**2)*dE*father*np.sqrt((father-M)*(father+M))*steady
        H=np.sqrt(eos.state(T)[0]/3)/MPL
        s0=eos.state(T)[2]
        nn=r["stable_yields"][1]*s0
        rr=(np.array(r["stable_yields"])@np.array(r["redshifted_means"][:2]))*s0
        ss=nn*r["redshifted_means"][2]
        Eavg=float((father*par).sum()/par.sum()) if par.sum()>0 else M
        clock=width*M/Eavg
        mix.append(dict(TR=r["TR"],T=T,Q_phi=q,Q_chi=physical_qf,total=q+physical_qf,
                        snapshot_heat_over_chi_energy_per_H=q/(rr*H),
                        relative_heat_quadrature_refinement=abs(refine/q-1),
                        relative_heat_coarse_scalar_projection=abs(coarse/q-1),
                        snapshot_heat_over_phi_energy_per_H=q/(ss*H),
                        intermediate_density_relative_to_scalar=par.sum()/nn,
                        intermediate_energy_relative_to_prescribed_chi=float((par*father).sum()/rr),
                        normalized_heat_clock_over_estimated_parent_decay=q/(rr*clock)))
    return dict(thermal=output,independent_event_energy=E,nonthermal_snapshot=mix,
                scope="Lab-event scalar momentum transfer and opposite fermion ledger for prompt on-shell redistribution only; fixed-bath snapshots, not a full plasma or abundance.")

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument('eos');p.add_argument('spectra')
    a=p.parse_args()
    print(json.dumps(test(StandardModel(a.eos),json.loads(Path(a.spectra).read_text())),
            indent=2,allow_nan=False))
