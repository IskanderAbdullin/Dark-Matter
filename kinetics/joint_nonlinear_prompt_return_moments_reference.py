"""Joint scalar/fermion instantaneous number and energy for prompt resonant return."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre, expit
from factorized_prompt_return_graph_reference import C
from sm_dark_entropy_exchange_reference import StandardModel,MPL
from three_species_decay_onsager_reference import M,m,mp,P,Es

def geometry(p,q,Emax,order=6):
    p,q=np.asarray(p),np.asarray(q)
    u=np.hypot(p,mp)
    v=np.hypot(q,m)
    def changes(k,E,mi,E0):
        d=M*(E0*E0+k*k)/(E*E0+k*P)
        return np.r_[d,d+2*M*k*P/(mi*mi)]
    # Lab partition changes when either cell boundary is crossed by a
    # boosted endpoint, or when two constituent cell boundaries swap.
    cut=np.unique(np.clip(np.r_[M,Emax,
             changes(p,u,mp,Es),changes(q,v,m,M-Es),(u[:,None]+v).ravel(),
             M+np.geomspace(1e-9,Emax-M,46)],M,Emax))
    x,w=roots_legendre(order)
    l,r=np.log(cut[:-1]),np.log(cut[1:])
    E=np.exp((l[:,None]*(1-x)+r[:,None]*(1+x))/2).ravel()
    dW=((r-l)[:,None]*w/2).ravel()*E*C
    k=np.sqrt((E-M)*(E+M))
    up=(E*Es+k*P)/M
    lo=(Es*Es+(k*mp/M)**2)/up
    boundaries=np.concatenate(
         [np.broadcast_to(u,(len(E),len(u))),E[:,None]-v],axis=1)
    z=np.concatenate([np.clip(boundaries,lo[:,None],up[:,None]),lo[:,None],up[:,None]],axis=1)
    z.sort(axis=1)
    a,b=z[:,:-1],z[:,1:]
    mid=(a+b)/2
    xlen=b-a
    ip=np.clip(np.searchsorted(u,mid,side='right')-1,0,len(u)-1)
    iq=np.clip(np.searchsorted(v,E[:,None]-mid,side='right')-1,0,len(v)-1)
    return E,dW,mid,xlen,ip,iq

def eval_moments(g,F,G,*,thermal=None,stimulated=True):
    E,w,P,x,ip,iq=g
    if thermal is None:
        phi,chi=F[ip],G[iq]
    else:
        T,eta_phi,eta_chi=thermal
        phi=1/np.expm1(P/T-eta_phi)
        chi=expit(eta_chi-(E[:,None]-P)/T)
    up=x*phi*chi
    down=x*(1-chi)*(1+phi if stimulated else 1.)
    A=up.sum(axis=1)
    B=down.sum(axis=1)
    total=A+B
    h=np.divide(A,total,out=np.zeros_like(A),where=total>0)
    plus=w[:,None]*down*h[:,None]
    minus=w[:,None]*up*(1-h[:,None])
    net=plus-minus
    mag=np.sum(abs(plus)+abs(minus))
    scal=np.bincount(ip.ravel(),weights=net.ravel(),minlength=len(F))
    ferm=np.bincount(iq.ravel(),weights=net.ravel(),minlength=len(G))
    sq=np.bincount(ip.ravel(),weights=(net*P).ravel(),minlength=len(F))
    fq=np.bincount(iq.ravel(),weights=(net*(E[:,None]-P)).ravel(),minlength=len(G))
    assert max(abs(scal.sum()),abs(ferm.sum()))<3e-13*max(mag,1e-260)
    assert abs(sq.sum()+fq.sum())<3e-13*max(mag*max(E),1e-260)
    return scal,ferm,sq,fq,float(mag),float(h.max())

def rows(eos,data):
    out=[]
    for r in data['rows']:
        if r['TR'] not in (.02,.05,.1) or abs(r['TF']/r['TR']-.9)>1e-6:continue
        T=r['TF']
        ss=eos.state(T)[2]
        H=np.sqrt(eos.state(T)[0]/3)/MPL
        pp=np.array(r["scalar_differential_bins"])*T
        results=[]
        for step in (4,3):
            ia=np.unique(np.r_[np.arange(0,len(pp),step),len(pp)-1])
            u=pp[ia]
            Vol=np.diff(u**3)/(6*np.pi**2)
            nphi=ss*r['stable_yields'][1]
            nchi=ss*sum(r['stable_yields'])
            phi=nphi*np.add.reduceat(r["scalar_normalized_counts"],ia[:-1])
            chi=nchi*np.add.reduceat(r["stable_normalized_counts"],ia[:-1])
            F=np.r_[phi/Vol,0.]
            # Sum ep/anti and two spins; source references use this convention.
            G=np.r_[chi/(4*Vol),0.]
            assert max(G)<1.
            g=geometry(u,u,2*u[-1]+M)
            s,k,a,b,n,ph=eval_moments(g,F,G)
            z=eval_moments(g,F,G,thermal=(T,-.15,-6.))
            defect=(np.sum(abs(z[0]))+np.sum(abs(z[1])))/z[4]
            assert defect<2e-14
            void=F.copy();void[1::4]=0.
            v=G.copy();v[1::4]=0.;v[2::4]=1.
            sv,kv,*_=eval_moments(g,void,v)
            assert np.all(sv[1::4]>=0.) and np.all(kv[1::4]>=0.) and np.all(kv[2::4]<=0.)
            rphi=np.array([np.hypot(mom,mp) for mom in np.r_[np.cbrt((u[1:]**3+u[:-1]**3)/2),u[-1]]])
            rchi=np.hypot(np.r_[np.cbrt((u[1:]**3+u[:-1]**3)/2),u[-1]],m)
            # Deliberately report why advancing only the numbers with a
            # representative node energy is not conservative.
            defect_node=(rphi@s+rchi@k)/max(abs(a.sum()),1e-250)
            echi=ss*(np.array(r["stable_yields"])@np.array(r["redshifted_means"][:2]))
            results.append(dict(bins=len(F)-1,Q_phi=float(a.sum()),Q_chi=float(b.sum()),
                 Q_phi_over_chi_energy_per_H=float(a.sum()/(echi*H)),thermal_binwise_defect=defect,
                 normalized_number_only_energy_defect=float(defect_node),max_intermediate_parent_f=ph))
        err=abs(results[0]['Q_phi']/results[1]['Q_phi']-1)
        assert err<.5,err
        out.append(dict(TR=r['TR'],T=T,controls=results,relative_heat_mesh_shift=err))
    return dict(rows=out,scope="Instantaneous joint projected scalar/fermion prompt-resonant nonlinear operator with exact event energy; not a time evolution or physical abundance.")

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("eos");p.add_argument("spectra")
    a=p.parse_args()
    print(json.dumps(rows(StandardModel(a.eos),json.loads(Path(a.spectra).read_text())),indent=2,allow_nan=False))
