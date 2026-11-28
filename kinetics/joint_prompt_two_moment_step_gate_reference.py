"""Positive instantaneous number/energy step and reconstruction gate."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.special import expit, roots_legendre
from sm_dark_entropy_exchange_reference import StandardModel,MPL
from joint_nonlinear_prompt_return_moments_reference import geometry,eval_moments
from fermion_two_moment_closure_reference import infer as f_infer,mesh
from regular_bose_two_moment_closure_reference import infer as b_infer,occupation

from three_species_decay_onsager_reference import M,m,mp

def kinetic(p,m):
    return p*p/(np.hypot(p,m)+m)

def test(r,eos,step):
    T=r['TF'];S=eos.state(T)[2];H=np.sqrt(eos.state(T)[0]/3)/MPL
    ed=np.array(r['scalar_differential_bins'])*T
    ids=np.unique(np.r_[np.arange(0,len(ed),step),len(ed)-1])
    p=ed[ids]
    V=np.diff(p**3)/(6*np.pi**2)
    Ns=np.r_[S*r['stable_yields'][1]*np.add.reduceat(r['scalar_normalized_counts'],ids[:-1]),0.]
    Nf=np.r_[S*sum(r['stable_yields'])*np.add.reduceat(r['stable_normalized_counts'],ids[:-1]),0.]
    F=Ns/np.r_[V,1.]
    G=Nf/np.r_[4*V,1.]
    q=geometry(p,p,2*p[-1]+M)
    Ds,Df,Qs,Qf,_,_=eval_moments(q,F,G)
    Dks=Qs-mp*Ds
    Dkf=Qf-m*Df
    z,w=roots_legendre(44)
    x=(p[:-1,None]*(1-z)+p[1:,None]*(1+z))/2
    W=w*x*x*(p[1:]-p[:-1])[:,None]/2/(2*np.pi**2)
    K=[]; caps=[]
    for n,d,j,mi,g in ((Ns,Ds,Dks,mp,1),(Nf,Df,Dkf,m,4)):
        k=(W*kinetic(x,mi)).sum(axis=1)*n[:-1]/V
        K.append(np.r_[k,0.])
        low=kinetic(p[:-1],mi)
        high=kinetic(p[1:],mi)
        R=np.r_[n[:-1],k-low*n[:-1],high*n[:-1]-k]
        Rd=np.r_[d[:-1],j[:-1]-low*d[:-1],high*d[:-1]-j[:-1]]
        if g==4:
            R=np.r_[R,g*V-n[:-1]]
            Rd=np.r_[Rd,-d[:-1]]
        caps.extend((R[Rd<0]/-Rd[Rd<0]).tolist())
    dt=.12*min(a for a in caps if a>0)
    bad=0
    for a in range(12):
        trial=[]
        try:
            for n,d,e,j,mi,g,fn in ((Ns,Ds,K[0],Dks,mp,1,b_infer),(Nf,Df,K[1],Dkf,m,4,f_infer)):
                u=n+dt*d
                v=e+dt*j
                assert np.all(u>=-1e-22)
                for i,x0 in enumerate(u[:-1]):
                    if x0<=1e-240:continue
                    lo,hi=p[i:i+2]
                    t,w0,l,h=mesh(lo,hi,mi,96)
                    mean=(v[i]/x0-l)/(h-l)
                    fill=x0/(g*V[i])
                    if g==4:
                        mu=fill*(hi-lo)*(hi**2+lo*hi+lo**2)
                        yy,wg=roots_legendre(30)
                        lowq=np.cbrt(lo**3+mu*(1+yy)/2)
                        highq=np.cbrt(hi**3-mu*(1+yy)/2)
                        x=(wg@kinetic(lowq,mi))/2
                        y=(wg@kinetic(highq,mi))/2
                        e_now=v[i]/x0
                        margin=2e-11*max((h-l),1e-230)
                        if not x-margin<=e_now<=y+margin:
                            raise ValueError((i,g,"packing",x,y,e_now,fill))
                    try:
                        A,B=fn(lo,hi,mi,fill,mean)
                    except ValueError as exc:
                        raise ValueError((i,g,fill,mean,str(exc))) from exc
                    f=occupation(A,B,t) if g==1 else expit(A-B*t)
                    err=max(abs((w0@f)/fill-1),abs((w0@(t*f))/(fill*mean)-1))
                    if err>3e-7:raise ValueError('failed numerical moment recovery')
                    trial.append(err)
            break
        except (ValueError,AssertionError) as exc:
            blocker=str(exc);dt*=.25;bad+=1
    else:raise ArithmeticError(f'no admissible trial step passed the moment inversion: {blocker}')
    delta=dt*(Qs.sum()+Qf.sum())
    ratio=abs(delta)/(dt*max(abs(Qs.sum()),1e-260))
    assert ratio<2e-12
    return dict(TR=r['TR'],T=T,mesh_bins=len(p)-1, accepted_step_in_H_units=dt*H,
                backtracks=bad,maximum_moment_recovery_error=max(trial),
                physical_event_energy_relative_conservation=ratio,
                exterior_number_after_trial=(float(dt*Ds[-1]),float(dt*Df[-1])))
if __name__=="__main__":
    a=argparse.ArgumentParser();a.add_argument('eos');a.add_argument('spectra')
    b=a.parse_args();eos=StandardModel(b.eos)
    data=json.loads(Path(b.spectra).read_text())
    r=[z for z in data['rows'] if z['TR']==.05 and abs(z['TF']-.045)<1e-9][0]
    print(json.dumps(dict(rows=[test(r,eos,j) for j in (3,4)],
        scope="Single conservative snapshot trial; no long-time algorithm or tail reconstruction."),
        indent=2,allow_nan=False))
