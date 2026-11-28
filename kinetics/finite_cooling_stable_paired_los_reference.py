"""Cold-matched relativistic Vlasov response of the 20 MeV stable tracer."""
import argparse,json,re
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline,PchipInterpolator
from scipy.integrate import cumulative_simpson,simpson
from scipy.special import spherical_jn
from archived_collisionless_streaming_reference import Background
from finite_cooling_pair_split_tracer_reference import indices
from finite_cooling_matched_history_los_reference import sha
from three_species_decay_onsager_reference import m

def cumulative(y,t):
    return cumulative_simpson(y[...,::-1],x=-t[::-1],axis=-1,initial=0.)[...,::-1]


def bessel_deltas(x):
    q=x*x
    j0=np.where(abs(x)<.08,q*(-1./6+q*(1./120+q*(-1./5040+q/362880.))),spherical_jn(0,x)-1.)
    raw=np.divide(spherical_jn(1,x),x,out=np.full_like(x,1./3),where=x!=0.)
    j1=np.where(abs(x)<.08,q*(-1./30+q*(1./840+q*(-1./45360+q/3991680.))),raw-1./3)
    return j0,j1


def run(eos,tables,metrics,mass_factor=1.):
    assert mass_factor in (1.,2.)
    mass=m*mass_factor
    bg=Background(eos)
    raw=[np.loadtxt(p).T for p in tables]
    sp=[CubicSpline(q,f,extrapolate=False) for q,f in raw]
    lo=min(np.log(q[0]) for q,f in raw)
    hi=max(np.log(q[-1]) for q,f in raw)
    answer=[]
    for metric in metrics:
        k=float(re.search(r"k = (\S+)",metric.read_text().splitlines()[0])[1])
        I=indices(metric,"delta_cdm")
        M=np.loadtxt(metric,usecols=[I[x] for x in ("tau","a","psi","phi","delta_cdm","theta_cdm","delta_b")])
        _,j=np.unique(M[:,0],return_index=True)
        M=M[j]
        assert 0<k*M[0,0]<.1
        P=PchipInterpolator(M[:,0],M[:,1:],axis=0)
        tf=float(PchipInterpolator(np.log(M[:,1]),M[:,0])(np.log(1/21)))
        native=P(tf)[3]
        out=[]
        for nq,nt in ((460,9000),(850,18000)):
            for st in (M[0,0],.15/k):
                t=np.geomspace(st,tf,nt)
                a,psi,phi,cdm,theta,b=P(t).T
                phid=P.derivative()(t)[:,2]
                iv=cumulative_simpson((1./a)[::-1],x=-t[::-1],initial=0.)[::-1]
                cold=cdm[0]+3*(phi[-1]-phi[0])-theta[0]*a[0]*iv[0]-k*k*simpson(a*psi*iv,x=t)
                assert abs(cold/native-1.)<1e-4
                A=a*mass/bg.T0
                q=np.exp(np.linspace(lo,hi,nq))
                q0=np.hypot(q,A[0]); qf=np.hypot(q,A[-1])
                D=[];N=[];C=[];Fcoef=[];bd=[]
                for z in sp:
                    f=np.nan_to_num(z(q),nan=0.)
                    fp=q*np.nan_to_num(z(q,1),nan=0.)
                    w=q*q*f
                    ni=np.trapezoid(w*q0,q)
                    nf=np.trapezoid(w*qf,q)
                    di=-3.-np.trapezoid(w*q*q/q0,q)/ni
                    df=-3.-np.trapezoid(w*q*q/qf,q)/nf
                    edge=np.array([q[0]**3*qf[0]*f[0],q[-1]**3*qf[-1]*f[-1]])
                    bd.append(float(np.max(edge)/nf))
                    assert bd[-1]<2e-9
                    D.append(q*q*qf*fp)
                    N.append(nf)
                    C.append(cdm[0]/di)
                    Fcoef.append(df)
                G=[];Z=[]
                for s in range(0,nq,16):
                    xx=q[s:s+16,None]
                    eps=np.hypot(xx,A)
                    de=xx**2/(eps+A)
                    eta=cumulative(1./eps,t)
                    deta=cumulative(-de/(A*eps),t)
                    chi=k*xx*eta
                    Dj0,Dj1=bessel_deltas(chi)
                    Dj1v=k*((de*eta+A*deta)/3.+eps*eta*Dj1)
                    dif=simpson(-phid*Dj0+k*psi*Dj1v,x=t,axis=1)
                    G.append(dif+theta[0]/k*Dj1v[:,0])
                    Z.append(Dj0[:,0])
                dG=np.concatenate(G)
                dJ=np.concatenate(Z)
                G0=(cdm[0]-cold)/3.
                diffs=[]
                for w,n,c,df in zip(D,N,C,Fcoef):
                    # Analytic q integration of the cold part, plus only the
                    # paired residual from momentum/time quadrature.
                    delta=df*(c+G0)-cold
                    delta+=np.trapezoid(w*(dG+c*dJ),q)/n
                    diffs.append(float(delta/cold))
                out.append(dict(nodes=[nq,nt],start_k_tau=k*st,cold_over_native=cold/native,
                    kinetic_minus_cold_over_cold=diffs,cold_delta=cold,native_cdm_delta=native,
                    exact_baryon_delta=float(P(tf)[5]),mass_eV=mass*1e9,weighted_tail_boundary_defect=bd))
        differences=np.array([x["kinetic_minus_cold_over_cold"] for x in out])
        assert np.max(abs(differences[0]-differences[2]))<2e-10
        assert np.max(abs(differences[1]-differences[3]))<2e-10
        answer.append(dict(k_h_Mpc=k/(bg.H0/100),metric_sha256=sha(metric),data=out))
    return dict(results=answer,spectrum_sha256=[sha(p) for p in tables],eos_sha256=sha(eos),
        scope="Fixed Newtonian metric, stable test-particle density and velocity matched to CDM initially. The numerical cold comparator follows the same metric. No prehistory, gravitational feedback, absolute matter transfer, or likelihood.")

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument('eos',type=Path);p.add_argument('tables',type=Path,nargs=2);p.add_argument('metrics',type=Path,nargs='+')
    p.add_argument("--mass-factor",type=float,default=1.)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.tables,a.metrics,a.mass_factor),indent=2,allow_nan=False))
