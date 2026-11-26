"""Joint scalar/fermion equilibrium elastic Gram with explicit invariant s cut."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre,expit
from scipy.linalg import eigh,solve_triangular
from soft_scalar_elastic_graph_reference import spin
from three_species_decay_onsager_reference import M,m,mp
from joint_resonant_kinetic_linear_modes_reference import configuration,gram,metric
from sm_dark_entropy_exchange_reference import StandardModel

def setup(T,k,temps,ef,n_u,n_l,n_a,n_az,tail,pole,perturb=None,cache=None,cache_only=False):
    assert not cache_only or (cache is not None and perturb is None)
    n=len(k.c[0].vol);D=4*n
    el=[np.hypot(c.p,c.mass) for c in k.c]
    tau=temps
    g=np.zeros((D,D))
    nonlin=np.zeros((0 if perturb is None else len(perturb),D))
    entropy=np.zeros(len(nonlin))
    rs,ws=roots_legendre(n_u);rl,wl=roots_legendre(n_l);ra,wa=roots_legendre(n_a)
    s0=m+mp;s1=np.sqrt(pole)*M
    br=np.log(np.r_[1e-13,np.geomspace(1e-11,s1-s0,16)])
    h=(br[1:]-br[:-1])[:,None]/2
    E0=np.exp((br[1:]+br[:-1])[:,None]/2+h*rs)
    az=np.cos((np.arange(n_az)+.5)*2*np.pi/n_az)
    for p0,d in zip(E0.ravel(),(h*ws*E0).ravel()):
        u=s0+p0
        s=u*u
        pcm=np.sqrt(p0*(u+s0)*(u-m+mp)*(u+m-mp))/(2*u)
        E=np.hypot(mp,pcm)
        cap=u+T*tail
        gates=[]
        for j,c in enumerate(k.c):
            e0=E if j==0 else u-E
            cc=el[j]*e0+c.p*pcm
            low=u*(e0*e0+c.p*c.p)/cc
            high=u*cc/c.mass**2
            gates.extend(np.r_[low,high])
        l=np.unique(np.clip(np.r_[u,cap,gates,u+T*np.array([.03,.12,.4,.8,1.4,2.5,4.5,7.,11.,17.,25.])],u,cap))
        h=(l[1:]-l[:-1])[:,None]/2
        labs=(l[1:]+l[:-1])[:,None]/2+h*rl
        for L,dL in zip(labs.ravel(),(h*wl).ravel()):
            q=np.sqrt((L-u)*(L+u))
            B=q*pcm/u
            if B<=0:continue
            A=L*E/u
            a=(E*E+(q*mp/u)**2)/(A+B)
            b=A+B
            lo=max(a,el[0][0],L-el[1][-1])
            hi=min(b,el[0][-1],L-el[1][0])
            if lo>=hi:continue
            edge=np.unique(np.clip(np.r_[lo,hi,el[0],L-el[1]],lo,hi))
            lf=edge[:-1]
            delta=np.log1p((edge[1:]-lf)/lf)
            p=(lf[:,None]*np.exp(delta[:,None]*(1+ra)/2)).ravel()
            w=(delta[:,None]*wa/2).ravel()*p/B
            f=L-p
            i=np.clip(np.searchsorted(el[0],p,side='right')-1,0,n-1)
            j=np.clip(np.searchsorted(el[1],f,side='right')-1,0,n-1)
            v=np.clip((p-a)/B-1.,-1.,1.);z=np.sqrt(np.maximum(0.,1-v*v))
            cos=v[:,None]*v+z[:,None]*z*az[:,None,None]
            amp=spin(s,-2*pcm*pcm*(1-cos)).mean(axis=0)
            ang=w.copy() if perturb is not None or cache is not None else None

            F=np.zeros((len(p),D));aa=np.arange(len(p))
            F[aa,2*i]=1.;F[aa,2*n+2*j]=1.
            F[aa,2*i+1]=(p-mp-k.c[0].low[i])/k.c[0].dk[i]
            F[aa,2*n+2*j+1]=(f-m-k.c[1].low[j])/k.c[1].dk[j]
            coef=d*dL*q*pcm**2/u
            if cache is not None:
                cache.append((p,f,L,F,amp*ang[:,None]*ang*(4*coef/(256*np.pi**5))))
            if cache_only: continue
            w*=expit(f/tau-ef[1])/(-np.expm1(ef[0]-p/tau))
            WW=amp*w[:,None]*w*np.exp(ef[0]+ef[1]-L/tau)
            g+=coef*(F.T@(WW.sum(axis=1)[:,None]*F-WW@F))
            if perturb is not None:
                for z,field in enumerate(perturb):
                    a=F[:,:2*n]@field[:2*n]
                    b=F[:,2*n:]@field[2*n:]
                    aa=(p/tau-ef[0])+a
                    bb=(f/tau-ef[1])+b
                    if np.min(aa)<=0:raise ValueError("scalar normal BE domain")
                    w=ang*expit(bb)/(-np.expm1(-aa))
                    G=amp*w[:,None]*w*np.exp(ef[0]+ef[1]-L/tau)
                    x=np.expm1(-a-b)
                    rates=G@x-x*G.sum(axis=1)
                    nonlin[z]+=coef*F.T@rates
                    diff=x[:,None]-x
                    entropy[z]+=coef*.5*np.sum(G*diff*((a+b)[None,:]-(a+b)[:,None]))
    if cache_only: return None
    scale=4/(256*np.pi**5)
    mat=(g+g.T)*scale/2
    return mat if perturb is None else (mat,nonlin*scale,entropy*scale)

def eval_c(r,eos,order,pole,mesh_step=5):
    T,H,k,N,t,Q,th=configuration(r,eos,mesh_step,8.0)
    n=len(k.c[0].vol);d=4*n
    R=gram(k,th);C,L,z=metric(k,th)
    E=setup(T,k,t,[j[3] for j in Q],*order,29.,pole)
    nul=np.linalg.norm(E@z)/(np.linalg.norm(E)*np.linalg.norm(z))
    assert nul<2e-12
    def whiten(G):return solve_triangular(L,solve_triangular(L,G,lower=True).T,lower=True).T
    Br=whiten(R);Be=whiten(E);B=Br+Be
    ev,P=eigh((B+B.T)/2)
    ee=eigh(Be,eigvals_only=True)
    assert min(ee)>-2e-11*max(ee)
    V=np.linalg.qr(L.T@z)[0]
    er,Ur=eigh((Br+Br.T)/2)
    el,Ul=eigh((Be+Be.T)/2)
    def apply(z,values,vectors,h):
        return vectors@(np.exp(-np.maximum(values,0.)*h/H)*(vectors.T@z))
    def time_test(x):
        y=L.T@x;y-=V@(V.T@y);y/=np.linalg.norm(y)
        u=apply(y,ev,P,.02)
        errs=[]
        overlaps=[]
        for nseg in (2,4,8,16):
            h=.02/nseg
            a=y.copy()
            for _ in range(nseg):
                a=apply(a,er,Ur,.5*h)
                a=apply(a,el,Ul,h)
                a=apply(a,er,Ur,.5*h)
            errs.append(float(np.linalg.norm(a-u)))
            overlaps.append(float(y@a))
        return dict(resonant_only_initial_overlap=float(y@apply(y,er,Ur,.02)),
                    exact_combined_overlap_at_point_02_H_time=float(y@u),
                    splitting_absolute_errors=errs,splitting_overlaps=overlaps)
    def mode(x):
        k=L.T@x;k-=V@(V.T@k);k/=np.linalg.norm(k)
        g=(P.T@k)**2
        return dict(resonant_over_H=float(k@Br@k/H),elastic_over_H=float(k@Be@k/H),
                    weight_below_H=float(g[ev<H].sum()))
    ans={}
    for a in (.005,.05,.2,1.,3.):
        i=np.searchsorted(k.c[0].p,a*T,side="right")-1
        u=np.zeros(d);u[2*i]=1.;ans[str(a)]=mode(u)
    e=z[:,2].copy();e[2*n:]*=-1.
    i=np.searchsorted(k.c[0].p,.005*T,side="right")-1
    low=np.zeros(d);low[2*i]=1.
    split=time_test(low)
    return dict(pole_cut_fraction=pole,orders=order,bins=n,null_error=nul,
                thermal_contrast=mode(e),scalar_local=ans,
                cold_linear_split_time_test=split,mesh_step=mesh_step,
                closed_projected_eigenvalues_above_H=int((ev>H).sum()))

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("eos");p.add_argument("spectra")
    a=p.parse_args();eos=StandardModel(a.eos)
    r=next(r for r in json.loads(Path(a.spectra).read_text())['rows'] if r['TR']==.05 and abs(r['TF']-.045)<1e-8)
    print(json.dumps(dict(rows=[eval_c(r,eos,o,k,j) for o,k,j in (((4,3,2,10),.82,5),((6,5,3,16),.82,5),((6,5,3,16),.72,5),((5,4,3,14),.82,4))],
        scope="Closed near-equilibrium two-species kinetic projection; no real-pole double counting, full transport gap, or physical abundance claim."),
        indent=2,allow_nan=False))
