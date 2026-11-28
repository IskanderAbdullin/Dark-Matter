"""Entropy-metric finite-shell spectrum near joint prompt-resonant balance."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.linalg import block_diag,cholesky,solve_triangular,eigh
from scipy.special import expit
from joint_prompt_two_moment_evolution_reference import Collision,maxent
from sm_dark_entropy_exchange_reference import StandardModel,MPL

def configuration(row,eos,step,cap):
    T=row['TF'];ss=eos.state(T)[2];H=np.sqrt(eos.state(T)[0]/3)/MPL
    p=np.array(row['scalar_differential_bins'])*T
    last=np.searchsorted(p,cap*T,side='right')-1
    ed=np.unique(np.r_[np.arange(0,last+1,step),last])
    k=Collision(p[ed])
    N=[ss*y*np.add.reduceat(np.array(row[s])[:last],ed[:-1]) for y,s in
       ((row['stable_yields'][1],'scalar_normalized_counts'),(sum(row['stable_yields']),'stable_normalized_counts'))]
    E=sum(x@(c.mass+c.mean) for x,c in zip(N,k.c))
    t,Q=maxent(k.c,[x.sum() for x in N],E,T)
    theta=[np.array([(c.mass+c.low)/t-Q[j][3],c.dk/t]) for j,c in enumerate(k.c)]
    return T,H,k,N,t,Q,theta

def gram(k,theta):
    n=len(k.c[0].vol);d=4*n
    f=k.c[0].f(theta[0],k.ph,k.I);g=k.c[1].f(theta[1],k.ch,k.J)
    A=k.pre*f*g;B=k.pre*(1+f)*(1-g)
    a=np.bincount(k.y,weights=A,minlength=k.n);b=np.bincount(k.y,weights=B,minlength=k.n)
    w=k.fac*B*np.divide(a,a+b,out=np.zeros_like(a),where=(a+b)>0.)[k.y]
    e=(k.k1-k.c[0].low[k.I])/k.c[0].dk[k.I]; h=(k.k2-k.c[1].low[k.J])/k.c[1].dk[k.J]
    ids=np.stack([2*k.I,2*k.I+1,2*n+2*k.J,2*n+2*k.J+1],axis=1)
    val=np.stack([np.ones_like(e),e,np.ones_like(h),h],axis=1)
    mu=np.zeros((k.n,d));W=np.bincount(k.y,weights=w,minlength=k.n)
    for i in range(4):np.add.at(mu,(k.y,ids[:,i]),w*val[:,i])
    mu=np.divide(mu,W[:,None],out=np.zeros_like(mu),where=W[:,None]>0.)
    G=np.zeros((d,d))
    for j in range(0,len(w),3072):
        q=slice(j,j+3072)
        c=-mu[k.y[q]].copy()
        c[np.arange(len(c))[:,None],ids[q]]+=val[q]
        z=c*np.sqrt(w[q])[:,None];G+=z.T@z
    return G

def metric(k,theta):
    n=len(k.c[0].vol);blocks=[];z=np.zeros((4*n,3));C=[]
    for j,c in enumerate(k.c):
        s=theta[j][0][:,None]+theta[j][1][:,None]*c.u
        f=expit(-s) if j else np.exp(-s)/(-np.expm1(-s))
        w=c.w*c.vol[:,None]*f*(1+(-1 if j else 1)*f)
        E=c.mass+c.low[:,None]+c.dk[:,None]*c.u
        C.append(np.sum(w*(E-np.sum(w*E)/np.sum(w))**2))
        for i in range(n):
            u=np.array([np.ones_like(c.u[i]),c.u[i]])
            blocks.append((u*w[i])@u.T)
        z[2*j*n:2*(j+1)*n:2,j]=1
        z[2*j*n:2*(j+1)*n:2,2]=c.mass+c.low
        z[2*j*n+1:2*(j+1)*n:2,2]=c.dk
    return np.asarray(C),cholesky(block_diag(*blocks),lower=True),z

def run(r,eos,step,cap):
    T,H,k,N,t,Q,th=configuration(r,eos,step,cap)
    n=len(k.c[0].vol);d=4*n
    G=gram(k,th)
    Cv,L,z=metric(k,th)
    B=solve_triangular(L,solve_triangular(L,G,lower=True).T,lower=True).T
    ev,V=eigh((B+B.T)/2)
    ev/=H
    null=np.linalg.norm(G@z)/(np.linalg.norm(G)*np.linalg.norm(z))
    assert null<1e-13 and ev[0]>=-1e-11*max(ev)
    n0=np.linalg.qr(L.T@z)[0]
    def signal(h):
        x=L.T@h;x-=n0@(n0.T@x)
        return x/np.linalg.norm(x)
    def inspect(h):
        x=signal(h);s=(V.T@x)**2
        assert abs(s.sum()-1)<1e-12
        return dict(rayleigh_over_H=float(x@B@x/H),
             weight_below_H=float(s[ev<1.].sum()),weight_above_1e6_H=float(s[ev>=1e6].sum()))
    v=z[:,2].copy();v[2*n:]*=-1.
    mac=inspect(v)
    from dark_resonant_thermal_slow_manifold_reference import conductance
    expectation=conductance(t,Q[1][3],Q[0][3])*t**2*np.sum(1/Cv)/H
    local={}
    for p0 in (.005,.05,.2,.5,1.,3.):
        i=np.searchsorted(k.c[0].p,p0*T,side='right')-1
        h=np.zeros(d);h[2*i]=1.;local[str(p0)]=inspect(h)
    q=np.random.default_rng(8431).uniform(-1.,1.,d)
    rs=[q[j*2*n:(j+1)*2*n].reshape(n,2).T for j in (0,1)]
    def rhs(t):
        y=k.rate([s+t*x for s,x in zip(th,rs)])
        for j,c in enumerate(k.c):y[j,1]=(y[j,1]-c.low*y[j,0])/c.dk
        return y.transpose(0,2,1).ravel()
    dz=1e-5
    check=np.linalg.norm((rhs(dz)-rhs(-dz))/(2*dz)-G@q)/np.linalg.norm(G@q)
    assert check<3e-7,check
    mask=ev>max(ev)*1e-12
    return dict(bins=n,p_cap_over_visible_T=float(k.c[0].p[-1]/T),common_T=t,
        number_and_energy_null_error=null,independent_finite_difference_relative_error=check,
        thermal_contrast=mac,independent_untruncated_thermal_contrast_estimate=expectation,
        scalar_shell_modes=local,resolved_eigenvalue_quantiles_over_H=np.quantile(ev[mask],[0,.25,.5,.75,1]).tolist(),
        count_below_relative_1e_minus_12_cut=int((~mask).sum()))

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument('eos');p.add_argument('spectra')
    a=p.parse_args();eos=StandardModel(a.eos)
    r=next(r for r in json.loads(Path(a.spectra).read_text())['rows'] if r['TR']==.05 and abs(r['TF']-.045)<1e-8)
    print(json.dumps(dict(rows=[run(r,eos,j,k) for j,k in ((4,12),(3,12),(3,24))],
          scope="Finite-volume near-equilibrium kinetic diagnostics in one fixed support. Projected eigenvalues are not a physical full-spectrum gap."),indent=2,allow_nan=False))
