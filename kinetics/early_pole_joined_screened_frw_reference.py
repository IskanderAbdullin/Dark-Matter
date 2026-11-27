"""Joined screened FRW split with a precise finite comoving-energy difference."""
import argparse
import json
from pathlib import Path
import sys

from types import MethodType
import numpy as np
from scipy import __version__ as scipy_version
from scipy.special import roots_legendre,expit
from scipy.integrate import quad,solve_ivp

from early_pole_archived_cold_screened_elastic_reference import Screened
from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_archived_live_compton_reference import ComptonKernel
from early_pole_archived_live_flavor_conversion_reference import ConversionKernel
from early_pole_live_joint_hard_packets_frw_reference import Expanding
from early_pole_screened_elastic_channel_sum_reference import ElasticChannels
from early_pole_frozen_dynamic_vlasov_reference import load
from early_pole_screened_elastic_discrete_gradient_reference import step
from early_pole_cold_screened_split_reference import atomic_write,state_digest

def exterior(M,grid):
    d=load(f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz')
    z,w=roots_legendre(8)
    spec=[]
    for group,name,m in zip(M.edge_sets[1:],('parent','fermion'),(.05,.003),strict=True):
        e=np.array(d[f'endpoint_{name}_momentum_edges'])
        f=np.array(d[f'endpoint_{name}_occupations'])
        def batch(a,b):
            lo=np.maximum(e[:-1],a)
            hi=np.maximum(lo,np.minimum(e[1:],b))
            h=(hi-lo)[:,None]/2
            p=(hi+lo)[:,None]/2+h*z
            wt=h*w*p*p*f[:,None]*2/np.pi**2
            live=wt>0.
            return p[live],wt[live]
        p0,n0=batch(0.,group[0])
        p1,n1=batch(group[-1],np.inf)
        s0,c0=batch(.45,group[0])
        spec.append((np.r_[p0,p1],np.r_[n0,n1],np.r_[s0,p1],np.r_[c0,n1],m))
    lo,hi=M.edge_sets[0][0],M.edge_sets[0][-1]
    fun=lambda x:x**3/np.expm1(x)/np.pi**2 if x>0. else 0.
    ph=quad(fun,0,lo)[0]+quad(fun,hi,80.)[0]
    def at(eta):
        a2=np.exp(2*eta);r=ph;dr=0.;j=[M.background[0]];dj=[0.]
        for p,n,s,c,m in spec:
            ep=np.hypot(p,np.sqrt(a2)*m);es=np.hypot(s,np.sqrt(a2)*m)
            r+=n@ep;dr+=a2*m*m*(n@(1/ep))
            j.append((np.pi**2/2)*(c@(1/es)))
            dj.append(-a2*m*m*(np.pi**2/2)*(c@(1/es**3)))
        return np.array(j),np.array(dj),r,dr
    return at

def finite_energy(M,old,new):
    j=M.spectators(M.eta)[0].copy()
    dj=np.zeros(3);r=0.;a2=np.exp(2*M.eta)
    for flavor,span in enumerate(M.r):
        for k in range(span.start,span.stop):
            C=M.cells[k]
            if flavor==0:
                a,b,*_=C.fit(old[k]);c,d,*_=C.fit(new[k])
                f=C.f(a,b);g=C.f(c,d)
                r+=M.shape[k]@(new[k]-old[k])
                mult=C.W*C.inv
            else:
                a=C.fit(old[k])[0];b=C.fit(new[k])[0]
                f,g=expit(a@C.B),expit(b@C.B)
                bare=np.sqrt(C.energy*C.energy+(a2-1)*C.m*C.m)
                r+=np.sum(C.w*bare*(g-f))
                mult=C.w*np.pi**2/2/bare
            j[flavor]+=mult@f; dj[flavor]+=mult@(g-f)
    v=np.r_[sum(dj[1:]),dj[0]+dj[1:]]
    return float(r+.5*M.c*(2*j+dj)@v)

def run(grid,duration,steps,expands,full_H,projection="steps",checkpoint=None,resume=False,
        checkpoint_every=12,progress_every=0,stop_after=None):
    if duration<=0. or steps<1:
        raise ValueError("positive duration and step count required")
    if resume and not checkpoint:
        raise ValueError("Resume requires a checkpoint path.")
    if checkpoint_every<1 or progress_every<0:
        raise ValueError("Invalid progress or checkpoint interval.")
    if stop_after is not None and (not checkpoint or not 1<=stop_after<=steps):
        raise ValueError("Stopping at a macro boundary requires a path and a valid boundary.")
    if checkpoint and not resume and Path(checkpoint).exists():
        raise FileExistsError("Checkpoint exists; resume explicitly or choose another path.")
    requested=dict(format=1,grid=grid,duration=duration,steps=steps,expands=expands,full_H=full_H,
                   projection=projection,numpy=np.__version__,scipy=scipy_version)
    old=json.loads(Path(checkpoint).read_text()) if resume else None
    if old is not None and any(old.get('physics',{}).get(k)!=v for k,v in requested.items()):
        raise ValueError("FRW checkpoint physics does not match these controls.")
    C=Screened(grid,1,8,duration,(1,1,2),.08,1.,None,projection)
    M=C.M
    M.weights=MethodType(Expanding.weights,M)
    M.energy_increment=MethodType(finite_energy,M)
    M.eta=0.
    M.spectators=exterior(M,grid)
    assert max(abs(M.spectators(0.)[0]-M.background))<3e-12
    elastic=ElasticChannels(C,grid,1.)
    co=ComptonKernel(M)
    cv=ConversionKernel(M)
    facs=C.factor,elastic.soft.fac,co.fac,cv.fac
    mnorm=M.norm
    y0=np.r_[M.init.ravel(),C.ns,0.,0.,0.]
    vis=np.pi**2*106.75/30
    def totals(y,t):
        M.eta=t if expands else 0.
        _,_,_,e,s=M.weights(y[:C.j0].reshape(-1,2),True)
        _,r,q,*_=C.scalar(y[C.j0:C.j0+C.n])
        return e+r,s+q
    calls=0
    def advance(y,a,b):
        def f(t,zz):
            nonlocal calls
            calls+=1
            M.eta=t if expands else 0.
            M.norm=mnorm*np.exp(-M.eta)
            R=Cold.rhs(C,zz[:-1],True)
            E,F,G,rho,_=M.weights(zz[:C.j0].reshape(-1,2),True)
            rate=np.exp(M.eta)
            for kernel,amp in zip((co,cv),facs[2:],strict=True):
                kernel.fac=amp*rate
                D,S,W,*_=kernel.collision(E,F,G)
                R[:C.j0]+=D.ravel(); R[-2:]+=[S,W]
            if full_H:
                R/=np.sqrt(1+(rho+C.scalar(zz[C.j0:C.j0+C.n])[1])/vis)
            return np.r_[R,M.trace if expands else 0.]
        h=solve_ivp(f,(a,b),y,method='DOP853',rtol=2e-10,atol=3e-13)
        assert h.success
        return h.y[:,-1]
    dd=duration/steps
    y=y0.copy();e0,s0=totals(y,0.)
    an=[]
    history=[float(C.ns.sum())]
    begin=0
    end=steps if stop_after is None else stop_after
    physics=dict(requested,initial=state_digest(y0),edges=[state_digest(e) for e in M.edge_sets])
    if old is not None:
        if old['physics']!=physics:
            raise ValueError("FRW checkpoint physics does not match these controls.")
        begin=old['steps_done']
        y=np.asarray(old['state'])
        history=old['history']
        an=old['implicit']
        calls=old['slow_rhs_calls']
        if (not isinstance(begin,int) or not 0<=begin<=end or y.shape!=y0.shape
            or not np.all(np.isfinite(y)) or old['state_digest']!=state_digest(y)
            or len(history)!=begin+1 or not np.all(np.isfinite(history)) or len(an)!=begin
            or history[0]!=float(C.ns.sum()) or history[-1]!=float(y[C.j0:C.j0+C.n].sum())
            or not isinstance(calls,int) or calls<0):
            raise ValueError("Invalid FRW checkpoint state or history.")
    if progress_every:
        print(f"initialized fixed supports and exterior; steps {begin}/{steps}",file=sys.stderr,flush=True)
    for j in range(begin,end):
        t=dd*(j+.5)
        y=advance(y,dd*j,t)
        e,_=totals(y,t)
        scale=np.exp(M.eta)/(np.sqrt(1+e/vis) if full_H else 1.)
        C.factor,elastic.soft.fac=(k*scale for k in facs[:2])
        Y,info=step(elastic,y[:C.j0].reshape(-1,2),dd)
        y[:C.j0]=Y.ravel()
        y[-3]+=info['entropy_gain']
        an.append(info)
        y=advance(y,t,dd*(j+1))
        history.append(float(y[C.j0:C.j0+C.n].sum()))
        k=j+1
        if checkpoint and (k%checkpoint_every==0 or k==end):
            atomic_write(checkpoint,dict(physics=physics,steps_done=k,state=y.tolist(),state_digest=state_digest(y),
                                        history=history,implicit=an,slow_rhs_calls=calls))
        if progress_every and (k%progress_every==0 or k==end):
            print(f"steps {k}/{steps}; cold gain {history[-1]-history[0]:.9g}",file=sys.stderr,flush=True)
    if end<steps:
        return dict(checkpoint_only=True,steps_done=end,total_steps=steps,scalar_gain_so_far=history[-1]-history[0])
    e,s=totals(y,duration)
    err=e-e0-y[-1]
    assert abs(err)<3e-12 and abs(s-s0-y[-3])<3e-11,(err,s-s0-y[-3])
    return dict(grid=grid,expands=expands,full_H=full_H,steps=steps,duration=duration,projection=projection,
                slow_rhs_calls=calls,state_history=history,initial_dark_visible=e0/vis,
                scalar_gain=float((y-y0)[C.j0:C.j0+C.n].sum()),energy_ledger=err,
                entropy_ledger=s-s0-y[-3],trace_integral=y[-1],implicit=an)

def fixture():
    out=[]
    for g in ('fine','third'):
        M=Screened(g,1,8,2e-6,(1,1,2),.08,1.).M
        M.weights=MethodType(Expanding.weights,M)
        M.spectators=exterior(M,g)
        A=M.init.copy()
        rng=np.random.default_rng(91)
        B=A.copy()
        for span in M.r[1:]:
            ii=np.arange(span.start,span.stop)
            B[ii]*=(1.+rng.uniform(-.06,.06,len(ii)))[:,None]
        def energy(y,t):
            M.eta=t
            R=M.weights(y,True)[3]
            return R,M.trace
        for t in (0.,2e-6,.0002,.02):
            h=2e-5
            x=[energy(A,t+j*h)[0] for j in (-2,-1,1,2)]
            numeric=(x[0]-8*x[1]+8*x[2]-x[3])/(12*h)
            a,analytic=energy(A,t)
            b,_=energy(B,t)
            de=finite_energy(M,A,B)
            assert abs(numeric-analytic)<2e-10
            assert abs(de-(b-a))<2e-15
            out.append(dict(grid=g,eta=t,difference_error=float(de-(b-a)),redshift_trace_error=numeric-analytic))
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=('fine','third'))
    p.add_argument('--duration',type=float,default=2e-6)
    p.add_argument('--steps',type=int,default=12)
    p.add_argument('--expands',action='store_true')
    p.add_argument('--full-H',dest='full_H',action='store_true')
    p.add_argument("--projection",choices=("steps","cell-conservative"),default="steps")
    p.add_argument("--fixture",action="store_true")
    p.add_argument("--checkpoint")
    p.add_argument("--resume",action="store_true")
    p.add_argument("--checkpoint-every",type=int,default=12)
    p.add_argument("--progress-every",type=int,default=0)
    p.add_argument("--stop-after",type=int)
    args=vars(p.parse_args())
    result=fixture() if args.pop("fixture") else run(**args)
    print(json.dumps(result,indent=2,allow_nan=False))
