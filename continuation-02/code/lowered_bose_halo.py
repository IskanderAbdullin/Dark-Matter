"""Positive truncated Bose (King-like) halo; requirements numpy scipy."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.integrate import solve_ivp
from numpy.polynomial.legendre import leggauss

def integral_I(psi,eta,nquad=96):
    p=np.asarray(psi,dtype=float)
    if np.any(p<0) or eta<=np.max(p): raise ValueError("0<=psi<eta required")
    nd,ww=leggauss(nquad)
    u=(nd+1)/2
    ww=ww/2
    a=eta-p[...,None]*(1-u*u)
    delta=p[...,None]*(1-u*u)
    diff=np.exp(a)*np.expm1(delta)/(np.expm1(a)*np.expm1(eta))
    out=p**1.5*np.sum(ww*u*u*diff,axis=-1)
    return float(out) if p.ndim==0 else out

def solve_halo(eta=4.,psi0=3.8,grid_size=600):
    if not 0<psi0<eta: raise ValueError("Require 0<psi0<eta")
    ps=np.linspace(0,psi0,grid_size)
    vals=integral_I(ps,eta)
    assert vals[0]==0 and np.all(np.diff(vals)>0)
    eos=PchipInterpolator(ps,vals/vals[-1],extrapolate=False)
    eps=1.e-5
    def rhs(r,w):
        psi,Q=w
        return [-Q/r**2,r*r*float(eos(max(0,min(psi,psi0))))]
    def edge(r,w):return w[0]
    edge.terminal=True
    edge.direction=-1
    res=solve_ivp(rhs,[eps,200],[psi0-eps**2/6,eps**3/3],
        rtol=2e-9,atol=1e-11,events=edge,max_step=.04,dense_output=True)
    if not res.success or not len(res.t_events[0]): raise RuntimeError("No boundary")
    rt=float(res.t_events[0][0])
    return res,eos,rt,float(res.sol(rt)[1])

def run():
    tests=0
    for eta,psi in ((4.,3.8),(4.,1.),(8.,6.)):
        assert integral_I(psi,eta)>0
        tests+=1
    asym={}
    for eta in (3.,4.,8.):
        ps=np.array([1e-6,1e-5,1e-4])
        vals=integral_I(ps,eta)
        sl=float(np.polyfit(np.log(ps),np.log(vals),1)[0])
        coef=2/15*np.exp(eta)/np.expm1(eta)**2
        assert abs(sl-2.5)<1e-3
        asym[str(eta)]=dict(exponent=sl,ratio=float(vals[0]/(coef*ps[0]**2.5)))
        tests+=1
    sols=[]
    for psi0 in (1.,2.,3.,3.8):
        res,eos,rt,Q=solve_halo(eta=4.,psi0=psi0)
        assert rt>0 and Q>0
        radii=np.array([.1,.25,.5,.75,.9,.99])
        values=eos(res.sol(rt*radii)[0])
        sols.append(dict(psi0=psi0,rt_over_R=rt,total_mass_Q=Q,
         density_by_fraction_rt={str(f):float(v) for f,v in zip(radii,values)}))
        tests+=1
    return dict(tests_passed=tests,eta=4.,solutions=sols,edge_asymptotic=asym,
      disclaimer="No stability, galactic fit or kinetic formation proof")

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--out",default="lowered-bose-results.json")
    args=parser.parse_args()
    data=run()
    Path(args.out).write_text(json.dumps(data,indent=2),encoding="utf-8")
    print(json.dumps(data,indent=2))
