"""Local-column form of the elastic panel flux and analytic Jacobian."""
import numpy as np
from scipy.special import expit
from joint_combined_centered_implicit_reference import Combined


def compile_panels(panels,n):
    out=[]
    for p,c,e,P,raw in panels:
        ix=np.flatnonzero(np.any(P!=0.,axis=0))
        local=np.ascontiguousarray(P[:,ix])
        sf=ix<2*n
        out.append((p,c,e,ix,local,sf,raw))
    return out


def compact_elastic(z,t,eta,H,compiled,num,pairwise=False):
    value=np.zeros(num)
    cache=[]
    for p,c,e,ix,P,sf,raw in compiled:
        local=z[ix]
        z1=P@(local*sf)
        z2=P@(local*(~sf))
        aa=p/t-eta[0]+z1
        bb=c/t-eta[1]+z2
        if np.min(aa)<=0.: raise ValueError("BE exponent")
        f=np.exp(-aa)/(-np.expm1(-aa))
        g=expit(-bb)
        block=(1+f)*(1-g)
        K=raw*(block[:,None]*block)*np.exp(sum(eta)-e/t)/H
        K1=K.sum(axis=1)
        lr=-z1-z2
        ratio=np.exp(lr)
        if pairwise:
            diff=lr[None,:]-lr[:,None]
            delta=np.sign(diff)*np.exp(np.maximum(lr[:,None],lr[None,:]))*(-np.expm1(-abs(diff)))
            net=K*delta
            v=net.sum(axis=1)
            shift=None
        else:
            shift=np.expm1(lr)
            v=K@shift-shift*K1
            net=None
        value[ix]+=P.T@v
        cache.append((ix,P,sf,f,g,K,K1,ratio,shift,net,v))
    jac=None
    def derivative(d):
        nonlocal jac
        if jac is None:
            jac=np.zeros((num,num))
            for ix,P,sf,f,g,K,K1,ratio,u,net,v in cache:
                C=P*np.where(sf[None,:],-f[:,None],g[:,None])
                D=-ratio[:,None]*P
                if pairwise:
                    term=v[:,None]*C+net@C+K@D-K1[:,None]*D
                else:
                    term=v[:,None]*C+K@(u[:,None]*C+D)-u[:,None]*(K@C)-K1[:,None]*D
                jac[np.ix_(ix,ix)]+=P.T@term
        return jac@d
    return value,derivative


class CompactCombined(Combined):
    def __init__(self,k,H,panels,pairwise=True):
        super().__init__(k,H,panels)
        self.compiled=compile_panels(panels,self.n)
        self.pairwise=pairwise

    def elastic(self,z):
        return compact_elastic(z,self.t,self.eta,self.H,self.compiled,self.d,self.pairwise)


def test(eos,row):
    from time import perf_counter
    from joint_below_pole_elastic_linear_modes_reference import setup
    from joint_resonant_kinetic_linear_modes_reference import configuration,metric
    from joint_resonant_centered_stiff_reference import advance
    T,H,k,N,t,Q,sh=configuration(row,eos,8,6.5)
    z0=np.concatenate([x.T.ravel() for x in sh])
    start=z0+.07*np.cos(np.arange(len(z0))*.37+.1)
    old=Combined(k,H,[])
    ref,R,t=old.reference(start,T)
    p=[]
    tm=perf_counter()
    setup(T,k,t,old.eta,4,3,2,10,30.,.82,cache=p)
    dense_build=perf_counter()-tm
    q=[]
    tm=perf_counter()
    setup(T,k,t,old.eta,4,3,2,10,30.,.82,cache=q,cache_only=True)
    cache_build=perf_counter()-tm
    assert len(p)==len(q)>0
    for a,b in zip(p,q):
        for x,y in zip(a,b): assert np.array_equal(x,y)
    old.panels=p
    x=start-ref
    d=np.cos(np.arange(len(x))*.27)
    tm=perf_counter(); f,j=old.elastic(x); j0=j(d); dense=perf_counter()-tm
    checks=[]
    for pair in (False,True):
        new=CompactCombined(k,H,q,pairwise=pair)
        new.reference(start,T)
        tm=perf_counter(); g,kd=new.elastic(x); j1=kd(d); speed=perf_counter()-tm
        err=np.linalg.norm(g-f)/np.linalg.norm(f)
        je=np.linalg.norm(j0-j1)/np.linalg.norm(j0)
        assert max(err,je)<8e-13
        eps=4e-6
        fd=(new.elastic(x+eps*d)[0]-new.elastic(x-eps*d)[0])/(2*eps)
        fe=np.linalg.norm(fd-j1)/np.linalg.norm(j1)
        assert fe<5e-8
        _,_,inv=metric(k,sh)
        inic=new.moment(start)[0]
        fin,stat=advance(new,start,.02,T)
        defect=max(abs((inv.T@new.moment(fin)[0])/(inv.T@inic)-1))
        assert defect<5e-11
        checks.append(dict(pairwise=pair,flux_error=float(err),JVP_error=float(je),
                           finite_difference_error=float(fe),invariant_defect=float(defect),
                           measured_value_and_JVP_seconds=speed,
                           converged_residual=stat["residual"]))
    sizes=np.array([len(entry[3]) for entry in new.compiled])
    return dict(panels=len(p),global_dimension=len(x),mean_active_local_columns=float(sizes.mean()),
                measured_geometry_dense_seconds=dense_build,measured_geometry_cache_only_seconds=cache_build,
                measured_dense_value_and_JVP_seconds=dense,checks=checks)


if __name__=="__main__":
    import argparse,json
    from pathlib import Path
    from sm_dark_entropy_exchange_reference import StandardModel
    p=argparse.ArgumentParser();p.add_argument("eos");p.add_argument("spectra")
    a=p.parse_args()
    r=next(r for r in json.loads(Path(a.spectra).read_text())["rows"] if r["TR"]==.05 and abs(r["TF"]-.045)<1e-8)
    print(json.dumps(test(StandardModel(a.eos),r),indent=2,allow_nan=False))
