"""Fixed injection spectrum test, no scalar stimulation or bath dynamics."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
from scipy.linalg import solve
from scipy.optimize import brentq
from soft_scalar_elastic_graph_reference import graph,independent_walk
from decay_scalar_ir_relaxation_reference import limits
from three_species_decay_onsager_reference import M,m,mp,P,width
from sm_dark_entropy_exchange_reference import StandardModel,MPL

class Bath:
    def __init__(self,row,eos,smooth=True):
        T=row["TF"]
        bins=np.array(row["scalar_differential_bins"])*T
        counts=np.array(row["stable_normalized_counts"])
        n=sum(row["stable_yields"])*eos.state(T)[2]
        f=n*counts/(4*np.diff(bins**3)/(6*np.pi**2))
        # Interpolate finite positive supported means without extrapolating
        # the sparse high-momentum tail; retain a distinct cell-means control.
        z=np.sqrt(bins[:-1]*bins[1:])
        ok=(f>1e-300)&(z>0.)
        p=z[ok]
        self.pp=p
        self.f0=f
        self.bins=bins
        self.smooth=smooth
        self.normalization=1.
        if smooth:
            self.curve=PchipInterpolator(np.log(p),np.log(f[ok]),extrapolate=False)
        a=np.geomspace(max(p[0]*.1,1e-10),bins[-1],2200)
        q=np.sqrt(a*a+m*m)
        v=self(q)
        nn=np.trapezoid(2*a*a*v/np.pi**2,a)
        assert abs(nn/n-1)<.18,(n,nn)
        self.relative_renormalization=n/nn
        self.normalization=n/nn
        self.interp_integral=-cumulative_trapezoid((v*self.normalization)[::-1],q[::-1],initial=0.)[::-1]
        self.q=q

    def __call__(self,E):
        p=np.sqrt(np.maximum(0.,(E-m)*(E+m)))
        if self.smooth:
            x=self.curve(np.log(np.maximum(p,1e-300)))
            v=np.exp(np.where(np.isfinite(x),x,-800.))
        else:
            i=np.searchsorted(self.bins,p,side="right")-1
            v=np.where(i<len(self.f0),self.f0[np.clip(i,0,len(self.f0)-1)],0.)
        ans=self.normalization*v
        assert np.min(ans)>=0. and np.max(ans)<.08
        return ans

    def rate(self,p):
        # Negligible chi1 by prompt-decay control. Other inverse damping
        # from the injected chi2 spectrum, with no tail extrapolation.
        e,L,d=limits(p)
        if p==0.:return 0.
        a,b=np.interp([L-e,L-e+d],self.q,self.interp_integral)
        return 4*width*M*M/(2*P*p*e)*max(0.,a-b)

def evaluate(row,eos,smooth=True,n=(6,5,4,14),pole=.82):
    T=row["TF"]
    H=np.sqrt(eos.state(T)[0]/3)/MPL
    b=Bath(row,eos,smooth)
    p_s=np.geomspace(1e-7*T,3*T,180)
    fast=np.array([b.rate(p)/H for p in p_s])
    cr=[(p,q) for p,q,x,y in zip(p_s[:-1],p_s[1:],fast[:-1],fast[1:]) if x<1<=y]
    if not cr:return dict(TR=row["TR"],TF=T,has_fast_decay=False,max_decay_over_H=fast.max())
    q=brentq(lambda p:b.rate(p)/H-1,*cr[0],xtol=2e-15)
    support=np.r_[0.,np.geomspace(.5*mp,q,15)]
    A=graph(T,support,n=n,pole=pole,tail=32.,tagged_bath=b)
    mesh=np.diff(support**3)/(6*np.pi**2)
    off=A[:-1,:-1].copy();np.fill_diagonal(off,0.)
    outgoing=off.sum(axis=1)+A[:-1,-1]
    assert np.all(outgoing>0.)
    L=np.diag(outgoing)-off
    tau=solve(L,mesh)*H
    assert np.all(tau>0.)
    c=solve(L,A[:-1,-1])
    assert np.max(abs(c-1))<2e-8
    tm,se,j=independent_walk(off/outgoing[:,None],outgoing,mesh,H,n=18000)
    assert abs(tm-tau[0])<5*se
    return dict(TR=row["TR"],TF=T,has_fast_decay=True,fast_threshold_p_over_T=q/T,
                max_decay_over_H=fast.max(),occupation_interpolation_scale=b.relative_renormalization,
                first_passage_H_time=tau[0],first_cell_direct_H_rate=A[0,-1]/mesh[0]/H,
                mc_mean=tm,mc_se=se,mean_jump_count=j)
def run(eosfile,datafile):
    d=json.loads(Path(datafile).read_text())
    eos=StandardModel(eosfile)
    rows=[r for r in d["rows"] if r["TR"] in (.02,.05,.1) and abs(r["TF"]/r["TR"]-.9)<1e-6]
    out=[]
    for r in rows:
        a=evaluate(r,eos)
        b=evaluate(r,eos,smooth=True,n=(7,5,5,16))
        c=evaluate(r,eos,smooth=False,n=(7,5,5,16))
        if all(z["has_fast_decay"] for z in (a,b,c)):
            d=max(abs(a["first_passage_H_time"]/b["first_passage_H_time"]-1),
                  abs(c["first_passage_H_time"]/b["first_passage_H_time"]-1))
            assert d<.18,d
        out.append(dict(smooth=a,refined=b,piecewise_constant=c))
    return dict(rows=out,scope="Tag-particle diagnostic in the fixed collisionless injection spectrum; no scalar Bose effects, no evolving bath, no thermalization or actual late abundance claim.")
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument('eos');p.add_argument('source')
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.source),indent=2,allow_nan=False))
