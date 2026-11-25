"""Published Ghisoiu--Laine high-temperature sterile-neutrino rate screen.

Data: http://www.laine.itp.unibe.ch/production-highT/mesh_imPiR.dat
Paper: https://arxiv.org/abs/1411.1765, eq. (2.1), physical column 11.
This maps the per-occupation Majorana rate to the existing width convention;
it is NOT a completed Dirac spin/chirality model or a dark-sector thermal
dispersion calculation. We keep existing g=4 for a controlled rate comparison.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.special import expit

URL='http://www.laine.itp.unibe.ch/production-highT/mesh_imPiR.dat'
SHA='722e54d77007f3c591bf167ee8d19efb9fd856d6a35883b467328c4e49eca6b7'

class PortalTable:
    def __init__(self,path):
        p=Path(path)
        self.sha=hashlib.sha256(p.read_bytes()).hexdigest()
        if self.sha != SHA:
            raise ValueError('Unexpected source hash; inspect table provenance')
        a=np.loadtxt(p)
        axes=tuple(np.unique(a[:,j]) for j in (0,2,3))
        assert a.shape==(66185,12) and tuple(map(len,axes))==(35,61,31)
        i,j,k=(np.searchsorted(t,a[:,c]) for t,c in zip(axes,(0,2,3)))
        values=np.full(tuple(map(len,axes)),np.nan)
        values[i,j,k]=a[:,10]
        assert np.isfinite(values).all()
        self.axes=axes
        self.nonpositive=int(np.count_nonzero(values<=0))
        self.interp=RegularGridInterpolator(tuple(np.log(t) for t in axes),values,
                                            bounds_error=True)
    def pi(self,T,z,q):
        q=np.asarray(q,dtype=float)
        if not 160<=T<=1e15 or not z>0 or np.any(q<0):
            raise ValueError('Outside high-T domain')
        # Authors flag a spurious soft-momentum dip at M/T << 1.
        # Use M/T=.2 in the soft corner, taper correction over q=.15-.30.
        # Never extrapolate the temperature dependence.
        qclip=np.clip(q,.1,100.)
        soft=np.clip((.30-qclip)/.15,0.,1.)
        soft=soft*soft*(3-2*soft)
        zm=np.clip(z,.02,20.)
        zz=zm+soft*max(.2-zm,0.)
        p=self.interp(np.stack((np.full_like(q,np.log(T)),np.log(zz),np.log(qclip)),axis=-1))
        if np.any(p<=0):
            raise ValueError('Nonpositive physical tabulated rate after soft repair')
        # For z>20, match the tabulated ratio to the vacuum width at z=40.
        # No claim to know the uncomputed nonrelativistic tail.
        if z>20.:
            u=np.clip(np.log(z/20.)/np.log(2.),0.,1.)
            u=u*u*(3-2*u)
            p=((1-u)*p/(20.**2)+u/(8*np.pi))*z*z
        return p
    def factor(self,T,z,q):
        return 8*np.pi*self.pi(T,z,q)/z**2
    def rate(self,x,theta,Hf,parent_GeV,mratio,kappa,q,dilation):
        z=x*mratio/theta
        T=parent_GeV*theta/x
        return kappa*x*x*Hf*dilation*self.factor(T,z,np.asarray(q)/theta)
    def rate_with_derivative(self,x,theta,Hf,parent_GeV,mratio,kappa,q,dilation):
        r=self.rate(x,theta,Hf,parent_GeV,mratio,kappa,q,dilation)
        h=1e-5
        # Derivative wrt ln(theta), holding Hf and x fixed, as base Jacobian expects.
        a=self.rate(x,theta*np.exp(h),Hf,parent_GeV,mratio,kappa,q,dilation)
        b=self.rate(x,theta*np.exp(-h),Hf,parent_GeV,mratio,kappa,q,dilation)
        return r,(a-b)/(2*h)

def make_models():
    from yukawa_rate_bath import YukawaRateBath
    from resolved_parent_bath import ResolvedParent
    class TableClosure(YukawaRateBath):
        def __init__(self,table,parent_GeV=86000.,**kwargs):
            super().__init__(**kwargs)
            self.table=table
            self.parent_GeV=parent_GeV
        def bath_rates(self,x,th):
            theta=th['theta'];Hf=th['Hfactor'];b=x*self.fermion_ratio
            E=th['E'];p=np.sqrt(np.maximum(E*E-x*x,0.))
            eq=expit(-E/theta);eqT=eq*(1-eq)*E/theta
            ne=self.parent_w@eq;net=self.parent_w@eqT
            rate,rateT=self.table.rate_with_derivative(x,theta,Hf,self.parent_GeV,1.,self.kappa_sm1,p,x/E)
            if ne>1e-200:
                lam1=float(self.parent_w@(rate*eq)/ne)
                lam1T=float((self.parent_w@(rateT*eq+rate*eqT)-lam1*net)/ne)
            else:
                lam1=float(rate[0]);lam1T=float(rateT[0])
            E2=self.energy_chi;q=np.sqrt(np.maximum((E2-b)*(E2+b),0.))
            dil2=self.trace_weights[self.n:]/(self.g_fermion*self.vol*b)
            r2,t2=self.table.rate_with_derivative(x,theta,Hf,self.parent_GeV,self.fermion_ratio,self.kappa2,q,dil2)
            return r2,t2,lam1,lam1T,True
        def run(self,*args,**kwargs):
            out=super().run(*args,**kwargs)
            out['rate_model']='Ghisoiu-Laine 2014 high-T per-occupation rate, interpolated, effective spin convention'
            out['rate_source_sha256']=self.table.sha
            return out
    class TableResolved(ResolvedParent):
        def __init__(self,table,parent_GeV=86000.,**kwargs):
            super().__init__(**kwargs)
            self.table=table;self.parent_GeV=parent_GeV
        def rates(self,x,theta,Hf):
            b=x*self.fermion_ratio
            p=np.sqrt(np.maximum((self.e1-x)*(self.e1+x),0.))
            q=np.sqrt(np.maximum((self.energy_chi-b)*(self.energy_chi+b),0.))
            r1,t1=self.table.rate_with_derivative(x,theta,Hf,self.parent_GeV,1.,self.kappa_sm1,p,self.t1/x)
            r2,t2=self.table.rate_with_derivative(x,theta,Hf,self.parent_GeV,self.fermion_ratio,self.kappa2,q,self.trace_weights[self.n:]/(self.v2*b))
            return np.r_[np.zeros(self.n),r2,r1],np.r_[np.zeros(self.n),t2,t1]
    return TableClosure,TableResolved

def screen(table,parent=86000.):
    rows=[]
    for x in [.1,.2,.3,.5,1.,2.,5.,10.,20.,40.,80.]:
        q=np.array([.1,.2,.5,1.,3.,10.])
        # Compare to the previously used vacuum + tapered gamma=.01 rate.
        E=np.hypot(x,q)
        old=x/E+8*np.pi*.01/(x*(1+x*x))
        new=x/E*table.factor(parent/x,x,q)
        rows.append(dict(x=x,T_GeV=parent/x,q=q.tolist(),published_over_vacuum=table.factor(parent/x,x,q).tolist(),
                         published_over_old=(new/old).tolist()))
    return dict(source=URL,sha256=table.sha,nonpositive_grid_values=table.nonpositive,
                convention='Gamma=|y|^2 ImPi_R/E; Gamma_vac=|y|^2 M/(8pi); g=4 held fixed',rows=rows)

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--table',default='.runtime/mesh_imPiR.dat');p.add_argument('--output',required=True)
    p.add_argument('--mode',choices=['screen','closure','resolved'],default='screen')
    p.add_argument('--n',type=int,default=120);p.add_argument('--nparent',type=int,default=120);p.add_argument('--g-fermion',type=int,choices=[2,4],default=4)
    p.add_argument('--kappa',type=float,default=.005278180162600994)
    p.add_argument('--kappa-sm1',type=float,default=47.);p.add_argument('--kappa2',type=float,default=.03)
    p.add_argument('--fermion-ratio',type=float,default=.9)
    p.add_argument('--x0',type=float,default=.1);p.add_argument('--xend',type=float,default=80.)
    p.add_argument('--wall-limit',type=float,default=0.);p.add_argument('--rtol',type=float,default=2e-6)
    p.add_argument('--fermion-atol',type=float,default=1e-15);p.add_argument('--parent-floor',type=float,default=1e-15)
    p.add_argument('--method',choices=['BDF','Radau'],default='BDF');p.add_argument('--auto-handoff',action='store_true')
    a=p.parse_args();table=PortalTable(a.table)
    if a.mode=='screen':
        result=screen(table)
    else:
        C,R=make_models()
        kwargs=dict(table=table,n=a.n,g_fermion=a.g_fermion,scalar_ratio=1e-12,fermion_ratio=a.fermion_ratio,kappa=a.kappa,kappa2=a.kappa2,
                    kappa_sm1=a.kappa_sm1,qmin=1e-8)
        if a.mode=='closure':
            m=C(**kwargs);r=m.run(a.x0,a.xend,a.rtol,1e-15,True,86000.,'empty')
        else:
            m=R(nparent=a.nparent,parent_floor=a.parent_floor,**kwargs);r=m.evolve(a.x0,a.xend,a.rtol,'empty',True,86000.,a.wall_limit,a.fermion_atol,a.method,a.auto_handoff)
        r['rate_source_sha256']=table.sha
        result=[r]
        print('final',json.dumps(r['final']),flush=True)
    Path(a.output).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')

if __name__=='__main__':main()
