"""Two finite daughter masses with closed background energy (constant bath g).
An isotropic vacuum two-body decay and phenomenological parent chemistry.
All masses/momenta in Tref units; actual channel width defines kappa.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre, expit
from energy_closed_bath import EnergyClosedBath
from two_mass_geometry import cell_area


def kinetic(p, m):
    return p*p/(np.hypot(p,m)+m)


def kinetic_bounds(p,m,b,K,D):
    """Bounds on OTHER daughter's kinetic energy for fixed |p|."""
    if m==0:
        return (b*p-K)**2/(2*K*p), np.full_like(p,np.inf)
    E=np.hypot(p,m)
    den=K*p+D*E
    # Rationalized signed lower momentum: avoids subtraction near degeneracy.
    pminus=(b*b*p*p-D*D)/den
    pplus=den/(m*m)
    return kinetic(pminus,b),kinetic(pplus,b)


def invariants(x,m,b):
    if not x>m+b or min(m,b)<0:
        raise ValueError("channel closed or negative mass")
    # Avoid subtracting near threshold.
    c=.5*(x-m-b)*(x+m+b)
    K=c+m*b
    D=np.sqrt(c*(c+2*m*b))
    return K,D


def rectangle_areas(edges,other,x,m,b,order=10):
    """Fast positive FV measure, integrating only boundary-crossed cells.

    Energy widths use exact rationalized differences. Integration in p rather
    than E resolves near-rest cells; at most four kinematic crossings split
    every partial cell into smooth pieces. Full/empty cells are analytic.
    """
    K,D=invariants(x,m,b)
    pl,ph=edges[:-1],edges[1:]
    sl,sh=other[:-1],other[1:]
    vl,vh=kinetic(sl,b),kinetic(sh,b)
    de=(ph-pl)*(ph+pl)/(np.hypot(ph,m)+np.hypot(pl,m))
    dr=(sh-sl)*(sh+sl)/(np.hypot(sh,b)+np.hypot(sl,b))
    a,c=kinetic_bounds(pl,m,b,K,D)
    d,e=kinetic_bounds(ph,m,b,K,D)
    mn=np.minimum(a,d)
    if b>0:
        pzero=D/b
        mn=np.where((pl<=pzero)&(pzero<=ph),0.,mn)
    mx=np.maximum(a,d)
    full=(mx[:,None]<=vl[None,:])&(c[:,None]>=vh[None,:])
    empty=(mn[:,None]>=vh[None,:])|(e[:,None]<=vl[None,:])
    area=de[:,None]*dr[None,:]*full
    i,j=np.nonzero(~(full|empty))
    if len(i)==0:return area
    left,right=pl[i],ph[i]
    low,high=vl[j],vh[j]
    knots=[left,right]
    for ss in (sl[j],sh[j]):
        u,v=kinetic_bounds(ss,b,m,K,D)
        for vv in (u,v):
            # Clip before multiplying to avoid inf / overflow for tiny masses.
            vv=np.minimum(vv,kinetic(right,m))
            pp=np.sqrt(vv*(vv+2*m))
            knots.append(np.clip(pp,left,right))
    knots=np.sort(np.stack(knots),axis=0)
    result=np.zeros(len(i))
    node,weight=roots_legendre(order)
    for L,H in zip(knots[:-1],knots[1:]):
        width=H-L
        active=width>0
        if not np.any(active):continue
        lo,hi=low[active],high[active]
        p=L[active,None]+width[active,None]*(node+1)/2
        vmin,vmax=kinetic_bounds(p,m,b,K,D)
        overlap=np.maximum(0.,np.minimum(hi[:,None],vmax)-np.maximum(lo[:,None],vmin))
        result[active]+=width[active]/2*((p/np.hypot(p,m)*overlap)@weight)
    area[i,j]=np.clip(result,0.,de[i]*dr[j])
    return area


class TwoMassBath(EnergyClosedBath):
    def __init__(self,n=200,scalar_ratio=1e-9,fermion_ratio=.5,order=10,**kwargs):
        super().__init__(n=n,**kwargs)
        invariants(1.,scalar_ratio,fermion_ratio)
        self.scalar_ratio=scalar_ratio;self.fermion_ratio=fermion_ratio;self.order=order
        l,h=self.edges[:-1],self.edges[1:];u,w=roots_legendre(8)
        self.ps=(h[:,None]+l[:,None])/2+(h-l)[:,None]/2*u
        self.ws=(h-l)[:,None]/2*w*self.ps**2
        self.ws/=self.ws.sum(axis=1)[:,None]
        self.energy_chi=self.mean

    def geometry(self,x):
        if x==self.geom_x:return
        m=x*self.scalar_ratio;b=x*self.fermion_ratio
        if m==b==0:
            super().geometry(x);return
        ep=np.hypot(self.ps,m);ef=np.hypot(self.ps,b)
        self.energy_phi=np.sum(self.ws*ep,axis=1)
        self.energy_chi=np.sum(self.ws*ef,axis=1)
        trp=np.sum(self.ws*m*m/ep,axis=1)
        trf=np.sum(self.ws*b*b/ef,axis=1)
        self.ew=np.r_[self.vol*self.energy_phi,self.g_fermion*self.vol*self.energy_chi]
        self.trace_weights=np.r_[self.vol*trp,self.g_fermion*self.vol*trf]
        self.pair_E=self.energy_phi[:,None]+self.energy_chi[None,:]
        area=rectangle_areas(self.edges,self.edges,x,m,b,self.order)
        K,D=invariants(x,m,b)
        eta=2*D/(x*x)
        self.R0=self.g_fermion*x/(2*np.pi**2)*area/eta*(self.kappa*x*x)
        self.geom_x=x

    def test_balance_jac(self):
        x=.2;z=np.log(x);self.geometry(x)
        N=self.thermo(x,0.,0.)['Neq']
        y=np.r_[1/np.expm1(self.energy_phi),expit(-self.energy_chi),N,0.,0.]
        self.set_state(x,y)
        dy=self.collision_rhs(y[:2*self.n])
        bal=(np.abs(dy[:self.n])@self.vol+np.abs(dy[self.n:])@(self.g_fermion*self.vol))/max(self.A.sum(),1e-30)
        assert bal<1e-9,bal
        y=np.r_[np.full(self.n,.1),np.full(self.n,.05),N*.8,.01,0.]
        rng=np.random.default_rng(17);d=rng.normal(size=len(y));d/=np.linalg.norm(d)
        eps=1e-6;j=self.jac_log(z,y)@d
        fd=(self.rhs_log(z,y+eps*d)-self.rhs_log(z,y-eps*d))/(2*eps)
        err=np.linalg.norm(j-fd)/max(np.linalg.norm(j),1e-30)
        assert err<1e-5,err
        return dict(balance=float(bal),jacobian=float(err))

    def run_massive(self,x0,xend,rtol,keep_psd,parent_gev):
        check=self.test_balance_jac()
        out=self.evolve_closed(x0,xend,rtol,keep_psd)
        out.update(scalar_ratio=self.scalar_ratio,fermion_ratio=self.fermion_ratio,parent_GeV=parent_gev,checks=check)
        for s in out['history']:
            # The inherited delta Neff is the MASSLESS moment only. A massive
            # stable fermion is a relic, not a constant radiation contribution.
            if self.fermion_ratio>0:
                s['delta_Neff_chi2_massless_momentum_proxy']=s.pop('delta_Neff_chi2_if_stable')
                s['omega_chi2_h2_if_stable_NR']=s['omega_phi_h2_m86']*self.fermion_ratio*parent_gev*1e9/86.
        return out


def test_geometry():
    from massive_daughter_bath import allowed_area
    rng=np.random.default_rng(7321)
    max_ref=max_swap=max_old=0.;count=0;nonzero=0
    for x,m,b in [(1.,0.,0.),(1.,.03,0.),(1.,0.,.9),(1.,.03,.2),(1.,.49,.49),(1.,1e-9,.999),(1.,.1,.89999)]:
        p=np.geomspace(1e-5,12.,100);r=np.geomspace(2e-5,14.,110)
        a=rectangle_areas(p,r,x,m,b)
        swapped=rectangle_areas(r,p,x,b,m).T
        ep=np.hypot(p,m);er=np.hypot(r,b)
        de=(p[1:]-p[:-1])*(p[1:]+p[:-1])/(ep[1:]+ep[:-1])
        dr=(r[1:]-r[:-1])*(r[1:]+r[:-1])/(er[1:]+er[:-1])
        scale=de[:,None]*dr[None,:]
        # Compare all cells in both orientations, including very thin cells.
        max_swap=max(max_swap,float(np.max(np.abs(a-swapped)/scale)))
        assert np.min(a)>=0 and np.max(a-scale)<1e-12
        if b==0:
            old=allowed_area(ep[:-1,None],ep[1:,None],er[None,:-1],er[None,1:],x,m)
            max_old=max(max_old,float(np.max(np.abs(a-old)/scale)))
        partial=np.argwhere((a>0)&(a<.999999*scale))
        idx=rng.integers([0,0],[len(p)-1,len(r)-1],size=(100,2))
        if len(partial):
            idx=np.vstack((idx,partial[rng.integers(len(partial),size=160)]))
        for i,j in idx:
            ref=cell_area(ep[i],ep[i+1],er[j],er[j+1],x,m,b)
            rel=abs(a[i,j]-ref)/scale[i,j]
            max_ref=max(max_ref,float(rel));count+=1;nonzero+=bool(a[i,j]>0)
    result=dict(sampled_cells=count,nonzero=nonzero,max_reference_scaled_error=max_ref,max_swap_scaled_error=max_swap,max_previous_scaled_error=max_old)
    print(json.dumps(result),flush=True)
    assert max_ref<5e-5 and max_swap<5e-5 and max_old<5e-5,result
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--test-only',action='store_true');p.add_argument('--n',type=int,default=160)
    p.add_argument('--scalar-ratio',type=float,default=1e-9)
    p.add_argument('--fermion-ratios',nargs='+',type=float,default=[0.,.5,.9,.99])
    p.add_argument('--parent-gev',type=float,default=86.)
    p.add_argument('--x0',type=float,default=.1);p.add_argument('--xend',type=float,default=40.)
    p.add_argument('--qmin',type=float,default=1e-7);p.add_argument('--kappa',type=float,default=.1383)
    p.add_argument('--relaxation',type=float,default=100.);p.add_argument('--rtol',type=float,default=2e-6)
    p.add_argument('--order',type=int,default=10);p.add_argument('--keep-psd',action='store_true')
    p.add_argument('--output');a=p.parse_args()
    if a.test_only:
        test_geometry();return
    assert a.output
    out=[]
    for b in a.fermion_ratios:
        m=TwoMassBath(a.n,scalar_ratio=a.scalar_ratio,fermion_ratio=b,order=a.order,qmin=a.qmin,kappa=a.kappa,relaxation=a.relaxation)
        row=m.run_massive(a.x0,a.xend,a.rtol,a.keep_psd,a.parent_gev)
        out.append(row)
        print(json.dumps(dict(n=a.n,x0=a.x0,kappa=a.kappa,fermion_ratio=b,checks=row['checks'],**row['final'])),flush=True)
        # Save after EACH completed case (including its history and optional PSD).
        Path(a.output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
