"""Independent parent-angle decay quadrature with conservative energy projection.

Use one representative parent energy, integrate the uniform scalar energy
interval exactly against piecewise-linear nodal basis functions, then project
the remaining fermion energy between adjacent representatives. This guarantees
discrete energy and number balance and eliminates projection of a finite
daughter rectangle onto an anomalously small parent momentum volume.

This is an independent finite-volume/nodal approximation, not an exact
collision operator. Nodal interpolation can spread support slightly beyond
physical angular bounds; contributions outside the represented daughter
energy domain are dropped and their equilibrium-weighted measure is recorded.
The fermion and scalar grids can be distinct.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre,expit
from resolved_parent_bath import ResolvedParent
from two_mass_bath import invariants,rectangle_areas
from laine_portal_rates import PortalTable,make_models


def positive_hat_integrals(nodes,lower,upper):
    """Mass-lumped first-order basis integrals, with truncated support."""
    nodes=np.asarray(nodes);lower=np.asarray(lower);upper=np.asarray(upper)
    l=np.maximum(nodes[:-1,None],lower[None,:])
    h=np.minimum(nodes[1:,None],upper[None,:])
    width=np.maximum(h-l,0.)
    center=(l+h)*.5
    t=np.clip((center-nodes[:-1,None])/np.diff(nodes)[:,None],0.,1.)
    a=width*(1-t);b=width*t
    out=np.zeros((len(nodes),len(lower)))
    out[:-1]+=a;out[1:]+=b
    return out


class ParentAngle(ResolvedParent):
    def __init__(self,fermion_grid='balanced',**kwargs):
        super().__init__(**kwargs)
        self.fermion_grid=fermion_grid
        if fermion_grid=='shared':
            self.fermion_edges=self.edges.copy()
        elif fermion_grid=='balanced':
            nl=self.n//4;nm=5*self.n//8;nh=self.n-nl-nm
            self.fermion_edges=np.r_[np.linspace(0,2,nl+1),
                np.linspace(2,24,nm+1)[1:],np.geomspace(24,self.edges[-1],nh+1)[1:]]
        elif fermion_grid=='uniform':
            nm=7*self.n//8;nh=self.n-nm
            self.fermion_edges=np.r_[np.linspace(0,24,nm+1),
                np.geomspace(24,self.edges[-1],nh+1)[1:]]
        else:raise ValueError(fermion_grid)
        l,h=self.fermion_edges[:-1],self.fermion_edges[1:]
        u,w=roots_legendre(8)
        self.p2=(l[:,None]+h[:,None])/2+(h-l)[:,None]*u/2
        self.w2=(h-l)[:,None]*w[None,:]*self.p2**2/2
        self.w2/=self.w2.sum(axis=1)[:,None]
        self.v2=self.g_fermion*(h-l)*(h*h+h*l+l*l)/(6*np.pi**2)
        self.vols=np.r_[self.vol,self.v2,self.v1]
        self.source2=np.r_[np.zeros(self.n),self.v2,np.zeros(self.n1)]
        self.max_shell_excess=0.;self.max_weighted_shell_loss=0.

    def geometry(self,x):
        if x==self.geom_x:return
        n=self.n;m=x*self.scalar_ratio;b=x*self.fermion_ratio
        K,D=invariants(x,m,b)
        ep=np.hypot(self.ps,m);ef=np.hypot(self.p2,b);e1=np.hypot(self.p1,x)
        self.energy_phi=np.sum(self.ws*ep,axis=1)
        self.energy_chi=np.sum(self.w2*ef,axis=1)
        self.e1=np.sum(self.w1*e1,axis=1)
        self.t1=np.sum(self.w1*x*x/e1,axis=1)
        self.trace_weights=np.r_[self.vol*np.sum(self.ws*m*m/ep,axis=1),
            self.v2*np.sum(self.w2*b*b/ef,axis=1)]
        self.energies=np.r_[self.energy_phi,self.energy_chi,self.e1]
        self.energy_weights=self.vols*self.energies
        self.traces=np.r_[self.trace_weights,self.v1*self.t1]
        p1=np.sqrt(np.maximum((self.e1-x)*(self.e1+x),0.))
        estar=(m*m+K)/x;pstar=D/x
        # Stable backward energy, including near-degenerate nonzero scalar mass.
        low=x*estar/(self.e1+p1)+p1*m*m/(x*(estar+pstar))
        high=(self.e1*estar+p1*pstar)/x
        width=2*p1*pstar/x
        assert np.all(width>0)
        parent_rate=self.v1*self.kappa*x*self.t1
        weights=positive_hat_integrals(self.energy_phi,low,high)/width[None,:]
        ii,jj=np.nonzero(weights>0)
        E2=self.e1[jj]-self.energy_phi[ii]
        R=weights[ii,jj]*parent_rate[jj]
        valid=(E2>=self.energy_chi[0])&(E2<=self.energy_chi[-1])
        ii,jj,E2,R=ii[valid],jj[valid],E2[valid],R[valid]
        hi=np.clip(np.searchsorted(self.energy_chi,E2,side='right'),1,n-1)
        lo=hi-1
        wh=np.clip((E2-self.energy_chi[lo])/(self.energy_chi[hi]-self.energy_chi[lo]),0,1)
        wl=1-wh
        load=np.bincount(jj,weights=R,minlength=self.n1)
        ratio=load/parent_rate
        self.max_shell_excess=max(self.max_shell_excess,float(np.max(ratio-1,initial=0)))
        eq=expit(-self.e1)
        loss=float(np.dot(parent_rate-load,eq)/max(np.dot(parent_rate,eq),1e-300))
        self.max_weighted_shell_loss=max(self.max_weighted_shell_loss,loss)
        self.max_dropped_equilibrium_fraction=max(self.max_dropped_equilibrium_fraction,loss)
        err=np.abs(self.energy_phi[ii]+wl*self.energy_chi[lo]+wh*self.energy_chi[hi]-self.e1[jj])/self.e1[jj]
        self.max_projection_energy_error=max(self.max_projection_energy_error,float(np.max(err,initial=0)))
        self.ii,self.jp,self.lo,self.hi=ii,jj,lo,hi
        self.wl,self.wh,self.event_rate0=wl,wh,R
        self.shell_rate_ratio=ratio
        self.event_idx=np.vstack((ii,n+lo,n+hi,2*n+jj))
        self.stoich=np.vstack((1/self.vol[ii],wl/self.v2[lo],wh/self.v2[hi],-1/self.v1[jj]))
        self.jidx=(self.event_idx[:,None,:]*self.dim+self.event_idx[None,:,:]).ravel()
        self.geom_x=x

    def dark(self,yd,Hf=1.,jac=False):
        n=self.n;phi=yd[self.ii];f=yd[n+self.lo];h=yd[n+self.hi];g=yd[2*n+self.jp]
        if jac:
            F,Fl,Fh=self.geometric(f,h,self.wl,self.wh,True)
            B,Bl,Bh=self.geometric(1-f,1-h,self.wl,self.wh,True)
        else:
            F=self.geometric(f,h,self.wl,self.wh)
            B=self.geometric(1-f,1-h,self.wl,self.wh)
        R=self.event_rate0*Hf
        net=R*(g*B+phi*(g*B-(1-g)*F))
        dy=np.bincount(self.event_idx.ravel(),weights=(self.stoich*net).ravel(),minlength=self.dim)
        if not jac:return dy
        deriv=np.vstack((R*(g*B-(1-g)*F),
            -R*(g*(1+phi)*Bl+(1-g)*phi*Fl),
            -R*(g*(1+phi)*Bh+(1-g)*phi*Fh),
            R*((1+phi)*B+phi*F)))
        J=np.bincount(self.jidx,weights=(self.stoich[:,None,:]*deriv[None,:,:]).ravel(),
                      minlength=self.dim*self.dim).reshape(self.dim,self.dim)
        return dy,J

    def validate(self):
        rng=np.random.default_rng(192)
        rows=[]
        for x in [.1,.2,1.,5.,20.]:
            self.geometry(x)
            eq=expit(-self.energies);eq[:self.n]=1/np.expm1(self.energy_phi)
            d=self.dark(eq)
            bal=float(np.abs(d)@self.vols/max(self.event_rate0@expit(-self.e1[self.jp]),1e-300))
            assert bal<1e-8,bal
            g=expit(-self.e1);p=np.r_[np.zeros(2*self.n),g]
            c=self.dark(p)
            expect=self.kappa*x*(self.v1*self.t1)@g
            ratio=float(self.vol@c[:self.n]/expect)
            assert .85<ratio<=1+1e-10,ratio
            assert np.max(self.shell_rate_ratio)<=1+1e-10
            rows.append(dict(x=x,balance=bal,spontaneous_ratio=ratio,minimum_shell_ratio=float(self.shell_rate_ratio.min()),
                maximum_shell_ratio=float(self.shell_rate_ratio.max())))
        self.geometry(.2)
        y=np.r_[rng.uniform(.05,.3,self.n),rng.uniform(.03,.2,self.n+self.n1),.01,0.,0.,0.]
        c=self.dark(y[:self.dim]);denom=np.abs(c)@self.vols
        Np=self.vol@c[:self.n];N2=self.v2@c[self.n:2*self.n];N1=self.v1@c[2*self.n:]
        en=float(abs(c@self.energy_weights)/max(np.abs(c)@self.energy_weights,1e-300))
        num=float(max(abs(Np-N2),abs(N1+Np))/max(denom,1e-300))
        dr=rng.normal(size=len(y));dr/=np.linalg.norm(dr);eps=1e-6
        fd=(self.rhs(np.log(.2),y+eps*dr)-self.rhs(np.log(.2),y-eps*dr))/(2*eps)
        jd=self.evaluate(np.log(.2),y,True)@dr
        jerr=float(np.linalg.norm(fd-jd)/max(np.linalg.norm(jd),1e-300))
        assert en<1e-12 and num<1e-12 and jerr<1e-5,(en,num,jerr)
        for bound in [0.,1.]:
            q=y[:self.dim].copy();q[self.n::3]=bound
            mask=np.r_[np.zeros(self.n,dtype=bool),q[self.n:]==bound]
            dc=self.dark(q)
            assert np.min(dc[mask]*(1-2*bound),initial=0.)>=-1e-10
        return dict(balance=max(r['balance'] for r in rows),energy=en,number=num,jacobian=jerr,
            spontaneous_rate=rows,max_shell_excess=self.max_shell_excess,
            max_dropped_equilibrium_fraction=self.max_dropped_equilibrium_fraction,
            max_projection_energy_error=self.max_projection_energy_error)

    def evolve(self,*args,**kwargs):
        r=super().evolve(*args,**kwargs)
        r.update(kernel='parent-angle',fermion_grid=self.fermion_grid,
            max_shell_excess=self.max_shell_excess,max_weighted_shell_loss=self.max_weighted_shell_loss)
        if 'frozen_cell_spectrum' in r:r['frozen_cell_spectrum']['chi2_edges']=self.fermion_edges.tolist()
        return r


class SeparatedRectangle(ParentAngle):
    """Original daughter-rectangle quadrature, independent heavy grids.

    Unlike ParentAngle, the scalar cell has a piecewise-constant finite-volume
    source and parent is fractionally projected. Comparison on matching grids
    separates angular interpolation effects from fermion mesh effects.
    """
    def geometry(self,x):
        if x==self.geom_x:return
        n=self.n;m=x*self.scalar_ratio;b=x*self.fermion_ratio
        ep=np.hypot(self.ps,m);ef=np.hypot(self.p2,b);e1=np.hypot(self.p1,x)
        self.energy_phi=np.sum(self.ws*ep,axis=1)
        self.energy_chi=np.sum(self.w2*ef,axis=1)
        self.e1=np.sum(self.w1*e1,axis=1)
        self.t1=np.sum(self.w1*x*x/e1,axis=1)
        self.trace_weights=np.r_[self.vol*np.sum(self.ws*m*m/ep,axis=1),
            self.v2*np.sum(self.w2*b*b/ef,axis=1)]
        self.energies=np.r_[self.energy_phi,self.energy_chi,self.e1]
        self.energy_weights=self.vols*self.energies
        self.traces=np.r_[self.trace_weights,self.v1*self.t1]
        self.pair_E=self.energy_phi[:,None]+self.energy_chi[None,:]
        _,D=invariants(x,m,b)
        eta=2*D/(x*x)
        area=rectangle_areas(self.edges,self.fermion_edges,x,m,b,self.order)
        self.R0=self.g_fermion*x/(2*np.pi**2)*area/eta*(self.kappa*x*x)
        i,j=np.nonzero(self.R0>0)
        E=self.pair_E[i,j];R=self.R0[i,j]
        valid=(E>=self.e1[0])&(E<=self.e1[-1])
        eqweighted=R*expit(-E)
        loss=float(eqweighted[~valid].sum()/max(eqweighted.sum(),1e-300))
        self.max_dropped_equilibrium_fraction=max(self.max_dropped_equilibrium_fraction,loss)
        i,j,E,R=i[valid],j[valid],E[valid],R[valid]
        hi=np.clip(np.searchsorted(self.e1,E,side='right'),1,self.n1-1);lo=hi-1
        wh=np.clip((E-self.e1[lo])/(self.e1[hi]-self.e1[lo]),0.,1.);wl=1-wh
        self.max_projection_energy_error=max(self.max_projection_energy_error,
            float(np.max(np.abs(wl*self.e1[lo]+wh*self.e1[hi]-E)/E,initial=0.)))
        self.ii,self.jj,self.lo,self.hi=i,j,lo,hi
        self.wl,self.wh,self.event_rate0=wl,wh,R
        self.event_idx=np.vstack((i,n+j,2*n+lo,2*n+hi))
        self.stoich=np.vstack((1/self.vol[i],1/self.v2[j],-wl/self.v1[lo],-wh/self.v1[hi]))
        self.jidx=(self.event_idx[:,None,:]*self.dim+self.event_idx[None,:,:]).ravel()
        self.geom_x=x

    dark=ResolvedParent.dark
    validate=ResolvedParent.validate

    def evolve(self,*args,**kwargs):
        r=ResolvedParent.evolve(self,*args,**kwargs)
        r.update(kernel='separated-rectangle',fermion_grid=self.fermion_grid)
        if 'frozen_cell_spectrum' in r:r['frozen_cell_spectrum']['chi2_edges']=self.fermion_edges.tolist()
        return r

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--table',default='.runtime/mesh_imPiR.dat')
    p.add_argument('--kernel',choices=['angle','rectangle'],default='angle')
    p.add_argument('--n',type=int,default=120);p.add_argument('--nparent',type=int,default=120)
    p.add_argument('--kappa',type=float,default=.006);p.add_argument('--kappa-sm1',type=float,default=47.)
    p.add_argument('--kappa2',type=float,default=.03);p.add_argument('--fermion-ratio',type=float,default=.9)
    p.add_argument('--g-fermion',type=int,choices=[2,4],default=4)
    p.add_argument('--parent-grid',choices=['uniform','balanced','legacy'],default='balanced')
    p.add_argument('--fermion-grid',choices=['uniform','balanced','shared'],default='balanced')
    p.add_argument('--parent-floor',type=float,default=1e-18)
    p.add_argument('--fermion-atol',type=float,default=1e-11);p.add_argument('--rtol',type=float,default=1e-5)
    p.add_argument('--x0',type=float,default=.1);p.add_argument('--xend',type=float,default=80.)
    p.add_argument('--method',choices=['BDF','Radau'],default='BDF');p.add_argument('--local-sparse-jac',action='store_true');p.add_argument('--history-psd',action='store_true')
    p.add_argument('--wall-limit',type=float,default=480.);p.add_argument('--history-dense',action='store_true')
    p.add_argument('--auto-handoff',action='store_true');p.add_argument('--test-only',action='store_true');p.add_argument('--old-rates',action='store_true')
    p.add_argument('--output');a=p.parse_args()
    kw=dict(n=a.n,nparent=a.nparent,scalar_ratio=1e-12,fermion_ratio=a.fermion_ratio,
        kappa=a.kappa,kappa_sm1=a.kappa_sm1,kappa2=a.kappa2,g_fermion=a.g_fermion,
        parent_grid=a.parent_grid,fermion_grid=a.fermion_grid,parent_floor=a.parent_floor,local_sparse_jac=a.local_sparse_jac,record_scalar_history=a.history_psd,qmin=1e-8)
    K=ParentAngle if a.kernel=='angle' else SeparatedRectangle
    if a.old_rates:m=K(**kw);table=None
    else:
        _,R=make_models()
        class TableParent(K,R):pass
        table=PortalTable(a.table);m=TableParent(table=table,**kw)
    if a.test_only:
        print(json.dumps(m.validate(),indent=2));return
    if not a.output:raise ValueError('output required')
    r=m.evolve(a.x0,a.xend,a.rtol,'empty',True,86000.,a.wall_limit,a.fermion_atol,a.method,a.auto_handoff,True,history_sample_x=(np.r_[np.linspace(.13,.22,37),np.linspace(1.,2.5,31)] if a.history_dense else ()))
    if table:r['rate_source_sha256']=table.sha
    Path(a.output).write_text(json.dumps([r],indent=2,allow_nan=False)+'\n')
    print('final',json.dumps({k:r[k] for k in ['success','message','n','nparent','wall_seconds','final']}),flush=True)
    if not r['success'] and not r['handoff_reached']:raise RuntimeError(r['message'])

if __name__=='__main__':main()
