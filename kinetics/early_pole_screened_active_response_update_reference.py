"""Active-cell retarded polarization increments with a frozen archived exterior."""
import argparse
import json

import numpy as np
from scipy.special import expit, roots_legendre

from early_pole_archived_cold_screened_elastic_reference import Screened
from early_pole_frozen_dynamic_vlasov_reference import C, real_kernels
from early_pole_frozen_elastic_direct_q_reference import traces
from early_pole_screened_elastic_discrete_gradient_reference import step


class ActiveResponse:
    def __init__(self,M,smear=.12,mesh=1,radial_order=None):
        if not 0.<=smear<=.35 or mesh<1:
            raise ValueError("smear must stay within adjacent cell widths")
        if radial_order is not None and radial_order<8:
            raise ValueError("response radial quadrature must have at least eight nodes per cell")
        if radial_order is not None:
            rz,rw=roots_legendre(radial_order)
        self.model=M
        self.M=M.M
        self.base=M.polarization
        self.g2=.005**2
        self.sample=np.unique(np.r_[np.linspace(0,.99,501*mesh),
                                    1-np.geomspace(1e-8,.01,500*mesh)])
        zz,ww=roots_legendre(6)
        v=[]
        self.cells=[]
        for k in self.M.r[1:]:
            cells=self.M.cells[k]
            widths=np.array([c.hi-c.lo for c in cells])
            h=np.r_[widths[0],np.minimum(widths[:-1],widths[1:]),widths[-1]]
            for j,c in enumerate(cells,k.start):
                left,right=np.array([c.lo,c.hi])
                if radial_order is None:
                    energy,measure,u=c.energy,c.w,c.u
                    q2=energy*energy-c.m*c.m
                else:
                    half=(right-left)/2
                    p=(left+right)/2+half*rz
                    q2=p*p
                    energy=np.hypot(p,c.m)
                    measure=C*half*rw*q2
                    u=c.basis(p)
                v.extend(np.sqrt(q2)/energy)
                for b,radius in zip((left,right),smear*h[j-k.start:j-k.start+2],strict=True):
                    p=b+radius*zz
                    v.extend(p/np.hypot(p,c.m))
                left_factor=-C*left*np.hypot(left,c.m)*ww/2
                right_factor=C*right*np.hypot(right,c.m)*ww/2
                interp_gain=measure*(energy+q2/energy)/q2
                self.cells.append((j,u,measure,c.de,left_factor,right_factor,interp_gain))
        phase=self.sample[:,None]
        velocity=np.asarray(v)
        a,b,inside,z=real_kernels(phase,velocity)
        v2=velocity**2
        self.kernel=np.stack([a,b*v2,np.pi/2*np.where(inside,z,0.),
                       np.pi/4*np.where(inside,z*(1-z*z),0.)*v2],axis=0)
        self.F0=self.M.weights(self.M.init,True)[1]
        self.coeff0,self.direct0=self.occupation_weights(self.F0)
        self.initial_graph=M.graph
        self.soft=M.graph[-2]<M.cut
        ix,u,b,q,W=M.graph
        p,p1,_=M.original_geometry
        ix=self.soft
        self.a,self.b,self.q=p[ix],b[ix],q[ix]
        self.x=(p1[ix]-p[ix])/q[ix]
        self.old_trace=traces(self.a,self.b,self.q,self.x,1.,0.,self.base)[1]
        assert min(self.old_trace)>0.

    def occupation_weights(self,F):
        weights=[]
        stat=0.
        for j,u,measure,scale,left,right,static in self.cells:
            a,b=F[j]
            f=expit(a+b*u)
            # -C*p*E*f'(p)*dp = -c.w*b*f(1-f)/Delta_E.
            weights.extend(-measure*b/scale*f*(1-f))
            weights.extend(left*expit(a-b))
            weights.extend(right*expit(a+b))
            stat+=static@f
        return np.asarray(weights),stat

    def delta(self,Y):
        F=self.M.weights(Y,True)[1]
        weight,stat=self.occupation_weights(F)
        update=(self.kernel@(weight-self.coeff0)).T*self.g2
        exact=(stat-self.direct0)*self.g2
        assert abs(update[0,0]-exact)<1e-15
        return update,exact

    def refresh(self,Y,mode):
        update,ds=self.delta(Y)
        if mode not in ('frozen','debye','full'):
            raise ValueError(mode)
        base=self.base
        mD2=self.model.mD2
        if mode=='full':
            def pol(x):
                return base(x)+np.stack([np.interp(x,self.sample,w) for w in update.T],axis=-1)
        elif mode=='debye':
            def pol(x):
                return base(x)*(1.+ds/mD2)
        else:
            def pol(x):
                return base(x)
        new=traces(self.a,self.b,self.q,self.x,1.,0.,pol)[1]
        ratio=new/self.old_trace
        assert np.all(np.isfinite(ratio)) and min(ratio)>0.
        W=self.initial_graph[-1].copy()
        W[self.soft]*=ratio
        self.model.graph=(*self.initial_graph[:-1],W)
        rel=ds/mD2
        return dict(mode=mode,fractional_mD2_change=rel,
                    max_response_change_over_initial_mD2=float(np.max(abs(update))/mD2),
                    response_phases=[0.,.2,.5,.8,.95,.99,.999,.9999],
                    response_delta_over_initial_mD2=[
                        (np.array([np.interp(x,self.sample,w) for w in update.T])/mD2).tolist()
                        for x in [0.,.2,.5,.8,.95,.99,.999,.9999]],
                    soft_packet_mean_ratio=float(W[self.soft].sum()/self.initial_graph[-1][self.soft].sum()),
                    soft_packet_ratio_quantiles=np.quantile(ratio,[0.,.01,.5,.99,1.]).tolist())


def run(grid,duration,smear=.12,resolution=1):
    M=Screened(grid,1,8,duration,(1,1,2),.08,1.)
    S=ActiveResponse(M,smear,resolution)
    a=S.refresh(M.M.init,'full')
    assert np.allclose(a['soft_packet_ratio_quantiles'],1.,rtol=0.,atol=2e-14)
    Y,diag=step(M,M.M.init,duration,'backward')
    rows=[]
    for mode in ('frozen','debye','full'):
        r=S.refresh(Y,mode)
        D,ent,E,heat=M.collision(*M.M.weights(Y,True)[:3])
        rows.append(dict(response=r,daughter_heat=heat,entropy=ent,work=E))
    return dict(grid=grid,duration=duration,smear=smear,resolution=resolution,elastic=diag,rows=rows,
                scope='Active-occupation increment to a frozen exterior, one short elastic step.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=('fine','third'))
    p.add_argument('--duration',type=float,default=8e-8)
    p.add_argument("--smear",type=float,default=.12)
    p.add_argument("--resolution",type=int,choices=(1,2),default=1)
    a=p.parse_args()
    print(json.dumps(run(**vars(a)),indent=2,allow_nan=False))
