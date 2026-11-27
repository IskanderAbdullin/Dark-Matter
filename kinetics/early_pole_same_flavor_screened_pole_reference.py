"""Same-flavor soft t-pole redistribution on a disjoint finite exchange window."""
import argparse
import json
import numpy as np
from scipy.special import roots_legendre
from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_frozen_elastic_direct_q_reference import traces
from early_pole_frozen_screened_elastic_amplitude_reference import tabulated_polarization
from early_pole_live_joint_hard_pair_packets_reference import gauss_each
from early_pole_same_flavor_hard_redistribution_reference import SameFlavorKernel


class SelfSoft(SameFlavorKernel):
    def __init__(self,M,grid,order,cut):
        self.M=M
        self.strength=1.
        self.cut=cut
        self.reject=[0.,0.]
        self.polarization,self.mD2=tabulated_polarization(grid)
        self.graph=[self.packets(k,order) for k in (1,2)]
        self.fac=.005**4/(np.sqrt(np.pi**2*106.75/90)*(86/.05)/2.435e18)/(32*np.pi**5)

    def packets(self,k,quadrature):
        nr,nq,nz=quadrature
        M=self.M
        edge=[np.asarray(M.edge_sets[k])]*2
        hi=[e[-1] for e in edge]
        lo=[e[0] for e in edge]
        masses=[M.cells[M.r[k].start].m]*2
        zq,dz=roots_legendre(nz)
        nodes=np.array([.0001,.001,.0025,.005,.01,.02,.04,.08,self.cut])
        nodes=np.unique(nodes[nodes<=self.cut])
        zz,qw=gauss_each(np.log(nodes),nq)
        columns=[[] for _ in range(5)]
        spans=[M.r[k]]*4
        if len(M.cells)>np.iinfo(np.int32).max:
            raise ValueError("too many cells for the packed graph")
        for q,wq in zip(np.exp(zz),qw,strict=True):
            group=[[] for _ in range(5)]
            def panel(e):
                x=np.r_[e,e-q,e+q] if q<.08 else e
                return gauss_each(np.unique(x[(x>=e[0])&(x<=e[-1])]),nr)
            incoming=[panel(e) for e in edge]
            pb,wb=incoming[1]
            for p,wp in zip(*incoming[0],strict=True):
                xl=np.maximum(np.maximum((abs(p-q)-p)/q,(lo[0]-p)/q),
                              np.maximum(-1.,(pb-hi[1])/q))
                xh=np.minimum(np.minimum(1.,(hi[0]-p)/q),
                              np.minimum((pb-abs(pb-q))/q,(pb-lo[1])/q))
                candidates=np.column_stack([np.broadcast_to((edge[0]-p)/q,(len(pb),len(edge[0]))),
                                            (pb[:,None]-edge[1])/q])
                if q<.01:
                    offs=np.r_[0.,.5*q*q/self.mD2,-.5*q*q/self.mD2]
                    candidates=np.column_stack([candidates,np.broadcast_to(offs,(len(pb),3))])
                X=np.sort(np.column_stack([xl,xh,np.clip(candidates,xl[:,None],xh[:,None])]),axis=1)
                hx=.5*np.diff(X,axis=1)
                i,j=np.nonzero((hx>1e-14)&(xh>xl)[:,None])
                if i.size==0:
                    continue
                X=.5*(X[i,j]+X[i,j+1])[:,None]+hx[i,j,None]*zq
                P=p+q*X
                B=pb[i,None]-q*X
                W=wp*wb[i,None]*wq*q*q*hx[i,j,None]*dz
                _,screen=traces(p,pb[i,None],q,X,1.,0.,self.polarization)
                W*=screen*.5
                keep=W>0.
                P,B,W=P[keep],B[keep],W[keep]
                pa=np.full(P.shape,p)
                pp=np.broadcast_to(pb[i,None],keep.shape)[keep]
                ps=(pa,pp,P,B)
                ids=np.stack([np.searchsorted(e[1:],v)+span.start for v,e,span in
                              zip(ps,edge+edge,spans,strict=True)]).astype(np.int32)
                moment=np.hypot(np.stack(ps),np.asarray(masses*2)[:,None])
                xv=(moment-M.shape[ids,0])/M.shape[ids,1]
                assert np.max(abs(xv),initial=0.)<1.+2e-12
                for column,arr in zip(group,(ids[:3].copy(),xv[:3].copy(),pp,np.full(W.shape,q),W),strict=True):
                    column.append(arr)
            for column,batch in zip(columns,group,strict=True):
                if batch:
                    column.append(np.concatenate(batch,axis=-1))
                    batch.clear()
        out=[]
        # Do not retain six raw momenta, or unused fourth-leg storage,
        # while assembling the packed physical graph for large apertures.
        for column in columns:
            out.append(np.concatenate(column,axis=-1))
            column.clear()
        assert np.min(out[-1])>0.
        return tuple(out)

def run(grid,level,cut):
    M=Cold(grid,1,0,8,True)
    q=[(1,1,2),(1,2,3)][level]
    K=SelfSoft(M.M,grid,q,cut)
    e,F,E,*_=M.M.weights(M.M.init,True)
    D,S,W=K.collision(e,F,E)
    Y=np.r_[M.M.init.ravel(),M.ns,0.,0.]
    select=slice(M.j0,M.j0+M.n)
    base=M.rhs(Y,False)[select]
    tangent={}
    for k in ('both','parent','daughter'):
        R=D.copy()
        if k=='parent':
            R[M.M.r[2]]=0.
        if k=='daughter':
            R[M.M.r[1]]=0.
        v=np.r_[R.ravel(),np.zeros(M.n+2)]
        tangent[k]=[((M.rhs(Y+h*v,False)[select]-base)/h).tolist()
                    for h in (1e-14,5e-14,2e-13)]
    for j,span in enumerate(M.M.r[1:]):
        F[span]=-E[span]*(1.+.7*j)
        F[span,0]+=[.2,-.4][j]
    x,S0,W0=K.collision(e,F,E)
    assert max(abs(z) for z in S0)<1e-18 and max(abs(x).ravel())<3e-9
    return dict(grid=grid,order=q,upper_Q=cut,entropy=S,species_energy_drift=W,
                initial_scalar_source_tangents=tangent,reject=K.reject,
                equilibrium_rhs_inf=float(np.max(abs(x))),max_separate_temperature_entropy=max(abs(z) for z in S0),
                packets=[x[-1].size for x in K.graph],
                scope='Frozen longitudinal/transverse screening of only the same-flavor t pole, soft Q window.')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=('fine','third'))
    p.add_argument('--level',type=int,choices=(0,1),default=0)
    p.add_argument('--cut',type=float,choices=(.04,.08),default=.08)
    a=p.parse_args()
    print(json.dumps(run(**vars(a)),indent=2,allow_nan=False))
