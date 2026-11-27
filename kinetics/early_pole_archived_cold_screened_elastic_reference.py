"""Reciprocal cold model in direct plasma-frame Q, with fixed archived screening."""
import argparse,json
import numpy as np
from scipy.special import roots_legendre, expit
from scipy.integrate import solve_ivp

from early_pole_archived_cold_elastic_packets_reference import Elastic
from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_frozen_screened_elastic_amplitude_reference import tabulated_polarization
from early_pole_frozen_elastic_direct_q_reference import traces
from early_pole_live_joint_hard_pair_packets_reference import gauss_each


class Screened(Elastic):
    def __init__(self,grid,wide,count,exposure,quadrature,cut,strength):
        Cold.__init__(self,grid,wide,0,count,True,exposure=exposure)
        self.cut=cut
        self.strength=strength
        self.reject_E=0.
        self.polarization,self.mD2=tabulated_polarization(grid)
        self.graph=self.packets(quadrature)
        self.factor=.005**4/(np.sqrt(np.pi**2*106.75/90)*(86/.05)/2.435e18)/(32*np.pi**5)

    def packets(self,quadrature):
        nr,nq,nz=quadrature
        M=self.M
        edge=[np.asarray(e) for e in M.edge_sets[1:]]
        hi=[e[-1] for e in edge]
        lo=[e[0] for e in edge]
        masses=[M.cells[f.start].m for f in M.r[1:]]
        zq,dz=roots_legendre(nz)
        nodes=np.unique(np.r_[.0001,.001,.0025,.005,.01,.02,.04,.08,.16,.32,.64,1.28,2.56,5.12,
                              self.cut,2*max(hi)])
        nodes=nodes[nodes<=2*max(hi)]
        zz,qw=gauss_each(np.log(nodes),nq)
        packets=[]
        for q,wq in zip(np.exp(zz),qw,strict=True):
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
                vacc,sc=traces(p,pb[i,None],q,X,1.,0.,self.polarization)
                W*=sc if q<self.cut else vacc
                keep=W>0.
                packets.append(np.array([np.full(P.shape,p)[keep],
                     np.broadcast_to(pb[i,None],P.shape)[keep],P[keep],B[keep],
                     np.full(P.shape,q)[keep],W[keep]]))
        pa,pb,P,B,Q,W=np.concatenate(packets,axis=1)
        ids=np.stack([np.searchsorted(e[1:],x)+f.start for x,e,f in
                      zip((pa,pb,P,B),edge+edge,M.r[1:]*2,strict=True)])
        moment=np.hypot(np.stack([pa,pb,P,B]),np.asarray(masses*2)[:,None])
        xv=(moment-M.shape[ids,0])/M.shape[ids,1]
        assert max(np.ravel(abs(xv)))<1.+2e-12 and min(W)>0.
        self.original_geometry=(pa,P,B)
        return ids[:3],xv[:3],pb,Q,W

    def collision(self,e,f,grad):
        M=self.M
        ix,u,pb,Q,W=self.graph
        d=M.r[2]
        val=e[ix,0]+e[ix,1]*u
        ec=val[0]+val[1]-val[2]
        hi=e[d,0]+e[d,1]
        assert min(np.diff(hi))>0.
        j=d.start+np.searchsorted(hi,ec)
        bounded=j<d.stop
        j=np.minimum(j,d.stop-1)
        v=(ec-e[j,0])/e[j,1]
        bare=M.shape[j,0]+M.shape[j,1]*v
        q=M.cells[d.start].pfrombare(bare)
        ok=bounded&(abs(v)<=1.)&(abs(pb-q)<=Q)&(Q<=pb+q)
        self.reject_E=max(self.reject_E,float(W[~ok].sum()/W.sum()))
        v=np.where(ok,v,0.)
        a,b,c=f[ix,0]+f[ix,1]*u
        fn=f[j,0]+f[j,1]*v
        affinity=a+b-c-fn
        A,B,C,D=map(expit,(a,b,c,fn))
        # This symmetric odds form avoids subtracting almost equal rates at tiny Q.
        size=np.sqrt(A*B*(1-C)*(1-D)*C*D*(1-A)*(1-B))
        raw=np.where(abs(affinity)<10.,2*size*np.sinh(np.clip(affinity,-10.,10.)/2),
                     A*B*(1-C)*(1-D)-C*D*(1-A)*(1-B))
        R=self.factor*self.strength*W*ok*raw
        ent=float(R@affinity)
        out=np.zeros_like(e)
        for ids,basis,flux in zip([*ix,j],[*u,v],[-R,-R,R,R],strict=True):
            out[:,0]+=np.bincount(ids,weights=flux,minlength=len(out))
            out[:,1]+=np.bincount(ids,weights=flux*basis,minlength=len(out))
        work=np.sum(grad*out)
        norm=max(1.,np.sum(abs(out)))
        assert abs(work)<5e-12*norm and ent>=-1e-13
        assert max(abs(np.sum(out[k,0])) for k in M.r[1:])<5e-12*norm
        heat=float(np.sum(grad[d]*out[d]))
        leaf=R*(ec-val[1])
        assert abs(heat-leaf.sum())<5e-12*max(1.,abs(heat))
        self.last_split=[float(leaf[Q<self.cut].sum()),float(leaf[Q>=self.cut].sum())]
        return out,ent,work,heat

    def one(self,pair):
        z=np.r_[self.M.init.ravel(),self.ns,0.,0.]
        def totals(a):
            _,_,_,r,s=self.M.weights(a[:self.j0].reshape(-1,2),True)
            _,rs,ss,_,pole=self.scalar(a[self.j0:self.j0+self.n])
            return r+rs,s+ss,pole
        e0,s0,_=totals(z)
        sol=solve_ivp(lambda _,x:self.rhs(x,pair),(0.,self.exposure),z,method='DOP853',
            first_step=min(1e-11,self.exposure/25),max_step=self.exposure/8,
            rtol=2e-10,atol=3e-13,t_eval=np.linspace(0.,self.exposure,13))
        assert sol.success,sol.message
        a=sol.y[:,-1]
        e1,s1,pole=totals(a)
        assert abs(e1-e0)<2e-13 and abs(e1-e0-a[-1])<2e-13
        assert abs(s1-s0-a[-2])<3e-11
        d=np.sum((a[:self.j0]-z[:self.j0]).reshape(-1,2)[:,0])
        assert abs(d)<2e-13
        dn=(a-z)[self.j0:self.j0+self.n]
        tot=np.sum(sol.y[self.j0:self.j0+self.n],axis=0)
        return dict(pair=pair,absolute_scalar_gain=dn.sum(),scalar_gain_by_q=dn.tolist(),
                    net_gain_timeline=(tot-tot[0]).tolist(),
                    scalar_initial_bin_counts=self.ns.tolist(),final_scalar_pole=pole,
                    gauge_mixed_count_drift=d,energy_drift=e1-e0,entropy_gain=s1-s0,
                    entropy_ledger=a[-2],rhs_calls=sol.nfev)

    def geometry_audit(self):
        M=self.M
        pa,Pa,original=self.original_geometry
        ids,u,pb,Q,W=self.graph
        e=M.weights(M.init,True)[0]
        v=e[ids,0]+e[ids,1]*u
        goal=v[0]+v[1]-v[2]
        daughter=M.r[2]
        j=daughter.start+np.searchsorted(e[daughter,0]+e[daughter,1],goal)
        bounded=j<daughter.stop
        j=np.minimum(j,daughter.stop-1)
        x=(goal-e[j,0])/e[j,1]
        out=M.cells[daughter.start].pfrombare(M.shape[j,0]+M.shape[j,1]*x)
        ok=bounded&(abs(x)<=1.)&(abs(pb-out)<=Q)&(Q<=pb+out)
        ap=((Pa-pa)*(Pa+pa)-Q*Q)/(2*Q)
        def allowed(p_after):
            bp=((pb-p_after)*(pb+p_after)+Q*Q)/(2*Q)
            tr=np.sqrt(np.maximum(0.,(pa*pa-ap*ap)*(pb*pb-bp*bp)))
            c=np.clip((2*(pa*pb-ap*bp)-1)/np.maximum(2*tr,1e-300),-1.,1.)
            return 1-np.arccos(c)/np.pi
        d=abs(allowed(original)-np.where(ok,allowed(out),0.))
        # These monitor geometry of the remap, not a reconstructed in-medium amplitude.
        V=W*(original-pb)**2
        return dict(mean_change_allowed_s_fraction=float(W@d/W.sum()),
                    energy_transfer_squared_weighted_change=float(V@d/V.sum()))

    def static(self):
        M=self.M
        rows=[]
        for mode in ('archive','common_T','parent_hotter'):
            Y=M.init.copy()
            if mode!='archive':
                for _ in range(48):
                    old=Y.copy()
                    e=M.weights(Y,True)[2]
                    for j in range(M.r[1].start,len(Y)):
                        b=1. if j<M.r[2].start or mode=='common_T' else 1/.65
                        chem=.1 if j<M.r[2].start else -.4
                        C=M.cells[j]
                        Y[j]=C.I@expit(-b*(e[j,0]+e[j,1]*C.u)+chem)
                    if np.max(abs(Y-old))<3e-15:
                        break
            e,F,E,_,_=M.weights(Y,True)
            D,S,work,heat=self.collision(e,F,E)
            row=dict(mode=mode,entropy=S,energy=work,daughter_heat=heat,
                     rhs_inf=float(np.max(abs(D))))
            if mode=='common_T':
                exact=F.copy()
                for j,mu in zip(M.r[1:],(.1,-.4),strict=True):
                    exact[j]=-E[j]
                    exact[j,0]+=mu
                eq,ent,_,h=self.collision(e,exact,E)
                row.update(exact_affinity_heat=h,exact_affinity_rhs=float(np.max(abs(eq))),
                           exact_affinity_entropy=ent,
                           fitted_affinity_residual=float(np.max(abs(exact[M.r[1].start:]-F[M.r[1].start:]))))
            if mode=='parent_hotter' and self.strength:
                assert heat>0. and S>0.
            rows.append(row)
        return rows


def run(grid,wide,n,exposure,level,cut,strength,static):
    order=[(1,1,2),(1,2,3),(2,2,3)][level]
    M=Screened(grid,wide,n,exposure,order,cut,strength)
    rows=M.static()
    out=dict(scope='Live two-flavor elastic subset with a frozen polarization, not full gauge transport.',
             grid=grid,wide=wide,order=order,matching_Q=cut,strength=strength,
             packets=M.graph[-1].size,static=rows,max_static_remap_rejection=M.reject_E)
    M.reject_E=0.
    if not static:
        out['evolution']=M.run()
        out['max_live_remap_rejection']=M.reject_E
    M.collision(*M.M.weights(M.M.init,True)[:3])
    out['initial_daughter_heat_by_split']=M.last_split
    out['initial_geometry_audit']=M.geometry_audit()
    return out


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--grid',default='fine',choices=('fine','third'))
    p.add_argument('--wide',type=int,default=1,choices=(0,1,2))
    p.add_argument('--n',type=int,default=8)
    p.add_argument('--exposure',type=float,default=8e-8)
    p.add_argument('--level',type=int,default=0,choices=(0,1,2))
    p.add_argument('--cut',type=float,default=.08,choices=(.04,.08,.16))
    p.add_argument('--strength',type=float,default=1.)
    p.add_argument('--evolve',dest='static',action='store_false',help='run the short resolved initial transient')
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
