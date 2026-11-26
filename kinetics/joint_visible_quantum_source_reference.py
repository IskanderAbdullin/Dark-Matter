"""Signed two-species portal source on finite supports at a fixed state."""
import numpy as np
from scipy.special import roots_legendre
from empty_lepton_pair_spectrum_reference import Source
from three_species_decay_onsager_reference import m
from joint_visible_fd_source_reference import Visible
from joint_visible_heavy_adiabatic_reference import shell_projection,snapshot_functions

def signed_shells(T,theta,edges,order=(6,12,14),npolar=6):
    src=Source(order,dark_mass=m)
    q,p,v,N,A=src.event_moments(T)
    hs=(p*np.sqrt(v)[:,None]/2).ravel()
    q=q.ravel(); N=N.ravel(); A=A.ravel()
    es=np.hypot(edges,m)
    x,w=roots_legendre(npolar)
    out=np.zeros((2,2,len(es)+1))
    sums=np.zeros((2,2)); vac=np.zeros(2)
    for i in range(0,len(q),256):
        s=slice(i,i+256)
        center=q[s,None]*.5
        h=hs[s,None]
        cuts=np.sort(np.clip(np.concatenate([np.full_like(center,-1.),np.full_like(center,1.),
             (es[None,:]-center)/h,(q[s,None]-es[None,:]-center)/h],axis=1),-1.,1.),axis=1)
        mid=(cuts[:,1:]+cuts[:,:-1])*.5
        dh=np.diff(cuts,axis=1)*.5
        z=mid[:,:,None]+dh[:,:,None]*x
        z=z.reshape(len(center),-1)
        aw=(dh[:,:,None]*w).reshape(len(center),-1)
        E1=center+h*z
        E2=q[s,None]-E1
        W=(N[s,None]+A[s,None]*(3*z*z-1)*.5)*aw
        assert np.all(W>=0.)
        vac+=np.array([W.sum(),np.sum(W*q[s,None]*.5)])
        b,c=theta(E1),theta(E2)
        lf=-np.logaddexp(0.,-b)-np.logaddexp(0.,-c)
        aff=q[s,None]/T-b-c
        signed=-np.sign(aff)*np.exp(lf+np.maximum(0.,aff))*(-np.expm1(-abs(aff)))
        idx=np.searchsorted(es,E1,'right').ravel()
        for j in range(2):
            pos=W*np.maximum(signed*(-1)**j,0.)
            out[j,0]+=np.bincount(idx,pos.ravel(),minlength=len(es)+1)
            out[j,1]+=np.bincount(idx,(pos*(E1-m)).ravel(),minlength=len(es)+1)
            sums[j]+=(np.sum(pos),np.sum(pos*q[s,None]/2))
    assert np.max(abs(out[:,0].sum(axis=1)-sums[:,0]))/vac[0]<1e-12
    assert np.max(abs(out[:,1].sum(axis=1)+m*out[:,0].sum(axis=1)-sums[:,1]))/vac[1]<1e-12
    return out

class QuantumVisible(Visible):
    def __init__(self, order=(6,12,14), polar=16, direct_angle=5, decay_angle=4):
        super().__init__(order=order,polar=polar)
        self.direct_angle=direct_angle
        self.decay_angle=decay_angle

    def rates_at_state(self,T,q,a,th):
        x=dict(lab_momentum_edges_GeV=np.asarray(q)/a,
               entropy_variable_scalar=np.asarray(th[0]),
               entropy_variable_fermion=np.asarray(th[1]))
        fp,sg,fc,cg=snapshot_functions(x)
        def tc(E):
            f=fc(E)
            return np.where(f>0.,np.log1p(-f)-np.log(np.where(f>0.,f,1.)),1000.)
        sd=signed_shells(T,tc,np.asarray(q)/a,self.direct.orders,self.direct_angle)
        sh=shell_projection(T,fp,fc,np.asarray(q)/a,sg,cg,order=self.heavy.orders,polar=self.polar,decay=self.decay_angle)
        res=np.stack([sd,sh[1],sh[0]])
        assert np.all(res>=0.)
        return res

def selftest():
    from scipy.special import expit
    T=.045
    ed=.05*np.r_[0.,np.geomspace(.02,16.,20)]
    vac=signed_shells(T,lambda E:np.full_like(E,1000.),ed)
    ref=Visible((6,12,14)).rates(T,ed,1.)[0]
    order=(6,12,14)
    eq=signed_shells(T,lambda E:E/T,ed)
    direct_neg=signed_shells(T,lambda E:E/T-.3,ed)
    neg=shell_projection(T,lambda E:1/np.expm1(E/T+.1),
                         lambda E:expit(.3-E/T),ed,order=order,decay=4)
    l1=np.sum(abs(vac[0]-ref),axis=1)/ref.sum(axis=1)
    k=eq.sum(axis=(0,2))/ref.sum(axis=1)
    assert max(l1)<3e-12 and max(k)<3e-13
    assert direct_neg[0].sum()==0. and np.all(direct_neg[1].sum(axis=-1)>0.)
    assert neg[:,0].sum()==0. and np.all(neg[:,1].sum(axis=-1)>0.)
    return dict(vacuum_shell_L1=l1.tolist(),equilibrium_shell_L1=k.tolist(),
                synthetic_net_removal_had_only_negative_rates=True)

if __name__=="__main__":
    import json
    print(json.dumps(selftest(),indent=2,allow_nan=False))
