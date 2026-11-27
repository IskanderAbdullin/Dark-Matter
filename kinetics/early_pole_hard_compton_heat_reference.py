"""Frozen reciprocal hard Compton heat on BE photons and archived fermions."""
import argparse
import json

import numpy as np
from scipy.special import expit, roots_legendre

from early_pole_frozen_dynamic_vlasov_reference import load
from early_pole_live_joint_hard_pair_packets_reference import gauss_each


def compton(aa,bb,m):
    delta=4*m*m*(1/aa-1/bb)
    return aa/bb+bb/aa+delta+delta*delta/4


def check_spin():
    m=.05
    e,c=np.meshgrid(np.geomspace(.01,1e3,73),np.linspace(-1,1,53))
    out=e/(1+e/m*(1-c))
    expect=e/out+out/e-(1-c*c)
    assert np.max(abs(compton(2*m*e,2*m*out,m)/expect-1))<2e-13


def integrate(level, bounds, temperature, angular_floor, mixed_fixture, flavor=None):
    np_,ns,nt,naz=[(2,4,4,4),(3,6,6,5),(4,8,8,5),(1,3,4,3)][level]
    profiles=[]
    for g in ('fine','third'):
        d=load(f'results/early-pole-spectrum-{g}-m06-2026-11-26.json.gz')
        profiles.extend((np.asarray(d[f'endpoint_{s}_momentum_edges']),
                     np.r_[0.,d[f'endpoint_{s}_occupations'],0.]) for s in ('parent','fermion'))
    names=['fine_parent','fine_daughter','third_parent','third_daughter']
    m=0.
    if flavor is not None:
        m=.05 if flavor=='parent' else .003
        pairs=[(k,r) for k,r in zip(names,profiles,strict=True) if k.endswith('_'+flavor)]
        names,profiles=map(list,zip(*pairs,strict=True))
    child=len(names)
    low,high=bounds
    if not 0.<low<high or temperature<=0. or not 0.<angular_floor<1.:
        raise ValueError("Positive temperature, valid radial and angular cuts required")
    base=[.45,.55,.65,.8,1.,1.25,1.32,1.56,1.7,2.2,2.8,3.7,5.,6.5,8.,10.,14.]
    knots=lambda a:np.unique(np.r_[low,[p for p in a if low<p<high],high])
    ge=knots(base)
    fe=knots(base+[p for e,occ in profiles for p,j in zip(e,abs(np.diff(occ)),strict=True)
                           if 1.25<p<1.6 and j>.012])
    p,wp=gauss_each(ge,np_)
    mid,ww=roots_legendre(np_)
    z,w=roots_legendre(ns)
    u0=(1-m*m)/(2*p)
    left=np.maximum(fe[:-1][None,:],((u0*u0-m*m)/(2*u0))[:,None])
    h=np.maximum(0.,(fe[1:]-left)/2)
    B=left[...,None]+h[...,None]*(1+mid)
    G=np.broadcast_to(p[:,None,None],B.shape)
    X=(wp*p*p/np.pi**2)[:,None,None]*(2*h[...,None]*ww*B*B/np.pi**2)
    eb=np.hypot(B,m)
    s0=np.maximum(1.,m*m+2*G*(eb-B))
    ds=np.maximum(m*m+2*G*(eb+B)-s0,0.)
    u=s0[...,None]+ds[...,None]*(1+z)/2
    rate=(u-m*m)/(2*(G*eb)[...,None])
    # Isotropic incoming angle ds/(4 p_gamma p_f); quadrature and 1/2
    # reciprocal convention give ds_panel/(16 p_gamma p_f).
    W=(X*ds/(16*G*B))[...,None]*rate*w
    good=W>0
    s=u[good]
    G,B=[np.broadcast_to(k[...,None],W.shape)[good] for k in (G,B)]
    W=W[good]
    Ein=np.hypot(B,m)
    E=G+Ein
    P=np.sqrt(np.maximum(E*E-s,0.))
    R=(s-m*m)/(2*s)
    A=np.clip(np.divide(G/R-E,P,out=np.zeros_like(P),where=P>0.),-1.,1.)
    xa,xw=roots_legendre(naz)
    q=np.hypot(knots([p for e,_ in profiles for p in e]),m)
    transfer=np.unique(np.r_[angular_floor,[i for i in (.15,.25,.4,.7,1.) if i>angular_floor]])
    tt,tw=gauss_each(transfer,nt)
    bose=lambda p,T,chem=0.: 1/np.expm1(p/T-chem)
    fg=bose(G,temperature)
    fi=[x[1][np.searchsorted(x[0],B)] for x in profiles]
    out=np.zeros((len(transfer)-1,child+2 if mixed_fixture else child,3))
    integrated_phase=0.
    for jt,(t,width) in enumerate(zip(tt,tw,strict=True)):
        zone=jt//nt
        ct=2*t-1
        center=R*(E+P*A*ct)
        span=np.maximum(R*P*np.sqrt((1-ct*ct)*(1-A*A)),1e-30)
        aa=s-m*m
        bb=aa-(aa*aa/s)*(1-t)
        angular=2*width*compton(aa,bb,m)/(16*np.pi*s)
        integrated_phase+=np.sum(W*angular)
        for i0 in range(0,len(E),224):
            sl=slice(i0,i0+224)
            H=E[sl]
            C=center[sl]
            D=span[sl]
            # Photon has smooth Bose occupations; split only its box ends,
            # and all actual outgoing fermion occupation jumps.
            chops=np.c_[([low,high]-C[:,None])/D[:,None],
                        (H[:,None]-q-C[:,None])/D[:,None],
                        np.broadcast_to([-1.,1.],(len(H),2))]
            Pha=np.sort(np.arccos(np.clip(chops,-1.,1.)),axis=1)
            dz=np.diff(Pha,axis=1)/2
            i,j=np.nonzero(dz>1e-14)
            half=dz[i,j]
            az=(Pha[i,j]+Pha[i,j+1])/2
            go=C[i,None]+D[i,None]*np.cos(az[:,None]+half[:,None]*xa)
            Eout=H[i,None]-go
            fo=np.sqrt(np.maximum((Eout-m)*(Eout+m),0.))
            valid=(go>=low)&(go<=high)&(fo>=low)&(fo<=high)
            wt=(W*angular)[sl][i,None]*(half[:,None]*xw)/(np.pi)*valid
            dE=Eout-Ein[sl][i,None]
            for k in range(out.shape[1]):
                if k<child:
                    g0=fg[sl][i,None]; g1=bose(go,temperature)
                    f0=fi[k][sl][i,None]
                    ed,fc=profiles[k]; f1=fc[np.searchsorted(ed,fo)]
                else:
                    g0=bose(G[sl][i,None],1.,-.2); g1=bose(go,1.,-.2)
                    beta=1. if k==child else 1.8
                    f0=expit(.35-beta*Ein[sl][i,None]); f1=expit(.35-beta*Eout)
                net=g0*f0*(1+g1)*(1-f1)-g1*f1*(1+g0)*(1-f0)
                u=wt*net
                out[zone,k]+=np.r_[np.sum(u*dE),np.sum(abs(u)),np.sum(abs(u*dE))]
    # Sum cross section before final support, analytic antiderivative in tau=-u/s.
    primitive=lambda r:np.log(r)+.5*r*r
    if m==0.:
        norm_exact=np.sum(W/(8*np.pi*s))*(primitive(1.)-primitive(angular_floor))
    else:
        v,dw=gauss_each(np.array([angular_floor,1.]),72)
        aa=s-m*m
        norm_exact=sum(h*np.dot(W/(8*np.pi*s),compton(aa,aa-aa*aa/s*(1-x),m))
                       for x,h in zip(v,dw,strict=True))
    defect=abs(integrated_phase/norm_exact-1)
    assert defect<2e-4
    if mixed_fixture:
        assert max(abs(out[:,child,0]))<5e-16 and min(out[:,child+1,0])>0.
    H=np.sqrt(np.pi**2*106.75/90)*(86/.05)/2.435e18
    factor=.005**4/H
    if mixed_fixture:
        names+=['common_T','photons_hotter']
    return dict(scope='Strict vacuum hard gamma-f Compton; imposed photons at a fixed BE temperature, no soft or collinear completion.',
                flavor=flavor,fermion_mass=m,
                bounds=bounds,photon_temperature=temperature,order=[np_,ns,nt,naz],
                tau_edges=transfer.tolist(),angular_normalization_error=defect,
                rows={k:dict(total_heat_into_fermion=float(v[:,0].sum()*factor),
                             heat_by_angle_band=(v[:,0]*factor).tolist(),
                             absolute_heat=float(v[:,2].sum()*factor))
                      for k,v in zip(names,out.transpose(1,0,2),strict=True)})


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--level',type=int,choices=range(4),default=0)
    p.add_argument('--bounds',type=float,nargs=2,default=[.45,8.])
    p.add_argument('--temperature',type=float,default=1.)
    p.add_argument('--angular-floor',type=float,default=.15)
    p.add_argument('--without-fixtures',dest='mixed_fixture',action='store_false')
    p.add_argument("--flavor",choices=("parent","daughter"),help="include the chosen fermion bare mass")
    a=p.parse_args()
    check_spin()
    print(json.dumps(integrate(**vars(a)),indent=2,allow_nan=False))
