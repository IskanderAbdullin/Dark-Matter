"""Finite-support frozen real-photon gain/loss with a prescribed soft scalar."""
import argparse
import gzip
import json

import numpy as np
from scipy.special import roots_legendre

from early_pole_hard_yukawa_photoproduction_reference import amplitudes, com, invariant_total, dot
from early_pole_native_aperture_witness_reference import GRID, PROFILES, mode, occup, refine


def edges_between(e,lo,hi):
    return np.unique(np.r_[lo, e[(lo<e)&(e<hi)],hi])


def orientation(s,M,m,ms,z,cstar,st,az,trace=False):
    angle=z[:,None]*cstar+np.sqrt(np.maximum(0.,1-z*z))[:,None]*st*np.cos(az)
    angle=np.clip(angle,-1.,1.)
    if not trace:
        return invariant_total(s,angle,M,m,ms)
    p,k,b,_,_=com(s,angle,M,m,ms)
    return sum(amplitudes(p,k,b,M,m,v)[0] for v in
               (np.array([0,1,0,0]),np.array([0,0,1,0])))


def coordinate_fixture(M,m,ms,q):
    worst=0.
    Q=np.longdouble
    q,ms,M,m=map(Q,(q,ms,M,m))
    for p in map(Q,(.6,1.8,4.5)):
        E=np.hypot(p,m)
        Es=np.hypot(q,ms)
        for u in map(Q,(-.8,.6)):
            b=np.array([E,0.,0.,p],dtype=Q)
            phi=np.array([Es,q*np.sqrt(1-u*u),0.,q*u],dtype=Q)
            S=b+phi
            lab=np.linalg.norm(S[1:])
            direction=S[1:]/lab
            s=dot(S,S)
            if s<=M*M:
                continue
            root=np.sqrt(s)
            gamma=S[0]/root
            velocity=lab/S[0]
            par=b[1:]@direction
            axis=np.array([-direction[2],0.,direction[0]])
            normal=np.cross(direction,axis)
            bperp=b[1:]@axis
            bpar=gamma*(par-velocity*E)
            bc=np.hypot(bpar,bperp)
            for z in (-.6,.2,.9):
                a=.83
                n=z*direction+np.sqrt(1-z*z)*(np.cos(a)*axis+np.sin(a)*normal)
                kstar=(s-M*M)/(2*root)
                kin=np.r_[kstar,n*kstar]
                ray=np.r_[gamma*kstar*(1+velocity*z),
                          kin[1:]+kstar*((gamma-1)*z+gamma*velocity)*direction]
                angle=z*bpar/bc+np.sqrt(1-z*z)*np.cos(a)*bperp/bc
                pcom,kcom,bcom,_,_=com(s,np.array(angle),M,m,ms)
                # The equality is invariant; lab forward denominator can be
                # checked independently of the cos(theta) and boost inversion.
                errs=[abs(dot(S-ray,S-ray)-M*M),abs(dot(ray,ray)),
                      abs(dot(ray,b)-dot(kcom,bcom)),abs(dot(ray,S-ray)-dot(kcom,pcom))]
                worst=max(worst,*errs)
    assert worst<1e-12,worst
    return float(worst)


def frozen(q,ms,M,m,ie,fi,oe,fo,n,temperature=1.,chemical=False,pole_panels=True,trace=False):
    # Inverse-coordinate phase space: scalar(q)+fermion(m,p)
    # -> photon(k)+fermion(M,E'-k). All momentum cuts are plasma-frame.
    Es=np.hypot(ms,q)
    x,w=roots_legendre(n)
    a,aw=roots_legendre(n+2)
    az,aw=np.pi*(1+a),np.pi*aw
    d=edges_between(oe,.45,5.)
    h=np.diff(d)[:,None]/2
    pp=((d[1:]+d[:-1])[:,None]/2+h*x).ravel()
    wp=(h*w).ravel()
    Nx=np.hypot(ie,M)
    lines=np.zeros(3)
    for p,W in zip(pp,wp,strict=True):
        E=np.hypot(p,m)
        ET=Es+E
        f=occup(np.array(p),oe,fo) if not chemical else 1/(np.exp((E+.10)/temperature)+1)
        mx=(ms*ms+m*m+2*Es*E-M*M)/(2*p*q)
        if mx<=-1:
            continue
        mhi=min(mx,1.)
        mu=(mhi-1)/2+(mhi+1)/2*x
        wm=(mhi+1)/2*w
        for u,Wu in zip(mu,wm,strict=True):
            s=ms*ms+m*m+2*(Es*E-p*q*u)
            root=np.sqrt(s)
            kstar=(s-M*M)/(2*root)
            P=np.sqrt(p*p+q*q+2*p*q*u)
            gamma=ET/root
            width=kstar*P/root
            center=gamma*kstar
            lo=max(.5,center-width,ET-np.hypot(5.,M))
            hi=min(4.,center+width,ET-np.hypot(.45,M))
            if hi<=lo:
                continue
            parallel=p*(p+q*u)/P
            trans=p*q*np.sqrt(1-u*u)/P
            bst=gamma*(parallel-P/ET*E)
            ca,sa=np.array([bst,trans])/np.hypot(bst,trans)
            raw=ET-Nx
            if pole_panels:
                eb=(s+m*m-ms*ms)/(2*root)
                pb=np.sqrt(max(0.,(eb-m)*(eb+m)))
                if pb>0.:
                    rho=m*m/(pb*(eb+pb))
                    # Resolve the c≈1 virtual-fermion denominator; retain
                    # native final-population discontinuities independently.
                    raw=np.r_[raw,center+width*(1-rho*np.array([.0625,.25,1,4,16,64,256]))]
            edges=edges_between(raw,lo,hi)
            dh=np.diff(edges)[:,None]/2
            K=((edges[1:]+edges[:-1])[:,None]/2+dh*x).ravel()
            WK=(dh*w).ravel()
            z=(K-center)/width
            azint=orientation(s,M,m,ms,z,ca,sa,az,trace)@aw
            F=occup(np.sqrt(np.maximum(0.,(ET-K)**2-M*M)),ie,fi) if not chemical else 1/(np.exp((ET-K-.12)/temperature)+1)
            B=1/np.expm1(K/temperature)
            forward=F*B*(1-f)
            inverse=f*(1-F)*(1+B)
            if chemical:
                assert np.max(abs((inverse/forward)/np.exp((Es-.22)/temperature)-1))<2e-13
            pref=p*p*W*Wu/E*kstar/root/width/(256*np.pi**4*Es)
            lines+=pref*np.r_[WK@(azint*forward),WK@(azint*inverse),WK@(azint*(forward-inverse))]
    return lines*2  # charge conjugate; both spins and photon pols summed by the trace.


def parent_input(grid,a,witness_mix):
    if not 0<=witness_mix<=1:
        raise ValueError("mix must be in [0,1]")
    p=np.asarray(a['endpoint_parent_momentum_edges'])
    P=np.asarray(a['endpoint_parent_occupations'])
    if witness_mix:
        profile=next(x for x in json.loads(PROFILES.read_text()) if x['grid']==grid)
        p,P,change,_=refine(p,P,profile['parent'],a['x_end'])
        P=P+witness_mix*(change-P)
    return p,P


def case(g,n,q,pole_panels=True,trace=False,witness_mix=0.):
    with gzip.open(GRID/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz','rt') as f:
        a=json.load(f)
    ms=a['checkpoints'][-1]['scalar_pole_over_T']
    heavy,light=.05,.003
    p,P=parent_input(g,a,witness_mix)
    b=np.asarray(a['endpoint_fermion_momentum_edges']); B=np.asarray(a['endpoint_fermion_occupations'])
    ph=frozen(q,ms,heavy,light,p,P,b,B,n,pole_panels=pole_panels,trace=trace)
    bh=frozen(q,ms,light,heavy,b,B,p,P,n,pole_panels=pole_panels,trace=trace)
    E,I,_=mode(heavy,ms,light,q,p,P,P,b,B,0.)
    vacuum=((heavy+light)**2-ms**2)/(4*np.pi*q*E)*np.array([I[2],I[0]-I[1]])
    assert np.all(np.isfinite(np.r_[ph,bh,vacuum]))
    assert np.max(abs(np.r_[ph[0]-ph[1]-ph[2],bh[0]-bh[1]-bh[2]]))<1e-10
    thermal=frozen(q,ms,heavy,light,p,P,b,B,min(n,4),chemical=True)
    assert abs((thermal[1]/thermal[0])/np.exp(E-.22)-1)<5e-13
    return dict(grid=g,q=q,quadrature=n,
                max_independent_kinematic_invariant_defect=max(coordinate_fixture(heavy,light,ms,q),
                                                            coordinate_fixture(light,heavy,ms,q)),
                parent_plus_photon=ph.tolist(),daughter_plus_photon=bh.tolist(),
                direct_seed_and_stimulated_without_y2=vacuum.tolist(),
                associated_seed_and_stimulated_relative_at_g005=[
                    float(x*.005**2/y) if abs(y)>1e-30 else None for x,y in zip((ph+bh)[[0,2]],vacuum)],
                **(dict(witness_mix=witness_mix) if witness_mix else {}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=('fine','third'))
    p.add_argument('--q',type=float,default=4.1e-6)
    p.add_argument('--order',type=int,default=6)
    p.add_argument('--no-pole-panels',action='store_true')
    p.add_argument('--trace',action='store_true')
    p.add_argument('--witness-mix',type=float,default=0.)
    a=p.parse_args()
    print(json.dumps(case(a.grid,a.order,a.q,not a.no_pole_panels,a.trace,a.witness_mix),indent=2,allow_nan=False))
