"""Frozen, finite-momentum scalar scattering through a below-pole Yukawa channel."""
import argparse
import gzip
import json
import numpy as np
from scipy.optimize import brentq
from scipy.special import expit, roots_legendre

from early_pole_archive_scalar_rest_below_pole_reference import spin,fold
from early_pole_decay_mode_clock_reference import G,KAPPA
from early_pole_hard_yukawa_photoproduction_reference import dot,slash,adj,I
from early_pole_native_aperture_witness_reference import GRID,occup,mode


def fixture(M,m,ms):
    errors=[]
    for s in np.geomspace(1.05*(m+ms)**2,(.90*M)**2,13):
        u=np.sqrt(s)
        r=np.sqrt((s-(m+ms)**2)*(s-(m-ms)**2))/(2*u)
        ob,os=np.hypot(r,m),np.hypot(r,ms)
        f=np.array([ob,0.,0.,-r])
        pair=np.array([u,0.,0.,0.])
        for c in (-.99,-.4,.2,.8,.999):
            z=np.array([os,r*np.sqrt(1-c*c),0.,r*c])
            out=pair-z
            cross=f-z
            V=(slash(pair)+M*I)/(s-M*M)+(slash(cross)+M*I)/(dot(cross,cross)-M*M)
            d=float((np.trace((slash(out)+m*I)@V@(slash(f)+m*I)@adj(V))/2).real)
            reduced=spin(s,-2*r*r*(1-c),M,m,ms)
            errors.append(abs(d/reduced-1))
    assert max(errors)<4e-12,max(errors)
    return max(errors)


def integrate(a,q,order,fraction):
    M,m=.05,.003
    ms=a['checkpoints'][-1]['scalar_pole_over_T']
    E0=np.hypot(ms,q)
    pe=np.asarray(a['endpoint_fermion_momentum_edges'])
    ferm=np.asarray(a['endpoint_fermion_occupations'])
    qs=np.asarray(a['scalar_momentum_edges'])
    be=np.maximum(0.,np.asarray(a['endpoint_scalar_occupations']))
    ftag=float(occup(np.array(q),qs,be)) if q else float(be[0])
    assert ftag>0
    x,w=roots_legendre(order)
    aa,aw=roots_legendre(order+2)
    az=np.cos(np.pi*(1+aa))
    minE,maxE=np.hypot([.45,5.],m)
    cap=(fraction*M)**2
    def S(E,sgn):
        return m*m+ms*ms+2*(E*E0+sgn*q*np.sqrt((E-m)*(E+m)))
    # Resolve both radial thresholds where the angular range changes.
    cuts=np.hypot(pe,m)
    c=[minE,maxE,*cuts[(minE<cuts)&(cuts<maxE)]]
    for sign in (-1.,1.):
        if S(minE,sign)<cap<S(maxE,sign):
            c.append(brentq(lambda E:S(E,sign)-cap,minE,maxE,xtol=2e-15))
    # The radial archive already resolves all incoming shell steps.
    cuts=np.unique(c)
    he=np.diff(cuts)[:,None]/2
    energies=((cuts[1:]+cuts[:-1])[:,None]/2+he*x).ravel()
    energy_weights=(he*w).ravel()
    qsc=qs*qs/(np.hypot(qs,ms)+ms)
    foe=np.hypot(pe,m)
    F=np.zeros((2,4))
    max_thermal=0.
    for E,dE in zip(energies,energy_weights,strict=True):
        p=np.sqrt((E-m)*(E+m))
        if q:
            ml=max(-1.,(m*m+ms*ms+2*E*E0-cap)/(2*p*q))
            if ml>=1.:continue
        else:
            if S(E,0.)>=cap:continue
            ml=-1.
        fa=occup(np.array(p),pe,ferm)
        angles=(1+ml)/2+(1-ml)/2*x
        angle_weights=(1-ml)*w/4
        for mu,W in zip(angles,angle_weights,strict=True):
            invariant=E*E0-p*q*mu
            s=m*m+ms*ms+2*invariant
            lab=np.sqrt(p*p+q*q+2*p*q*mu)
            root=np.sqrt(s)
            r=np.sqrt((s-(m+ms)**2)*(s-(m-ms)**2))/(2*root)
            Ec=np.hypot(ms,r)
            ET=E+E0
            gamma=ET/root
            speed=lab/ET
            D=lab*r/root
            # Cancellation-free kinetic-energy lower turning point.
            t0=(gamma*(r-speed*Ec))**2/(gamma*(Ec-speed*r)+ms)
            t1=t0+2*D
            lo=max(t0,ET-ms-maxE)
            hi=min(t1,ET-ms-minE)
            if hi<=lo:continue
            par=(q*q+p*q*mu)/lab
            orth=p*q*np.sqrt(1-mu*mu)/lab
            cc=gamma*(par-speed*E0)
            c0,s0=np.array([cc,orth])/np.hypot(cc,orth)
            cuts=np.r_[qsc,ET-ms-foe]
            cuts=np.unique(np.r_[lo,hi,cuts[(lo<cuts)&(cuts<hi)]])
            h=np.diff(cuts)[:,None]/2
            K=((cuts[1:]+cuts[:-1])[:,None]/2+h*x).ravel()
            WK=(h*w).ravel()
            z=np.clip((K-t0)/D-1,-1.,1.)
            cb=np.clip(c0*z[:,None]+s0*np.sqrt(1-z*z)[:,None]*az,-1.,1.)
            matrix=spin(s,-2*r*r*(1-cb),M,m,ms)@(aw/2)
            assert np.min(matrix)>=-1e-12
            benergy=ET-ms-K
            qb=np.sqrt(np.maximum(0.,(benergy-m)*(benergy+m)))
            qs1=np.sqrt(K*(K+2*ms))
            fb=occup(qb,pe,ferm)
            b=occup(qs1,qs,be)
            em=fb*(1-fa)*b*(1+1/ftag)
            loss=fa*(1-fb)*(1+b)
            vac=fa*(1-fb)
            net=b*(fb-fa)+b*fb*(1-fa)/ftag-vac
            Fa,Fb=expit(-E),expit(-benergy)
            Ga,Gb=1/np.expm1(.003+E0-ms),1/np.expm1(.003+K)
            lhs=Fb*(1-Fa)*Gb*(1+1/Ga)
            rhs=Fa*(1-Fb)*(1+Gb)
            max_thermal=max(max_thermal,float(np.max(abs(lhs/rhs-1))))
            vel=np.sqrt(max(0.,invariant**2-m*m*ms*ms))/(E*E0)
            pref=4*p*E*dE*W*vel/(2*np.pi**2)/(32*np.pi*s*D)
            weight=pref*matrix*WK
            for k,tag in enumerate((qs1<1e-5,qs1>=1e-5)):
                F[k]+=np.array([weight[tag]@u[tag] for u in (em,loss,net,vac)])
    assert max_thermal<1e-13
    H=np.sqrt(np.pi**2*106.75/90)*(86./a['x_end'])/2.435e18
    fac=H*(16*np.pi*KAPPA*M)**2
    pn=np.asarray(a['endpoint_parent_momentum_edges'])
    pf=np.asarray(a['endpoint_parent_occupations'])
    if q:
        e,v,_=mode(M,ms,m,q,pn,pf,pf,pe,ferm,0.)
        direct=G*KAPPA*M*((M+m)**2-ms**2)/(e*q)*(v[0]-v[1]+v[2]/ftag)
    else:
        ref=fold(a,order,fraction)
        direct=ref['normalized_direct_scalar_rest_net_per_reference_H']
        expected=np.asarray(ref['cold_and_exterior_columns_emission_loss_net_vacuum_per_rest_occupation_without_y4'])
        np.testing.assert_allclose(F,expected,rtol=4e-10,atol=1e-14)
    return dict(q=q,order=order,pole_guard=fraction,tagged_bose_occupation=ftag,
                cold_and_exterior_columns_emission_loss_net_vacuum_without_y4=F.tolist(),
                net_by_band_per_reference_H=(fac*F[:,2]).tolist(),direct_per_reference_H=float(direct),
                maximum_thermal_balance_defect=max_thermal,independent_dirac_spin_error=fixture(M,m,ms))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=('fine','third'))
    p.add_argument('--order',type=int,default=8)
    p.add_argument('--q',type=float,nargs='+',default=[0.,4.1e-6])
    o=p.parse_args()
    with gzip.open(GRID/f'early-pole-spectrum-{o.grid}-m06-2026-11-26.json.gz','rt') as fh:
        a=json.load(fh)
    print(json.dumps(dict(grid=o.grid,cases=[integrate(a,q,o.order,c) for q in o.q for c in (.70,.82,.90)]),
                     indent=2,allow_nan=False))
