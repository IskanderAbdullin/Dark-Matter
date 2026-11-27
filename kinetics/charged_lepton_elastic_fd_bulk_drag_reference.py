"""Local common-temperature FD momentum and energy response of the elastic portal."""
import json

import numpy as np
from scipy.special import expit, roots_legendre

from boost_axis_elastic_heat_reference import integrate, leptons
from fd_number_heat_state_reference import evaluate
from timelike_portal_onsager_coverage_reference import V, mchi, gprod


def response(T,eta,ml,n,tail,drift=0.):
    """Independent relative-azimuth sum; integrate energy and invariant transfer together."""
    M,L=mchi/T,ml/T
    z,w=roots_legendre(n)
    theta=(np.arange(2*n)+.5)*np.pi/n
    st=np.sqrt(1-z*z)
    cosphi=np.cos(theta)
    cg=z[:,None,None]*z[None,:,None]+st[:,None,None]*st[None,:,None]*cosphi
    dz=z[None,:]-z[:,None]
    transverse=st[:,None,None]**2+st[None,:,None]**2-2*st[:,None,None]*st[None,:,None]*cosphi
    polar=w[:,None]*w[None,:]
    X=(1+z)*np.sqrt(tail)/2
    dX2=w*np.sqrt(tail)*X
    tot=np.zeros(4 + bool(drift))
    for u,du in zip(M+L+X*X,dX2,strict=True):
        s=u*u
        pcm2=(u-M-L)*(u+M+L)*(u+M-L)*(u-M+L)/(4*s)
        pcm=np.sqrt(pcm2)
        tm=2*pcm2*(cg-1)
        A=s-M*M-L*L
        amp=2*(gprod/V**2)**2*(A*A+(A+tm)**2+2*(M*M+L*L)*tm)/(1-tm*(T/V)**2)**2
        E=u+X*X
        pp=np.sqrt((E-u)*(E+u))
        ed=(E[:,None]*(s+M*M-L*L)/(2*u)+pp[:,None]*pcm*z)/u
        el=E[:,None]-ed
        fd,fs=expit(eta-ed),expit(-el)
        differential=(fd*fs)[:,:,None]*((1-fd)*(1-fs))[:,None,:]
        dE=ed[:,None,:]-ed[:,:,None]
        W=dX2*pp*pcm2/u
        base=polar*amp.mean(axis=-1)
        t_weight=polar*np.mean(amp*(-tm),axis=-1)
        # Independent boost reconstruction of |Delta p|^2, including relative
        # azimuth and longitudinal gamma; no use of invariant transfer t here.
        direct=pcm2*(E/u)[:,None,None]**2*dz*dz
        transverse_w=polar*np.mean(amp*transverse,axis=-1)*pcm2
        e=dE*dE
        for k,kernel in enumerate((e*base,t_weight,(e*base+t_weight)/3,
                                    (direct*base+transverse_w)/3)):
            tot[k]+=du*float(W@np.sum(differential*kernel,axis=(1,2)))
        if drift:
            # Fully displaced dark FD factors (both in and out), unlike
            # the linear moment above. Average flow along six directions
            # giving the rotationally invariant second moment as u->0.
            star=(s+M*M-L*L)/(2*u)
            pz=(pp[:,None]*star+E[:,None]*pcm*z)/u
            before=(pcm*st[:,None,None],0.,pz[:,:,None,None])
            after=(pcm*st[None,:,None]*cosphi,pcm*st[None,:,None]*np.sin(theta),pz[:,None,:,None])
            iE,fE=ed[:,:,None,None],ed[:,None,:,None]
            iS,fS=fs[:,:,None,None],fs[:,None,:,None]
            gamma=1/np.sqrt(1-drift*drift)
            fwd_visible=iS*(1-fS)
            rev_visible=fS*(1-iS)
            for pi,pf in zip(before,after,strict=True):
                kick=pf-pi
                odd=0.
                for vel in (-drift,drift):
                    iD=expit(eta-gamma*(iE-vel*pi))
                    fD=expit(eta-gamma*(fE-vel*pf))
                    F=fwd_visible*iD*(1-fD)-rev_visible*fD*(1-iD)
                    odd+=np.sign(vel)*F
                gk=np.mean(-odd*kick*amp,axis=-1)*polar/(6*drift)
                tot[4]+=du*float(W@np.sum(gk,axis=(1,2)))
    tot*=16*T**9/(512*np.pi**5)
    if min(tot)<=0. or abs(tot[2]/tot[3]-1)>5e-14:
        raise ValueError("Independent lab-momentum reconstruction failed.")
    return tot


def run():
    rows=[]
    for T,eta in ((.005,0.),(.015,0.),(.03,1.),(.05,0.),(.1,0.),(.3,0.),
                  (.6882487043,-1.),(.6882487043,0.),(.6882487043,2.)):
        a=np.array([response(T,eta,m,18,29.) for m in leptons])
        b=np.array([response(T,eta,m,26,36.) for m in leptons])
        x,y=a.sum(axis=0),b.sum(axis=0)
        old=sum(integrate(T,T,eta,m,n=26,cut=36.,linear=True) for m in leptons)
        test1=float(np.max(abs(y/x-1)))
        test2=abs(y[0]/old-1)
        if test1>7e-5 or test2>3e-13:
            raise ValueError(f"Unexpected heat or quadrature error: {(T,eta,test1,test2)}")
        dz=evaluate(T,eta)
        Kp,Kq=y[2],y[0]
        rows.append(dict(T=T,eta_dark=eta,momentum_per_velocity=Kp,heat_per_temperature_ratio=Kq,
                         drag_to_heat=Kp/Kq,momentum_fraction_by_lepton=(b[:,2]/Kp).tolist(),
                         heat_fraction_by_lepton=(b[:,0]/Kq).tolist(),invariant_spacelike_per_ratio=y[1],
                         momentum_to_heat_rate_under_common_dark_bulk=(Kp/(dz["rho"]+dz["pressure"])) /
                             (Kq/(T*dz["cv_fixed_number"])),
                         angular_and_cutoff_discrepancy=test1,previous_heat_integral_discrepancy=test2))
    finite=[]
    for T,eta in ((.05,0.),(.6882487043,2.)):
        # Equal quadrature for the linear and the independently displaced
        # weights; check decrease under halving the probe velocity.
        err=[]
        for v in (.001,.0005):
            s=sum((response(T,eta,m,15,30.,drift=v) for m in leptons),np.zeros(5))
            err.append(abs(s[4]/s[2]-1))
        if err[1]>6e-6 or err[0]<2.8*err[1]:
            raise ValueError(f"Finite drift did not approach the expected linear response: {T,eta,err}")
        finite.append(dict(T=T,eta_dark=eta,probe_velocities=[.001,.0005],relative_discrepancy=err))
    return dict(scope="Elastic charged-lepton vacuum channel, fixed common temperature, linear local drift. No kinetic decoupling or transfer.",
                rows=rows,independent_finite_velocity_checks=finite)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
