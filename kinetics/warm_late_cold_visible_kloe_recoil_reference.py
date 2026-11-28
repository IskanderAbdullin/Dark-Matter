"""Finite measured charged-pion window in the vacuum cold recoil fold."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.special import k1,roots_legendre
from archived_collisionless_streaming_reference import Background
from decay_scalar_ir_relaxation_reference import limits
from kloe_pion_band_portal_reference import table,POINT_MUON
from late_radiation_background_reference import TabulatedRadiation
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import Es,M,P,mp
from timelike_portal_onsager_coverage_reference import V,width,mpi,factor
from warm_late_cold_visible_vacuum_dose_reference import run as base


def rates(T,pc,grid,cs,sigma,n,*,spectral_profile=None):
    C,L,span=limits(pc)
    z,w=roots_legendre(n)
    x=(z+1)/2;w=w/2
    ss=(grid[:,None]+.01*(x-.5)).ravel()
    cs=(np.repeat(cs,n) if spectral_profile is None else np.asarray(spectral_profile(ss)))
    assert cs.shape==ss.shape and np.min(cs)>0.
    err=np.repeat(sigma,n)
    W=np.tile(.01*w,len(grid))
    R=np.sqrt(ss)
    v=(R-2*M)*(R+2*M)/ss
    beta=np.sqrt(v)
    Q0=np.where(R>=2*L,R,(4*L*L+v*ss)/(2*L+2*beta*np.sqrt((L-M)*(L+M))))
    assert np.all(Q0>=R*(1-1e-14))
    Q0=np.maximum(Q0,R)
    pref=2*factor*(ss+2*M*M)*ss*beta*W*ss/POINT_MUON/((ss-V*V)**2+(V*width)**2)
    vpi=(R-2*mpi)*(R+2*mpi)/ss
    assert np.all(vpi>0.)
    totals=np.zeros(3)
    bounds=(0,.25,1,3,6,10,18,31,49,72)
    for lo,hi in zip(bounds[:-1],bounds[1:],strict=True):
        Q=Q0[:,None]+T*(lo+(hi-lo)*x)
        dQ=T*(hi-lo)*w
        gP=np.sqrt((Q-R[:,None])*(Q+R[:,None]))
        h=gP*beta[:,None]/2
        Mn=np.maximum(L,Q/2-h);Mx=np.minimum(L+span,Q/2+h)
        D=np.maximum(Mx-Mn,0.)
        e=Mn[:,:,None]+D[:,:,None]*x
        mu=(e-Q[:,:,None]/2)/h[:,:,None]
        gamma=e/M
        g=np.sqrt((gamma-1)*(gamma+1))
        E0=Es/(gamma+g)+(mp*mp/(Es+P))*g
        angular=np.maximum(C-E0,0.)/(2*P*g)
        cos=gP[:,:,None]*np.sqrt(vpi)[:,None,None]*z/(2*T)
        u=Q[:,:,None]/(2*T)
        boltz=1./((-np.expm1(-u+cos))*(-np.expm1(-u-cos)))
        A=boltz@w
        B=-(2*v/(3-v))[:,None]*(boltz@(w*(3*z*z-1)/2))
        result=((A[:,:,None]+B[:,:,None]*(3*mu*mu-1)/2)*angular)@w
        assert np.min(result)>=0.
        WQ=pref[:,None]*np.exp(-Q/T)*gP*dQ*D/(2*h)
        mb=(WQ*(angular@w)).sum(axis=1)
        be=(WQ*result).sum(axis=1)
        totals+=np.array([cs@mb,cs@be,err@be])
    assert 0<totals[0]<=totals[1]*(1+3e-15)
    return totals


def independent_controls():
    # Explicit azimuth average of the spin-zero pion vector current.
    x,w=roots_legendre(24)
    phi=2*np.pi*np.arange(32)/32.
    cos=x[:,None,None]*x[None,:,None]+np.sqrt((1-x[:,None,None]**2)*(1-x[None,:,None]**2))*np.cos(phi)
    P2=(3*x*x-1)/2
    err=0.
    for v in (.01,.5,.99):
        raw=(1-v*cos*cos).mean(axis=2)/(1-v/3)
        folded=1-2*v/(3-v)*P2[:,None]*P2
        err=max(err,float(np.max(abs(raw-folded))))
    assert err<3e-15
    # Compare absolute normalization of the uncut MB q integral:
    # integral_R^infty sqrt(Q^2-R^2) exp(-Q/T)dQ = R*T*K1(R/T).
    T=.05
    d=(x+1)/2;wp=w/2
    check=[]
    t=np.array([0,.4,1.,2.,4.,7.,12.,20.,40.,72.])
    q=(t[:-1,None]+np.diff(t)[:,None]*d*d).ravel()
    dw=(np.diff(t)[:,None]*2*d*wp).ravel()
    for R in (.317,.45,.775,.965):
        Q=R+T*q
        numeric=T*(dw@(np.sqrt((Q-R)*(Q+R))*np.exp(-Q/T)))
        exact=R*T*k1(R/T)
        check.append(float(abs(numeric/exact-1.)))
    assert max(check)<3e-13
    return dict(max_direct_azimuth_vs_p2_error=err,max_uncut_MB_q_weight_vs_Bessel_error=max(check))


def run(eos,archive,data):
    grid,cs,e=table(data)
    bg=Background(eos)
    comparison=base(eos,archive)["comparison_archive_held_equality_cold_yield"]
    TR=.05
    rad=TabulatedRadiation(eos,TR)
    rows=[]
    for o,n in ((8,14),(12,20)):
        u,w=roots_legendre(n)
        l,h=np.log([.012,TR])
        vals=[]
        for T in np.exp((h+l)/2+(h-l)*u/2):
            a=bg.scale(T)
            _,rho=rad.state(a/bg.scale(TR))
            s=(2*np.pi**2/45)*np.exp(bg.gs(np.log(T)))*T**3
            d=(3+bg.gs.derivative()(np.log(T)))/3
            vals.append(rates(T,.1*bg.T0/a,grid,cs,e,o)*d/(s*np.sqrt(rho/3)/MPL))
        y=(h-l)/2*(w@np.asarray(vals))
        rows.append(dict(kinetic_order=o,time_order=n,MB=float(y[0]),BE=float(y[1]),
                         linear_published_marginal_variation=float(y[2]),BE_over_held_yield=float(y[1]/comparison)))
    assert abs(rows[0]["BE"]/rows[-1]["BE"]-1)<.006
    return dict(local_independent_controls=independent_controls(),data_sha256=hashlib.sha256(data.read_bytes()).hexdigest(),
                scope="Numerical 12–50 MeV, KLOE exclusive vacuum bin centers treated as constant cross sections. "
                      "Incoming thermal charged pions including BE; open unpolarized heavy decay; no initial heavy "
                      "and no unmeasured bins, scattering or full thermal hadronic response.",controls=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","archive","kloe"):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.archive,a.kloe),indent=2,allow_nan=False))
