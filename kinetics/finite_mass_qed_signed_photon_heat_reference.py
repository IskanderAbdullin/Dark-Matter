"""Real photon heat cancels its leading Bose piece in a different signed combination."""
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from finite_mass_qed_signed_support_bridge_reference import BridgeSource
from finite_mass_qed_integrated_signed_bridge_reference import seed


def spectrum(q,T,Fp,Rp,Fm,Rm):
    hs=1/(2*np.tanh(q/(2*T)))
    return hs*((Fp-Rp)-(Fm-Rm))+(Fp+Rp+Fm+Rm)/2


def artificial():
    x=np.array([2e-7,2e-5,.002],dtype=np.longdouble)
    T=np.longdouble('.67')
    p0=np.longdouble('1.2')
    a=np.longdouble('.35')
    Pp,Pm=p0+.12*x,p0-.12*x
    s=np.tanh(x/(2*T))
    thermal=spectrum(x,T,(1-s)*Pp/2,(1+s)*Pp/2,(1+s)*Pm/2,(1-s)*Pm/2)
    Np,Nm=.11+a*x,.11-a*x
    real=spectrum(x,T,(Pp+Np)/2,(Pp-Np)/2,(Pm+Nm)/2,(Pm-Nm)/2)
    expect=float(2*T*a+p0)
    if max(abs(thermal))>7e-12 or abs(real[0]-expect)>1e-7:
        raise ArithmeticError('Manufactured soft heat limit or equilibrium failed')
    return dict(thermal_absolute_defect=float(max(abs(thermal))),linear_polynomial_limit=expect,
                photon_heat_fixture=real.astype(float).tolist())


def control(S,fs,lo,hi,n,o):
    A0=seed(S)*[1+np.longdouble(1)/fs,1.]
    N0=A0[0]-A0[1]
    B=N0/S.B
    x,w=roots_legendre(o)
    a,b=np.log(np.longdouble(lo)),np.log(np.longdouble(hi))
    powers=np.exp((a+b)/2+(b-a)*x/2)
    a1=[]
    defect=0.
    for q in powers:
        Fp,Rp=S.signed(q,n,9,28)[0,:2]*[1+np.longdouble(1)/fs,1.]
        Fm,Rm=S.signed(-q,n,9,28)[0,:2]*[1+np.longdouble(1)/fs,1.]
        nB=1/np.expm1(q)
        Q=spectrum(q,1.,Fp,Rp,Fm,Rm)
        P=(1+nB)*Fp-nB*Rp
        A=nB*Fm-(1+nB)*Rm
        defect=max(defect,abs(P-A-Q))
        # Signed photon heat needs no virtual scalar counterterm: its plus
        # and minus real events carry +k and -k respectively. Scalar power
        # needs the common leading subtraction before an IR continuation.
        C=(P+A)-(1+2*nB)*N0
        np0=Fp-Rp;nm0=Fm-Rm
        Pp=Fp+Rp;Pm=Fm+Rm
        stable=(.5+nB)*(np0+nm0-2*N0)+(Pp-Pm)/2
        if abs(C-stable)>3e-15:
            raise ArithmeticError('Independent scalar local subtraction failed')
        a1.append([q*Q,stable, -(q*Q+S.Es*stable)])
    Gam,Signed,Charged=(b-a)/2*(w@np.asarray(a1))
    c=.005**2/(4*np.pi**2)
    return dict(edges=[lo,hi],photon_heat_in_unscaled_proxy_units=float(Gam),
                scalar_net_in_unscaled_proxy_units=float(Signed),
                charged_ledger_implied=float(Charged),
                photon_heat_over_Born_scalar_energy_g_point005=float(c*Gam/(B*S.Es)),
                scalar_number_relative_to_Born_g_point005=float(c*Signed/B),
                direct_cancellation_replay_abs_defect=float(defect))


def run(grids):
    tracks=[]
    for g in grids:
        root=Path(__file__).parent/'results'
        with gzip.open(root/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as f:d=json.load(f)
        S=BridgeSource(d)
        fs=d['endpoint_scalar_occupations'][np.searchsorted(d['scalar_momentum_edges'],float(S.q))-1]
        bounds=[1e-5,2e-5,4e-5,8e-5,1.6e-4,2e-4]
        a,b=([control(S,fs,l,h,n,o) for l,h in zip(bounds[:-1],bounds[1:],strict=True)]
             for n,o in ((28,4),(46,8)))
        drift=max(abs(x['photon_heat_over_Born_scalar_energy_g_point005']-y['photon_heat_over_Born_scalar_energy_g_point005'])
                  for x,y in zip(a,b,strict=True))
        if drift>1e-9:raise ArithmeticError(('Photon heat integration was not refined',g,drift))
        tracks.append(dict(grid=g,panels=b,max_two_mesh_abs_change_in_photon_ratio=drift,
                           total_photon_ratio=sum(x['photon_heat_over_Born_scalar_energy_g_point005'] for x in b),
                           total_scalar_ratio=sum(x['scalar_number_relative_to_Born_g_point005'] for x in b)))
    return dict(scope='Fixed cold scalar mode and frozen archived charged populations. Signed real photon heat in a finite band; no in-medium photon dressing.',manufactured=artificial(),result=tracks)


if __name__=='__main__':
    import sys
    print(json.dumps(run(tuple(sys.argv[1:]) or ('fine','third')),indent=2,allow_nan=False))
