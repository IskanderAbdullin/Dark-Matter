"""Noncommuting slow-chemistry and long-wavelength limits in a static control."""
import argparse
import json
from pathlib import Path

import numpy as np

from charged_lepton_portal_local_longitudinal_reference import model
from positive_fd_portal_quadrature_reference import PositivePortal


def run(receipt):
    portal=PositivePortal((13,36,26),tail=65.)
    out=[]
    for row in receipt['rows'][:2]:
        T=row['T']
        xi=10.
        get=lambda k,e:model(T,xi,portal,row,k,elastic_only=e)
        L,H,*_=get(0.,False)
        L0,H0,*_=get(0.,True)
        np.testing.assert_array_equal(H,H0)
        h,U=np.linalg.eigh(H)
        W=(U*np.sqrt(h))@U.T
        J=(U/np.sqrt(h))@U.T
        transform=lambda X:W@X@J
        B=-transform(get(1.,True)[0]-L0).imag
        C=-transform(L).real
        Cel=-transform(L0).real
        Cpair=C*(1+row['timelike_heat']/row['spacelike_heat'])-Cel
        for D in (Cel,Cpair):
            if np.max(abs(D-D.T))>3e-13 or min(np.linalg.eigvalsh((D+D.T)/2))<-3e-13:
                raise ValueError("Pair/elastic split lost positivity.")
        t0,U0=np.linalg.eigh(Cel)
        t1,U1=np.linalg.eigh(Cel+Cpair)
        Ne,N=U0[:,:3],U1[:,:2]
        u,S=np.linalg.eigh(Ne.T@Cpair@Ne)
        chem=Ne@S[:,-1]
        rate=float(u[-1])
        speed,VS=np.linalg.eigh(N.T@B@N)
        acoustic=N@VS[:,-1]
        acoustic_speed=float(speed[-1])
        frozen_speed=float(np.linalg.eigvalsh(Ne.T@B@Ne)[-1])
        singular=float((chem@B@acoustic)**2/rate)
        if not rate>0 or not singular>0 or not frozen_speed>acoustic_speed:
            raise ValueError("Inconsistent singular coefficient or propagation speeds.")
        sweep=[]
        for eps in (.001,.0001,.00001,.000001):
            Ct=Cel+eps*Cpair
            R,V=np.linalg.eigh(Ct)
            slow=R[2]
            F=V[:,2:]
            d=float(np.sum((F.T@B@acoustic)**2/R[2:]))
            path=[]
            for over in (.01,.05,.2,1.,4.,20.,100.):
                k=slow*over
                ev,v=np.linalg.eig(-Ct-1j*k*B)
                j=int(np.argmin(ev.imag))
                a=v[:,j]
                decay=float(R[2:]@abs(F.T@a)**2/np.vdot(a,a).real)
                path.append(dict(k_over_slow_chemical_rate=over,wave_number=k,phase_speed=float(-ev[j].imag/k),
                                 damping_over_k_squared=decay/k**2))
            sweep.append(dict(timelike_weight=eps,chemical_rate=slow,pair_scaled_secondary_damping=eps*d,frequency_path=path))
        endl=sweep[-1]
        if abs(endl['chemical_rate']/(endl['timelike_weight']*rate)-1)>3e-5 or abs(endl['pair_scaled_secondary_damping']/singular-1)>3e-5:
            raise ValueError("Frozen-to-equilibrated singular coefficient was not reproduced.")
        path=endl['frequency_path']
        gap=frozen_speed-acoustic_speed
        if abs(path[0]['phase_speed']-acoustic_speed)>0.03*gap or abs(path[-1]['phase_speed']-frozen_speed)>0.03*gap:
            raise ValueError(f"Acoustic crossover missing: {T}, {path[0]}, {path[-1]}, {acoustic_speed}, {frozen_speed}.")
        out.append(dict(T=T,enthalpy_ratio=xi,frozen_sound_speed=frozen_speed,full_equilibrium_sound_speed=acoustic_speed,
                        projected_chemical_rate_per_timelike_weight=rate,
                        projected_singular_attenuation_weight=singular,sweep=sweep))
    return dict(scope="Static local limit with ideal thermal species and scaled timelike kernel, not the physical decoupling transition.",
                rows=out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('receipt')
    a=p.parse_args()
    print(json.dumps(run(json.loads(Path(a.receipt).read_text())),indent=2,allow_nan=False))
