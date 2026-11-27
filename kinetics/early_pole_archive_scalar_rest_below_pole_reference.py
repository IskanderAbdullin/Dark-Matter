"""Below-parent-pole Yukawa elastic diagnostic for a tagged early scalar at rest."""
import argparse
import gzip
import json

import numpy as np
from scipy.special import expit, roots_legendre

from early_pole_native_aperture_witness_reference import GRID, occup
from early_pole_decay_mode_clock_reference import G, KAPPA


def spin(s,t,M,m,ms):
    u=2*(m*m+ms*ms)-s-t
    a,b=1/(s-M*M),1/(u-M*M)
    A=M*(a+b)
    d=.5*(a*(s+m*m-ms*ms)+b*(u+m*m-ms*ms))
    C=s*a*a+u*b*b+2*(m*m-ms*ms)*a*b
    return 2*A*A*(2*m*m-t/2)+8*A*m*d+4*d*d+t*C


def fold(a,order,fraction):
    M,m=.05,.003
    ms=a['checkpoints'][-1]['scalar_pole_over_T']
    qe=np.asarray(a['scalar_momentum_edges'])
    bose=np.maximum(0.,np.asarray(a['endpoint_scalar_occupations']))
    pe=np.asarray(a['endpoint_fermion_momentum_edges'])
    ferm=np.asarray(a['endpoint_fermion_occupations'])
    Eedges=np.hypot(pe,m)
    lo=np.hypot(.45,m)
    hi=min(np.hypot(5.,m),(fraction**2*M*M-m*m-ms*ms)/(2*ms))
    outer=np.unique(np.r_[lo,hi,Eedges[(lo<Eedges)&(Eedges<hi)]])
    w,x=roots_legendre(order)
    h=np.diff(outer)[:,None]/2
    X=((outer[1:]+outer[:-1])[:,None]/2+h*w).ravel()
    W=(h*x).ravel()
    u=qe*qe/(np.hypot(ms,qe)+ms)
    bands=np.zeros((2,4))
    defect=0.
    for E,dE in zip(X,W,strict=True):
        p=np.sqrt((E-m)*(E+m))
        s=m*m+ms*ms+2*ms*E
        d=ms*p*p/s
        smax=min(2*d,E-lo)
        if smax<=0:
            continue
        bnd=np.r_[u,E-Eedges]
        bnd=np.unique(np.r_[0.,smax,bnd[(0<bnd)&(bnd<smax)]])
        h=np.diff(bnd)[:,None]/2
        k=((bnd[1:]+bnd[:-1])[:,None]/2+h*w).ravel()
        WK=(h*x).ravel()
        q=np.sqrt(k*(k+2*ms))
        t=-2*(ms*p/np.sqrt(s))**2*k/d
        f0=occup(np.array(p),pe,ferm)
        f1=occup(np.sqrt((E-k-m)*(E-k+m)),pe,ferm)
        g=occup(q,qe,bose)
        matrix=spin(s,t,M,m,ms)
        assert np.min(matrix)>=-1e-12
        rate=4*p*p*dE/(2*np.pi**2)*matrix*WK/(32*np.pi*s*d)
        em=f1*(1-f0)*g*(1+1/bose[0])
        loss=f0*(1-f1)*(1+g)
        vac=f0*(1-f1)
        net=g*(f1-f0)+g*f1*(1-f0)/bose[0]-vac
        FB,FA=expit(-E),expit(-(E-k))
        G0,G1=1/np.expm1(.003),1/np.expm1(.003+k)
        therm_emit=FA*(1-FB)*G1*(1+1/G0)
        therm_loss=FB*(1-FA)*(1+G1)
        defect=max(defect,float(np.max(abs(therm_emit/therm_loss-1))))
        for i,mask in enumerate((q<1e-5,q>=1e-5)):
            bands[i]+=np.array([rate[mask]@v[mask] for v in (em,loss,net,vac)])
    assert defect<1e-13
    # Gamma0/H = KAPPA M^2 and Gamma0/T = y^2 M/(16 pi),
    # with the same massless vacuum-width convention as the direct solver.
    H=np.sqrt(np.pi**2*106.75/90)*(86./a['x_end'])/2.435e18
    match=H*(16*np.pi*KAPPA*M)**2
    Ep=(M*M+ms*ms-m*m)/(2*ms)
    Eb=Ep-ms
    P=occup(np.sqrt((Ep-M)*(Ep+M)),np.asarray(a['endpoint_parent_momentum_edges']),
            np.asarray(a['endpoint_parent_occupations']))
    B=occup(np.sqrt((Eb-m)*(Eb+m)),pe,ferm)
    D=np.sqrt((M*M-(ms+m)**2)*(M*M-(ms-m)**2))/2
    direct=G*KAPPA*M*((M+m)**2-ms*ms)*2*D/ms**3*(P-B+P*(1-B)/bose[0])
    return dict(max_sqrt_s_over_parent_mass=fraction,terminal_incoming_energy=hi,
                cold_and_exterior_columns_emission_loss_net_vacuum_per_rest_occupation_without_y4=bands.tolist(),
                normalized_below_pole_net_by_band_per_reference_H=(match*bands[:,2]).tolist(),
                normalized_direct_scalar_rest_net_per_reference_H=float(direct),
                used_f_rest=float(bose[0]),thermal_balance_defect=defect)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=('fine','third'))
    p.add_argument('--order',type=int,default=8)
    a=p.parse_args()
    with gzip.open(GRID/f'early-pole-spectrum-{a.grid}-m06-2026-11-26.json.gz','rt') as f:
        data=json.load(f)
    print(json.dumps(dict(grid=a.grid,order=a.order,tests=[fold(data,a.order,r) for r in (.70,.82,.90)]),
                     indent=2,allow_nan=False))
