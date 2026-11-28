"""Conditional equilibrated scalar-only capacity with conserving Hartree energy.

This cannot infer whether elastic collisions act quickly or preserve scalar
particle count in the actual expanding, driven system.
"""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import roots_legendre


def scalar_capacity(n,rho,m0,lam,order):
    zz,ww=roots_legendre(order)
    y=26*(zz+1)
    weights=26*ww*y*y/(2*np.pi**2)
    a=lam/2.
    assert n>0 and rho>0 and m0>0 and lam>0

    def thermal(th,m,offset=0.):
        p=th*y
        E=np.hypot(p,m)
        k_over_T=(th*y*y)/(E+m)
        w=th**3*weights/np.expm1(k_over_T+offset)
        return w.sum(), w@(1/E), w@E

    def regular(th,offset=0.):
        upper=m0+20*np.sqrt(lam)*th+1e-14
        gap=lambda m:m*m-m0*m0-a*thermal(th,m,offset)[1]
        m=brentq(gap,m0,upper,xtol=1e-15)
        nt,inv,E=thermal(th,m,offset)
        U=(m*m-m0*m0)**2/(2*lam)
        return nt,E-U,m

    def condensed(th):
        upper=m0+1.001*(a*n)**(1/3)
        def gap(m):
            nt,inv,E=thermal(th,m)
            return m*m-m0*m0-a*((n-nt)/m+inv)
        m=brentq(gap,m0,upper,xtol=1e-15)
        nt,inv,E=thermal(th,m)
        U=(m*m-m0*m0)**2/(2*lam)
        return (n-nt)*m+E-U,n-nt,m

    # Saturated regular distribution at the input total number.
    upper=(n*np.pi**2/1.2020569031595943)**(1/3)*1.1
    satT=brentq(lambda t:regular(t)[0]-n,upper*1e-7,upper,xtol=1e-14)
    nt,rt,ms=regular(satT)
    normal=[]
    # Verify fixed-N regular solutions for nonsaturated fugacities.
    for v in (.0001,.001,.01,.1,.3,1.,2.,4.,8.):
        th=brentq(lambda t:regular(t,v)[0]-n,upper*1e-7,upper*np.exp(v/3)*1.2,xtol=1e-13)
        normal.append(regular(th,v)[1]/rt)
    assert normal[0]>1 and min(np.diff(normal))>0.
    # Any monotonicity cited is only sampled on the represented branch.
    er=[regular(t)[1] for t in np.geomspace(.02*satT,satT,17)]
    assert min(np.diff(er))>0
    if rho>=rt:
        return dict(saturated_regular_temperature_over_T=satT,
                    saturated_regular_energy_over_T4=rt,ratio_to_threshold=rho/rt,
                    lowest_tested_nonsaturated_regular_energy_over_saturated=float(normal[0]),
                    condensed_branch_needed=False)

    l=1e-8*satT
    e0,n0,mzero=condensed(l)
    if rho<e0:
        raise ValueError("Input energy below the Hartree zero momentum floor")
    T=brentq(lambda t:condensed(t)[0]-rho,l,satT,xtol=1e-14)
    energy,n0,m=condensed(T)
    assert 0<n0<n and abs(energy/rho-1)<1e-9
    return dict(saturated_regular_temperature_over_T=satT,
                saturated_regular_energy_over_T4=rt,ratio_to_threshold=rho/rt,
                condensed_branch_needed=True,
                lowest_tested_nonsaturated_regular_energy_over_saturated=float(normal[0]),
                zero_temperature_energy_floor_over_T4=e0,
                conditional_zero_mode_number_fraction=n0/n,excitation_temperature_over_T=T,
                equilibrium_scalar_pole_over_T=m,
                equilibrium_energy_relative_residual=energy/rho-1)


def assess(spec):
    r=spec['checkpoints'][-1]
    n=float(np.sum(spec['endpoint_scalar_signed_shell_number_over_T3']))
    rho=r['integrated_scalar_energy_over_T4']
    m0=1e-9*r['x']
    cases=[scalar_capacity(n,rho,m0,7.65175e-8,k) for k in (96,192)]
    if cases[-1]['condensed_branch_needed']:
        assert abs(cases[0]['conditional_zero_mode_number_fraction']-
                   cases[-1]['conditional_zero_mode_number_fraction'])<3e-7
    q=np.array(spec['scalar_momentum_nodes'])
    atoms=np.maximum(np.array(spec['endpoint_scalar_signed_shell_number_over_T3']),0.)
    hot=np.array(spec['scalar_momentum_edges'])[:-1]>=.1
    return dict(input_x=r['x'],input_number_over_T3=n,input_scalar_energy_over_T4=rho,
                original_q_above_point1_number_fraction=float(atoms[hot].sum()/atoms.sum()),
                original_q_above_point1_share_of_second_momentum_moment=float(atoms[hot]@(q[hot]**2)/(atoms@(q*q))),
                coarse_quadrature=cases[0],refined_quadrature=cases[-1])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('spectra',nargs='*')
    p.add_argument('--ideal-control',action='store_true')
    args=p.parse_args()
    if args.ideal_control:
        n=.015
        th=(n*np.pi**2/1.2020569031595943)**(1/3)
        reference=np.pi**2/30*th**4
        check=scalar_capacity(n,2*reference,1e-10,1e-18,192)
        assert abs(check['saturated_regular_temperature_over_T']/th-1)<3e-8
        assert abs(check['saturated_regular_energy_over_T4']/reference-1)<1e-7
        assert not check['condensed_branch_needed']
        print(json.dumps(check,indent=2,allow_nan=False))
        raise SystemExit(0)
    if not args.spectra:
        p.error('provide at least one spectrum or --ideal-control')
    result={}
    for f in args.spectra:
        file=Path(f)
        raw=gzip.decompress(file.read_bytes()) if file.suffix=='.gz' else file.read_bytes()
        result[file.name]=assess(json.loads(raw))
    print(json.dumps(result,indent=2,allow_nan=False))
