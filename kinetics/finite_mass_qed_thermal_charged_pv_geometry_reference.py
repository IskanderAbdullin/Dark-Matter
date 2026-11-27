"""Geometry and angular principal value of the positive bare charged cuts."""
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad

from early_pole_decay_mode_clock_reference import intervals,Q
from early_pole_hard_yukawa_photoproduction_reference import dot
from finite_mass_qed_thermal_photon_vertex_cut_reference import trace


def chart(M,m,s,E,cut):
    Es=np.hypot(s,Q)
    P=np.sqrt((E-M)*(E+M))
    ca=(E*Es-(M*M+s*s-m*m)/2)/(P*Q)
    z=P*ca
    perp=P*np.sqrt((1-ca)*(1+ca))
    if cut=='parent':
        mass,other,sig,e,up=M,m,1.,E,z
    else:
        mass,other,sig,e,up=m,M,-1.,E-Es,z-Q
    p=np.sqrt((e-mass)*(e+mass))
    c=up/p
    d=(M*M+m*m-s*s)/2
    U=np.array([e,perp,0.,up])
    W=sig*np.array([Es,0.,0.,Q])
    V=U-W
    C=(M+m)**2/(d+M*m)
    D=4-2*(M+m)*mass/(d+M*m)

    def evaluate(rog):
        a=rog if cut=='parent' else rog-Es
        r=np.sqrt((a-mass)*(a+mass))
        den=2*sig*Q*r
        dc=mass*mass+s*s-other*other-2*sig*Es*a
        root=-dc/den
        eta=np.arcsinh(p/mass)
        tau=np.arcsinh(r/mass)
        base=2*mass*mass*np.sinh((eta-tau)/2)**2
        th=np.arccos(c)
        def photon(u):
            x=np.arccos(u)
            diff=base+2*p*r*np.sin((x-th)/2)**2
            sum_=base+2*p*r*np.sin((x+th)/2)**2
            return -1/(2*np.sqrt(diff*sum_))

        def numerator(u):return C+4*d*photon(u)

        rad=np.sqrt((2*base)/(p*r))
        guides=[v for v in ([c]+[c+k*rad for k in (-64,-16,-4,-1,1,4,16,64)]) if -1<v<1]
        if -1<root<1:
            # Piecewise quadrature prevents the sharp near-collinear photon
            # factor from being confused with the simple fermion pole.
            side=sorted({-1.,*guides,1.})
            pv=0.
            for l,h in zip(side[:-1],side[1:],strict=True):
                if l<root<h:
                    pv+=quad(numerator,l,h,weight='cauchy',wvar=root,
                             epsabs=4e-11,epsrel=5e-11,limit=250)[0]
                else:
                    pv+=quad(lambda u:numerator(u)/(u-root),l,h,
                             epsabs=4e-11,epsrel=5e-11,limit=250)[0]
            n0=numerator(root)
            def diff(u):
                if abs(u-root)<1e-12:
                    delta=1e-8
                    return (numerator(root+delta)-numerator(root-delta))/(2*delta)
                return (numerator(u)-n0)/(u-root)
            alt=quad(diff,-1.,1.,points=sorted(set([*guides,root])),
                     epsabs=4e-11,epsrel=5e-11,limit=250)[0]+n0*np.log((1-root)/(1+root))
            check=abs(pv-alt)/max(abs(alt),1.)
            if check>2e-8:
                raise ValueError(f'Independent principal-value methods disagreed: {check:g}')
            inverse_at_root=float(photon(root))
        else:
            pv=quad(lambda u:numerator(u)/(u-root),-1.,1.,points=sorted(set(guides)),
                    epsabs=4e-11,epsrel=5e-11,limit=250)[0]
            inverse_at_root=None
            check=None
        reg=quad(photon,-1.,1.,points=sorted(set(guides)),
                 epsabs=4e-11,epsrel=5e-11,limit=250)[0]*D/2
        # Verify the pole split from an independent gamma trace at off-axis
        # directions, avoiding the intentionally excluded simultaneous pole.
        t=.31 if root<.06 else -.38
        psi=.74
        trf=r*np.sqrt((1-t)*(1+t))
        k=U-np.array([a,trf*np.cos(psi),trf*np.sin(psi),r*t])
        K2=dot(k,k)
        qk=dot(W,k)
        n0=trace(U,V,k,mass,other)
        ideal=K2*C+4*d+D*qk
        spin=abs(n0-ideal)/max(abs(ideal),1e-6)
        phase=abs((a-Es*sig)**2-((r*t-Q*sig)**2+trf**2)-other**2-2*qk)
        if spin>1e-10 or phase>3e-15:
            raise ValueError('On-shell scalar spin/denominator identities did not close.')
        return dict(independent_spin_trace_defect=float(spin),uncut_denominator_identity_error=float(phase),
                    internal_energy=float(a),angle_root=float(root),support=bool(abs(root)<1),
                    principal_value_coefficient=float(pv/den+reg),
                    pv_algorithm_defect=None if check is None else float(check),
                    azimuth_average_inverse_photon_at_root=inverse_at_root,
                    radial_delta_Jacobian_times_dE=float(r/(2*abs(den))))
    return evaluate


def run():
    rows=[]
    for grid in ('fine','third'):
        path=Path(__file__).parent/f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz'
        with gzip.open(path) as f:d=json.load(f)
        M=d['x_end'];m=M*d['mass_ratio'];s=d['checkpoints'][-1]['scalar_pole_over_T']
        Es,lo,hi,*_=intervals(M,s,m,Q)
        e=(lo+hi)/2
        p,b=[chart(M,m,s,e,k) for k in ('parent','daughter')]
        data=[]
        for t in (-.12,.08,.24,.76,.92,1.12):
            x=lo+t*(hi-lo)
            a,c=p(x),b(x)
            if a['support']!= (0<t<1) or c['support']!=a['support']:
                raise ValueError('Positive cut roots did not map to the same aperture.')
            ar,cr=[z['azimuth_average_inverse_photon_at_root'] for z in (a,c)]
            witness=abs(ar/cr-1) if ar is not None else None
            if witness is not None and witness>3e-8:
                raise ValueError(f'On the simultaneous cut the uncut photon did not agree: {witness:g}.')
            ds=max(abs(a['radial_delta_Jacobian_times_dE']-1/(4*Q)),
                   abs(c['radial_delta_Jacobian_times_dE']-1/(4*Q)))
            if ds>3e-11:
                raise ValueError('The radial energy Jacobian did not agree.')
            data.append(dict(parent_window_fraction=t,parent=a,daughter=c,paired_inverse_photon_relative_defect=witness))
        rows.append(dict(grid=grid,external_q=Q,external_parent_energy=e,internal_parent_Born_edges=[lo,hi],rows=data))
    return dict(scope='Bare Feynman-gauge positive-energy charged cut angular PV geometry, separated from soft coincident radial point and without occupation or real-time cut weights.',rows=rows)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
