"""IR-focused positive energy hats for physical Bose Yukawa moments.

Synthetic target: exact thermal vs curved IR Bose log-fugacity weak sources.
No physical early-universe relic, zero-mode condensate or full χ-SM plasma.
"""
from __future__ import annotations
import math
import numpy as np
from kinetics.phi_bose_moment_galerkin import Mode,yukawa_moment_source,physical_moments
from kinetics.phi_bose_implicit_moment_step import physically_weighted_yukawa_events
from kinetics.phi_bose_ir_shape_gate import eta_phi_analytic,reference


def nodes(mass,upper,n,method='uniform',*,split=.25):
    if not(mass>=0 and upper>mass and isinstance(n,int) and n>=8):
        raise ValueError('positive mass/support and nodes required')
    x=np.linspace(0.,1.,n+1)
    if method=='uniform':
        E=mass+(upper-mass)*x
    elif method=='power2':
        E=mass+(upper-mass)*x*x
    elif method=='power3':
        E=mass+(upper-mass)*x**3
    elif method=='two_zone':
        if not(mass<split<upper):raise ValueError('split in domain')
        low=n//2
        E=np.r_[np.linspace(mass,split,low+1)[:-1],
                np.linspace(split,upper,n-low+1)]
    else:raise ValueError('unknown energy-hat scheme')
    if abs(E[0]-mass)>1e-12 or abs(E[-1]-upper)>1e-12 or np.any(np.diff(E)<=0):
        raise ArithmeticError('nonmonotone energy grid')
    return E


def moment_weak_rate(shape='infrared_curved',bins=32,scheme='uniform',
                     *,pmax=2.,nr=24,nangle=48):
    m1,m2,mphi=.1,.02,86e-9
    upper=math.hypot(pmax,m1)+.2
    modes=[Mode(m,nodes(m,upper,bins,'uniform'),fermion=True,degeneracy=4)
           for m in (m1,m2)]
    modes.append(Mode(mphi,nodes(mphi,upper,bins,scheme),degeneracy=1))
    eta=[modes[0].nodes/.9,modes[1].nodes,eta_phi_analytic(modes[2].nodes,shape)]
    events,_=physically_weighted_yukawa_events(m1=m1,m2=m2,mphi=mphi,
                        y=.001,pmax=pmax,nr=nr,nangle=nangle,g1=4)
    out=yukawa_moment_source(modes,eta,events)
    return {'shape':shape,'scheme':scheme,'bins':bins,
            'Rphi':out['rate_phi'],'H':out['entropy_from_event'],
            'Q_residual':out['total_energy_residual'],
            'Nphi_thermal_if_thermal':(physical_moments(modes[2],eta[2],order=36)['N']
                if shape=='thermal' else None),
            'min_phi_spacing':float(np.min(np.diff(modes[2].nodes)))}


def audit(*,pmax=2.,nrad=24,nangle=48,binlist=(16,32,64,128,256)):
    rows=[]
    for shape in ('thermal','broader','infrared_curved'):
        ref=reference(shape,nrad=nrad,nangle=nangle,pmax=pmax)
        for n in binlist:
            for scheme in ('uniform','power2','power3','two_zone'):
                r=moment_weak_rate(shape,n,scheme,pmax=pmax,nr=nrad,nangle=nangle)
                r['rate_relative_error']=r['Rphi']/ref['R_phi']-1
                r['H_relative_error']=r['H']/ref['entropy']-1
                rows.append(r)
    trials=[implicit_trial(n,scheme,dt=1e3) for n,scheme in ((16,'uniform'),(16,'power2'),(32,'power3'),(64,'uniform'),(64,'power3'))]
    return {'schema_version':1,'implicit_near_critical_trials':trials,'pmax':pmax,'nr':nrad,'nangle':nangle,
            'models':'synthetic fixed Pauli/Bose inverse fugacity, no physical FRW or BEC',
            'rows':rows,'class_ready':False}


def implicit_trial(bins=16,scheme='uniform',*,pmax=2.,dt=1e7):
    """Actual physical Γrest implicit source; failure is a retained negative.

    Nonuniform phi hats have a positive-p radial measure but their
    moment inversion near ηφ(0)=mφ/T can become infeasible under a
    sparse event mesh. A failure is NOT proof of a real BEC.
    """
    from kinetics.phi_bose_implicit_moment_step import implicit_step
    m1,m2,mphi=.1,.02,86e-9
    upper=math.hypot(m1,pmax)+.2
    modes=[Mode(m,nodes(m,upper,bins),degeneracy=4,fermion=True)
           for m in (m1,m2)]
    modes.append(Mode(mphi,nodes(mphi,upper,bins,scheme)))
    eta=[modes[0].nodes/.9,modes[1].nodes,modes[2].nodes]
    events,_=physically_weighted_yukawa_events(pmax=pmax,nr=12,nangle=12)
    try:
        out=implicit_step(modes,eta,events,dt,quad_order=24)
        return {'scheme':scheme,'bins':bins,'dt':dt,'accepted':True,
                'Nphi_change':out['increments']['N_phi'],
                'H_gain':out['increments']['S_total'],
                'min_phi_eta':out['min_phi_eta'],
                'N_E_defects':[out['increments'][name] for name in
                   ('N1_plus_N2','N1_plus_Nphi','U_total')]}
    except ArithmeticError as exc:
        return {'scheme':scheme,'bins':bins,'dt':dt,'accepted':False,
                'failure':str(exc),
                'not_a_real_condensate':True}


if __name__=='__main__':
    import json
    print(json.dumps(audit(),indent=2))
