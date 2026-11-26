"""Fixed physical 1<->2 Yukawa rate: independent Bose IR mode shape resolution.

Synthetic source only; mphi=86eV benchmark, no FRW/SM and no inference of BEC.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.special import expit
from kinetics.phi_bose_moment_galerkin import Mode,yukawa_moment_source
from kinetics.phi_bose_implicit_moment_step import physically_weighted_yukawa_events
from kinetics.moving_yukawa_bose_pauli_weak import weak_moments


def eta_phi_analytic(E,shape):
    E=np.asarray(E,dtype=float)
    if shape == 'thermal':
        return E
    if shape == 'infrared_curved':
        return E + .40*np.exp(-((E-.04)/.02)**2)
    if shape == 'broader':
        return E + .30*np.exp(-((E-.10)/.06)**2)
    raise ValueError(shape)


def rate(shape='infrared_curved',bins=64,nrad=24,nangle=48,pmax=2.):
    m1,m2,mphi=.1,.02,86e-9
    upper=math.hypot(m1,pmax)+.2
    modes=[Mode(m,np.linspace(m,upper,bins+1),degeneracy=4 if i<2 else 1,
                fermion=i<2) for i,m in enumerate((m1,m2,mphi))]
    eta=[modes[0].nodes/.9,modes[1].nodes,eta_phi_analytic(modes[2].nodes,shape)]
    events,width=physically_weighted_yukawa_events(m1=m1,m2=m2,mphi=mphi,
                       y=.001,pmax=pmax,nr=nrad,nangle=nangle,g1=4)
    x=yukawa_moment_source(modes,eta,events)
    return {'shape':shape,'bins':bins,'nr':nrad,'nangle':nangle,
            'Rphi':x['rate_phi'],'H':x['entropy_from_event'],
            'Q_residual':x['total_energy_residual'],
            'gamma_rest':width}


def reference(shape='infrared_curved',nrad=24,nangle=48,pmax=2.):
    return weak_moments(lambda E: expit(-np.asarray(E)/.9),
        lambda E: expit(-np.asarray(E)),
        lambda E: 1/np.expm1(eta_phi_analytic(E,shape)),
        m1=.1,m2=.02,mphi=86e-9,y=.001,pmax=pmax,nr=nrad,nangle=nangle,g1=4)


def audit():
    data={}
    for shape in ('thermal','broader','infrared_curved'):
        old=reference(shape)
        rows=[rate(shape,bins=n) for n in (16,32,64,128,256,512)]
        for r in rows:
            r['Q_abs_rel_error']=r['Rphi']/old['R_phi']-1
            r['H_abs_rel_error']=r['H']/old['entropy']-1
        data[shape]={'physical_parent_rate_reference':old,'rows':rows}
    data['parent_cutoff_ir_curved_nrad72_nangle144']=cutoff_audit()
    data['parent_resolution_control']={k:[{'nr':n,'na':a,**reference(k,nrad=n,nangle=a)} for n,a in ((24,48),(48,96),(72,144))] for k in ('broader','infrared_curved')}
    data['scope']='synthetic physical vacuum Yukawa p/c weak integral, not q/FRW/SM convergence'
    data['class_ready']=False
    return data


def cutoff_audit(shape='infrared_curved', cutoffs=(2.,4.,6.,8.,12.,16.),
                 nrad=72,nangle=144):
    """Physical parent p cutoff, independent of any φ energy-hat mesh."""
    return [{'parent_pmax_GeV':float(p),
             'Rphi':reference(shape,nrad=nrad,nangle=nangle,pmax=p)['R_phi'],
             'entropy':reference(shape,nrad=nrad,nangle=nangle,pmax=p)['entropy']}
            for p in cutoffs]


if __name__=='__main__':
    import json
    print(json.dumps(audit(),indent=2))
