"""Couple Pauli χ1ℓ,χ2ℓ 2→2 to real-φ Yukawa χ1↔χ2+φ.

All FOUR node populations evolve on shared per-species physical E grids;
there is ONE shared fermionic lepton bath and one dynamical real Bose φ.
Two collision source operators independently use identical χ grids,
capacities and current distributions. Stoichiometric nullspace:
Nχ1+Nχ2, Nχ1+Nφ, Nℓ, Etot. Entropy derivative ≥0 by
quantum F/R at every on-shell event. Adaptive TWO half-Euler steps are
FIRST ORDER and do not certify cosmological rates, BEC or 86-eV relic.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.special import expit,logit
from kinetics.chi1_chi2_shared_lepton_pauli import (
    system as elastic_system, rhs as elastic_rhs, ledger as elastic_ledger
)
from kinetics.yukawa_bose_pauli_bounded_time import (
    Config as YukawaConfig, density_cells, rhs as yukawa_rhs,
    fermion_entropy,boson_entropy
)


def configuration(*,bins=24,masses=(.1,.02),mphi=.001,mlepton=0.,
                  elastic_cut=6.,yukawa_pmax=12.,y=.001,
                  elastic_G=.01,elastic_mv=.1,
                  elastic_nr=4,elastic_nmu=4,elastic_nc=8,elastic_nphi=8,
                  yukawa_nr=12,yukawa_nangle=12,
                  elastic_radial_scheme='uniform_p',elastic_segments=1,
                  elastic_integration_cut=None):
    if (elastic_integration_cut is not None and
        not (math.isfinite(elastic_integration_cut) and
             0<elastic_integration_cut<=elastic_cut)):
        raise ValueError('elastic integration cutoff outside common mode grid')
    elastic=elastic_system(masses=masses,mlepton=mlepton,bins=bins,cut=elastic_cut)
    yuk=YukawaConfig(m1=masses[0],m2=masses[1],mphi=mphi,y=y,
                     pmax=yukawa_pmax,bins=bins,nr=yukawa_nr,
                     nangle=yukawa_nangle)
    upper=elastic['lg'][-1]
    if math.hypot(masses[0],yukawa_pmax)>upper:
        raise ValueError('shared physical grids must contain moving parent and daughters')
    # EXACT same chi nodes and physical capacities in BOTH sources.
    yuk.e1=elastic['chi_grid'][0]
    yuk.e2=elastic['chi_grid'][1]
    yuk.w1=elastic['chi_w'][0]
    yuk.w2=elastic['chi_w'][1]
    yuk.ephi=np.linspace(mphi,upper,bins+1)
    yuk.wphi=density_cells(yuk.ephi,mphi,yuk.gphi)
    return {'elastic':elastic,'yukawa':yuk,'G':elastic_G,'mv':elastic_mv,
            'elastic_nr':elastic_nr,'elastic_nmu':elastic_nmu,
            'elastic_nc':elastic_nc,'elastic_nphi':elastic_nphi,
            'shared_energy_upper':upper,
            'elastic_radial_scheme':elastic_radial_scheme,
            'elastic_segments':elastic_segments,
            'elastic_integration_cut':elastic_integration_cut}


def initial(cfg,*,T1=1.1,T2=1.,Tlep=1.,Tphi=1.,mu_phi=0.):
    eg=cfg['elastic']; y=cfg['yukawa']
    if not(T1>0 and T2>0 and Tlep>0 and Tphi>0 and mu_phi<y.mphi):
        raise ValueError('physical initial Fermi/Bose thermal parameters')
    return (expit(-eg['chi_grid'][0]/T1),
            expit(-eg['chi_grid'][1]/T2),
            expit(-eg['lg']/Tlep),
            1/np.expm1((y.ephi-mu_phi)/Tphi))


def ledger(cfg,state):
    e=cfg['elastic'];y=cfg['yukawa']
    f1,f2,fl,b=state
    n1=float(np.dot(y.w1,f1));n2=float(np.dot(y.w2,f2))
    nl=float(np.dot(e['lw'],fl));nb=float(np.dot(y.wphi,b))
    energy=float(np.dot(y.w1*y.e1,f1)+np.dot(y.w2*y.e2,f2)+
                 np.dot(e['lw']*e['lg'],fl)+np.dot(y.wphi*y.ephi,b))
    entropy=sum((fermion_entropy(f1,y.w1),fermion_entropy(f2,y.w2),
                 fermion_entropy(fl,e['lw']),boson_entropy(b,y.wphi)))
    return {'N_dark_fermions':n1+n2,'N1_plus_phi':n1+nb,
            'Nlepton':nl,'energy':energy,'entropy':entropy,
            'N1':n1,'N2':n2,'Nphi':nb}


def rhs(cfg,state):
    f1,f2,fl,b=state
    e=cfg['elastic'];y=cfg['yukawa']
    de,dl,ele=elastic_rhs(e,[f1,f2],fl,nr=cfg['elastic_nr'],
            nmu=cfg['elastic_nmu'],nc=cfg['elastic_nc'],
            nphi=cfg['elastic_nphi'],G=cfg['G'],mv=cfg['mv'],
            radial_scheme=cfg['elastic_radial_scheme'],
            segments=cfg['elastic_segments'],
            incoming_cut=cfg['elastic_integration_cut'])
    d1,d2,db,yuk=yukawa_rhs(y,f1,f2,b)
    ds=(de[0]+d1,de[1]+d2,dl,db)
    sgrad=(-float(np.dot(y.w1*ds[0],logit(f1)))
           -float(np.dot(y.w2*ds[1],logit(f2)))
           -float(np.dot(e['lw']*ds[2],logit(fl)))
           +float(np.dot(y.wphi*ds[3],np.log1p(1/b))))
    sed=ele['entropy_direct']+yuk['entropy_event']
    charges=[float(np.dot(y.w1,ds[0])+np.dot(y.w2,ds[1])),
             float(np.dot(y.w1,ds[0])+np.dot(y.wphi,ds[3])),
             float(np.dot(e['lw'],ds[2]))]
    dE=float(np.dot(y.w1*y.e1,ds[0])+np.dot(y.w2*y.e2,ds[1])+
            np.dot(e['lw']*e['lg'],ds[2])+np.dot(y.wphi*y.ephi,ds[3]))
    return ds,{'entropy_event':sed,'entropy_node':sgrad,
               'entropy_identity_abs_error':abs(sed-sgrad),
               'elastic_entropy':ele['entropy_direct'],
               'yukawa_entropy':yuk['entropy_event'],
               'energy_residual':dE,'charge_residuals':charges,
               'elastic_Qchi':ele['Qchi'],'elastic_Qlep':ele['Qlep'],
               'Yukawa_Rphi':yuk['Rphi'],
               'Yukawa_Qtotal':yuk['heat_residual']}


def guarded_step(cfg,state,dt):
    derivative,_=rhs(cfg,state)
    after=tuple(f+dt*d for f,d in zip(state,derivative))
    if any(np.any(x<=0) for x in after) or any(np.any(x>=1) for x in after[:3]):
        return None
    if ledger(cfg,after)['entropy']+1e-12<ledger(cfg,state)['entropy']:
        return None
    return after


def evolve(cfg,state0,*,end_time=2e7,initial_dt=1e7,atol=2e-5,max_attempts=600):
    if not(end_time>0 and initial_dt>0 and atol>0):raise ValueError
    x=tuple(np.array(f,copy=True) for f in state0)
    first=ledger(cfg,x)
    t=0.;dt=initial_dt;accept=reject=0;minS=math.inf;maxerr=0.
    for _ in range(max_attempts):
        if t>=end_time*(1-1e-14):break
        dt=min(dt,end_time-t)
        whole=guarded_step(cfg,x,dt)
        h=guarded_step(cfg,x,dt/2)
        hh=guarded_step(cfg,h,dt/2) if h is not None else None
        if whole is not None and hh is not None:
            err=max(*(float(np.max(abs(a-b))) for a,b in zip(whole[:3],hh[:3])),
                    float(np.max(abs(np.log1p(whole[3])-np.log1p(hh[3])))))
            if err<=atol:
                before=ledger(cfg,x)['entropy']
                x=tuple(f.copy() for f in hh)
                minS=min(minS,ledger(cfg,x)['entropy']-before)
                maxerr=max(maxerr,err)
                accept+=1;t+=dt
                dt*=min(2.,max(1.05,.9*math.sqrt(atol/max(err,1e-30))))
                continue
        reject+=1;dt*=.5
        if dt<1e-16:raise ArithmeticError('unable to preserve Pauli/Bose/entropy')
    else:raise RuntimeError('max adaptive attempts exhausted')
    last=ledger(cfg,x)
    # Explicitly detect actual energy exchange with ONE common lepton bath.
    egrid=cfg['elastic']['lg']; capacity=cfg['elastic']['lw']
    delta_lepton_E=float(np.dot(egrid*capacity,x[2]-state0[2]))
    return {'time':t,'accepted':accept,'rejected':reject,
            'max_local_spectral_error':maxerr,'min_step_entropy_gain':minS,
            'dNdark':last['N_dark_fermions']-first['N_dark_fermions'],
            'dN1phi':last['N1_plus_phi']-first['N1_plus_phi'],
            'dNlep':last['Nlepton']-first['Nlepton'],
            'dEnergy':last['energy']-first['energy'],
            'dElepton':delta_lepton_E,
            'dEntropy':last['entropy']-first['entropy'],
            'Nphi_created':last['Nphi']-first['Nphi'],
            'N2_created':last['N2']-first['N2'],
            'N1_destroyed':first['N1']-last['N1'],
            'Fermi_min':min(float(np.min(f)) for f in x[:3]),
            'Fermi_max':max(float(np.max(f)) for f in x[:3]),
            'Bose_min':float(np.min(x[3])),'Bose_max':float(np.max(x[3])),
            'not_physical_cosmological_gamma_or_86eV':True,
            'first_order_method':True}

if __name__=='__main__':
    import json
    c=configuration(bins=16)
    x=initial(c,T1=1.1,T2=1.,Tlep=1.)
    _,r=rhs(c,x)
    print(json.dumps({'instant':r,'time':evolve(c,x)},indent=2))
