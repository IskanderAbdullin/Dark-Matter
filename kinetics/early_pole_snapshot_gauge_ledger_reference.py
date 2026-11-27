"""Frozen reciprocal gauge-energy scale audit on archived snapshots."""
import json
import numpy as np
from scipy.integrate import quad
from scipy.special import expit

from early_pole_frozen_dynamic_vlasov_reference import load
from early_pole_snapshot_hard_loop_screening_reference import fd


def cell_J(edges,f,m,cut):
    p=np.asarray(edges)
    ans=0.
    for l,h,o in zip(p[:-1],p[1:],f):
        if h>max(l,cut):
            ans+=o*quad(lambda q:q*q/np.hypot(q,m),max(l,cut),h)[0]
    return ans


def photon(T,alpha,cut):
    return quad(lambda p:p/np.expm1(p/T-alpha),cut,80.)[0]


def moment(m,N,T,cut):
    y=fd(m,N,T,0.)
    a=y['alpha']
    return quad(lambda p:p*p/np.hypot(p,m)*expit(a-np.hypot(p,m)/T),cut,80.)[0]


def ledger(c,x,g,J,scalar_energy):
    I=g*g/np.pi**4*(.5*(J[1]**2+J[2]**2)+J[0]*(J[1]+J[2]))
    vis=np.pi**2*106.75/30
    return dict(case=c,x=x,g=g,field_energy_over_T4=I,field_energy_over_visible=I/vis,
                field_energy_over_scalar_energy=I/scalar_energy if scalar_energy else None)


def run():
    d=load('results/early-pole-spectrum-third-m06-2026-11-26.json.gz')
    S=load('results/early-pole-separate-fermion-heat-2026-11-26.json.gz')
    out=[]
    x=d['x_end'];r=d['mass_ratio']
    for cut in (0.,.05,.5):
        jp,jb=(cell_J(d[f'endpoint_{lab}_momentum_edges'],d[f'endpoint_{lab}_occupations'],m,cut)
               for lab,m in [('parent',x),('fermion',r*x)])
        out.append(dict(cut=cut,**ledger('initial',x,.005,[photon(1.,0.,cut),jp,jb],
                                        d['checkpoints'][-1]['integrated_scalar_energy_over_T4'])))
        for lab,s in S.items():
            if s['g']==0:continue
            xe=s['xend']
            jp=moment(xe,s['final_parent_count_over_T3'],s['final_parent_T'],cut)
            jb=moment(r*xe,s['final_daughter_count_over_T3'],s['final_daughter_T'],cut)
            jg=photon(s['final_photon_T'],s['final_photon_chemical_alpha'],cut)
            out.append(dict(cut=cut,**ledger(lab,xe,s['g'],[jg,jp,jb],None)))
    return dict(scope='Frozen energy size only; no momentum- or time-local augmented collision operator.',rows=out)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
