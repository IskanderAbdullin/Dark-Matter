"""Coherent bare numerator identity on the coincident sequential-pole locus."""
import gzip,json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from early_pole_hard_yukawa_photoproduction_reference import I,G,dot,slash,spinors
from finite_mass_qed_nonnull_photon_ward_poles_reference import setup
from finite_mass_qed_frozen_transverse_plasmon_gate_reference import Q,Species
from finite_mass_qed_frozen_transverse_cut_reference import response


def geometry_root(geom,K,arm,az):
    ok=[]
    for E in np.linspace(.3,3.4,220):
        try:
            p,b=geom(E,az)
        except ValueError:
            continue
        ok.append(E)
    def fun(E):
        P,b=geom(E,az)
        return dot(K,K)+(2*dot(P,K) if arm==0 else -2*dot(b,K))
    for a,b in zip(ok[:-1],ok[1:],strict=True):
        if fun(a)*fun(b)<0:return brentq(fun,a,b,xtol=2e-14)
    raise ArithmeticError('Coalescence locus lost one pole')


def sample(d,k,c):
    m0=d['x_end'];m1=m0*d['mass_ratio']
    es=np.hypot(Q,d['checkpoints'][-1]['scalar_pole_over_T'])
    wp=Q*k*c/es
    M,m,q,K,_,geo=setup(d,k,wp,c)
    ch=[Species(mass,d[f'endpoint_{x}_momentum_edges'],d[f'endpoint_{x}_occupations'])
        for mass,x in ((m0,'parent'),(m1,'fermion'))]
    L=float(response(ch,wp/k)[1])
    if L<=0 or dot(K,K)>=0 or abs(dot(q,K))>1e-22:
        raise ArithmeticError('Not an occupied spacelike coalescence locus')
    rows=[]
    for az in (0.,1.07,2.4,3.1):
        E0,E1=(geometry_root(geo,K,i,az) for i in (0,1))
        if abs(E0-E1)>4e-12:
            raise ArithmeticError('Both propagators did not coalesce')
        p,b=geo(E0,az)
        P,B,S=map(slash,(p,b,K))
        up=spinors(p,M)
        ub=spinors(b,m).conj().T@G[0]
        num=[]
        for pol in np.eye(4)[1:3]:
            V=slash(pol)
            A=ub@(P+S+M*I)@V@up
            C=ub@V@(B-S+m*I)@up
            target=2*dot(q,pol)*ub@up
            norm=np.sum(abs(A)**2)+np.sum(abs(C)**2)
            n=np.sum(abs(A-C)**2)
            identity=np.sqrt(np.sum(abs(A-C-target)**2)/norm)
            if identity>5e-9:
                raise ArithmeticError(('Exact external current identity',identity))
            num.append([n,norm,identity])
        coherent,diag=np.array(num)[:,:2].sum(axis=0)
        occ=[]
        for v,s in ((p,'parent'),(b,'fermion')):
            i=np.searchsorted(d[f'endpoint_{s}_momentum_edges'],np.linalg.norm(v[1:]))-1
            occ.append(d[f'endpoint_{s}_occupations'][i])
        if not all(0<x<1 for x in occ):
            raise ArithmeticError('Coalescent roots do not overlap occupied bins')
        rows.append(dict(azimuth=az,coincident_parent_energy=E0,pole_root_abs_difference=abs(E0-E1),
                         coherent_residue_over_incoherent_diagonals=float(coherent/diag),
                         amplitude_identity_relative_defect=max(x[2] for x in num),
                         actual_external_occupations=occ))
    return dict(k=k,cosine=c,coalescent_photon_frequency=wp,minus_Im_retarded_T=L,rows=rows)


def run():
    r=[]
    for g in ('fine','third'):
        with gzip.open(Path(__file__).parent/'results'/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as f:d=json.load(f)
        r.append(dict(grid=g,checks=[sample(d,k,c) for k in (.0008,.00128,.00256) for c in (.47,.91)]))
    return dict(scope='Bare coherent residue on q·K=0 only; not an arbitrary width assumption or a prescription at a merged pole.',data=r)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
