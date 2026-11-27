"""Split initialized response curvature by fixing packet geometry or daughter occupations."""
import argparse
import gzip
import json
from math import fsum

import numpy as np
from scipy.special import expit

from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_live_pauli_gauge_moment_decay_reference import Tool
from early_pole_decay_mode_clock_reference import G,KAPPA


def initial(bundle):
    with gzip.open(bundle,'rt') as f:
        z=json.load(f)['branches']
    a,b=(z[k] for k in ('steps','cell-conservative'))
    n=np.asarray(a['snapshots'][0]['cold_bins'])
    dep=a.get('aperture_subdivision',1)
    if dep!=b.get('aperture_subdivision',1) or dep not in (1,4):
        raise ValueError('Use matched validated grids.')
    np.testing.assert_array_equal(n,b['snapshots'][0]['cold_bins'])
    c=Cold(a['grid'],1,0,len(n),True,exposure=a['duration'],
           reference_aperture_grid=a['aperture_grid'],archive_projection='steps',
           aperture_subdivision=dep)
    M=c.M
    for j,x,y in zip(M.edge_sets,a['hard_edge_sets'],b['hard_edge_sets'],strict=True):
        np.testing.assert_array_equal(j,x)
        np.testing.assert_array_equal(j,y)
    y0,y1=(np.asarray(v['snapshots'][0]['hard_moments']) for v in (a,b))
    np.testing.assert_array_equal(y0,M.init)
    np.testing.assert_array_equal(c.ns,n)
    edge=np.asarray(M.edge_sets[2])
    lo,hi=1.4138888888888888,1.4666666666666666
    j=np.flatnonzero((edge[:-1]>=lo-2e-14)&(edge[1:]<=hi+2e-14))
    if len(j)!=2*dep or np.any(np.diff(j)!=1):
        raise ValueError('Unmatched physical daughter interval.')
    np.testing.assert_allclose(edge[[j[0],j[-1]+1]],(lo,hi),rtol=0.,atol=2e-14)
    ii=M.r[2].start+j
    ym=y0.copy()
    ym[ii]+=.5*(y1[ii]-y0[ii])
    yd=y0.copy()
    yd[ii]=y1[ii]
    return c,n,(y0,ym,yd),a['snapshots'][0]['instantaneous_source_by_q'],dep


def compare(bundle):
    c,n,Y,saved,depth=initial(bundle)
    E,_,_,scalar_f,m=c.scalar(n)
    x=.05*np.exp(getattr(c.M,'eta',0.))
    coeff=c.sw*G*KAPPA*x*((1.06*x)**2-m*m)/(c.q*E)
    states=[c.M.weights(q,True)[:2] for q in Y]
    geom,logit=zip(*states,strict=True)
    for F in logit[1:]:
        np.testing.assert_array_equal(F[c.M.r[1]],logit[0][c.M.r[1]])
    def rate(e,f,mix=None):
        out=np.zeros(len(n))
        for j,es in enumerate(E):
            for i,k,u,v,w in Tool.events(c.lambdas[j],dict(coeff=e[c.M.start:],es=es)):
                i+=c.M.start
                k+=c.M.start
                Pr=expit(f[i,0]+f[i,1]*u)
                if mix is None:
                    Dr=expit(f[k,0]+f[k,1]*v)
                else:
                    Dr=(expit(logit[0][k,0]+logit[0][k,1]*v)+
                        expit(logit[2][k,0]+logit[2][k,1]*v))/2
                out[j]+=(w*(Pr*(1-Dr)+scalar_f[j]*(Pr-Dr))*coeff[j]).sum()
        return out
    replay=np.asarray([Cold.rhs(c,np.r_[y.ravel(),n,0.,0.],False)[c.j0:-2] for y in Y])
    physical=np.asarray([rate(e,f) for e,f in states])
    np.testing.assert_array_equal(replay[0],saved)
    np.testing.assert_array_equal(physical,replay)
    fixed=np.asarray([rate(geom[0],f) for f in logit])
    geometry_only=np.asarray([rate(e,logit[0]) for e in geom])
    linear=rate(geom[0],logit[0],True)
    if np.max(abs(linear-(fixed[0]+fixed[-1])/2))>8*np.finfo(float).eps*np.max(abs(physical)):
        raise ValueError('Pointwise affine occupation control did not close.')
    curve=lambda v:v[1]-(v[0]+v[2])/2
    p,f,g=map(curve,(physical,fixed,geometry_only))
    interaction=p-f-g
    return dict(scope='Frozen diagnostic controls at native packet geometry or logits, not independent physical sources or propagated states.',
                hard_aperture_depth=depth,cold_modes=len(n),
                total_midpoint_curvature=fsum(p),fixed_geometry_curvature=fsum(f),
                fixed_logits_curvature=fsum(g),mixed_diagnostic_curvature=fsum(interaction),
                linear_daughter_occupation_control_max_abs=float(max(abs(linear-(fixed[0]+fixed[-1])/2))),
                total_midpoint_curvature_by_mode=p.tolist(),fixed_geometry_curvature_by_mode=f.tolist(),
                fixed_logits_curvature_by_mode=g.tolist(),mixed_diagnostic_curvature_by_mode=interaction.tolist())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('archives',nargs='+')
    a=p.parse_args()
    print(json.dumps([compare(f) for f in a.archives],indent=2,allow_nan=False))
