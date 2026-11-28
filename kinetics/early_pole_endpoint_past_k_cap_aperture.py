"""Past-spacelike parent-cap geometry as a function of held photon momentum."""
import argparse
from functools import lru_cache
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from early_pole_cell_endpoint_collision_pilot import CellChart
from early_pole_cell_endpoint_transverse_audit import charge_field,continuum_moments,root
from early_pole_endpoint_cold_threshold_enclosure import energy
from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_stencil_remap_pilot import SPECS,X,load


def aperture(k,q,w,scalar,parent,daughter,Ecap):
    u=w-np.hypot(scalar,q)
    if u<=0:
        raise ArithmeticError('This is not the past-time-transfer test')
    p=np.sqrt((Ecap-parent)*(Ecap+parent))
    a=parent*parent-daughter*daughter+u*u+2*u*Ecap
    if a<=0:
        raise ArithmeticError('Not the positive past-transfer threshold')
    vcap=a/(np.sqrt(p*p+a)+p)
    def delta(v):
        return v*p-u*Ecap-((parent-daughter)*(parent+daughter)+(u-v)*(u+v))/2
    if abs(delta(vcap))>5e-14*max(vcap*p,1e-20) or not vcap>u:
        raise ArithmeticError('Cap inversion or branch diagnostics failed')
    cosine=((k-vcap)*(k+vcap)+q*q)/(2*k*q)
    fraction=(np.clip(cosine,-1.,1.)+1)/2
    # Verify every external classification geometrically, or the root and cap
    # against the independent parent-energy inversion when the root lies inside.
    if -1<cosine<1:
        v=np.sqrt((k-q)**2+2*k*q*(1-cosine))
        if abs(energy(parent,daughter,u,v)/Ecap-1)>8e-12:
            raise ArithmeticError('The parent-energy solution disagrees at the clipped boundary')
    elif (cosine<=-1 and delta(k+q)>3e-15) or (cosine>=1 and delta(abs(k-q))<-3e-15):
        raise ArithmeticError('The full-angle cap classification does not close')
    test_n=4096
    nodes=-1.+(np.arange(test_n)+.5)*2/test_n
    v=np.sqrt((k-q)**2+2*k*q*(1-nodes))
    approx=np.count_nonzero((v>u)&(delta(v)>=0))/test_n
    if abs(approx-fraction)>.51/test_n+2e-14:
        raise ArithmeticError('Direct finite cosine sweep did not reproduce the cap measure')
    return dict(cosine_upper_before_clipping=cosine,allowed_cosine_measure_fraction=fraction,
                vcap=vcap,maximum_v_minus_cap_requirement=(k+q)-vcap,
                independent_4096_cosine_defect=abs(approx-fraction))


def run(c):
    g=c['grid']
    if g not in ('fine','third') or c['x_end']!=.05 or c['subcells_per_archived_cell'] not in (1,2,4,8) or c['reference_x']!=X:
        raise ArithmeticError('Expected a saved native charge-only endpoint')
    x=c['x_end']
    d=CellChart(load(g)[0],SPECS[g],c['subcells_per_archived_cell'],reference_x=X)
    y=np.asarray(c['endpoint_occupation'])
    if y.shape!=(d.n,) or c['counts']!=[d.ns,d.nf,d.np] or min(y)<0 or not np.all(np.isfinite(y)) or max(y[d.ns:])>1:
        raise ArithmeticError('Archived state did not match its chart')
    m=d.state(x,y)[0]
    mf=x*d.ratio
    Ecap=float(np.hypot(d.p_edges[-1],x))
    P=[]
    for mass,e,n in ((x,d.p_edges,d.wp*y[d.iP:]),(mf,d.d_edges,d.wd*y[d.ns:d.iP])):
        p,w,_,_,_,_=continuum_moments(e,mass,n,36)
        P.append((p,w,mass))
    @lru_cache(maxsize=100)
    def solve(k):
        _,w,R=root(lambda u:sum(charge_field(*z,u,36) for z in P),k)
        return w,R
    out=[]
    ks=np.linspace(.00135,.00165,25)
    qs=(1e-6,3.5e-6,1.5e-5)
    for k in ks:
        w,R=solve(k)
        out.append(dict(k=k,omega=w,residue=R,scalars=[dict(q=q,**aperture(k,q,w,m,x,mf,Ecap)) for q in qs]))
    crossings=[]
    for i,q in enumerate(qs):
        cs=np.array([v['scalars'][i]['cosine_upper_before_clipping'] for v in out])
        for boundary in (-1.,1.):
            def objective(k):
                w,_=solve(k)
                return aperture(k,q,w,m,x,mf,Ecap)['cosine_upper_before_clipping']-boundary
            for j in np.where((cs[:-1]-boundary)*(cs[1:]-boundary)<0)[0]:
                k=brentq(objective,ks[j],ks[j+1],xtol=5e-16)
                if abs(objective(k))>2e-8:
                    raise ArithmeticError('Optical–geometric crossing failed a replay')
                crossings.append(dict(q=q,target_cosine=boundary,bracket=ks[j:j+2].tolist(),sample_crossing_k=k))
    return dict(grid=g,x=x,charge_subcells=c['subcells_per_archived_cell'],scalar_mass=m,parent_energy_cap=Ecap,
                sample_rows=out,sign_change_crossings=crossings,
                scope='Parent-cap-only past-spacelike geometric permission with held fitted charges. No occupation guarantee, inverse/forward weighting or full angular support claim.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('archives',type=Path,nargs='+')
    a=p.parse_args()
    print(json.dumps([run(read(f)) for f in a.archives],indent=2,allow_nan=False))
