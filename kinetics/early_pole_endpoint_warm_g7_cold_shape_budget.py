"""Audit a uniform cold-shape release fraction for the actual G6/G7 contrast."""
import argparse
from functools import lru_cache
import json
from pathlib import Path

import numpy as np

from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack
from early_pole_endpoint_warm_g7_continuous_mass_sign import model,tighter_error
from early_pole_endpoint_warm_stress_conditional_hull import measure,values


def release_loss(n,N,lower,upper):
    """Difference between common-shape and independent lower bounds."""
    return np.minimum(n,N)*(upper-lower)


def run(native, states, prefix, covered):
    z,hs=read(native)
    lhs=prefix['rows'][-1]
    d,y=unpack(z,False,.02)
    ref=d.state(.02,y)[0]
    if ((z['x_end'],z['steps'],hs,ref)!=(.02,2048,lhs['common_native_sha256'],lhs['mass'])
            or (covered['x'],covered['reference_mass'])!=(.02,ref)
            or covered['uncertified_intervals'] or covered['proposed_observation_range']!=[0.,.1]):
        raise ArithmeticError('Require the matched .02 archive and continuous sign cover')
    charts=[]
    for k,p in enumerate(states):
        a,h=read(p)
        if h!=prefix['audit_receipts'][k][-1]['uncompressed_state_sha256']:
            raise ArithmeticError('Changed actual chart behind the cold-release control')
        charts.append(model(a,d,ref))
    cold=[[a for a in p if not a[-1]] for p in charts]
    e,f=[np.array([u[:2] for u in x])**2 for x in cold]
    if not np.array_equal(e,f):
        raise ArithmeticError('Shared cold support changed')
    n,N=[np.array([r[2] for r in x])[:,None] for x in cold]
    def available(a,b):
        A=values(e[:,0],b)
        B=values(e[:,1],a)
        # Let A,B enclose one kernel on the cell.  The worst lower
        # shared contribution is min((N-n)*A,(N-n)*B).  Under
        # independent release it becomes N*A-n*B.  Their exact
        # difference is min(n,N)*(B-A), not (n+N)*(B-A).
        S=np.sum(release_loss(n,N,A,B),axis=0)
        if np.min(S)<=0 or not np.all(np.isfinite(S)):
            raise ArithmeticError('Invalid monotone cold support width')
        return S,np.sum(n*A,axis=0)
    out=[]
    records=covered['certified_intervals']
    for k,r in enumerate(records):
        a,b=r['mass']
        if (a if k==0 else (a-records[k-1]['mass'][1]))!=0. or not a<b:
            raise ArithmeticError('Unexpected mass-cover topology')
        # Give G6 a rigorous denominator floor from all q supports.
        loss,_=available(a,b)
        floor=np.sum([r[2]*values(r[0]**2,b) for r in charts[0]],axis=0)
        rel=np.asarray(r['relative_pressure_and_pseudopressure_lower'])
        if min(rel)<=0. or np.min(floor)<=0.:
            raise ArithmeticError('The accepted sign margin is missing')
        budget=rel*floor/loss
        out.append(dict(mass=[a,b],cold_release_fraction_floor=budget.tolist()))
    if records[-1]['mass'][1]!=.1:
        raise ArithmeticError('Cover narrowed')
    # At the actual reference mass, recompute a tighter fixed-mass margin
    # independently of the partition; use live count and mean constraints.
    def moment(c):
        warm=[r for r in c if r[-1]]
        _,_,_,inner=measure(warm,ref,ref,96)
        err=np.sum([q[2]*tighter_error(q[0],q[1],ref,ref,96) for q in warm],axis=0)
        return inner+np.stack((-err,err))
    X,Y=map(moment,charts)
    C=np.sum(np.minimum((N-n)*values(e[:,0],ref),(N-n)*values(e[:,1],ref)),axis=0)
    fixed=Y[0]-X[1]+C
    loss,_=available(ref,ref)
    budget=fixed/loss
    if min(budget)<=0:
        raise ArithmeticError('No usable anchor cold-release margin')
    checks=[]
    for mass in (0.,1e-6,2e-6,3e-6,4e-6,ref,.1):
        q0=values(e[:,0],mass)
        q1=values(e[:,1],mass)
        g6,g7=[np.sum([a[5]@values(a[4]**2,mass) for a in u if a[-1]],axis=0) for u in charts]
        nodal=g7-g6
        shared=np.sum(np.minimum((N-n)*q0,(N-n)*q1),axis=0)
        free=np.sum(N*q0-n*q1,axis=0)
        penalty,_=available(mass,mass)
        if np.max(abs(shared-free-penalty))>3e-15*np.max(g6+g7):
            raise ArithmeticError('Explicit endpoint construction does not reproduce independent loss')
        checks.append(dict(formal_mass=mass,
                compatible_shared_shape_endpoint_witness_difference=(nodal+shared).tolist(),
                compatible_independent_cold_endpoint_witness_difference=(nodal+free).tolist(),
                independent_release_threshold_for_this_witness=((nodal+shared)/penalty).tolist()))
    @lru_cache(None)
    def warmed(mass):
        out=[]
        for c in charts:
            warm=[r for r in c if r[-1]]
            _,_,_,inner=measure(warm,mass,ref,96)
            err=np.sum([q[2]*tighter_error(q[0],q[1],mass,ref,96) for q in warm],axis=0)
            guard=256*np.finfo(float).eps*inner[1]
            out.append(inner+np.stack((-err-guard,err+guard)))
        return out

    def certify(eps):
        cells=[]
        failed=[]
        def visit(a,b,depth):
            aa,bb=warmed(a),warmed(b)
            low=values(e[:,0],b)
            high=values(e[:,1],a)
            shared=np.sum(np.minimum((N-n)*low,(N-n)*high),axis=0)
            unc=np.sum(N*low-n*high,axis=0)
            l=bb[1][0]-aa[0][1]+(1-eps)*shared+eps*unc
            if np.min(l)>0. or depth==25:
                dest=cells if np.min(l)>0. else failed
                dest.append(([a,b],l.tolist()))
            else:
                m=(a+b)/2
                visit(a,m,depth+1)
                visit(m,b,depth+1)
        visit(0.,.1,0)
        cuts=sorted([u[0] for u in cells+failed])
        if cuts[0][0]!=0. or cuts[-1][1]!=.1 or any(x[1]!=y[0] for x,y in zip(cuts,cuts[1:])):
            raise ArithmeticError('Cold-release interval cover has gap or overlap')
        return dict(release_fraction=eps,interval_count=len(cells),unresolved=failed,
                    worst_absolute_difference_lower=(np.min([u[1] for u in cells],axis=0).tolist() if cells else None),
                    corroborating_partition=cells)

    cert=certify(.035)
    if cert["unresolved"]:
        raise ArithmeticError("Uniform 3.5% cold release bound has an unresolved mass interval")
    # Give the adverse witness strictly interior q locations so that it
    # does not rely on which finite-volume cell owns a shared endpoint.
    t=1e-5
    left,right=np.sqrt(e).T
    selected=2e-6
    lo=values((left+t*(right-left))**2,selected)
    hi=values((right-t*(right-left))**2,selected)
    extra=.038
    g6,g7=[np.sum([a[5]@values(a[4]**2,selected) for a in u if a[-1]],axis=0) for u in charts]
    witness=g7-g6+(1-extra)*np.sum(np.minimum((N-n)*lo,(N-n)*hi),axis=0)+extra*np.sum(N*lo-n*hi,axis=0)
    if witness[1]>=-1e-16:
        raise ArithmeticError('The opposite pseudo-pressure witness is not strictly resolved')


    return dict(x=.02,range_mass=[0.,.1],reference_mass=ref,
                cover_fraction_bound=np.min([r['cold_release_fraction_floor'] for r in out],axis=0).tolist(),
                reference_fraction_bound=budget.tolist(),reference_warm_plus_shared_cold_numerator_lower=fixed.tolist(),
                reference_full_cold_release_penalty_upper=loss.tolist(),point_constructive_endpoint_witnesses=checks,intervals=out,
                direct_adaptive_release_cover=cert,
                explicit_contrast_reversal_inside_cells=dict(formal_mass=selected,release_fraction=extra,
                    inward_cell_width_fraction=t,relative_order_numerator=witness.tolist()),
                scope='Same-cell perturbation of a common cold shape independently in each chart, given their counts. Fractions strictly below bound preserve sign; a small sufficient budget does not establish ambiguity beyond it.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('native','prefix','cover'):
        p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--states',type=Path,nargs=2,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.native,a.states,read(a.prefix)[0],read(a.cover)[0]),indent=2,allow_nan=False))
