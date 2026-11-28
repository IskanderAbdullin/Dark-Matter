"""Execute the independent metric checks after direct midpoint child-rate integration."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_continuum_cold_band_pilot import refine_panels
from early_pole_endpoint_continuum_cold_band_receipt import run as receipt
from early_pole_endpoint_cold_non_nested_cdf_contrast import audit, overlap_tv
from early_pole_endpoint_cold_normalized_profile_envelope import geometry
from early_pole_endpoint_cold_panel_moment_bounds import bounds
from early_pole_endpoint_cold_panel_transport_upper import most
from early_pole_endpoint_cold_profile_interval_transport import fixed,transport


CHANNELS = (('forward_number',0,0),('inverse_number',0,1),('forward_excess',2,0))


def probability(q):
    return np.asarray(q) / fsum(q)


def project(e,parent,f,child):
    """Conserve the original panel masses using only the measured child ratios."""
    idx = np.searchsorted(e,f[:-1],side='right')-1
    if (np.any(idx < 0) or np.any(idx >= len(parent)) or f[-1] != e[-1]
            or not np.all(np.isin(e,f)) or np.any(f[1:] > e[idx+1])
            or not np.all(np.asarray(parent)>0) or not np.all(np.asarray(child)>0)):
        raise ArithmeticError('Child projection requires positive exact subdivisions')
    subtotal = np.bincount(idx,weights=child,minlength=len(parent))
    z = child*(parent/subtotal)[idx]
    if np.max(abs(np.bincount(idx,weights=z)/parent-1))>2e-15:
        raise ArithmeticError('Projected within-parent conservation check failed')
    dist = .5*fsum(abs(probability(child)-probability(z)))
    delta = np.max(abs(subtotal/parent-1.))
    return z,dist,float(delta)


def metrics(e,f,a,b,m):
    a,b=map(probability,(a,b))
    lo=transport(fixed(dict(base=a)),fixed(dict(base=b)),e,f,'fixed',independent_bracket=True)
    hi=most(e,f,a,b)
    tv=overlap_tv(e,f,a,b)
    return dict(allowance_q=fsum(a*np.diff(e))+fsum(b*np.diff(f)),
                minimum_q_bracket=lo['independent_cost_lower_upper_q'],sharp_upper_q=hi['sharp_upper_q'],
                minimum_TV_bracket=tv['minimum_TV_bounds'],
                mean_rms_speed_shift=(bounds(f,b,m)-bounds(e,a,m)[:,::-1]).tolist())


def run(root,archive,design):
    d='2026-11-28'
    feed=[(g,[root/f'{g}-d2-endpoint-continuum-cold-{i}.jsonl.gz' for i in (1,2)],
           root/f'early-pole-cell-endpoint-{g}-d2-x05-state-{d}.json.gz') for g in ('fine','third')]
    original=[audit(v,s,archive,g) for g,v,s in feed]
    p=design['plans'][0]
    if [o[1] for o in original]!=design['states_sha256_uncompressed'] or p['original_panels_to_split']!=8:
        raise ArithmeticError('Need the locked eight-target design and the original exact fine/third states')
    refined=[]
    checks=[]
    for (g,v,s),o,split in zip(feed,original,p['split_plans_by_grid'],strict=True):
        if not split:
            new=v
        else:
            new=[root/f'{g}-d2-placement-eight-endpoint-continuum-cold-{i}.jsonl.gz' for i in (1,2)]
        check=receipt(new,archive)
        if check['full_pole_and_snapshot']!=o[2]['full_pole_and_snapshot'] or np.max(abs(np.asarray(
            check['raw_number_scalar_energy_and_excess'])/o[2]['raw_number_scalar_energy_and_excess']-1.))>6e-5:
            raise ArithmeticError('Direct certification changed the original optical or gross reference')
        q=[geometry(f) for f in new]
        for (now,_),(before,_) in zip(q,o[0],strict=True):
            exact=refine_panels(before,split)[0] if split else before
            if not np.array_equal(now,exact):
                raise ArithmeticError('Refined support does not implement the locked exact split')
        refined.append(q)
        checks.append(check)
    m=original[0][2]['full_pole_and_snapshot']['scalar_mass']
    row = []
    for paired,orig_pair in zip(zip(*refined,strict=True),zip(original[0][0],original[1][0],strict=True),strict=True):
        (e,A),(f,B) = paired
        (g,G),(h,H) = orig_pair
        ans = {}
        for label,i,j in CHANNELS:
            a,b = A[:,i,j],B[:,i,j]
            za,tva,da = project(g,G[:,i,j],e,a)
            zb,tvb,db = project(h,H[:,i,j],f,b)
            empirical=metrics(e,f,a,b,m)
            conditional=metrics(e,f,za,zb,m)
            baseline=metrics(g,h,G[:,i,j],H[:,i,j],m)
            if label=='forward_number' and abs(empirical['allowance_q']-p['post_refinement_placement_allowance_q'])>2e-12:
                raise ArithmeticError('Actual child totals do not reproduce the equal-midpoint width allowance')
            # Changing normalized panel totals by TV can change any W1 extremum
            # by at most band diameter times TV, since common panel mass
            # can preserve its locations and the rest can move anywhere in the band.
            diameter=e[-1]-e[0]
            if [e[0],e[-1]]!=[f[0],f[-1]] or not diameter>0:
                raise ArithmeticError('The two profiles have different retained bands')
            cap=diameter*(tva+tvb)
            guard=2e-13
            for key,tol in (('minimum_q_bracket',cap+guard),('sharp_upper_q',cap+guard),
                            ('minimum_TV_bracket',tva+tvb+2e-9)):
                eb=np.atleast_1d(empirical[key]); cb=np.atleast_1d(conditional[key])
                if max(abs(eb-cb))>tol:
                    raise ArithmeticError('Parent-total normalization changes a metric past its support bound')
            if (conditional['minimum_q_bracket'][1]+guard<baseline['minimum_q_bracket'][0]
                    or conditional['sharp_upper_q']>baseline['sharp_upper_q']+guard
                    or conditional['minimum_TV_bracket'][1]+2e-9<baseline['minimum_TV_bracket'][0]):
                raise ArithmeticError('Exact parent conditioning violates nested support implications')
            ans[label]=dict(original=baseline,empirical=empirical,parent_conservative=conditional,
                            child_acceptance=dict(maximum_relative_parent_receipt_defect=[da,db],
                                normalized_panel_TV_correction=[tva,tvb],
                                max_W1_perturbation_q_after_exact_parent_projection=cap))
        row.append(ans)
    for label,_,_ in CHANNELS:
        for method in ('original','empirical','parent_conservative'):
            for key in row[0][label][method]:
                if np.max(abs(np.asarray(row[0][label][method][key])-row[1][label][method][key]))>(2e-8 if key=='minimum_TV_bracket' else 2e-12):
                    raise ArithmeticError('Two q meshes do not independently reproduce refined diagnostics')
    return dict(x=.05,original_grid_labels=['fine','third'],design=design,refined_pair_receipts=checks,
                forward_number=row[1]['forward_number']['empirical'],compared_measures=row[1],
                scope='Direct eight-parent subpanel integration at fixed native states; exact-parent projection preserves the old gross panel ledger. No remaining intra-panel shape inferred.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','archive','design'):
        p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.root,a.archive,json.loads(a.design.read_text())),indent=2,allow_nan=False))
