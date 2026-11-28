"""Bind the varied-k raw rates, check parent and daughter caps, and compare fixed q."""
import argparse
import json
from pathlib import Path

import numpy as np

from early_pole_endpoint_past_k_cap_aperture import aperture
from early_pole_endpoint_peer_d8_k_rate_points import K,Q,panels,RULES
from early_pole_endpoint_scalar_subshell_chain_receipt import read


def run(points,states,caps,tail):
    final=(.00150,Q[0])
    u={}
    for s in points:
        key=s['grid'],s['k'],s['q']
        if key in u:
            raise ArithmeticError('Repeated k/q/control')
        u[key]=s
    if len(u)!=2*len(K)*len(Q):
        raise ArithmeticError('Need the complete declared k/q grid from both charts')
    out=[]
    seen=set()
    for c,sha in states:
        g=c['grid']
        if g not in ('fine','third') or g in seen:
            raise ArithmeticError('Need an independent chart pair')
        seen.add(g)
        cap=[m for m in caps if m['grid']==g and m['charge_subcells']==8]
        if len(cap)!=1:
            raise ArithmeticError('Missing held cap geometry')
        r=[]
        for p,ch,cb,E in panels(c,cap[0],sha):
            k,q=p['k'],p['q']
            row=u[(g,k,q)]
            if any(row.get(j)!=v for j,v in p.items()) or not 0<=row['joint_relative_defect']<=2e-6:
                raise ArithmeticError('Raw point provenance changed')
            ar=np.array(row['coefficients'])
            if ar.shape!=(3,) or not np.all(np.isfinite(ar)) or min(ar)<0 or (
                    abs(ar[0]-ar[1]-ar[2])>6e-13*max(max(ar),1e-80)):
                raise ArithmeticError('Invalid gross affine ledger')
            if (None if ar[2]==0 else ar[1]/ar[2])!=row['neutral_occupation']:
                raise ArithmeticError('Incorrect conditional neutrality')
            P,F=ch
            effective_f_cap=F.E[-1]+E-p['w']
            fc=aperture(k,q,p['w'],p['m'],P.m,F.m,effective_f_cap)['allowed_cosine_measure_fraction']
            if (row['cap_fraction']==0 and ar[0]!=0) or (fc==0 and ar[1]!=0):
                raise ArithmeticError('Nonzero rate in a forbidden cap projection')
            r.append(dict(k=k,q=q,forward_parent_cap_cosine_fraction=row['cap_fraction'],
                          reverse_daughter_cap_cosine_fraction=fc,coefficients=row['coefficients'],
                          inverse_gross_resolved_zero=bool(ar[1]==0),conditional_neutral=row['neutral_occupation'],
                          joint_mesh_defect=row['joint_relative_defect']))
        audit=[x for x in tail if x['grid']==g]
        if len(audit)!=1:
            raise ArithmeticError('Need an independent fixed-k parent-tail baseline')
        ref_tail=audit[0]
        kfinal=next(z['k'] for z in r if abs(z['k']-final[0])<1e-15 and z['q']==final[1])
        for s in r:
            base=np.array(next(v['coefficients'] for v in r if (v['k'],v['q'])==(kfinal,s['q'])))
            if min(base)<=0:
                raise ArithmeticError('The baseline k must be nonzero in every tested coefficient')
            s['ratios_to_same_q_at_k_point0015']=(np.array(s['coefficients'])/base).tolist()
        controls={}
        for q in Q:
            uu=[x for x in ref_tail['points'] if x['q']==q and x['cosine'] is None]
            a=u[(g,kfinal,q)]
            if len(uu)!=1 or ref_tail['omega']!=a['w'] or ref_tail['state_scalar_mass']!=a['m']:
                raise ArithmeticError('Independent baseline has a different fitted state')
            A,B,D=np.sum(uu[0]['group_affine'],axis=0)
            f=a['scalar_f']
            ref=q*q/np.hypot(q,a['m'])*np.array([(1+f)*A,f*B,A+f*D])
            fmax=np.max(abs(np.array(a['coefficients'])/ref-1))
            if fmax>3e-10:
                raise ArithmeticError('k=.0015 point differs from the independent angular tail sum')
            controls[str(q)]=fmax
        out.append(dict(grid=g,rows=r,independent_base_point_max_relative_defects=controls))
    return dict(x=.05,charge_subcells=8,modes='past-spacelike transverse',columns=('gross_forward','gross_inverse','signed_slope'),
                rules=RULES,grids=out,scope='Only twelve discrete (grid,k,q) diagnostics on the retained spectra. No interpolation, k/q integration, or physical net source assumed.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('points',type=Path)
    p.add_argument('caps',type=Path)
    p.add_argument('states',type=Path,nargs=2)
    p.add_argument('--tail',type=Path,required=True)
    a=p.parse_args()
    points=[json.loads(x) for x in a.points.read_text().splitlines()]
    print(json.dumps(run(points,[read(x) for x in a.states],json.loads(a.caps.read_text()),json.loads(a.tail.read_text())),indent=2,allow_nan=False))
