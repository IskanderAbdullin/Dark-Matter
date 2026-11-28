"""Evaluate actual retained cold-point coefficients across the D8 optical cap brackets."""
import argparse
import json
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_transverse_audit import charge_field,continuum_moments,model_callback,piecewise_logit,root
from early_pole_endpoint_peer_charge_shape_transport import unpack
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from finite_mass_qed_timelike_pole_angular_reference import angle

K=(.00140,.00145,.00150)
Q=(3.5e-6,1.5e-5)
RULES=((4,10,18),(8,16,24))


def panels(card, cap, signature):
    g=card['grid']
    d,y,m=unpack(card,g,8)
    x=.05
    if (cap['grid'],cap['x'],cap['charge_subcells'],cap['scalar_mass'])!=(g,x,8,m):
        raise ArithmeticError('Wrong cap comparison')
    mf=x*d.ratio
    F=[]; ch=[]; cb=[]
    for mass,edges,n in ((x,d.p_edges,d.wp*y[d.iP:]),(mf,d.d_edges,d.wd*y[d.ns:d.iP])):
        p,w,_,_,a,b=continuum_moments(edges,mass,n,36)
        F.append((p,w,mass))
        cb.append(model_callback(mass,edges,a,b))
        ch.append(piecewise_logit(mass,edges,a,b,n,1))
    for k in K:
        c=[z for z in cap['sample_rows'] if abs(z['k']-k)<1e-15]
        if len(c)!=1:
            raise ArithmeticError('Missing or ambiguous cap benchmark k')
        c=c[0]
        k=c['k']
        _,w,R=root(lambda beta:sum(charge_field(*a,beta,36) for a in F),k)
        if (w,R)!=(c['omega'],c['residue']):
            raise ArithmeticError('The cap archive does not use the recovered optical pole')
        for q in Q:
            b=[z for z in c['scalars'] if z['q']==q]
            if len(b)!=1:
                raise ArithmeticError('Missing cap benchmark q')
            f=y[np.searchsorted(d.ds.ed,q,side='right')-1]
            E=np.hypot(q,m)
            yield dict(grid=g,x=x,k=k,q=q,m=m,w=w,R=R,state_sha256_uncompressed=signature,
                       cap_fraction=b[0]['allowed_cosine_measure_fraction'],scalar_f=f),ch,cb,E


def point(tag,ch,cb,E):
    k,q,omega,R,f=(tag[v] for v in ('k','q','w','R','scalar_f'))
    z=[]
    nn=[]; dr=[]
    for i,j,l in RULES:
        a,n,panels,p=angle(ch,k,omega,R,E,'T',i,(j,l),q,occupation_callbacks=cb,mirror_azimuth=True)
        if np.min(a)<0 or abs(a[0]-a[1]-a[2])>4e-13*max(np.max(a),1e-80):
            raise ArithmeticError('The retained affine ledger is not positive')
        z.append((q*q/E)*np.array([(1+f)*a[0],f*a[1],a[0]+f*a[2]]))
        nn.append(panels)
        dr.append(p)
    if nn[0]!=nn[1] or max(dr)>3e-12:
        raise ArithmeticError('The paired point aperture geometries do not agree')
    if np.any(np.sign(z[0])!=np.sign(z[1])) or (tag['cap_fraction']==0 and np.any(z[1]!=0)):
        raise ArithmeticError(('Paired point support does not agree with cap',tag,z))
    nonzero=z[1]>0
    e=float(np.max(abs(z[0][nonzero]/z[1][nonzero]-1))) if np.any(nonzero) else 0.
    if e>2e-6:
        raise ArithmeticError(('Paired retained cold point not resolved',tag,e))
    a,b,s=z[1]
    return dict(tag,coefficients=z[1].tolist(),joint_relative_defect=e,angular_panels=nn[1],
                neutral_occupation=b/s if s>0 else None)


def run(states, cap_path, output):
    caps=json.loads(cap_path.read_text())
    s={}
    if output.exists():
        for raw in output.read_text().splitlines():
            r=json.loads(raw)
            key=r['grid'],r['k'],r['q']
            if key in s:
                raise ArithmeticError('Point table contains repeats')
            s[key]=r
    keys=set()
    with output.open('a',encoding='utf-8',buffering=1) as f:
        for path in states:
            c,sha=read(path)
            refs=[a for a in caps if a['grid']==c['grid'] and a['charge_subcells']==8]
            if len(refs)!=1:
                raise ArithmeticError('Need one cold cap per D8 chart')
            for a,ch,cb,E in panels(c,refs[0],sha):
                key=a['grid'],a['k'],a['q']
                if key in keys:
                    raise ArithmeticError('The input states are not independent')
                keys.add(key)
                if key in s:
                    if any(s[key].get(n)!=v for n,v in a.items()):
                        raise ArithmeticError('Cannot reuse point data for an altered state')
                    continue
                data=point(a,ch,cb,E)
                f.write(json.dumps(data,allow_nan=False)+'\n')
    if not set(s)<=keys:
        raise ArithmeticError('Unexpected data outside the declared point set')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('fine',type=Path)
    p.add_argument('third',type=Path)
    p.add_argument('caps',type=Path)
    p.add_argument('output',type=Path)
    a=p.parse_args()
    run((a.fine,a.third),a.caps,a.output)
