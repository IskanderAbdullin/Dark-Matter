"""Verify the retained T-mode scalar shoulder and the cold-band redistribution."""
import argparse
import gzip
import json
from math import fsum
from pathlib import Path

import numpy as np

def open_rows(path):
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'rt',encoding='utf-8') as f:
        return [json.loads(line) for line in f]


def rows_at(root,file,grid,region,support,key):
    rows=open_rows(root/file)
    total=rows[0]['total']
    if not total==len(rows) or [x['panel'] for x in rows]!=list(range(total)):
        raise ArithmeticError(('Incomplete panel sequence',file,total,len(rows)))
    if any(a['limits'][1]!=b['limits'][0] for a,b in zip(rows[:-1],rows[1:],strict=True)):
        raise ArithmeticError(('Band is not contiguous',file))
    if (rows[0]['limits'][0],rows[-1]['limits'][1])!=support:
        raise ArithmeticError(('Band support changed',file,support))
    if any(x['grid']!=grid or x['region']!=region or x['mesh']!=rows[0]['mesh'] for x in rows):
        raise ArithmeticError(('Grid/method changed within receipt',file))
    a=np.array([r[key] for r in rows])
    if not np.all(np.isfinite(a)) or np.min(a[:,:,:2])<0:
        raise ArithmeticError('Invalid gross collision coefficient')
    scale=np.maximum(np.max(abs(a),axis=-1),np.finfo(float).tiny)
    if np.max(abs((a[:,:,0]-a[:,:,1]-a[:,:,2])/scale))>3e-13:
        raise ArithmeticError(('Unseparated gain and loss',file))
    return a,[r['limits'] for r in rows]


def integral(a):
    return np.array([[fsum(a[:,i,j]) for j in range(a.shape[2])] for i in range(a.shape[1])])


def cold_redistribution(root,grid,interval,entry,pack):
    chosen=[];panels=[]
    for i in (1,2):
        a,span=rows_at(root,f'{grid}-cold-center-{i}.jsonl'+pack,grid,'cold',interval,
                       'number_scalar_energy_and_centered_coefficients')
        chosen.append(a)
        panels.append(span)
    if panels[0]!=panels[1]:
        raise ArithmeticError('Cold grid changed support on refinement')
    previous,high=chosen
    b=integral(high)
    if not np.array_equal(b,np.array(entry['integrated_number_scalar_energy_and_excess_coefficients'])):
        raise ArithmeticError('Archived cold ledger differs')
    Ns=[]
    for i,n in enumerate(entry['number_and_scalar_energy_threshold']):
        x=high[:,:,2]*n-high[:,:,1]
        p=previous[:,:,2]*n-previous[:,:,1]
        kill=abs(p[:,0]-x[:,0])
        signs=x[:,0]>0
        npos=int(np.sum(signs))
        if npos==0 or npos==len(signs) or tuple(signs)!=(True,)*npos+(False,)*(len(signs)-npos):
            raise ArithmeticError(('No single *panel-wise* cold crossover',grid,i,signs))
        err=kill+4*np.spacing(high[:,0,1])+4*np.spacing(high[:,0,2]*n)
        if np.any(abs(x[:,0])<20*err):
            raise ArithmeticError(('Unresolved sign of a cold panel',grid,i))
        gross=np.array([fsum(x[signs,0]),-fsum(x[~signs,0])])
        heat=np.array([fsum(x[signs,1]),fsum(x[~signs,1])])
        if i==0 and not abs(fsum(x[:,0]))<1e-11*min(gross):
            raise ArithmeticError('Number closure lost')
        if i==1 and not abs(fsum(x[:,1]))<1e-11*min(abs(heat)):
            raise ArithmeticError('Scalar energy closure lost')
        Ns.append(dict(trial=n,signed_number_panels=int(len(signs)),positive_count=npos,
                       first_negative_lower_q=panels[0][npos][0],
                       positive_and_negative_absolute_number=gross.tolist(),
                       signed_scalar_energy_groups=heat.tolist(),number=fsum(x[:,0]),
                       scalar_energy=fsum(x[:,1]),min_panel_number_to_mesh_uncertainty=float(np.min(abs(x[:,0])/err)),
                       conservative_turnover_fraction_of_inverse=float(min(gross)/b[0,1])))
    return b,Ns


def run(root,pack='.gz',cold_root=None):
    cold_root=root if cold_root is None else cold_root
    base=json.loads((cold_root/'finite-mass-qed-fixed-pole-cold-scalar-band-2026-11-28.json').read_text())
    result=[]
    for z in base:
        g=z['grid']
        a=[];edges=[]
        for mesh in (0,1):
            raw,ed=rows_at(root,f'{g}-shoulder-{mesh}.jsonl'+pack,g,'shoulder',
                           (4e-5,3e-4),'number_and_scalar_energy_coefficients')
            edges.append(ed)
            a.append(integral(raw))
        if edges[0]!=edges[1]:
            raise ArithmeticError('Refining the shoulder mesh changed exact support')
        b=a[-1]
        dr=float(max(abs((a[0]/b-1).ravel())))
        sc=float(b[0,1]/b[0,2])
        se=float(b[1,1]/b[1,2])
        cold=z['number_and_scalar_energy_threshold'][0]
        test=1e-8
        Cold,redistribution=cold_redistribution(cold_root,g,(1e-6,1.5e-5),z,
                                               '.gz' if cold_root!=root else pack)
        C=Cold[0]
        if dr>6e-5 or not 0.<cold<test<min(sc,a[0][0,1]/a[0][0,2]):
            raise ArithmeticError(('Unresolved common-rate interval',g,dr,cold,sc))
        two=np.array([Cold[:2],b])
        merged=Cold[:2]+b
        n=merged[0,1]/merged[0,2]
        local=two[:,0,2]*n-two[:,0,1]
        if not (cold<n<sc and local[0]>0>local[1] and
                abs(sum(local))<1e-11*min(abs(local))):
            raise ArithmeticError('Retained-band tradeoff lost')
        # Distinguish scalar energy changes due to number transfer between
        # the two disjoint intervals from the unequal mean energies of
        # forward and inverse contributions *within* each interval.
        wide=two.astype(np.longdouble)
        mean_D=wide[:,1,2]/wide[:,0,2]
        mean_B=wide[:,1,1]/wide[:,0,1]
        group=wide.sum(axis=0)
        n_wide=group[0,1]/group[0,2]
        between=wide[:,0,2]*n_wide-wide[:,0,1]
        by_transfer=between[0]*(mean_D[0]-mean_D[1])
        within=wide[:,0,1]*(mean_D-mean_B)
        derived=by_transfer+sum(within)
        direct=group[1,2]*n_wide-group[1,1]
        if abs(derived-direct)>1e-13*max(abs(derived),abs(direct)):
            raise ArithmeticError(('Two-interval energy audit',g,derived,direct))
        result.append(dict(grid=g,shoulder_panels=len(edges[0]),last_whole_coefficient_drift=dr,
                           raw_number_and_scalar_energy=b.tolist(),cold_number_neutral=cold,
                           shoulder_number_neutral=sc,shoulder_scalar_energy_neutral=se,
                           conditional_net_number_at_n_1e_8=[C[2]*test-C[1],b[0,2]*test-b[0,1]],
                           ratio_of_shoulder_to_cold_forward_coefficients=b[0,0]/C[0],
                           cold_redistribution_at_number_and_energy_neutral=redistribution,
                           two_band_number_neutral=n,two_band_number_by_band=local.tolist(),
                           two_band_scalar_energy_neutral=merged[1,1]/merged[1,2],
                           two_band_scalar_energy_at_number_neutral=float(direct),
                           two_band_scalar_energy_transfer_between_bands=float(by_transfer),
                           two_band_scalar_energy_mean_D=mean_D.astype(float).tolist(),
                           two_band_scalar_energy_mean_B=mean_B.astype(float).tolist(),
                           two_band_scalar_energy_asymmetry_within_bands=within.astype(float).tolist()))
    if {x['grid'] for x in result}!={'fine','third'}:
        raise ArithmeticError('Both real archives are required')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root',nargs='?',type=Path,default=Path(__file__).parent/'results')
    parser.add_argument('--raw-panels',action='store_true')
    parser.add_argument('--cold-root',type=Path)
    args=parser.parse_args()
    print(json.dumps(run(args.root,'' if args.raw_panels else '.gz',args.cold_root),
                     indent=2,allow_nan=False))
