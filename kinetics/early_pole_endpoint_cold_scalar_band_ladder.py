"""Summarize completed fixed-charge cold-band scalar refinements without interpolating panels."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from early_pole_endpoint_continuum_cold_band_pilot import run as topology
from early_pole_endpoint_continuum_cold_band_receipt import run as receipt,read as panels
from early_pole_endpoint_scalar_subshell_chain_receipt import read


def verify_pair(files,card_path,archive,scale):
    c,s=read(card_path)
    if (c['grid'],c['x_end'],c['steps'])!=('fine',.05,2048):
        raise ArithmeticError('Use the evolved fine endpoint on the same clock')
    u=receipt(files,archive)
    t=u['full_pole_and_snapshot']
    if scale==1:
        yes=c.get('subcells_per_archived_cell')==2 and 'scalar_div' not in t and 'scalar_div' not in c
    else:
        sig=hashlib.sha256(json.dumps(c,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
        yes=(c['charge_div'],c['scalar_div'],c['cut'],t.get('state_sha256'),t.get('mirror_azimuth'))==(2,scale,1e-5,sig,True)
    if not yes:
        raise ArithmeticError('State and cold-panel scalar level do not match')
    calc=topology(c,2,None,topology_only=True,mirror_azimuth=t.get('mirror_azimuth',False))
    actual=panels(files[1])
    if any(calc[k]!=v for k,v in t.items()) or calc['edges']!=[actual[0]['bounds'][0]]+[r['bounds'][1] for r in actual]:
        raise ArithmeticError('The held scalar state does not reproduce its integrated support')
    return dict(level=scale,state_sha256_uncompressed=s,diagnostic=u)


def run(pairs,states,archive):
    levels=(1,2,4)
    entries=[verify_pair(f,p,archive,i) for f,p,i in zip(pairs,states,levels,strict=True)]
    out=[]
    for a,b in zip(entries[:-1],entries[1:],strict=True):
        u,v=[c['diagnostic'] for c in (a,b)]
        p,q=[r['raw_number_scalar_energy_and_excess'] for r in (u,v)]
        U,W=[r['full_pole_and_snapshot'] for r in (u,v)]
        nn=np.array([u['conditional_number_energy_neutral'],v['conditional_number_energy_neutral']])
        common=float(np.sqrt(nn[0,0]*nn[1,0]))
        net=(common/nn[:,0]-1).tolist()
        out.append(dict(levels=[a['level'],b['level']],relative_all_coefficients=(np.asarray(q)/np.asarray(p)-1).tolist(),
                        relative_pole_parameters=[W[k]/U[k]-1 for k in ('scalar_mass','omega','residue')],
                        relative_neutrality_roots=(nn[1]/nn[0]-1).tolist(),
                        constructed_common_occupation=common,net_over_inverse_on_common_occupation=net))
    return dict(x=.05,grid='fine',charge_subcells=2,band=[1e-6,1.5e-5],entries=entries,adjacent_contrasts=out,
                scope='Same charge partition and retained band, three independently evolved scalar charts and their own optical fits. No physical photon occupation, net rate or asymptotic error claimed.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('s1','s2','s4'):
        p.add_argument('--'+k,nargs=2,type=Path,required=True)
    p.add_argument('--states',nargs=3,type=Path,required=True)
    p.add_argument('--archive',type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run((a.s1,a.s2,a.s4),a.states,a.archive),indent=2,allow_nan=False))
