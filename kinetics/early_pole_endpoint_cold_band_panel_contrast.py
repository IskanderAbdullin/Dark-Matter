"""Compare frozen cold-band histograms without interpolating between q panels."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_continuum_cold_band_pilot import run as topology
from early_pole_endpoint_continuum_cold_band_receipt import read,run as validate
from early_pole_endpoint_scalar_subshell_chain_receipt import read as read_state

CUTS=(3e-6,5e-6,8e-6,1e-5,1.25e-5)
PROBS=(.1,.5,.9,.99)


def bound(edges,w):
    if (len(edges)!=len(w)+1 or not np.all(np.isfinite(w)) or np.min(w)<0
            or not np.all(np.diff(edges)>0) or fsum(w)<=0):
        raise ArithmeticError('Nonpositive distribution, missing edges, or invalid panel totals')
    v=w/fsum(w)
    c=np.cumsum(v)
    if abs(c[-1]-1)>2e-14:
        raise ArithmeticError('Panel distribution failed normalization')
    # Quantiles are only enclosed by a panel. No intrapanel reconstruction.
    j=np.searchsorted(c,PROBS)
    q=np.stack((edges[j],edges[j+1]),axis=1)
    a=[]
    for p in CUTS:
        lo=fsum(v[edges[:-1]>=p])
        hi=fsum(v[edges[1:]>p])
        if not 0<=lo<=hi<=1+2e-14:
            raise ArithmeticError('Tail enclosure did not order')
        a.append([lo,hi])
    return q.tolist(),a


def pair(meshes,state,archive,d):
    state,sig=read_state(state)
    if (state['grid'],state['x_end'],state['subcells_per_archived_cell'],state['steps'])!=('fine',.05,d,2048):
        raise ArithmeticError('The supplied fine endpoint does not match the declared charge level')
    r=validate(meshes,archive)
    t=r['full_pole_and_snapshot']
    if t['grid']!='fine' or t['time_x']!=.05 or 'scalar_div' in t:
        raise ArithmeticError('Expected unrefined scalar chart with changed charge')
    j=read(meshes[1])
    edges=np.array([j[0]['bounds'][0]]+[u['bounds'][1] for u in j])
    reproduced=topology(state,2,None,topology_only=True,mirror_azimuth=t.get('mirror_azimuth',False))
    if not np.array_equal(edges,reproduced['edges']) or any(reproduced[k]!=v for k,v in t.items()):
        raise ArithmeticError('Archived states do not reproduce declared pole and cold support')
    a=np.array([u['coeff'] for u in j])
    results={}
    for name,i,k in (('gross_forward_number',0,0),('inverse_number',0,1),
                     ('signed_number_slope',0,2),('gross_forward_excess',2,0),('inverse_excess',2,1)):
        q,p=bound(edges,a[:,i,k])
        results[name]=dict(quantile_q_enclosures=q,tail_fraction_enclosures=p)
    return dict(charge_subcells=d,state_sha256_uncompressed=sig,panels=r['panels'],pole=t,
                number_energy_and_excess=r['raw_number_scalar_energy_and_excess'],
                number_neutral=r['conditional_number_energy_neutral'][0],
                energy_neutral=r['conditional_number_energy_neutral'][1],
                relative_neutral_split=r['relative_neutrality_split'],
                mesh_relative_difference=r['maximum_componentwise_mesh_difference'],
                forward_number_relative_to_original_archive=r['number_A_B_D_new_to_old'][0],
                support=results)


def run(d1,d2,states,archive):
    a,b=[pair(u,s,archive,d) for u,s,d in zip((d1,d2),states,(1,2),strict=True)]
    if a['pole']['scalar_mass']==b['pole']['scalar_mass']:
        raise ArithmeticError('The two declared charge cards did not produce distinct gap states')
    u,v=[np.asarray(z['number_energy_and_excess']) for z in (a,b)]
    n1,n2=[z['number_neutral'] for z in (a,b)]
    n=np.sqrt(n1*n2)
    net=n*np.array([u[0,2],v[0,2]])-np.array([u[0,1],v[0,1]])
    return dict(band_q=[1e-6,1.5e-5],probabilities=PROBS,tail_cut_q=CUTS,cards=[a,b],
                d2_over_d1_same_band=(v/u).tolist(),
                common_diagnostic_photon_occupation=n,
                relative_net_to_inverse_gross_on_common_occupation=(net/np.array([u[0,1],v[0,1]])).tolist(),
                scope='Two evolved scalar/charge states and two reconstructed poles at the same fixed k=.0015. No interpolation or physical photon population specified.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--d1',nargs=2,type=Path,required=True)
    p.add_argument('--d2',nargs=2,type=Path,required=True)
    p.add_argument('--states',nargs=2,type=Path,required=True)
    p.add_argument('--archive',type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.d1,a.d2,a.states,a.archive),indent=2,allow_nan=False))
