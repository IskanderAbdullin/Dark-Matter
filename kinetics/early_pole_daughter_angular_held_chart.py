"""Optional held unsplit chart with the validated positive angular scalar quadrature."""
import argparse
import gzip
import json
from hashlib import sha256
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np

import early_pole_joined_simultaneous_chart_reference as chart
from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_daughter_angular_mode_initializer import regrid
from early_pole_frozen_dynamic_vlasov_reference import load


def run(bundle,mesh,local=False,tight=False,initialize_only=False):
    with gzip.open(bundle,'rt') as stream:
        original=json.load(stream)
    volume=json.loads(Path(mesh).read_text())['shared_scalar_volume_edges']
    inputid=sha256(np.asarray(volume,dtype='<f8').tobytes()).hexdigest()
    grid=original['branches']['steps']['grid']
    old=load(f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz')
    out={}
    moved={}
    original_ref=chart.Screened

    class Angular(original_ref):
        def __init__(self,grid,wide,count,exposure,quadrature,cut,strength,*a,**k):
            if count!=2*(len(volume)-1):
                raise ValueError('The declared angular quadrature does not have the expected number of nodes.')
            super().__init__(grid,wide,8,exposure,quadrature,cut,strength,*a,**k)
            regrid(self,old,volume,2)

    items={}
    for pr,v in original['branches'].items():
        if pr not in ('steps','cell-conservative'):
            continue
        n=v['snapshots'][0]
        depth=v.get('aperture_subdivision',1)
        model=Cold(grid,1,0,8,True,exposure=v['duration'],reference_aperture_grid=v['aperture_grid'],
                   archive_projection=pr,aperture_subdivision=depth)
        np.testing.assert_array_equal(model.ns,n['cold_bins'])
        np.testing.assert_array_equal(model.M.init,n['hard_moments'])
        regrid(model,old,volume,2)
        source=Cold.rhs(model,np.r_[model.M.init.ravel(),model.ns,0.,0.],False)[model.j0:-2]
        p=dict(v)
        p['snapshots']=[dict(n,cold_bins=model.ns.tolist(),instantaneous_source_by_q=source.tolist())]
        items[pr]=p
        moved[pr]=dict(seed=model.ns.tolist(),nodes=model.q.tolist(),weights=model.sw.tolist(),
                       initial_rate=float(source.sum()))
    for k in ('seed','nodes','weights'):
        np.testing.assert_array_equal(moved['steps'][k],moved['cell-conservative'][k])
    if initialize_only:
        return dict(initialization_only=True,angular_mesh_sha256=inputid,receipts=moved)
    with TemporaryDirectory() as tmp:
        p=Path(tmp)/'angular-archive.json.gz'
        with gzip.open(p,'wt') as f:
            json.dump(dict(branches=items),f,allow_nan=False)
        with patch.object(chart,'Screened',Angular):
            out=chart.run(str(p),'steps',3.3333333333333335e-7,16,
                          3e-10 if tight else 3e-9,2e-15 if tight else 2e-14,
                          local_from_paired_archive=local,record_mode_history=True)
    out['angular_mesh_sha256']=inputid
    out['positive_gauss_nodes_per_panel']=2
    out['initialized_scalar_volume_weights']=moved['steps']['weights']
    np.testing.assert_array_equal(out['cold_momentum_nodes_over_T'],moved['steps']['nodes'])
    np.testing.assert_array_equal(out['cold_mode_count_history'][0],moved['steps']['seed'])
    return out


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('bundle')
    p.add_argument('mesh')
    p.add_argument('--local',action='store_true')
    p.add_argument('--tight',action='store_true')
    p.add_argument('--initialize-only',action='store_true')
    v=p.parse_args()
    print(json.dumps(run(**vars(v)),indent=2,allow_nan=False))
