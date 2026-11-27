"""Refine the retained scalar mode integral of the frozen pointwise archive control."""
import argparse
import gzip
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from early_pole_daughter_pointwise_hard_ladder import build
from early_pole_daughter_archive_pointwise_frozen_packets import compare


def run(archive):
    with gzip.open(archive,'rt') as f:
        raw=json.load(f)
    base=np.asarray(raw['branches']['steps']['snapshots'][0]['cold_bins'])
    out={}
    with TemporaryDirectory() as d:
        p=Path(d)/'initialized.json.gz'
        for depth in (1,4):
            part=[]
            reference=None
            for q in (8,16,32,64,128,256):
                pack=build(raw,depth,q)
                seed=np.asarray(pack['branches']['steps']['snapshots'][0]['cold_bins'])
                hard=np.asarray(pack['branches']['steps']['snapshots'][0]['hard_moments'])
                if reference is None:
                    reference=hard
                np.testing.assert_array_equal(hard,reference)
                group=seed.reshape(8,q//8).sum(axis=1)
                np.testing.assert_allclose(group,base,atol=6e-18,rtol=0.)
                edges=np.geomspace(1e-6,1e-5,q+1)
                np.testing.assert_array_equal(edges[::q//8],np.geomspace(1e-6,1e-5,9))
                with gzip.open(p,'wt') as f:
                    json.dump(pack,f,allow_nan=False)
                value=compare(p,quadrature_order=12)
                rows=np.asarray(value['aperture_terms_by_mode_and_cell'])
                rates=rows.sum(axis=1)
                grouped_rates=rates.reshape(8,q//8,3).sum(axis=1)
                part.append(dict(cold_modes=q,
                                 max_grouped_seed_defect=float(max(abs(group-base))),
                                 direct_pair=value['pointwise_archive_paired_total'],
                                 fitted_pair=value['reconstructed_paired_total'],
                                 fixed_packet_curvature=value['reconstructed_curvature_fixed_packet_geometry'],
                                 sample_counter=value['sample_counter'],
                                 largest_positive_direct_paired_mode=float(max(rates[:,0])),
                                 grouped_rate_by_diagnostic=grouped_rates.tolist()))
            for v in part:
                v['delta_direct_vs_256']=v['direct_pair']-part[-1]['direct_pair']
                v['delta_fitted_vs_256']=v['fitted_pair']-part[-1]['fitted_pair']
            out[str(depth)]=part
    return dict(scope='Frozen occupied cold window; volume-centroid scalar mode ladders with no reassigned hard occupations or evolved state.',
                hard_depths=out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('archive')
    q=p.parse_args()
    print(json.dumps(run(q.archive),indent=2,allow_nan=False))
