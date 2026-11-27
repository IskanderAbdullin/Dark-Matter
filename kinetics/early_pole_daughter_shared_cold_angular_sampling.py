"""Test a shared event-aligned scalar mesh with one or two points per volume interval."""
import argparse
import json
from pathlib import Path

import numpy as np

from early_pole_daughter_direct_cold_volume_quadrature import angular_edge_modes,direct_cold
from early_pole_daughter_local_fraction_curvature_attribution import initial


def run(original,refined,reference):
    g={}
    files=(original,refined)
    for f in files:
        C,seed,states,_,dep=initial(f)
        e=C.M.weights(states[0],True)[0]
        g[dep]=angular_edge_modes(C.M,e,C.scalar(seed)[-1],257)
    if sorted(g)!=[1,4]:
        raise ValueError('Expected both validated hard geometries.')
    # Exact duplicate knots only are removed; preserve actually distinct
    # calculated faces, even if the new intervals are narrow.
    shared=np.unique(np.concatenate(list(g.values())))
    running={}
    global_edges=None
    for o in (1,2,4):
        for file in files:
            r=direct_cold(file,o,1,shared_angular_edges=shared)
            edges=r.pop('shared_scalar_volume_edges')
            if global_edges is None:
                global_edges=edges
            if edges!=global_edges:
                raise ValueError('Changed cold integration nodes or physical support.')
            k=str(r['hard_depth'])
            base=reference['hard_depths'][k]['reference_paired_rate']
            g0=reference['hard_depths'][k]['geometric_quadratures'][-1]['grouped_rate']
            r['relative_abs_error_vs_converged_volume']=abs((r['direct_scalar_volume_rate']-base)/base)
            r['max_abs_coarse_group_error_vs_converged_volume']=float(max(abs(np.asarray(r['grouped_rate'])-g0)))
            running.setdefault(k,[]).append(r)
    return dict(scope='Shared frozen event-aligned static scalar quadrature; no dynamics and no online remeshing.',
                two_mesh_knot_union=shared.tolist(),shared_scalar_volume_edges=global_edges,
                per_hard_depth=running)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('original','refined','volume_reference'):
        p.add_argument(k)
    v=p.parse_args()
    print(json.dumps(run(v.original,v.refined,json.loads(Path(v.volume_reference).read_text())),
                     indent=2,allow_nan=False))
