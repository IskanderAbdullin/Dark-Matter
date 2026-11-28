"""Verify that hard same-flavor geometry batches retain the exact packet graph."""
import argparse
import hashlib
import json
import resource
from pathlib import Path

import numpy as np

from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_same_flavor_hard_redistribution_reference import SameFlavorKernel


def fingerprint(items):
    return [hashlib.sha256(np.asarray(x,dtype='<i8' if i==0 else '<f8').tobytes()).hexdigest()
            for i,x in enumerate(items)]


def run(grid,order):
    box=Cold(grid,1,0,8,True)
    M=box.M
    K=SameFlavorKernel(M,order,'both',.15)
    memory=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    physical=M.weights(M.init,True)[:3]
    d,entropy,energy=K.collision(*physical)
    return dict(grid=grid,order=order,peak_rss_kib=memory,
                graphs=[dict(packets=len(pack[-1]),bytes_visible=sum(a.nbytes for a in pack),sha=fingerprint(pack))
                        for pack in K.graph],
                initial_entropy=entropy,initial_energy_residual=energy,
                collision_checksum=hashlib.sha256(np.asarray(d,dtype='<f8').tobytes()).hexdigest())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=('fine','third'))
    p.add_argument('--order',type=int,nargs=4,default=[1,3,4,1])
    p.add_argument('--expected',type=Path)
    a=p.parse_args()
    out=run(a.grid,a.order)
    if a.expected:
        calc=json.loads(a.expected.read_text())
        expected=next(x for x in calc['cases'] if x['grid']==a.grid and x['before']['order']==a.order)
        for graphs in (expected['before']['graphs'],expected['after']['graphs']):
            assert [x['sha'] for x in graphs]==[x['sha'] for x in out['graphs']]
    print(json.dumps(out,indent=2,allow_nan=False))
