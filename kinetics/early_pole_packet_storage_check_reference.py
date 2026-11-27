"""Fingerprint screened elastic packets independently of integer storage width."""
import argparse
import hashlib
import json
import resource
from pathlib import Path

import numpy as np

from early_pole_archived_cold_screened_elastic_reference import Screened
from early_pole_same_flavor_screened_pole_reference import SelfSoft


def digest(graph):
    h=[]
    for i,x in enumerate(graph):
        a=np.asarray(x,dtype='<i8' if i==0 else '<f8')
        h.append(hashlib.sha256(a.tobytes()).hexdigest())
    return h


def run(grid,n):
    M=Screened(grid,1,8,8e-8,(1,1,2),.08,1.,subdivision=n)
    k=SelfSoft(M.M,grid,(1,1,2),.08)
    # Measure construction separately, before the hash can allocate arrays.
    rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    data=[M.graph+M.original_geometry,*k.graph]
    lanes=[dict(packets=len(data[i][4]),bytes_visible=sum(v.nbytes for v in x),sha=digest(x))
           for i,x in enumerate(data)]
    return dict(grid=grid,subdivision=n,construction_peak_rss_kib=rss,
                graphs=dict(zip(('unlike','self_parent','self_daughter'),lanes,strict=True)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=('fine','third'))
    p.add_argument('--subdivision',type=int,default=1)
    p.add_argument('--expected',type=Path)
    a=p.parse_args()
    out=run(a.grid,a.subdivision)
    if a.expected:
        src=json.loads(a.expected.read_text())
        case=next(x for x in src['exact_graphs'] if x['grid']==a.grid and x['subdivision']==a.subdivision)
        for kind,v in out['graphs'].items():
            assert v['sha']==case['before'][kind]['sha']==case['after'][kind]['sha'],kind
    print(json.dumps(out,indent=2,allow_nan=False))
