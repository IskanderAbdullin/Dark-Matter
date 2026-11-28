"""Bound signed CDF displacements for original-display counts with shared bin shapes."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np


def band(edges,u):
    e=np.asarray(edges,dtype=float)
    v=np.asarray(u,dtype=float)
    if (len(e)!=len(v)+1 or np.min(np.diff(e))<=0 or not np.all(np.isfinite(e))
            or abs(fsum(v))>2e-14 or not np.all(np.isfinite(v))):
        raise ArithmeticError('A signed probability contrast must close on the shared support')
    cuts=np.r_[0.,np.cumsum(v)]
    cuts[-1]=0.
    a,b=cuts[:-1],cuts[1:]
    lo=np.where(a*b<=0.,0.,np.minimum(abs(a),abs(b)))
    hi=np.maximum(abs(a),abs(b))
    w=np.diff(e)
    # Sharp on a common conditional shape: each bin can realize
    # any constant CDF between its two endpoint values by mixing endpoint
    # atoms; left/right endpoint choices realize the upper extremes.
    i=int(np.argmax(abs(cuts)))
    return dict(flow_cdf_L1_q_lower_upper=[fsum(w*lo),fsum(w*hi)],
                fixed_common_shape_TV=fsum(abs(v))/2,
                signed_profile_by_bin=v.tolist(),signed_CDF_at_cuts=cuts.tolist(),
                absolute_max_CDF_cut_and_signed_value=[float(e[i]),float(cuts[i])])


def run(paths):
    out=[]
    edges=None
    for path, pair in zip(paths,([2,4],[4,8],[8,16]),strict=True):
        raw=path.read_bytes()
        x=json.loads(raw)
        if x['levels']!=pair or len(x['common_q_edges'])!=18 or x['maximum_paired_profile_defect']>2e-9:
            raise ArithmeticError('Need the three accepted direct two-order rectangles')
        if edges is None:
            edges=x['common_q_edges']
        if x['common_q_edges']!=edges:
            raise ArithmeticError('The common display must not change')
        r={}
        for name,z in x['observables'].items():
            a,b,c,d=np.asarray(z['normalized_profile_rectangle']).reshape(4,-1)
            u=(b-a+d-c)/2
            v=(c-a+d-b)/2
            if np.max(abs(u+v-(d-a)))>2e-14:
                raise ArithmeticError('Signed two-order composition failed')
            all_=[band(edges,w) for w in (d-a,u,v,c-a,d-b,b-a,d-c)]
            if (abs(all_[0]['fixed_common_shape_TV']-z['native_profile_TV'])>2e-14 or
                np.max(abs(u-z['symmetric_scalar_profile_contrast']))>2e-14):
                raise ArithmeticError('Displayed effects differ from the accepted integrated rectangle')
            # Triangle and reverse triangle hold at ANY one common shape,
            # not only at independently optimized endpoint witnesses.
            dl,du=all_[0]['flow_cdf_L1_q_lower_upper']
            pl,pu=all_[1]['flow_cdf_L1_q_lower_upper']
            kl,ku=all_[2]['flow_cdf_L1_q_lower_upper']
            if dl+2e-18<pl-ku or du>pu+ku+2e-18:
                raise ArithmeticError('CDF intervals violate the signed triangle inequality')
            # A difference of two normalized corner vectors, and each
            # symmetric contrast, may change each bin by at most 2d
            # when a corner changes by d. Add an explicit float guard.
            wire_bound=2*(len(edges)-1)*(edges[-1]-edges[0])*(
                         x['maximum_paired_profile_defect']+64*np.finfo(float).eps)
            r[name]=dict(native=all_[0],explicit_scalar=all_[1],combined_kernel=all_[2],
                         bare_kernel_at_each_population=all_[3:5],population_at_each_kernel=all_[5:],
                         paired_mesh_cdf_L1_change_cap_q=wire_bound,
                         population_to_kernel_min_ratio_before_wire_guard=(pl/ku if ku>0 else None),
                         guarded_population_to_kernel_min_ratio=max(0.,pl-wire_bound)/(ku+wire_bound))
        out.append(dict(levels=pair,source_sha256=hashlib.sha256(raw).hexdigest(),observables=r))
    return dict(common_q_edges=edges,pairs=out,
                scope='One normalized in-bin conditional law is shared by all four counterfactual charts. Ranges are sharp for each separate signed contrast under that condition, but need not be attained at the same conditional law; no actual shape or dynamics inferred.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--rectangles',nargs=3,type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.rectangles),indent=2,allow_nan=False))
