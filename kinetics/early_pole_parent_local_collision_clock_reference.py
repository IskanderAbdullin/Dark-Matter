"""Compare momentum-local drain clocks with the illustrative parent mix rates."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np

from early_pole_massive_daughter_energy_stencil_reference import MassiveFermion
from early_pole_live_parent_energy_reference import state_moments


def clocks(spec):
    ns,nf,np_=spec['cells']
    # Main spectral cohorts use these explicit deterministic grids.
    opts={(67,112,122):(65,90,72),(97,156,176):(95,125,108)}
    d=MassiveFermion(*opts[ns,nf,np_],spec['mass_ratio'],0.,True)
    x=spec['x_end']
    y=np.r_[spec['endpoint_scalar_occupations'],spec['endpoint_fermion_occupations'],
            spec['endpoint_parent_occupations']]
    assert np.array_equal(d.pe,spec['endpoint_parent_momentum_edges'])
    E,_D,_p,m=state_moments(d,x,y)
    p=y[d.iP:]
    v=d.at(x,m,y)[d.iP:]
    npop=d.pw*p
    drain=-d.pw*v
    pos=np.maximum(drain,0.)
    crossing=pos.sum()
    keep=(p>1e-10)&(npop>2e-10)&(pos>crossing*1e-4)
    assert keep.any()
    tau=np.zeros(len(p))
    tau[keep]=pos[keep]/npop[keep]
    i=np.argmax(tau)
    isolated=d.elastic_target(x,p)
    off=isolated-p
    restore=off>0
    eqrate=np.divide(-v,off,out=np.zeros_like(p),where=keep&restore)
    W=pos[keep]
    rate=tau[keep]
    ix=np.argsort(rate)
    cum=np.cumsum(W[ix])/W.sum()
    weighted=np.interp([.5,.9,.99],cum,rate[ix])
    return dict(x=x,mass_ratio=d.ratio,snapshot_elastic_rate=spec['elastic_rate'],cells=spec['cells'],
                active_parent_cells=int(keep.sum()), fraction_positive_drain_kept=float(pos[keep].sum()/crossing),
                maximum_local_fractional_drain_per_reference_efold=float(tau[i]),
                maximum_local_drain_cell=d.pe[i:i+2].tolist(),
                gross_positive_drain_weighted_quantiles_per_reference_efold=weighted.tolist(),
                max_local_restoration_rate_to_cancel_drain_where_BGK_target_higher=float(max(eqrate)),
                fraction_of_gross_drain_with_BGK_target_higher=float(pos[keep&restore].sum()/crossing))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('spectra',nargs='+')
    args=p.parse_args()
    result={}
    for f in args.spectra:
        file=Path(f)
        raw=gzip.decompress(file.read_bytes()) if file.suffix=='.gz' else file.read_bytes()
        result[file.name]=clocks(json.loads(raw))
    print(json.dumps(result,indent=2,allow_nan=False))
