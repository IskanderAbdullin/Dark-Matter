"""Formal frozen-snapshot vacuum Bose flux; no kinetic-validity claim."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np

from early_pole_scalar_elastic_stencil_reference import Elastic


def one(spec,order):
    q=np.array(spec['scalar_momentum_nodes'])
    edges=np.array(spec['scalar_momentum_edges'])
    raw=np.array(spec['endpoint_scalar_occupations'])
    f=np.maximum(raw,0.)
    r=spec['checkpoints'][-1]
    lam=7.65175e-8
    d=Elastic(q,r['scalar_pole_over_T'],edges,order,lam)
    n=d.vol*f
    assert abs(d.vol@np.minimum(raw,0.)) < 1e-10*n.sum()
    v,H=d.flux(f)
    Q=n.sum()
    ir=edges[1:]<=1e-5
    total=np.sum(abs(v))
    pk=int(np.argmax(f))
    # Href/T is tied to the explicit 86 GeV illustrative parent.
    h=np.sqrt(np.pi**2*106.75/90)*(86./r['x'])/2.435e18
    assert total>0 and H>0
    nr=abs(v.sum())/total
    er=abs(d.K@v)/sum(d.K*abs(v))
    # The choice of kinetic energy removes a large m*N cancellation.
    assert max(nr,er)<5e-12
    tau=total/Q/h
    decay=r['net_log_growth_rate']
    return dict(x=r['x'],mass_ratio=spec['mass_ratio'],angle_order=order,cells=d.p.size,
                projected_events=d.S.shape[0],unweighted_step_escape_fraction=d.raw_escape,
                number_rate_closure_relative=float(nr),kinetic_energy_rate_closure_relative=float(er),
                formal_absolute_node_redistribution_per_particle_per_reference_H=float(tau),
                formal_peak_signed_fractional_redistribution_per_reference_H=float(v[pk]/n[pk]/h),
                formal_IR_signed_fractional_change_per_reference_H=float(v[ir].sum()/n[ir].sum()/h),
                absolute_flux_within_initial_IR_fraction=float(abs(v[ir]).sum()/total),
                formal_redistribution_to_early_net_source_ratio=float(tau/abs(decay)),
                entropy_is_finite=bool(np.isfinite(H)),
                entropy_production_over_T=float(H) if np.isfinite(H) else None)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('spectra',nargs='+')
    p.add_argument('--angles',type=int,nargs='+',default=[6,10])
    a=p.parse_args()
    output={}
    for f in a.spectra:
        path=Path(f)
        raw=gzip.decompress(path.read_bytes()) if path.suffix=='.gz' else path.read_bytes()
        s=json.loads(raw)
        output[path.name]=[one(s,order) for order in a.angles]
    print(json.dumps(output,indent=2,allow_nan=False))
