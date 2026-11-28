"""Retain exact modes for a finite occupation ladder, with matched pure CDM."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np

from scalar_cold_packing_native_class_reference import MODES
from scalar_cold_packing_native_power_audit_reference import power


def run(base,ladder,out):
    k=np.array(MODES)
    res=[]
    for zid,z in ((1,5),(2,20)):
        ref={}
        ans=[]
        for tier,grid in (('compact',(100,1000,10000)),('dense',(5000,7000,10000))):
            cdm,cs=power(base/f'{tier}-cdm_00_z{zid}_pk.dat',k)
            ref[tier]=dict(power=cdm,sha256=cs)
            for cap in grid:
                p,s=power(ladder/f'{tier}-F{cap}_00_z{zid}_pk.dat',k)
                r=np.asarray(p)/cdm
                ans.append(dict(tier=tier,cap=cap,exact_power=p,sha256=s,ratio=r.tolist(),
                                 all_sampled_k_pass_reference_point95=bool(np.min(r)>=.95)))
        da=[x for x in ans if x['tier']=='dense']
        ca=[x for x in ans if x['tier']=='compact']
        assert all(np.all(np.diff(np.array([j['ratio'][i] for j in pp]))>0)
                   for i in range(1,len(k)) for pp in (da,ca))
        drift=float(np.max(abs(np.asarray(da[-1]['ratio'])/np.asarray(ca[-1]['ratio'])-1.)))
        assert drift<.003
        res.append(dict(z=z,matched_cdm=ref,cases=ans,dense_and_compact_F10000_maximum_relative_drift=drift))
    for tier in ('compact','dense'):
        data={}
        for path in sorted(ladder.glob(f'{tier}-F*_pk.dat')):
            raw=path.read_bytes()
            data[path.name]=dict(sha256=hashlib.sha256(raw).hexdigest(),raw_ascii=raw.decode('ascii'))
        assert len(data)==6
        packed=gzip.compress(json.dumps(data,sort_keys=True,separators=(',',':')).encode(),mtime=0)
        if out:
            (out/f'scalar-cold-packing-cap-ladder-{tier}-raw-2026-11-26.json.gz').write_bytes(packed)
    return dict(scope="Counterfactual fixed-abundance family at six exact modes, not an observational "
                      "cut, no continuous occupation or k inference. The reference .95 power ratio "
                      "is a chosen numerical target; cases use a smoothed filled-ball PSD.",
                modes_h_per_Mpc=k.tolist(),comparisons=res)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('base','ladder'):
        p.add_argument(k,type=Path)
    p.add_argument('--raw-out',type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.base,a.ladder,a.raw_out),indent=2,allow_nan=False))
