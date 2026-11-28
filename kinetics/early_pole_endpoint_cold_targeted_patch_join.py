from pathlib import Path
import argparse
import json
import numpy as np
from early_pole_endpoint_continuum_cold_band_receipt import read, reduce


def merge(base, patch, plan, output):
    old=read(base); new=read(patch)
    edges=plan['edges']
    replace=set(plan['recalculated_panel_indices'])
    if not old or not new or {r['panel'] for r in new}!=replace or len(new)!=len(replace):
        raise ArithmeticError('Incomplete targeted patch')
    meta={k:v for k,v in plan.items() if k not in ('edges','recalculated_panel_indices','total')}
    if any(r.get(k)!=v for r in old+new for k,v in meta.items()):
        raise ArithmeticError('Not the same state and quadrature')
    old_by_edge={tuple(r['bounds']):r for r in old}
    patch_by_index={r['panel']:r for r in new}
    replaced_original=[]
    uncovered=set(old_by_edge)
    out=[]
    for i,(a,b) in enumerate(zip(edges[:-1],edges[1:],strict=True)):
        if i in replace:
            r=patch_by_index[i]
            if r['bounds']!=[a,b] or r['total']!=len(edges)-1:
                raise ArithmeticError('Wrong patch support')
        else:
            r=old_by_edge.pop((a,b))
            uncovered.remove((a,b))
        out.append(dict(r,panel=i,total=len(edges)-1))
    for a,b in sorted(uncovered):
        section=[r for r in out if a<=r['bounds'][0]<r['bounds'][1]<=b]
        if len(section)!=2 or section[0]['bounds'][0]!=a or section[-1]['bounds'][1]!=b:
            raise ArithmeticError('Unexpected uncovered original panel')
        direct=np.sum([r['coeff'] for r in section],axis=0)
        before=np.array(old_by_edge[(a,b)]['coeff'])
        scale=np.maximum(np.max(abs(before),axis=1),1e-90)[:,None]
        err=float(np.max(abs(direct-before)/scale))
        if err>6e-5:
            raise ArithmeticError('Split panel original aggregate fails comparison')
        replaced_original.append([a,b,err])
    if 2*len(replaced_original)!=len(replace):
        raise ArithmeticError('Reintegration receipt has a consistency defect')
    if output.exists():
        raise FileExistsError(output)
    output.write_text(''.join(json.dumps(r,allow_nan=False)+'\n' for r in out))
    reduce(output)
    return dict(reintegrated_old_panels=replaced_original,result_panels=len(out))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('base',type=Path); p.add_argument('patch',type=Path)
    p.add_argument('plan',type=Path); p.add_argument('output',type=Path)
    a=p.parse_args()
    print(json.dumps(merge(a.base,a.patch,json.loads(a.plan.read_text()),a.output),indent=2,allow_nan=False))
