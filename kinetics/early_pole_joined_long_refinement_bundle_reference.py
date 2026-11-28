"""Validate recovered 432-step ledgers before comparing temporal and reconstruction effects."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np

from early_pole_joined_history_memory_reference import compare as memory
from early_pole_joined_split_checkpoint_ledger_reference import score
from early_pole_joined_temporal_refinement_reference import compare, fixture, select


def read(p):
    return json.loads(Path(p).read_text())


def same_width_control(old,new):
    p,q=old['physics'],new['physics']
    except_keys={'graph_code','elastic_energy_budget'}
    common=lambda d:{k:v for k,v in d.items() if k not in except_keys}
    if common(p)!=common(q) or 'elastic_energy_budget' in p or not q.get('elastic_energy_budget'):
        raise ValueError('Control may only differ in budget and implementation hash.')
    if not old['steps_done']==new['steps_done']==p['count']:
        raise ValueError('Expected complete, independently audited checkpoints.')
    oa,na=score(old),score(new)
    if oa['pass_gate'] or not na['pass_gate']:
        raise ValueError('Expected a rejected older ledger and an accepted replacement.')
    a,b=np.asarray(old['states']),np.asarray(new['states'])
    assert a.shape==b.shape==(p['count']+1,) and a[0]==b[0]
    return dict(rejected_energy_defect=oa['final_energy_defect'],
                accepted_energy_defect=na['final_energy_defect'],
                rejected_ledger_use='Same-mesh numerical sensitivity only; not a valid physical comparison.',
                maximum_absolute_count_history_change_from_tightening=float(np.max(abs(b-a))),
                final_count_change_from_tightening=float(b[-1]-a[-1]),
                accepted_checkpoint_audit=na)


def bundle(coarse,rows):
    fixture()
    if set(rows)!={'steps','cell-conservative'}:
        raise ValueError('Both hard inputs are required.')
    accepted={}
    control={}
    temporal={}
    for projection,(raw,corrected,original) in rows.items():
        a=select(coarse,projection)
        b=select(raw,projection)
        new=read(corrected)
        old=read(original)
        x=same_width_control(old,new)
        audit=x['accepted_checkpoint_audit']
        assert b['steps']==432 and new['physics']['projection']==projection
        assert b['state_history']==new['states']
        assert b['energy_drift']==audit['final_energy_defect']
        assert b['entropy_closure']==audit['final_entropy_defect']
        temporal[projection]=compare(a,b)
        control[projection]=x
        accepted[projection]=b
    return dict(scope='Finite fixed-scale conditional five-channel time refinement only.',
                temporal=temporal,budget_only_same_width=control,
                reconstruction=memory(accepted['steps'],accepted['cell-conservative']),
                branches=accepted)


def write(output,d):
    data=(json.dumps(d,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode()
    target=Path(output)
    target.write_bytes(gzip.compress(data,mtime=0))
    assert gzip.decompress(target.read_bytes())==data
    assert json.loads(data)==json.loads(gzip.decompress(target.read_bytes()))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('coarse')
    p.add_argument('--steps',nargs=3,required=True,metavar=('ACCEPTED_RAW','ACCEPTED_CKPT','REJECTED_CKPT'))
    p.add_argument('--cell-conservative',nargs=3,required=True,metavar=('ACCEPTED_RAW','ACCEPTED_CKPT','REJECTED_CKPT'))
    p.add_argument('--output',required=True)
    a=p.parse_args()
    write(a.output,bundle(a.coarse,{'steps':a.steps,'cell-conservative':a.cell_conservative}))
