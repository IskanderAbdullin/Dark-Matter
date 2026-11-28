"""Replay a q64 nine-cell update using a verified one-cell state only as Newton seed."""
import argparse
from pathlib import Path
import json, numpy as np
import coupled_warm_pair_continuation_reference as kernel
from joint_resonant_centered_stiff_reference import advance
from joint_resonant_frw_split_reference import pack
from warm_cold_light_pair_partition_reference import split_rate
from warm_coupled_serial_restart_reference import atomic,sha
from warm_history_lower_matched_audit_reference import measures
from warm_light_occupation_stage_audit_reference import digest, read
parser=argparse.ArgumentParser(description=__doc__)
for key in ("eos","native_start","native_end","one_cell_guide","output"):
    parser.add_argument(key,type=Path)
a=parser.parse_args()
C,N=[read(a.native_start),read(a.native_end)]
S,A=C['snapshot'],N['snapshot']
W=read(a.one_cell_guide)
for t in (C,N,W):
 assert digest(t['snapshot'])==t['snapshot_sha256']
G=W['snapshot']
assert C['snapshot_sha256']==N['previous']
b=N['plan']['step_expansion']
P=np.array(S['lab_momentum_edges_GeV'])
np.testing.assert_allclose(P/b,G['lab_momentum_edges_GeV'],rtol=8e-14)
eos=a.eos
assert sha(eos)==N['plan']['EOS_sha256']
seed=pack([np.array(G['entropy_variable_'+k]) for k in ('scalar','fermion')])
stats={}
rr=[]
old_r,old_s=kernel.rate,kernel.collide

def rate(p,theta,order):
 R,d,ratio,m=split_rate(p,theta,order,(4,8,6,8,8),cut_cells=9)
 rr.append(m)
 return R,d,ratio

def solver(model,u,h,T):
 def obs(i,r,k): status.append((i,r,k))
 status=[]
 try:
  U,m=advance(model,u,h,T,reference_guess=False,vacuum_scaled=True,
              initial_theta=seed,balanced_merit=True,observer=obs)
 except (RuntimeError,FloatingPointError,RuntimeWarning,ValueError,AssertionError) as e:
  stats.update(used_guide=False, reason=str(e), observation=status)
  return old_s(model,u,h,T,continuation=True,balanced=True)
 stats.update(used_guide=True, observation=status, **m)
 return U,m['residual'],1,0

out={}
try:
 kernel.rate, kernel.collide=rate,solver
 report=kernel.run(eos,S,1,b,checkpoint=out,track_exterior=True,local_prompt_jacobian=True,scheduled_single_block=True)
finally:
 kernel.rate,kernel.collide=old_r,old_s
assert abs(out['T_visible_GeV']/A['T_visible_GeV']-1)<3e-14
assert abs(report['energy_ledger'])<5e-13 and abs(report['count_ledger'])<5e-13
def f(s):return measures(s)[0,1,0]/(4*s['lab_momentum_edges_GeV'][1]**3/(6*np.pi**2))
R=dict(snapshot=out,snapshot_sha256=digest(out),source_sha256=C['snapshot_sha256'],
       original_output_sha256=N['snapshot_sha256'], guide_sha256=W['snapshot_sha256'],
       solver=stats, corrections=rr, ledgers=[report['energy_ledger'],report['count_ledger']],
       slopes=[(f(x)-f(S))/np.log(b) for x in (A,G,out)],
       L1_vs_native=(abs(measures(out)-measures(A)).sum(axis=-1)/measures(A).sum(axis=-1)).tolist())
atomic(a.output,R)
print(json.dumps({k:v for k,v in R.items() if k!='snapshot'},indent=2))
