"""Check a multi-channel split restart against an uninterrupted trajectory."""
import hashlib,json
from tempfile import TemporaryDirectory
from pathlib import Path
from early_pole_cold_screened_split_reference import run
kw=dict(grid='fine',duration=2e-7,count=2,pair=True,compton_strength=1.,
        conversion_strength=1.,same_soft_strength=1.,screening_update='full',snapshots=True)
def fingerprint(d):
 return hashlib.sha256(json.dumps(d,sort_keys=True,allow_nan=False).encode()).hexdigest()
whole=run(**kw)
with TemporaryDirectory() as cache:
 f=Path(cache)/'restart.json'
 first=run(**kw,checkpoint=str(f),stop_after=1,progress_every=1)
 z=dict(first)
 replay=run(**kw,checkpoint=str(f),resume=True,progress_every=1)
 h=[fingerprint(x) for x in (whole,replay)]
 assert h[0]==h[1], [k for k in whole if whole[k]!=replay[k]]
 saved=json.loads(f.read_text())
 saved['state'][0]+=.0001
 f.write_text(json.dumps(saved))
 try: run(**kw,checkpoint=str(f),resume=True)
 except ValueError as e: z['tampered_state_error']=str(e)
 else: raise AssertionError("tampered state was accepted")
print(json.dumps(dict(inputs=kw,first=z,identical_output_hashes=h,
                      identical_final_gain=whole['scalar_gain'],energy=whole['energy_drift'],
                      entropy=whole['entropy_closure']),indent=2))
