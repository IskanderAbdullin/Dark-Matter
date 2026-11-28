"""Independent parent-first witness at the same fixed cold mode."""
import gzip
import json
from pathlib import Path
from finite_mass_qed_parent_signed_soft_reverse_check import independent
from finite_mass_qed_signed_cold_witness_reference import build

p=Path(__file__).parent/'results'
with gzip.open(p/'early-pole-spectrum-fine-m06-2026-11-26.json.gz') as r:
    s=build('fine',json.load(r),.85)
for k in (1e-5,-1e-5):
    a=s.signed(k,96,18,36)[0,2]
    b=independent(s,k)
    err=float(abs(b-a))
    print(json.dumps(dict(mix=.85,k=k,photon_first=float(a),parent_first=float(b),absolute_defect=err)),flush=True)
    if err>2e-13:
        raise ArithmeticError(('Reversed witness shell',k,err))
