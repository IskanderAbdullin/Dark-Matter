"""Compare two separately integrated coordinate choices on the same numerical contour."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def run(i,j):
    x,y=[a.read_bytes() for a in (i,j)]
    a,b=json.loads(x),json.loads(y)
    if (a['table_input_sha256']!=b['table_input_sha256']
            or a['previous_segment_sha256']!=b['previous_archive_sha256']):
        raise ArithmeticError('Independent contour inputs must be exactly the same')
    tab=[]
    for p,q in zip(a['rows'],b['rows'],strict=True):
        if p['operator']!=q['operator']:
            raise ArithmeticError('Do not swap distinct physical models')
        z=np.array(p['points'][-len(q['rows']):])
        v=np.array([[k['log_a'],*k['state'],k['relative_entropy_added'],k['normalized_work'],k['comoving_number_ratio']-1.] for k in q['rows']])
        diff=np.max(abs(z[:,[0,1,2,3,5,6,7]]-v),axis=0)
        if np.max(diff)>4e-8:
            raise ArithmeticError('Independent extensive and intensive coordinates do not resolve common dynamics')
        tab.append(dict(operator=p['operator'],maximum_coordinate_absolute_defects=diff.tolist(),
                        maximum_temperature_state_defect=np.max(diff[1:4])))
    return dict(rows=tab,fields=['log_a','T_visible','eta_dark','log_T_dark','entropy','pressure_work','number'],
                input_sha256=[hashlib.sha256(z).hexdigest() for z in (x,y)],
                scope='Same finite restricted FD dynamics in two numerical coordinate systems only.')


if __name__=='__main__':
    c=argparse.ArgumentParser(description=__doc__)
    c.add_argument('intensive',type=Path)
    c.add_argument('extensive',type=Path)
    o=c.parse_args()
    print(json.dumps(run(o.intensive,o.extensive),indent=2,allow_nan=False))
