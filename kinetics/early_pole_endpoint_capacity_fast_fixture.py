"""Bitwise fast-path controls on frozen and short evolving endpoint cards."""
import json
from statistics import median
from time import process_time

import numpy as np

from early_pole_cell_endpoint_collision_pilot import CellChart
from early_pole_cell_endpoint_startup_pilot import drive
from early_pole_endpoint_capacity_fast_step import reference,step
from early_pole_endpoint_stencil_remap_pilot import SPECS,X,load


def run(g,n):
    b=load(g)[0]
    d=CellChart(b,SPECS[g],n,reference_x=X)
    x=.05
    y=d.initial
    m,_,_,H=d.state(x,y)
    ch,_=d.chart(x,m,3)

    caps=[]
    for z in (1e-6,1e-4,.1,100.):
        dt=z*H
        rate=reference(ch,d,y,dt/2)[0]
        a,qa=reference(ch,d,y,dt,rate_state=rate)
        r,qr=step(ch,d,y,dt,rate_state=rate)
        if not np.array_equal(a,r) or qa!=qr:
            raise ArithmeticError((g,n,z,qa,qr,'Fast stage is not bitwise identical'))
        caps.append(qr)
    if caps[0]!=0 or caps[-1]<.5:
        raise ArithmeticError('Need an inactive and a strongly active limiter fixture')

    timings=([],[])
    for j in (0,1,1,0,0,1):
        f=(reference,step)[j]
        t=process_time()
        for _ in range(24):
            f(ch,d,y,.00001*H)
        timings[j].append((process_time()-t)/24)

    stages=[]
    for f in (reference,step):
        t=process_time()
        row,state=drive(d,.05,16,x_start=.049997,initial=y,stepper=f)
        stages.append((row,state,process_time()-t))
    if stages[0][0]!=stages[1][0] or not np.array_equal(stages[0][1],stages[1][1]):
        raise ArithmeticError('Fast short continuation changed state or reported budgets')
    timings=[median(t) for t in timings]
    return dict(grid=g,charge_div=n,quadrature_events=len(ch[-1]),bitwise_stage_and_short_chain=True,
                trial_dt_without_clock=[1e-6,1e-4,.1,100.],limited_extent_fraction=caps,
                median_cpu_seconds_per_step=[float(v) for v in timings],step_speedup=timings[0]/timings[1],
                illustrative_short_chain_speedup=stages[0][2]/stages[1][2])


if __name__=='__main__':
    print(json.dumps([run(g,n) for g,n in (('fine',2),('fine',4),('third',4))],indent=2,allow_nan=False))
