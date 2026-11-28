"""Inspect CLASS's manual radial trapezoid against a positive continuous input."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import simpson


def run(psd):
    q,f=np.loadtxt(psd).T
    interp=CubicSpline(q,f)
    powers=np.array([1,2,3,4,6])
    table=simpson(q[None,:]**powers[:,None]*f,x=q,axis=1)
    vals=[]
    for n,K in ((44,29*44/70),(70,29.),(96,29.),(104,29*104/168),(128,29.),(168,29.),(336,29.)):
        step=K/n
        x=step*np.arange(1,n+1)
        v=interp(x)
        assert min(v)>0.
        w=np.full(n,step);w[-1]*=.5
        val=(x[None,:]**powers[:,None]*v)@w
        vals.append(dict(q_nodes=n,q_cap=K,step=step,positive_raw_number_defect=float(val[1]/table[1]-1),
             normalized_relative_errors_qminus1_q1_q2_q4=(np.delete(val/val[1],1)/np.delete(table/table[1],1)-1).tolist()))
    return dict(scope="Manual linear trapezoid as in CLASS, including zero-weight q=0 and half-weight outer edge. "
                      "Moment errors of the input do not certify its perturbation response.",
                PSD_sha256=hashlib.sha256(psd.read_bytes()).hexdigest(),rows=vals)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("PSD",type=Path)
    a=p.parse_args();print(json.dumps(run(a.PSD),indent=2,allow_nan=False))
