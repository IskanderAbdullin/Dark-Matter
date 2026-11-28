"""Remove an outer q band while preserving the exact low-q linear mesh spacing."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline


def prepare(wide,output):
    output.mkdir(parents=True,exist_ok=True)
    result=[]
    for case,bins,big in (("step70",44,70),("step168",104,168)):
        cap=29.*bins/big
        for arm in ("trace","physical"):
            q={}
            for l in (wide/("compact-"+arm+".ini")).read_text().splitlines():
                k,v=l.split("=",1);q[k.strip()]=v.strip()
            F=np.loadtxt(q["ncdm_psd_filenames"])
            f=CubicSpline(F[:,0],F[:,1])
            x=np.r_[F[F[:,0]<cap,0],cap]
            ratios=[float(np.trapezoid(x**k*f(x),x)/np.trapezoid(F[:,0]**k*F[:,1],F[:,0])) for k in (2,4)]
            assert 0.99999<min(ratios)<=max(ratios)<=1.+1e-14
            name=case+"-"+arm
            q["ncdm_N_momentum_bins"]=bins;q["ncdm_maximum_q"]=format(cap,".17g")
            q["root"]=str(output.resolve()/(name+"_"))
            (output/(name+".ini")).write_text("".join(f"{k} = {v}\n" for k,v in q.items()))
            result.append(dict(tag=name,cut_q=cap,q_count=bins,matching_untruncated_q_count=big,mesh_step=cap/bins,
                               retained_number_and_q2_fraction=ratios))
    return dict(scope="Keep exact spacing in the retained radial range while dropping only a high-q tail. "
                      "Native abundance remains normalized to the same density.",cases=result)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("wide",type=Path);p.add_argument("output",type=Path)
    a=p.parse_args();print(json.dumps(prepare(a.wide,a.output),indent=2,allow_nan=False))
