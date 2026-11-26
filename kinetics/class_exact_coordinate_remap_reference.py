"""Reparameterize native transformed momentum nodes without changing the PSD."""
import argparse
import json
from pathlib import Path

import numpy as np

from archived_class_massive_moment_dilation_reference import continuum,moments,native,spline_tail

X=np.array([.01,.3,2.,8.,40.,1e7])


def scan(path,n):
    ks,f=spline_tail(path)
    p,v=continuum(ks,None)
    ref=[]
    for x in X:
        a,b,k=moments(p,v,f,x)
        ref.append(a)
    ref=np.asarray(ref)
    q,w=native(n,None)
    def score(S):
        rows=[moments(q*S,w*S,f,x) for x in X]
        a=np.asarray([r[0] for r in rows]); d=np.asarray([r[1] for r in rows])
        dr=-np.asarray([[3*(A+B),5*B-r[2]] for (A,B),r in zip(a,rows)])
        e=np.max(abs(a/ref-1),axis=0)
        z=np.max(abs(d/dr-1),axis=0)
        return dict(scale=S,maximum_residual=float(max(e.max(),z.max())),
                    energy_pressure_residual=e.tolist(),dilation_residual=z.tolist())
    ans=min((score(s) for s in np.geomspace(.18,12,950)),key=lambda y:y["maximum_residual"])
    S=ans["scale"]
    def equivalent(z):
        g,d=f(S*z)
        return g,S*d
    err=[]
    for x in X:
        A,D,K=moments(p,v,f,x)
        B,E,J=moments(p/S,v/S,equivalent,x/S)
        err.append(max(np.max(abs(A/B/S**4-1)),np.max(abs(D/E/S**4-1)),abs(K/J/S**4-1)))
    assert max(err)<1e-12
    return dict(nodes=n,at_native_scale_1=score(1.),best=ans,max_exact_continuum_coordinate_defect=float(max(err)))


def setup(path,out,template,nodes):
    out.mkdir(parents=True,exist_ok=True)
    source=np.loadtxt(path)
    rows=[scan(path,n) for n in nodes]
    for d in rows:
        n=d["nodes"]
        S=d["best"]["scale"]
        name=f"psd_coordinate_remap{n}_sync"
        tab=out/(name+".dat")
        pp=source.copy()
        pp[:,0]/=S
        np.savetxt(tab,pp,fmt="%.16e")
        s=template.read_text()
        def put(field,value):
            import re
            return re.sub(r"(?m)^"+re.escape(field)+r" = .+$",f"{field} = {value}",s)
        for field,val in (("root",out/(name+"_")),("T_ncdm",f"{S:.17g}"),("ncdm_psd_filenames",tab),
                          ("ncdm_quadrature_strategy","2"),("ncdm_N_momentum_bins",n)):
            assert f"{field} = " in s
            s=put(field,val)
        (out/(name+".ini")).write_text(s)
    return rows


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("stable_table",type=Path)
    p.add_argument("template",type=Path)
    p.add_argument("out",type=Path)
    p.add_argument("--nodes",nargs="+",type=int,default=[128,224,320])
    a=p.parse_args()
    print(json.dumps(dict(momentum_mass_ratios_for_selection=X.tolist(),
        controls=setup(a.stable_table.resolve(),a.out.resolve(),a.template.resolve(),a.nodes),
        scope="Background and derivative variants are mathematically the same continuum physical shape under q'=q/S and T_ncdm'=S, chosen only by external moment residuals. No perturbation or physical transfer result here."),
        indent=2,allow_nan=False))
