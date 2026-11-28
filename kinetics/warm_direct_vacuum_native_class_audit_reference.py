"""Audit signed small-fraction and amplified native CLASS power differences."""
import argparse,hashlib,json,re
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator,CubicSpline


def power(root,name,index):
    f=root/f"{name}_00_z{index}_pk.dat"
    raw=f.read_bytes()
    z=float(re.search(rb"redshift z=([^\s]+)",raw.splitlines()[0])[1])
    a=np.loadtxt(f)
    assert a.shape[1]==2 and np.all(a>0) and np.all(np.diff(a[:,0])>0)
    return z,np.log(a[:,0]),np.log(a[:,1]),hashlib.sha256(raw).hexdigest()


def run(root,binary,partial=False):
    spec=json.loads((root/"spec.json").read_text())
    f=lambda m:m["rest_fraction"]
    zmode=np.array([.002,.01,.03,.1,.3,1.,2.])
    log=np.log(zmode)
    modes={}
    for entry in spec["cases"]:
        name=entry["label"]
        logname=root/(name+".log")
        if not logname.exists() or not (root/(name+"_00_z2_pk.dat")).exists():
            if partial:continue
            raise FileNotFoundError(name)
        tier,label=name.split("-",1)
        source=power(root,name,1),power(root,name,2)
        assert [t[0] for t in source]==[5.,20.]
        ier=[]
        for z,x,y,sha in source:
            v=PchipInterpolator(x,y)(log)
            alt=CubicSpline(x,y)(log)
            ier.append((z,v,alt,sha,x,y))
        o=logname.read_text()
        stored=[float(x) for x in re.findall(r"Non-Cold Species Nr\..+?omega = ([^\s]+)",o)]
        if f(entry)>0.0:
            assert len(stored)==1 and abs(stored[0]/(.12*f(entry))-1)<6e-6
        else:
            assert not stored
        assert "Error" not in o and "switch on ncdm fluid" not in o
        modes.setdefault(tier,{})[label]=(entry,ier)
    sample=[]
    for tier,m in modes.items():
        if "base" not in m: continue
        base=m["base"][1]
        val=[]
        for name,(e,r) in m.items():
            shifts=[]
            for a,b in zip(r,base,strict=True):
                za,y,alt,sha,x,raw=a
                assert b[0]==za
                effect=y-b[1]
                interp=alt-b[2]
                exact_grid=np.array_equal(x,b[4])
                if exact_grid:
                    direct=PchipInterpolator(x,raw-b[5])(log)
                    assert max(abs(effect-direct))<5e-4
                else:
                    direct=None
                shifts.append(dict(z=za,source_sha256=sha,base_sha256=b[3],log_power_ratio=effect.tolist(),
                   comparison_grid_identical=bool(exact_grid),
                   interpolation_difference=float(max(abs(effect-interp))),
                   matched_native_log_ratio_interpolated=None if direct is None else direct.tolist()))
            val.append(dict(case=name,scalar_fraction=f(e),shape=e["kind"],shifts=shifts))
        sample.append(dict(tier=tier,cases=val))
    find=lambda tier,label:next(x for j in sample if j["tier"]==tier for x in j["cases"] if x["case"]==label)
    matrix=lambda x:np.asarray([r["log_power_ratio"] for r in x["shifts"]])
    trace=[]
    for tier in ("std","fine"):
        cases=modes.get(tier,{})
        if not {"trace","physical","amp4"}<=cases.keys():continue
        a=find(tier,"trace");b=find(tier,"physical")
        p=matrix(a);ex=matrix(b)-p
        c=find(tier,"amp4")
        tangent=(matrix(c)-p)*(b["scalar_fraction"]-a["scalar_fraction"])/(c["scalar_fraction"]-a["scalar_fraction"])
        row=dict(tier=tier,trace_presence_log_grid_shift=p.tolist(),physical_minus_trace=ex.tolist(),
                 amplified_prediction_centered_on_trace=tangent.tolist(),
                 max_difference_from_centered_amp4_tangent=float(np.max(abs(ex-tangent))))
        if {"amp3","amp2","held3"}<=cases.keys():
            native_slope=np.array([(matrix(find(tier,x))-p)/(find(tier,x)["scalar_fraction"]-a["scalar_fraction"])
                         for x in ("amp4","amp3","amp2")])
            row["per_fraction_centered_slopes"]=native_slope.tolist()
            row["held_minus_direct_equal_abundance"]= (matrix(find(tier,"held3"))-matrix(find(tier,"amp3"))).tolist()
        if {"amp3","thermal3"}<=cases.keys():
            row["matched_FD_minus_direct_equal_abundance"]=(matrix(find(tier,"thermal3"))-matrix(find(tier,"amp3"))).tolist()
        trace.append(row)
    derived=[]
    if all(k in modes and {"physical","amp4","amp3","amp2","held3"}<=modes[k].keys() for k in ("std","fine")):
        for tier in ("std","fine"):
            p=find(tier,"physical");ap=[find(tier,k) for k in ("amp4","amp3","amp2")]
            slope=[matrix(x)/x["scalar_fraction"] for x in ap]
            physical=matrix(p)
            h=matrix(find(tier,"held3"))
            check=physical-slope[0]*p["scalar_fraction"]
            derived.append(dict(tier=tier,physical_log_power_ratio=physical.tolist(),
              from_amp4_tangent= (slope[0]*p["scalar_fraction"]).tolist(),
              max_physical_minus_tangent=float(np.max(abs(check))),
              amplifier_log_power_per_scalar_fraction=[s.tolist() for s in slope],
              held_minus_direct_at_equal_amp3=(h-matrix(ap[1])).tolist()))
        sz=matrix(find("fine","physical"))-matrix(find("std","physical"))
        diff=float(max(abs(sz).ravel()))
    else:diff=None
    return dict(scope="Matched post-production native CLASS linear transfers for the specified scalar PSD "
                      "and stable/CDM proxy, not evolution of the portal-production model or measured cosmology.",
                k_h_per_Mpc=zmode.tolist(),binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                cases=sample,trace_centered_audit=trace,amplification_audit=derived,physical_tier_refinement_max_abs_log_ratio_change=diff)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("root",type=Path);p.add_argument("binary",type=Path)
    p.add_argument("--partial",action="store_true")
    a=p.parse_args();print(json.dumps(run(a.root,a.binary,a.partial),indent=2,allow_nan=False))
