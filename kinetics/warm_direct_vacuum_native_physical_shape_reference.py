"""Measure the equal-physical-fraction shape change with two trace controls."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from warm_direct_vacuum_native_focus_reference import K
from warm_light_occupation_stage_audit_reference import read


def run(files,amplified):
    ref=read(amplified)
    old=np.asarray(ref["held_minus_direct_at_equal_fraction_sig07_and_sig11"][0])
    f=next(x["omega_ncdm"]/.12 for x in read(files/"spec.json")["rows"] if x["tag"]=="physical")
    case=[]
    for j,z in enumerate((5,20),1):
        vals=[];starts=[];sha={}
        for label in ("trace","physical","heldtrace","heldphysical"):
            p=files/f"{label}_00_z{j}_pk.dat"
            data=np.loadtxt(p)
            i=np.searchsorted(data[:,0],K)
            assert np.array_equal(data[i,0],K)
            sha[label]=hashlib.sha256(p.read_bytes()).hexdigest()
            vals.append(data[i,1])
            at=[]
            for i in range(len(K)):
                file=files/f"{label}_00_perturbations_k{i}_s.dat"
                with file.open() as handle:
                    next(handle);next(handle)
                    at.append(float(next(handle).split()[1]))
            starts.append(at)
        assert np.max(abs(np.array(starts)/np.asarray(starts)[0]-1))<1e-6
        t,p,T,P=vals
        base=np.log(T/t)
        ds=np.log(P/p)-base
        verify=np.log(P/T)-np.log(p/t)
        np.testing.assert_allclose(ds,verify,atol=5e-16,rtol=0)
        primary=np.log(p/t)
        assert ds[-1]<0 and abs(base[-1])<.01*abs(ds[-1])
        hint=float(old[j-1,-2]*(f-1e-10)/(1e-3-1e-10))
        case.append(dict(z=z,sha256=sha,initial_scale_factors=starts,
           held_minus_direct_trace=base.tolist(),held_minus_direct_raw_at_physical_fraction=np.log(P/p).tolist(),
           trace_centered_held_minus_direct=ds.tolist(),direct_minus_own_trace=primary.tolist(),
           absolute_shape_change_as_share_of_direct_effect_at_k1=float(ds[-1]/primary[-1]),
           amplified_standard_tangent_at_k1=hint))
    return dict(scope="Equal physical fraction of the scalar in two distinct PSDs, each centered on its own "
                      "near-zero trace. Early exact-k RK; not the difference of two physical production models.",
                k_h_per_Mpc=list(K),cases=case)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("native",type=Path);p.add_argument("amplified",type=Path)
    a=p.parse_args();print(json.dumps(run(a.native,a.amplified),indent=2,allow_nan=False))
