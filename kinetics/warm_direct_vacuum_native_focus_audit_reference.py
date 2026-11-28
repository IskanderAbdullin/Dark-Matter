"""Compare matching native exact-k modes and a separate native RK hierarchy."""
import argparse,hashlib,json,re
from pathlib import Path
import numpy as np
from warm_direct_vacuum_native_focus_reference import K
from warm_light_occupation_stage_audit_reference import read


def run(d,std):
    seed=read(d/"spec.json")
    f={t["tag"]:t["omega_ncdm"]/.12 for t in seed["rows"]}
    old=read(std)
    pilot=next(x for x in old["trace_centered_audit"] if x["tier"]=="std")
    stdk=np.asarray(old["k_h_per_Mpc"])
    modes=np.asarray(K)
    i=np.searchsorted(stdk,modes)
    assert np.array_equal(stdk[i],modes)
    first=np.asarray(pilot["physical_minus_trace"])[:,i]
    case=[]
    for j,z in enumerate((5,20),1):
        spectra={};hashes={}
        for n in ("base","trace","physical","amp4"):
            source=d/f"{n}_00_z{j}_pk.dat"
            a=np.loadtxt(source)
            assert f"redshift z={z}" in source.read_text().splitlines()[0]
            sought=np.searchsorted(a[:,0],modes)
            assert np.array_equal(a[sought,0],modes)
            spectra[n]=a[sought,1]
            assert np.all(spectra[n]>0)
            hashes[n]=hashlib.sha256(source.read_bytes()).hexdigest()
            l=(d/(n+".log")).read_text()
            assert "Error" not in l and "switch on ncdm fluid" not in l
            # Verify exact requested fractional composition; budget
            # printouts have fewer digits than the input specification.
            if n!="base":
                N=float(re.search(r"Non-Cold Species Nr\..+?omega = ([^\s]+)",l)[1])
                assert abs(N/(f[n]*.12)-1)<6e-6
        p=np.log(spectra["physical"]/spectra["trace"])
        t=np.log(spectra["trace"]/spectra["base"])
        s=np.log(spectra["amp4"]/spectra["trace"])*(f["physical"]-f["trace"])/(f["amp4"]-f["trace"])
        case.append(dict(z=z,raw_table_sha256=hashes,exact_k_physical_minus_trace=p.tolist(),
            exact_k_trace_presence_offset=t.tolist(),amp4_centered_prediction=s.tolist(),
            max_abs_native_centered_tangent_error=float(max(abs(p-s))),
            difference_from_standard_interpolated_mode=(p-first[j-1]).tolist()))
    starts=[]
    for n in ("base","trace","physical","amp4"):
        vals=[]
        for l,k in enumerate(modes):
            path=d/f"{n}_00_perturbations_k{l}_s.dat"
            with path.open() as handle:
                h=next(handle);next(handle)
                row=[float(a) for a in next(handle).split()[:2]]
            k0=float(re.search(r"mode k = ([^ ]+)",h)[1])
            assert abs(k0/(k*.6736)-1)<1e-12
            vals.append(row)
        starts.append(dict(case=n,first_tau_Mpc_and_scale_at_exact_modes=vals))
    return dict(initial_times=starts,scope="Narrowed alternate CLASS RK native hierarchy at exact requested output modes, compared "
                       "with the previous NDF pilot. Same conditional scalar composition and PSD. "
                       "Not complete high-k or cosmological convergence.",
                native_k_h_per_Mpc=list(K),cases=case)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("native",type=Path);p.add_argument("std",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.native,a.std),indent=2,allow_nan=False))
