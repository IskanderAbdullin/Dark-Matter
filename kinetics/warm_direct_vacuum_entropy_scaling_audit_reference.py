"""Check the kinematic meaning of algebraic dilution when applied as entropy release."""
import argparse,json
from pathlib import Path
import numpy as np
from archived_collisionless_streaming_reference import Background
from three_species_decay_onsager_reference import mp
from warm_light_occupation_stage_audit_reference import read


def run(eos,psd,N):
    bg=Background(eos)
    p=read(psd)
    a=read(N)
    d=a["dilution_factoring_old_fraction_to_reference_total"]
    g=d**(-1/3)
    q=np.array([float(x) for x in p["conserved_momentum_edges_q_over_T0"][:-1]])
    W=np.asarray(p["controls"][-1]["shell_yields"]).sum(axis=0)
    w=W/W.sum()
    # Compare after the archive in the same prescribed late background;
    # a putative entropy event must already have finished for this test.
    ai=bg.scale(.0073940494651976)
    cl=min(bg.or0,bg.ogamma/2*np.exp(np.min(bg.gr(bg.u))+4*(np.log(bg.gs0)-np.max(bg.gs(bg.u)))/3))
    limit=bg.D0*(bg.aeq-ai)/np.sqrt(cl)
    rows=[]
    for eta in (1.,g):
        x=np.maximum(q,1e-30)*eta*bg.T0
        D=bg.streams(x,mp,ai,bg.aeq,6)
        D[0]=0.
        E=bg.streams(x,mp,ai,bg.aeq,9);E[0]=0.
        np.testing.assert_allclose(D,E,atol=4e-10,rtol=5e-11)
        assert np.all(np.diff(D)>0) and limit>D[-1]
        rows.append(dict(momentum_scaling=eta,mean_post_archive_path_interval_Mpc=[float(w@D),float(w@np.r_[D[1:],limit])]))
    assert rows[1]["mean_post_archive_path_interval_Mpc"][1]<rows[0]["mean_post_archive_path_interval_Mpc"][0]
    return dict(scope="Only kinematics after a matched time, not a physical entropy source. A hypothetical uniform "
                      "visible entropy increase cools decoupled momenta relative to photons even at the same diluted "
                      "rest fraction; keeping the original spectrum is a distinct mathematical normalization.",
                dilution_factor=d,momentum_ratio=g,instantaneous_fixed_degrees_radiation_energy_ratio=d**(4/3),rows=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","psd","normalization"):
        p.add_argument(k,type=Path)
    a=p.parse_args();print(json.dumps(run(a.eos,a.psd,a.normalization),indent=2,allow_nan=False))
