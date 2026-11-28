"""Propagate a prescribed collision-only common-bath tangent in fixed coordinates."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from warm_common_pair_thermal_relaxation_reference import modes
from warm_light_occupation_stage_audit_reference import digest, read


def run(eos,tapes,wide):
    J=read(tapes)["cases"]
    indices=(0,6,12,18,24)
    assert len(J)==len(indices)
    raw=[read(wide/"guide-start.json"),*[read(wide/"q64"/f"step-{i}.json") for i in indices[1:]]]
    anchor=raw[0]["snapshot"]["scale_factor"]
    matrix=[]
    loc=[]
    piv=.0075
    for i,R,B in zip(indices,raw,J,strict=True):
        S=R["snapshot"]
        assert digest(S)==R["snapshot_sha256"]==B["snapshot_sha256"] and S["T_visible_GeV"]==B["T_GeV"]
        assert abs(S["scale_factor"]/anchor-1.10**(i/24))<2e-13
        fit=B["controls"][-1]
        Q=fit["equilibrium_FD_BE_response"]
        thermal=modes(Q,eos,18)
        T,mu=Q["T_GeV"],Q["mu_GeV"]
        physics=np.array([[1.,0.],[-mu/T,piv/T]])
        susceptibility=sum(np.linalg.inv(thermal[k]) for k in ("scalar_susceptibility","fermion_susceptibility"))
        F=np.array(Q["extrapolated"]["richardson_response"])*Q["equilibrium_forward_event_number_energy"][0]
        F=(F+F.T)/2
        norm=np.array(fit["synthetic_pair_tangent_rates_over_actual_anchor_H"])
        H=thermal["decay_rates_GeV"][0]/norm[0]
        A=np.linalg.solve(physics,F@susceptibility@physics)/H
        np.testing.assert_allclose(np.sort(np.linalg.eigvals(A)),norm,rtol=2e-12)
        matrix.append(A)
        loc.append(dict(index=i,snapshot_sha256=B["snapshot_sha256"],T_visible_GeV=B["T_GeV"],
                        relative_lna=np.log(S["scale_factor"]/anchor),
                        fit_T_GeV=T,fit_mu_GeV=mu,collision_decay_rates_over_H=norm.tolist(),
                        independent_relative_eigenvalue_defects=thermal["independent_relative_eigenvalue_defects"]))
    X=np.array([r["relative_lna"] for r in loc])

    def propagate(n):
        U=np.eye(2)
        for i,dt in enumerate(np.diff(X)):
            for j in range(n):
                q=(j+.5)/n
                A=(1-q)*matrix[i]+q*matrix[i+1]
                U=expm(-A*(dt/n)) @ U
        return U
    ideal_trace=0.
    for a,b,dx in zip(matrix[:-1],matrix[1:],np.diff(X),strict=True):
        ideal_trace+=.5*np.trace(a+b)*dx
    report=[]
    for n in (1,4,16,64):
        P=propagate(n)
        sign,det=np.linalg.slogdet(P)
        assert sign==1 and abs(det+ideal_trace)<2e-12
        eig=np.sort(np.linalg.eigvals(P))
        assert np.isreal(eig).all() and np.min(eig)>0.
        report.append(dict(midpoint_exponentials_per_anchor_segment=n, propagator=P.tolist(),
                           eigen_amplitude_factors=eig.tolist(), trace_logdet_defect=det+ideal_trace,
                           singular_amplitude_factors=np.linalg.svd(P,compute_uv=False).tolist()))
    return dict(scope="Piecewise-linear interpolation of fitted thermal collision matrices, in fixed coordinates "
                      "dN and dE divided by a constant pivot, with scalar+fermion summed invariants frozen for "
                      "the perturbation. Retained pair only; no expansion pressure work, velocity-space distortion, "
                      "external source, or changing equilibrium driving in this prescribed tangent.",
                energy_pivot_GeV=piv,local_fits=loc,propagator_checks=report)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for x in ("eos","tangent_tape","wide"):
        p.add_argument(x,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.tangent_tape,a.wide),indent=2,allow_nan=False))
