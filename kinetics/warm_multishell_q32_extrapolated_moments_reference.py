"""Compare admissible extrapolated moments on the corrected 12/24/48 clock."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import expit

from frw_exterior_restart_reference import import_restart
from joint_prompt_two_moment_evolution_reference import Cells, kin, rule
from three_species_decay_onsager_reference import m, mp
from warm_light_occupation_stage_audit_reference import read
from warm_multishell_q32_time_audit_reference import series
from warm_prompt_reentry_shape_audit_reference import get_mom


def run(eos, seed, folders):
    source = read(seed)[-1]["snapshot"]
    p0 = np.array(source["lab_momentum_edges_GeV"])
    mass = np.array([mp, m])
    scale = source["scale_factor"]
    old = import_restart(source["positive_exterior_restart"]).array().shape[1]

    def fetch(k, i):
        if i == 0:
            return source
        return read(folders[k]/f"step-{i*2**k}.json")["snapshot"]

    def moments(s):
        p = np.asarray(s["lab_momentum_edges_GeV"])
        a = s["scale_factor"]
        n, k, _ = get_mom(p, [s["entropy_variable_scalar"], s["entropy_variable_fermion"]])
        ext = import_restart(s["positive_exterior_restart"]).array()[:, old:]
        tagged = [ext[5,ext[1]==d].sum() for d in mass]
        return np.stack((n, k))*a**3, np.asarray(tagged)*a**3

    q = [series(eos, seed, D, 12*2**k, 12*2**k)[0] for k,D in enumerate(folders)]
    for k in (1,2):
        np.testing.assert_allclose(q[0]["T_GeV"],q[k]["T_GeV"][::2**k],rtol=1.5e-13)
    base = moments(source)[0]
    norm = base.sum(axis=-1)
    rows = []
    initial_pressure = None
    for i in range(13):
        snaps = [fetch(k,i) for k in range(3)]
        a = scale*1.10**(i/12)
        p = p0*scale/a
        X, R = [], []
        for s in snaps:
            np.testing.assert_allclose(p,s["lab_momentum_edges_GeV"],rtol=1e-13)
            u, r = moments(s)
            X.append(u)
            R.append(r)
        X, R = np.asarray(X), np.asarray(R)
        estimated = np.array([2*X[1]-X[0], 2*X[2]-X[1]])
        receipts = np.array([2*R[1]-R[0], 2*R[2]-R[1]])
        assert np.min(receipts)>=0.
        z,w = rule(192)
        t = p[:-1,None]+np.diff(p)[:,None]*(1+z)/2
        dnu = w*np.diff(p)[:,None]*t**2/(4*np.pi**2)
        pressures = []
        fits = []
        independent_defects = []
        for state in np.vstack((X,estimated)):
            errors = []
            shell_pressure = []
            independent = []
            for j, (d,g) in enumerate(zip(mass,(1,4),strict=True)):
                mesh = Cells(p,d,g)
                (A,B), defect, _ = mesh.reconstruct(state[:,j]/a**3)
                K = kin(t,d)
                s = A[:,None]+B[:,None]*(K-mesh.low[:,None])/mesh.dk[:,None]
                f = expit(-s) if g==4 else np.exp(-s)/(-np.expm1(-s))
                weight = g*dnu*f*a**3
                shell_pressure.append(np.sum(weight*t**2/(3*np.hypot(t,d)),axis=1))
                check = np.stack((weight.sum(axis=1), (weight*K).sum(axis=1)))
                independent.append(float(np.max(abs(check/state[:,j]-1.))))
                errors.append(defect)
            fits.append(errors)
            independent_defects.append(independent)
            pressures.append(shell_pressure)
        count = estimated[:,0].sum(axis=(1,2))+receipts.sum(axis=1)
        closure = count/base[0].sum()-1.
        assert np.max(abs(closure))<2e-12
        raw = abs(X[-1]-X[-2]).sum(axis=-1)/norm
        rich = abs(estimated[1]-estimated[0]).sum(axis=-1)/norm
        pp = np.asarray(pressures)
        if initial_pressure is None:
            initial_pressure = pp[0].sum(axis=-1)
        dp = abs(pp[2]-pp[1]).sum(axis=-1)/initial_pressure
        dr = abs(pp[4]-pp[3]).sum(axis=-1)/initial_pressure
        rows.append(dict(index=i,T_GeV=q[0]["T_GeV"][i],raw_NK_L1=raw.tolist(),
                         extrapolated_NK_L1=rich.tolist(),raw_recovered_P_L1=dp.tolist(),
                         extrapolated_recovered_P_L1=dr.tolist(),reconstruction_error=fits,
                         independent_order192_moment_error=independent_defects,
                         extrapolated_inclusive_count_defect=closure.tolist(),
                         extrapolated_retained_count=(estimated[:,0].sum(axis=-1)/base[0].sum(axis=-1)).tolist(),
                         extrapolated_new_receipts=(receipts/base[0].sum(axis=-1)).tolist(),
                         cold_cell_occupations=(np.vstack((X,estimated))[:,0,1,0]/(a**3*4*p[1]**3/(6*np.pi**2))).tolist()))
    return dict(scope="Extrapolate only comoving number and kinetic energy; reconstruct admissible phases "
                     "and recompute pressure. Exterior receipts are aggregate counts, not phase. "
                     "No extrapolant is an evolved state or a rigorous error bound.", frames=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for x in ("eos","seed","coarse","medium","fine"):
        p.add_argument(x,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.seed,(a.coarse,a.medium,a.fine)),indent=2,allow_nan=False))
