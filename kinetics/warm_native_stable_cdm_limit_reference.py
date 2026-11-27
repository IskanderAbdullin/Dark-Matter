"""Compare the archived massive stable fermion with its exactly cold limit on a fixed metric."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.integrate import cumulative_simpson,simpson
from scipy.special import spherical_jn
from archived_collisionless_streaming_reference import Background
from archived_scalar_pair_reaction_reference import points
from finite_cooling_relativistic_los_reference import trajectory
from joint_visible_heavy_adiabatic_reference import snapshot_functions
from three_species_decay_onsager_reference import m
from warm_light_occupation_stage_audit_reference import digest,read


def driven(metric,bg,q,weights,n,initial):
    k,b,t1,native,tab,_=metric
    time=np.geomspace(b[0],t1,n)
    a,psi,phi=PchipInterpolator(tab[:,0],tab[:,1:4],axis=0)(time).T
    T=bg.T0
    def B(x):
        return cumulative_simpson(x[...,::-1],x=-time[::-1],axis=-1,initial=0.)[...,::-1]
    # Same unit initial density contrast and zero initial velocity for
    # both comparisons. A dilational perturbation has F_ini=-1/3 for a
    # cold species (dln f/dln q integrates to -3).
    source_e=np.hypot(q,a[0]*m/T)
    a_init=-initial/(3+((q/source_e)**2)@(weights*source_e)/(weights@source_e))
    zeroF=-initial/3-phi[-1]+phi[0]+k*k/3*simpson(psi*a*B(1/a),x=time)
    Cold=-3*zeroF
    e=np.hypot(q,a[-1]*m/T)
    V=[]
    for j in range(0,len(q),24):
        p=q[j:j+24,None]
        v=p/np.hypot(p,a*m/T)
        u=v*(1-v*v)
        arg=k*B(v)
        darg=k*B(u)
        s0=spherical_jn(0,arg[:,0])
        s1=spherical_jn(1,arg[:,0])
        J=spherical_jn(1,arg)
        Jd=spherical_jn(1,arg,derivative=True)
        g=v*phi+psi/v
        dg=u*(phi-psi/v**2)
        F=(a_init+phi[0])*s0-phi[-1]+simpson(k*J*g,x=time,axis=-1)
        dF=(a_init+phi[0])*(-s1)*darg[:,0]+simpson(k*(Jd*darg*g+J*dg),x=time,axis=-1)
        R=-(3+(p[:,0]/e[j:j+24])**2)*F-dF
        V.extend(R.tolist())
    return e,np.asarray(V),float(Cold),a_init


def run(eos,archive,drivers,cdm=False):
    b=Background(eos)
    z=read(archive);S=z["snapshot"]
    assert z["snapshot_sha256"]==digest(S)
    ah=b.scale(S["T_visible_GeV"])
    e=np.asarray(S["lab_momentum_edges_GeV"])
    _,_,f,_=snapshot_functions(S)
    results=[]
    for path in drivers:
        raw=trajectory(path,cdm_only=cdm)
        assert raw[1][1]>ah
        alt=[]
        for nq,nt in ((12,6500),(24,12000)):
            p,_,w=points(e,m,f,nq)
            w/=w.sum()
            q=p*ah/b.T0
            D=b.streams(q*b.T0,m,ah,b.aeq,7)
            # Exact spherical-j0 bound for initial-only equal number
            # perturbation; not a bound on the forced response.
            ball=(raw[0]**2/6)*float(w@(D*D))
            ep,R,C,a_init=driven(raw,b,q,w,nt,1.)
            direct=float((w*ep)@(R-C)/((w@ep)*C))
            r=float((w*ep)@R/(w@ep))
            assert abs((r/C-1)-direct)<5e-15
            alt.append(dict(nodes_per_shell=nq,time_nodes=nt,fermion_response_in_initial_amplitude_units=r,
                cold_limit_same_fixture_and_initial_data=C,fermion_initial_derivative_coefficient_for_unit_energy_contrast=a_init,
                fractional_difference_from_cold_limit=direct,
                mean_unforced_path_to_equality_Mpc=float(w@D),
                rms_unforced_path_Mpc=float(np.sqrt(w@(D*D))),
                upper_ballistic_monopole_departure_at_fixtures_k_and_equality=ball))
        residual=abs(alt[0]["fractional_difference_from_cold_limit"]-alt[1]["fractional_difference_from_cold_limit"])
        assert residual<3e-10
        ep,R,C,_=driven(raw,b,q,w,12000,0.)
        normcheck=float((w*ep)@(R-C)/((w@ep)*C))
        assert abs(normcheck-alt[-1]["fractional_difference_from_cold_limit"])<1e-9
        results.append(dict(k_per_Mpc=raw[0],metric_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                            convergence_difference=residual,forcing_only_fractional_difference=normcheck,controls=alt))
    return dict(scope="Conditional stable-fermion stress on an identical fixed metric. Unit initial density and "
                      "zero peculiar velocity for both tracers; no independent physical adiabatic production, "
                      "native self-gravity, or correction to the full CLASS matter power.",
                archive=z["snapshot_sha256"],z=20,rows=results)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("eos",type=Path);p.add_argument("archive",type=Path)
    p.add_argument("metrics",nargs="+",type=Path)
    p.add_argument('--cdm-driver',action='store_true')
    a=p.parse_args();print(json.dumps(run(a.eos,a.archive,a.metrics,a.cdm_driver),indent=2,allow_nan=False))
