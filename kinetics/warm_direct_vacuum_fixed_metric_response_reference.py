"""Derive the fixed-metric scalar response directly from positive shell counts."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.integrate import cumulative_simpson,simpson
from scipy.special import roots_legendre,spherical_jn
from archived_collisionless_streaming_reference import Background
from archived_scalar_pair_reaction_reference import points
from finite_cooling_relativistic_los_reference import trajectory
from joint_visible_heavy_adiabatic_reference import snapshot_functions
from three_species_decay_onsager_reference import m as light_mass,mp
from warm_prompt_reentry_shape_audit_reference import get_mom
from warm_light_occupation_stage_audit_reference import digest,read


def kernel(bg,metric,q,nt,*,free=False):
    k,I,tf,native,M,_=metric
    theta=np.array([I[4]/4.,I[5]/(3*k),I[6]/2])
    tau=np.geomspace(I[0],tf,nt)
    a,psi,phi=PchipInterpolator(M[:,0],M[:,1:4],axis=0)(tau).T
    am=a*(mp/bg.T0)
    ep=np.hypot(q,am[-1])
    F=[];D=[];IF=[];ID=[]
    for start in range(0,len(q),18):
        z=q[start:start+18,None]
        v=z/np.hypot(z,am)
        w=v*(1-v*v)
        def back(x):
            return cumulative_simpson(x[:,::-1],x=-tau[::-1],axis=1,initial=0.)[:,::-1]
        A=k*back(v); dA=k*back(w)
        j0=spherical_jn(0,A[:,0]);j1=spherical_jn(1,A[:,0]);j2=spherical_jn(2,A[:,0])
        G=(-theta[0]+phi[0])*j0+3*theta[1]*j1-5*theta[2]*j2-phi[-1]
        dG=((-theta[0]+phi[0])*(-j1)+3*theta[1]*spherical_jn(1,A[:,0],derivative=True)
            -5*theta[2]*spherical_jn(2,A[:,0],derivative=True))*dA[:,0]
        form=v*phi+psi/v
        dform=w*(phi-psi/v**2)
        J=spherical_jn(1,A)
        dq=simpson(k*(spherical_jn(1,A,derivative=True)*dA*form+J*dform),x=tau,axis=1)
        if free:
            U=-theta[0]*j0+3*theta[1]*j1-5*theta[2]*j2
            dU=(theta[0]*j1+3*theta[1]*spherical_jn(1,A[:,0],derivative=True)-
                     5*theta[2]*spherical_jn(2,A[:,0],derivative=True))*dA[:,0]
            IF.append(U);ID.append(dU)
        F.append(G+simpson(k*J*form,x=tau,axis=1))
        D.append(dG+dq)
    F=np.concatenate(F);D=np.concatenate(D)
    # Numerator integral q^3 eps f'(q) F dq. Move derivative onto
    # the smooth metric kernel; this avoids differentiating shell noise.
    R=-(3*F+(q/ep)**2*F+D)/native[0]
    if free:
        base=-(3*np.concatenate(IF)+(q/ep)**2*np.concatenate(IF)+np.concatenate(ID))/native[0]
        return ep,F,D,R,base
    return ep,F,D,R


def layout(edges,r,shape):
    x,w=roots_legendre(r)
    w=w/2
    lo,hi=edges[:-1,None],edges[1:,None]
    if shape=="q":
        q=lo+(hi-lo)*(x+1)/2
    else:
        safe=np.maximum(lo,1e-33)
        q=safe*np.exp(np.log(hi/safe)*(x+1)/2)
        q[0]=hi[0]*(x+1)/2
    return q.ravel(),w


def run(eos,archive,psd,metrics):
    bg=Background(eos)
    z=read(archive)
    S=z["snapshot"]
    assert z["snapshot_sha256"]==digest(S)
    source=read(psd)
    ed=np.array([x if x is not None else 36. for x in source["conserved_momentum_edges_q_over_T0"]])
    u=np.array(source["controls"][-1]["normalized_shell_counts"])
    assert np.max(u[:,-1])<2e-9
    occupation,*_=snapshot_functions(S)
    p=np.asarray(S["lab_momentum_edges_GeV"])
    ah=bg.scale(S["T_visible_GeV"])
    N=get_mom(p,[S["entropy_variable_"+k] for k in ("scalar","fermion")])[0].sum(axis=1)
    T=S["T_visible_GeV"]
    visible=2*np.pi**2/45*T**3*np.exp(bg.gs(np.log(T)))
    fractions=[float(N[0]*mp/(N[0]*mp+N[1]*light_mass))]
    mix_y=np.asarray(source["controls"][-1]["conditional_total_scalar_yield"])
    new=visible*mix_y.sum()*mp
    fractions.append(float(new/(new+N[1]*light_mass)))
    results=[]
    for raw in metrics:
        m=trajectory(raw)
        assert m[1][1]>bg.scale(.012)>bg.scale(.05) and m[1][1]>ah
        tests=[]
        for nt,n,spacing in ((7000,8,"q"),(13000,12,"q"),(13000,12,"log_q")):
            q,w=layout(ed,n,spacing)
            p0,_,wH=points(p,mp,occupation,n)
            qH=p0*ah/bg.T0
            cq=np.r_[q,qH]
            Ep,F,_,K,Kfree=kernel(bg,m,cq,nt,free=True)
            numer=(Ep[:len(q)]*K[:len(q)]).reshape(len(ed)-1,n)@w
            denom=Ep[:len(q)].reshape(len(ed)-1,n)@w
            vals=(u@numer)/(u@denom)
            held=float(wH@(Ep[len(q):]*K[len(q):])/(wH@Ep[len(q):]))
            mix=float((mix_y*(u@numer)).sum()/(mix_y*(u@denom)).sum())
            fi=(Ep[:len(q)]*Kfree[:len(q)]).reshape(len(ed)-1,n)@w
            fmix=float((mix_y*(u@fi)).sum()/(mix_y*(u@denom)).sum())
            hfree=float(wH@(Ep[len(q):]*Kfree[len(q):])/(wH@Ep[len(q):]))
            tests.append(dict(time_nodes=nt,q_order=n,within_source_shell_measure=spacing,
                response_over_native_CDM=[*vals.tolist(),mix,held],
                initial_only_direct_mix_and_held_over_native_CDM=[fmix,hfree],
                prescribed_metric_component_of_mix_minus_held=float((mix-held)-(fmix-hfree)),
                total_amplitude_replacement_effect_if_other_matter_held_CDM_and_equal_scalar_rest_count=fractions[0]*(mix-held),
                total_amplitude_effect_if_swap_actual_direct_rest_yield_and_other_matter_is_CDM=
                    fractions[1]*(mix-1)-fractions[0]*(held-1)))
            if nt==13000 and spacing=="q":
                et,ff,_,_=kernel(bg,m,np.maximum(ed,1e-22),nt)
                jump=-np.diff(ed*et*ff)/np.diff(ed)
                within=-2*(Ep[:len(q)]*F[:len(q)]).reshape(len(ed)-1,n)@w
                v=(u@(jump+within))/(u@denom)/m[3][0]
                assert np.max(abs(v-vals))<5e-8
                tests[-1]["direct_jump_vs_kernel_IBP_max_difference"]=float(max(abs(v-vals)))
        # Compare analytic log-q derivative with independent finite differences.
        P=np.geomspace(.025,8.,10)
        eps=2e-5
        _,f,df,_=kernel(bg,m,P,12000)
        _,l,_,_=kernel(bg,m,P*np.exp(-eps),12000)
        _,r,_,_=kernel(bg,m,P*np.exp(eps),12000)
        err=float(max(abs((r-l)/(2*eps)-df))/max(max(abs(df)),1e-50))
        assert err<2e-5
        diff=max(abs(np.array(tests[0]["response_over_native_CDM"])-tests[1]["response_over_native_CDM"]))
        shape=max(abs(np.array(tests[1]["response_over_native_CDM"])-tests[2]["response_over_native_CDM"]))
        assert diff<.001 and shape<.005,(m[0],diff,shape)
        results.append(dict(k_per_Mpc=m[0],metric_sha256=hashlib.sha256(raw.read_bytes()).hexdigest(),
            driver_initial_scale=float(m[1][1]),analytic_log_q_kernel_gradient_relative_defect=err,
            maximum_time_momentum_refinement_shift=float(diff),maximum_between_shell_measure_shift=float(shape),controls=tests))
    return dict(scope="Post-production scalar tracer response on recorded Newtonian-gauge metrics with shared "
                      "approximate adiabatic initial monopole/dipole/quadrupole. No injection perturbations, "
                      "backreaction, inhomogeneous production or transfer prediction.",
                held_and_actual_direct_scalar_rest_fractions_with_same_archived_stable_count=fractions,
                z=20,source_T_selected_GeV=[.012,.05],response_columns=["direct leptons","measured pion window",
                       "yield-weighted direct mix","q64 held warm reference"],archive=z["snapshot_sha256"],rows=results)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","archive","psd"):
        p.add_argument(k,type=Path)
    p.add_argument("metrics",nargs="+",type=Path)
    a=p.parse_args();print(json.dumps(run(a.eos,a.archive,a.psd,a.metrics),indent=2,allow_nan=False))
