"""Scalar Yukawa chi-antichi <-> phi-phi on fixed two-moment archives."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from joint_visible_heavy_adiabatic_reference import snapshot_functions
from sm_dark_entropy_exchange_reference import MPL, StandardModel
from three_species_decay_onsager_reference import M, m, mp, y


def amplitude(s,x,light=m,heavy=M,scalar=mp):
    """Spin average over the two Dirac legs, y^4 included."""
    p=np.sqrt(np.maximum(s/4-light**2,0.))
    k=np.sqrt(np.maximum(s/4-scalar**2,0.))
    a=light*light+scalar*scalar-s/2-heavy*heavy
    d=2*p*k*x
    u=a*a-d*d
    # Equivalent t/u sum without direct cancellation at threshold.
    az=(2*p*a+2*k*x*d)/u
    ap=(2*k*d)*np.sqrt(np.maximum(1.-x*x,0.))/u
    b=2*heavy*a/u
    return 2*y**4*((s/4)*ap*ap+(light*az+p*b)**2)


def trace_check():
    z=np.zeros((2,2))
    I=np.eye(2)
    pauli=(np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),
           np.diag([1,-1]))
    ga=[np.block([[I,z],[z,-I]])]+[np.block([[z,v],[-v,z]]) for v in pauli]
    slash=lambda V: V[0]*ga[0]-sum(v*a for v,a in zip(V[1:],ga[1:]))
    one=np.eye(4)
    defects=[]
    for s in np.geomspace(4*m*m*1.0001,16.,10):
        p=np.sqrt(s/4-m*m);k=np.sqrt(s/4-mp*mp);E=np.sqrt(s)/2
        for x in (-.989,-.51,-.04,0.,.37,.92):
            p1=np.array([E,0.,0.,p]);p2=np.array([E,0.,0.,-p])
            k1=np.array([E,k*np.sqrt(1.-x*x),0.,k*x]);k2=2*np.array([E,0.,0.,0.])-k1
            def term(q):
                return (slash(q)+M*one)/(-q[1:]@q[1:]-M*M)
            A=term(p1-k1)+term(p1-k2)
            a=np.trace((slash(p2)-m*one)@A@(slash(p1)+m*one)@A).real*y**4/4
            b=amplitude(s,x)
            defects.append(abs(a/b-1.))
    assert max(defects)<2e-11
    q=np.geomspace(.0001,.001,25)
    # Verify the leading power independently of its coupling coefficient.
    sc=4*(m*m+q*q)
    vx=amplitude(sc,0.)/sc
    assert abs(vx[1]/vx[0]/(q[1]/q[0])**2-1.)<3e-5
    # Equal mediator and incoming mass, massless scalar: independent
    # leading expansion of the angular polynomial gives coefficient 3/128.
    p=m*1e-4;s=4*(m*m+p*p)
    w,x=roots_legendre(48)
    vrel=4*p/np.sqrt(s)
    sigma_v=np.sum(x*amplitude(s,w,heavy=m,scalar=0.))/(32*np.pi*s)
    limit=3*y**4*vrel*vrel/(128*np.pi*m*m)
    diagonal=abs(sigma_v/limit-1.)
    assert diagonal<3e-7
    return dict(max_spin_trace_relative_defect=max(defects),degenerate_internal_mass_pwave_relative_defect=diagonal)


def points(edges,mass,f,n):
    x,w=roots_legendre(n)
    a,b=edges[:-1,None],edges[1:,None]
    p=((b+a)+(b-a)*x)/2
    E=np.hypot(p,mass)
    d=(w*(b-a)/2)*p*p*f(E)/(2*np.pi**2)
    return p.ravel(),E.ravel(),d.ravel()


def gross(edges,inc,out,occupancy,final,reverse,radial=4,polar=9,azimuth=8,panel=3,split_final=True,project=False):
    """Ordered unit-spin incoming measures times pair-counting factors.

    Opposite-charge fermions distinguishable: 2 spin states each.
    Initial scalars identical: 1/2. Outgoing identical scalar factor
    lives in the 1/(32 pi s) sigma*v prefactor, not in this measure.
    """
    p,e,n=points(edges,inc,occupancy,radial)
    i,j=np.tril_indices(len(p))
    i,j=i[n[i]*n[j]>0.],j[n[i]*n[j]>0.]
    x,w=roots_legendre(polar)
    az=np.cos((np.arange(azimuth)+.5)*np.pi/azimuth)
    xf,wf=roots_legendre(panel)
    event=np.zeros(6)
    shell=np.zeros((2,2,len(edges)+1)) if project else None
    def record(dst, E, weights, energy_edges, mass):
        idx=np.searchsorted(energy_edges,E.ravel(),side="right")
        pure=weights.ravel()
        dst[0]+=np.bincount(idx,pure,minlength=len(edges)+1)
        dst[1]+=np.bincount(idx,pure*np.maximum(E.ravel()-mass,0.),minlength=len(edges)+1)
    in_e=np.hypot(edges,inc)
    out_e=np.hypot(edges,out) if split_final else np.array([])
    for start in range(0,len(i),72):
        I,J=i[start:start+72],j[start:start+72]
        E,F=e[I,None],e[J,None]
        prod=p[I,None]*p[J,None]
        cap=np.clip((E*F+inc*inc-2*out*out)/prod,-1.,1.)
        if not reverse:
            assert np.all(cap==1.)
        cos=(cap-1.)/2+(cap+1.)/2*x
        wind=(cap+1.)/2*w
        st=(2*inc**2+2*(E*F-prod*cos)).ravel()
        mask=(wind.ravel()>0.)&(st>4*out*out)&(st>4*inc*inc)
        if not np.any(mask):
            continue
        def flatten(a):
            return np.broadcast_to(a,wind.shape).ravel()[mask,None]
        total=flatten(E+F)
        s=st[mask,None]
        ci=np.sqrt(s/4-inc*inc)
        cf=np.sqrt(s/4-out*out)
        P=np.sqrt(np.maximum(total*total-s,0.))
        boost=P/np.sqrt(s)
        ci_cos=np.clip((flatten(E-F)/2)/np.maximum(boost*ci,1e-300),-1.,1.)
        cf_delta=boost*cf
        center=total/2
        lower=np.maximum((out_e-center)/np.maximum(cf_delta,1e-300),-1.)
        upper=np.minimum((center-out_e)/np.maximum(cf_delta,1e-300),1.)
        # Split outgoing lab-angle at both exact retained shell edges.
        cuts=np.sort(np.clip(np.concatenate([-np.ones_like(total),lower,upper,np.ones_like(total)],axis=1),-1.,1.),axis=1)
        d=np.diff(cuts,axis=1)/2
        z=(cuts[:,1:]+cuts[:,:-1])/2
        v=(z[:,:,None]+d[:,:,None]*xf).reshape(len(s),-1)
        wa=(d[:,:,None]*wf).reshape(len(s),-1)
        vcm=ci_cos*v
        vperp=np.sqrt(np.maximum(1.-ci_cos*ci_cos,0.))*np.sqrt(np.maximum(1.-v*v,0.))
        def scale(c):return np.clip(vcm[:,:,None]+vperp[:,:,None]*c,-1.,1.)
        # Dimensionless averages over outgoing azimuth, then dc_f.
        ma=amplitude(s[:,:,None],scale(az)).mean(axis=2)
        ef0=center+cf_delta*v
        ef1=center-cf_delta*v
        f0=final(ef0)
        f1=final(ef1)
        ff=(1-f0)*(1-f1) if reverse else (1+f0)*(1+f1)
        assert np.min(ff)>=0.
        av=np.sum(wa*ma,axis=1)
        a=wa*ma*ff
        aq=np.sum(a,axis=1)
        # The inverse has four final-spin states and no identical-final
        # division. Number of events (ordered incoming angular mean 1/2).
        weight=(wind.ravel()[mask]/2)*flatten((n[I]*n[J]*np.where(I==J,1.,2.))[:,None])[:,0]
        sig=weight*cf[:,0]/np.sqrt(s[:,0])/(64*np.pi*flatten(E*F)[:,0])
        # Forward: four counts of incoming spins. Inverse: times
        # eight for amplitude/final symmetry, then half for initial symmetry.
        sig *= 4.
        low=np.hypot(edges[0],out)
        high=np.hypot(edges[-1],out)
        n0=((ef0>=low)&(ef0<high))
        n1=((ef1>=low)&(ef1<high))
        kept_n=np.sum(a*(n0.astype(float)+n1)/2,axis=1)
        kept_energy=np.sum(a*(n0*ef0+n1*ef1),axis=1)
        min_scalar=np.sqrt(np.maximum(np.minimum(ef0,ef1)**2-out*out,0.))
        small=np.sum(a*(min_scalar<.00025)/2,axis=1)
        event += (sig@av,sig@aq,(sig*total[:,0])@aq,sig@kept_n,sig@kept_energy,sig@small)
        if project:
            for e_in in [flatten(E),flatten(F)]:
                record(shell[0],e_in,sig*aq,in_e,inc)
            rates=sig[:,None]*a
            record(shell[1],ef0,rates,out_e,out)
            record(shell[1],ef1,rates,out_e,out)
    if project and event[1]>0.:
        counts=shell[:,0].sum(axis=1)
        energy=shell[:,1].sum(axis=1)+np.array([inc,out])*counts
        assert np.max(abs(counts/(2*event[1])-1.))<2e-12
        assert np.max(abs(energy/event[2]-1.))<2e-12
    return (n.sum()*(1 if reverse else 4),event,shell)


def thermal_check():
    from scipy.special import expit
    T=.018
    mu=-.08
    F=lambda e:expit((mu-e)/T)
    B=lambda e:1./np.expm1((e-mu)/T)
    p=np.r_[0.,np.geomspace(.00015,1.1,34)]
    rows=[]
    for r,t,z,a in ((4,12,12,9),(6,18,16,12)):
        _,v,_=gross(p,m,mp,F,B,False,r,t,z,a,split_final=False)
        _,w,_=gross(p,mp,m,B,F,True,r,t,z,a,split_final=False)
        rows.append(dict(order=[r,t,z,a],inverse_over_forward_number=w[1]/v[1],
                         inverse_over_forward_energy=w[2]/v[2]))
    assert max(abs(rows[-1][k]-1.) for k in ("inverse_over_forward_number",
        "inverse_over_forward_energy"))<1e-5
    return rows


def cropped(snapshot,cut):
    p=np.asarray(snapshot["lab_momentum_edges_GeV"])
    cap=cut*snapshot["T_initial_GeV"]/snapshot["scale_factor"]
    k=np.argmin(abs(p/cap-1.))
    assert 0<k<len(p)-1 and abs(p[k]/cap-1.)<1e-12
    other=dict(snapshot)
    other["lab_momentum_edges_GeV"]=p[:k+1]
    for key in ("entropy_variable_scalar","entropy_variable_fermion"):
        other[key]=np.asarray(snapshot[key])[:,:k]
    fp,_,fc,_=snapshot_functions(other)
    other["represented_species_count_a3"]=(np.array([points(p[:k+1],mp,fp,24)[2].sum(),
        4*points(p[:k+1],m,fc,24)[2].sum()])*snapshot["scale_factor"]**3).tolist()
    return other


def run(eos,files,orders,crop_last=None):
    result=[]
    cases=[(Path(path).name,json.loads(Path(path).read_text())["snapshots"][-1]) for path in files]
    if crop_last is not None:
        label,last=cases[-1]
        cases.append((label+f": same-state crop q/T_R={crop_last:g}",cropped(last,crop_last)))
    for label,snap in cases:
        p=np.asarray(snap["lab_momentum_edges_GeV"])
        fp,_,fc,_=snapshot_functions(snap)
        T=snap["T_visible_GeV"]
        H=np.sqrt(eos.state(T)[0]/3)/MPL
        rows=[]
        for order in orders:
            project=order==orders[-1]
            n_ch,forward,fs=gross(p,m,mp,fc,fp,False,*order,project=project)
            n_phi,reverse,rs=gross(p,mp,m,fp,fc,True,*order,project=project)
            bound=y**4/(8*np.pi*M*M)
            rate=forward[0]/(n_ch/2)**2
            assert 0.<rate<bound
            assert reverse[1]<=reverse[0] and forward[1]>=forward[0]
            err=np.max(abs(np.array([n_phi,n_ch])*snap["scale_factor"]**3/
                snap["represented_species_count_a3"]-1.))
            rows.append(dict(order=order,relative_population_quadrature_defect=err,
                average_vacuum_sigma_v_GeV_minus2=rate,
                conditional_net_stable_number_loss_per_Hubble=2*(forward[1]-reverse[1])/(H*n_ch),
                conditional_net_scalar_number_gain_per_Hubble=2*(forward[1]-reverse[1])/(H*n_phi),
                mean_injected_scalar_energy_GeV=forward[2]/(2*forward[1]),
                stimulated_over_vacuum_forward=forward[1]/forward[0],
                reverse_over_stimulated_forward=reverse[1]/forward[1],
                retained_fraction_of_produced_scalar_number=forward[3]/forward[1],
                retained_fraction_of_produced_scalar_energy=forward[4]/forward[2],
                retained_fraction_of_produced_stable_number=reverse[3]/reverse[1] if reverse[1]>0. else None,
                retained_fraction_of_produced_stable_energy=reverse[4]/reverse[2] if reverse[2]>0. else None,
                produced_scalar_number_fraction_below_p_lab_0p25_MeV=forward[5]/forward[1],
                net_stable_number_loss_per_Hubble_in_retained_cells=2*(forward[1]-reverse[3])/(H*n_ch),
                net_scalar_number_gain_per_Hubble_in_retained_cells=2*(forward[3]-reverse[1])/(H*n_phi),
                net_scalar_energy_gain_rate_in_total_energy_units_GeV5=forward[2]-reverse[2]))
        net=np.array([fs[1]-rs[0],rs[1]-fs[0]])
        counts=net[:,0].sum(axis=1)
        energy=net[:,1].sum(axis=1)+np.array([mp,m])*counts
        err_n=abs(counts.sum())/(2*(forward[1]+reverse[1]))
        err_e=abs(energy.sum())/(forward[2]+reverse[2])
        assert max(err_n,err_e)<2e-12
        result.append(dict(snapshot=label,T_GeV=T,rows=rows,
            pair_number_ledger_defect=err_n,opposite_energy_ledger_defect=err_e,
            shell_lab_momentum_edges_GeV=p.tolist(),
            shell_signed_rates_phi_then_stable_N_then_K_lo_to_hi_including_tails=net.tolist()))
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos"); p.add_argument("snapshot",nargs="+")
    p.add_argument("--crop-last",type=float,help="also test the last archive with its identical inner cells up to this q/T_R edge")
    p.add_argument("--eos-lower-knot-limit",type=float,default=.009)
    p.add_argument("--smoke",action="store_true")
    a=p.parse_args()
    eos=StandardModel(a.eos,lower_knot_limit=a.eos_lower_knot_limit)
    if a.eos_lower_knot_limit<.009:
        from visible_frw_energy_ledger_reference import ConservativeEOS
        eos=ConservativeEOS(eos)
    result=dict(controls=trace_check(),thermal_balance=None if a.smoke else thermal_check(),archives=run(eos,a.snapshot,
                [(2,5,4,2)] if a.smoke else [(3,7,6,2),(5,10,8,3)],a.crop_last))
    print(json.dumps(result,indent=2,allow_nan=False))
