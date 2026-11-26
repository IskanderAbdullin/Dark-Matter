"""Heavy visible portal on an adiabatic parent and full daughter support.

The linear fast-decay response to an external parent source uses A+B (not
B alone). A and B are local inverse and decay blocking factors respectively.
The visible bath is prescribed; this is a fixed-snapshot response.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.special import expit, roots_legendre
from empty_lepton_pair_spectrum_reference import Source
from sm_dark_entropy_exchange_reference import StandardModel, MPL
from three_species_decay_onsager_reference import M,m,mp,P,Es,width


def angle_panels(parent, scalar_gates=(), fermion_gates=(), angle=14, soft=()):
    """Normalized lab angle, with exact cuts for both daughter supports."""
    dp=np.sqrt((parent-M)*(parent+M))
    hi=(parent*Es+dp*P)/M
    lo=(Es*Es+(mp*dp/M)**2)/hi
    width=2*dp*P/M
    assert min(width)>0
    gates=np.concatenate([lo[:,None],hi[:,None],
                          np.broadcast_to(np.r_[scalar_gates,soft],(len(parent),len(scalar_gates)+len(soft))),
                          parent[:,None]-np.asarray(fermion_gates)[None,:]],axis=1)
    e=np.sort(np.clip(gates,lo[:,None],hi[:,None]),axis=1)
    low=e[:,:-1]
    high=e[:,1:]
    z,w=roots_legendre(angle)
    ep=((low+high)[:,:,None]+(high-low)[:,:,None]*z)/2
    w=w*(high-low)[:,:,None]/(2*width[:,None,None])
    assert np.max(abs(w.sum(axis=(1,2))-1))<2e-10
    return ep,w


def parent_branch(parent,phi,chi,scalar_gates=(),fermion_gates=(),soft=(),angle=14,batch=256):
    x=parent.ravel()
    o=np.zeros((len(x),9+2*len(soft)))
    for i in range(0,len(x),batch):
        sl=slice(i,i+batch)
        p=x[sl]
        E,w=angle_panels(p,scalar_gates,fermion_gates,angle=angle,soft=soft)
        X=p[:,None,None]-E
        f=phi(E)
        c=chi(X)
        A=w*f*c
        B=w*(1+f)*(1-c)
        C=A+B
        A0=A.sum(axis=(1,2))
        B0=B.sum(axis=(1,2))
        C0=A0+B0
        # Internal eliminated-parent occupancy, normalized linear response
        # branching, and a separate decay-only angular comparator.
        o[sl,0]=A0/C0
        o[sl,1]=(C*E).sum(axis=(1,2))/C0
        o[sl,2]=(C*X).sum(axis=(1,2))/C0
        o[sl,3]=(B*E).sum(axis=(1,2))/B0
        o[sl,4]=(w*E).sum(axis=(1,2))
        for j,cx in enumerate(soft):
            o[sl,5+j]=(C*(E<cx)).sum(axis=(1,2))/C0
            o[sl,5+len(soft)+j]=(w*(E<cx)).sum(axis=(1,2))
        n=5+2*len(soft)
        if len(scalar_gates) and len(fermion_gates):
            outside_phi=(E<scalar_gates[0])|(E>scalar_gates[-1])
            outside_chi=(X<fermion_gates[0])|(X>fermion_gates[-1])
            o[sl,n]=(C*(outside_phi|outside_chi)).sum(axis=(1,2))/C0
            o[sl,n+1]=(C*E*outside_phi).sum(axis=(1,2))/C0
            o[sl,n+2]=(C*X*outside_chi).sum(axis=(1,2))/C0
        o[sl,-1]=C0
        if len(soft):
            part=o[sl,5:n]
            assert min(part.ravel())>=0 and max(part.ravel())<=1+2e-14
    assert np.max(abs((o[:,1]+o[:,2])/x-1))<5e-15
    return o.reshape(parent.shape+(-1,))


def portal(T,phi,chi,scalar_gates=(),fermion_gates=(),order=(6,14,14),polar=18,decay=10,H=None):
    src=Source(order,dark_mass=M)
    q,p,v,N,A=src.event_moments(T)
    x,w=roots_legendre(polar)
    E1=(q[:,:,None]+p[:,:,None]*np.sqrt(v)[:,None,None]*x)/2
    E2=q[:,:,None]-E1
    measure=(N[:,:,None]+A[:,:,None]*(3*x*x-1)/2)*w
    cuts=np.hypot(mp,T*np.array([.025,.05,.1,1.]))
    # E2 is the reflection across the symmetric angular grid, as are its
    # parent occupation and all angular decay averages.
    a=parent_branch(E1,phi,chi,scalar_gates,fermion_gates,cuts,angle=decay)
    f=a[...,0]
    g=f[:,:,::-1]
    block=(1-f)*(1-g)
    safe_f=np.where(f>0.,f,1.)
    log_f=np.where(f>0.,np.log(safe_f),-np.inf)
    log_rev=q[:,:,None]/T+log_f+log_f[:,:,::-1]
    log_block=np.log1p(-f)+np.log1p(-g)
    aff=log_rev-log_block
    back=np.exp(log_rev)
    net=-np.sign(aff)*np.exp(np.maximum(log_rev,log_block))*(-np.expm1(-abs(aff)))
    assert np.max(abs(E1[:,:,::-1]-E2)/E1)<3e-13
    assert np.max(abs(net-block+back))<5e-14*max(1.,np.max(abs(net)))
    vacuum=np.stack([np.ones_like(E1),E1],axis=-1)
    norm=np.sum(measure[...,None]*vacuum,axis=(0,1,2))
    def normalized(rows):
        return (np.sum(measure[...,None]*rows,axis=(0,1,2))/norm).tolist()
    moments={name:normalized(z[...,None]*vacuum) for name,z in (('net',net),('forward',block),('inverse',back))}
    en=lambda r: float(np.sum(measure*r)/norm[1])
    Ephi=a[...,1]; Echi=a[...,2]
    out=dict(portal_number_energy_over_vacuum=moments, largest_parent_occupation=float(f.max()),
             scalar_energy_over_vacuum_parent_energy=en(net*Ephi),
             fermion_energy_over_vacuum_parent_energy=en(net*Echi),
             scalar_decay_only_reference_over_vacuum_parent_energy=en(net*a[...,3]),
             isotropic_scalar_energy_over_vacuum_parent_energy=en(a[...,4]),
             net_scalar_count_below_p_over_T_over_vacuum_count=
               [float(np.sum(measure*net*a[...,j])/norm[0]) for j in range(5,5+len(cuts))],
             isotropic_scalar_count_below_p_over_T_over_vacuum_count=
               [float(np.sum(measure*a[...,j])/norm[0]) for j in range(5+len(cuts),5+2*len(cuts))],
             net_branch_outside_input_grid_number_over_vacuum_parent_count=
               float(np.sum(measure*net*a[...,-4])/norm[0]),
             net_branch_outside_input_grid_energy_over_vacuum_parent_energy=
               [en(net*a[...,-3]),en(net*a[...,-2])])
    if H is not None:
        lam=width*(M/E1)*a[...,-1]
        out["minimum_decay_response_rate_over_visible_H"]=float(np.min(lam)/H)
        out["energy_source_weighted_mean_visible_H_times_decay_response_time"]=en(H*E1/lam)
    assert abs(out['scalar_energy_over_vacuum_parent_energy']+out['fermion_energy_over_vacuum_parent_energy']-moments['net'][1])<2e-13
    return out


def from_bins(row,eos):
    T=row['TF']
    ed=T*np.array(row['scalar_differential_bins'])
    ss=eos.state(T)[2]
    # Scalar and total stable bin occupations; unknown exterior set to
    # vacuum but included in A/B angular normalization and final moments.
    def shape(mass,g,name,n):
        f=n*np.array(row[name])/(g*np.diff(ed**3)/(6*np.pi**2))
        edge=np.hypot(ed,mass)
        gate=edge[1:]
        def fun(E):
            ix=np.searchsorted(gate,E,side='right')
            v=f[np.minimum(ix,len(f)-1)]
            return np.where(ix<len(f),v,0.)
        return fun,edge,float(f.max())
    fp,ep,fcap=shape(mp,1,'scalar_normalized_counts',row['stable_yields'][1]*ss)
    fc,ec,ccap=shape(m,4,'stable_normalized_counts',sum(row['stable_yields'])*ss)
    assert ccap<1.
    return fp,fc,ep,ec,fcap,ccap


def shell_projection(T,phi,chi,p_mom,sg=(),cg=(),order=(6,12,14),polar=16,decay=1):
    """Resolved signed daughter rates and tails, seeded by the net portal.

    Output axes are species [phi, stable], sign [positive, negative],
    moment [number, kinetic energy], and [low, resolved shells, high].
    Negative values are stored as positive magnitudes in their own channel.
    """
    src=Source(order,dark_mass=M)
    q,p,v,N,A=src.event_moments(T)
    x,w=roots_legendre(polar)
    E=(q[:,:,None]+p[:,:,None]*np.sqrt(v)[:,None,None]*x)/2
    measure=(N[:,:,None]+A[:,:,None]*(3*x*x-1)/2)*w
    parents=parent_branch(E,phi,chi,sg,cg,angle=decay)
    f=parents[...,0]
    g=f[:,:,::-1]
    log_rev=q[:,:,None]/T+np.log(np.where(f>0.,f,1.))+np.log(np.where(g>0.,g,1.))
    log_rev=np.where((f>0)&(g>0.),log_rev,-np.inf)
    log_forward=np.log1p(-f)+np.log1p(-g)
    aff=log_rev-log_forward
    net=-np.sign(aff)*np.exp(np.maximum(log_rev,log_forward))*(-np.expm1(-abs(aff)))
    source=(measure*net).ravel()
    E=E.ravel()
    edges=[np.hypot(k,p_mom) for k in (mp,m)]
    sg=np.unique(np.r_[sg,edges[0]]); cg=np.unique(np.r_[cg,edges[1]])
    out=np.zeros((2,2,2,len(p_mom)+1))
    for j in range(0,len(E),180):
        sl=slice(j,j+180)
        ep,w=angle_panels(E[sl],sg,cg,angle=decay)
        ec=E[sl,None,None]-ep
        C=w*((1+phi(ep))*(1-chi(ec))+phi(ep)*chi(ec))
        C/=C.sum(axis=(1,2))[:,None,None]
        dst=[np.searchsorted(ed,ev,side="right").ravel() for ed,ev in zip(edges,(ep,ec))]
        for k in range(2):
            kinetic=((ep,ec)[k]-(mp,m)[k]).ravel()
            for sign in range(2):
                real=np.maximum(source[sl]*(-1)**sign,0.)[:,None,None]*C
                for moment in range(2):
                    row=(real.ravel() if moment==0 else real.ravel()*kinetic)
                    out[k,sign,moment]+=np.bincount(dst[k],weights=row,minlength=len(p_mom)+1)
    masses=np.array([mp,m])[:,None]
    counts=out[:,:,0].sum(axis=2)
    energy=out[:,:,1].sum(axis=2)+masses*counts
    expected=np.array([np.maximum(source,0.),np.maximum(-source,0.)])
    assert np.max(abs(counts-expected.sum(axis=1))/np.maximum(expected.sum(axis=1),1e-280))<3e-12
    assert np.max(abs(energy.sum(axis=0)-expected@E)/np.maximum(expected@E,1e-280))<3e-12
    return out


def shell_test(eos,row):
    from joint_visible_fd_source_reference import Visible
    T=row["TF"]
    p=T*np.r_[0.,np.geomspace(.005,16.,23)]
    order=(6,12,14); angle=16
    z=shell_projection(T,np.zeros_like,np.zeros_like,p,order=order,polar=angle)
    orig=Visible(order,polar=angle).rates(T,p,1.)[[2,1]]
    per_species=np.sum(abs(z[:,0]-orig),axis=-1)/orig.sum(axis=-1)
    assert np.max(per_species)<3e-12 and np.max(z[:,1])==0.
    fp,fc,sg,cg,_,_=from_bins(row,eos)
    out=shell_projection(T,fp,fc,p,sg,cg,order=order,polar=angle)
    ref=z[:,0].sum(axis=-1)
    diff=(out[:,0]-out[:,1]-z[:,0])/ref[:,:,None]
    return dict(T_GeV=T,lab_edges_over_T=(p/T).tolist(),
                vacuum_L1_error_number_energy=per_species.tolist(),
                quantum_shell_correction_L1_normalized=np.sum(abs(diff),axis=-1).tolist(),
                signed_shell_moments_over_vacuum_totals=((out[:,0]-out[:,1])/ref[:,:,None]).tolist(),
                negative_number_kinetic_over_vacuum_totals=(out[:,1].sum(axis=-1)/ref).tolist())

def snapshot_functions(snapshot, exterior="vacuum"):
    """Read represented entropy cells with an explicit exterior convention."""
    p=np.asarray(snapshot["lab_momentum_edges_GeV"])
    assert exterior=="vacuum"
    def cells(mass,key,BE):
        th=np.asarray(snapshot[key])
        edge=np.hypot(mass,p)
        lo= p[:-1]**2/(edge[:-1]+mass)
        dk= p[1:]**2/(edge[1:]+mass)-lo
        def fn(E):
            j=np.searchsorted(edge,E,"right")-1
            c=np.clip(j,0,len(lo)-1)
            u=np.clip(((E-mass)-lo[c])/dk[c],0.,1.)
            z=th[0,c]+th[1,c]*u
            if BE:
                assert np.min(z)>0.
                f=np.exp(-z)/(-np.expm1(-z))
            else:
                f=expit(-z)
            return np.where((j>=0)&(j<len(lo)), f, 0.)
        return fn,edge
    return (*cells(mp,"entropy_variable_scalar",True),*cells(m,"entropy_variable_fermion",False))


def evolved_selftest(eos,snapshots):
    result=[]
    for x in snapshots:
        fp,sg,fc,cg=snapshot_functions(x)
        T=x["T_visible_GeV"]
        rows=[]
        for order,polar,decay in (((6,12,14),16,4),((8,18,16),20,6)):
            r=portal(T,fp,fc,sg,cg,order=order,polar=polar,decay=decay,H=np.sqrt(eos.state(T)[0]/3)/MPL)
            rows.append(dict(order=order,polar=polar,decay=decay,**r))
        result.append(dict(T_GeV=T,steps=x['steps'],support_GeV=[x['lab_momentum_edges_GeV'][0],x['lab_momentum_edges_GeV'][-1]],controls=rows))
    return result

def forced_parent_test():
    """Two angular panels, fixed baths and an external parent source."""
    from scipy.integrate import solve_ivp
    A=np.array([.15,.72])
    B=np.array([.42,.81])
    initial=A.sum()/(A+B).sum()
    driver=.04
    def rhs(_,y):
        local=y[0]*B-(1-y[0])*A
        return np.r_[driver-local.sum(),local]
    t=11.
    base=np.r_[initial,0.,0.]
    forced=solve_ivp(rhs,(0.,t),base,rtol=3e-13,atol=1e-15).y[:,-1]
    steady=initial*B-(1-initial)*A
    diff=forced[1:]-t*steady
    observed=diff/diff.sum()
    predicted=(A+B)/(A+B).sum()
    assert max(abs(observed-predicted))<1e-12
    assert abs((forced[0]-initial)+diff.sum()-driver*t)<1e-12
    return dict(numerical_branch=observed.tolist(),response_branch=predicted.tolist(),
                decay_only_branch=(B/B.sum()).tolist())

def selftest(eos,rows):
    out=[]
    for row in rows:
        T=row['TF']
        H=np.sqrt(eos.state(T)[0]/3)/MPL
        phi,chi,sg,cg,fcap,ccap=from_bins(row,eos)
        runs=[]
        for order,polar in (((6,12,14),16),((8,18,16),20)):
            model=portal(T,phi,chi,sg,cg,order,polar,decay=1,H=H)
            runs.append(dict(order=order,polar=polar,**model))
        # At a common T and eta_chi+eta_phi=0, the angle-dependent detailed
        # balance ratio is exactly exp(-E_H/T) and yields zero portal net.
        eq=portal(T,lambda e:1/np.expm1(e/T+.1),lambda e:expit(.1-e/T))
        vac=portal(T,np.zeros_like,np.zeros_like)
        assert max(abs(v) for v in eq['portal_number_energy_over_vacuum']['net'])<2e-13
        assert max(abs(np.array(vac['portal_number_energy_over_vacuum']['net'])-1))<1e-14
        assert abs(vac['scalar_energy_over_vacuum_parent_energy']-Es/M)<5e-13
        out.append(dict(TR=row['TR'],TF=T,scalar_maximum_bin_occupation=fcap,fermion_maximum_bin_occupation=ccap,
                        controls=runs,equilibrium=eq,vacuum=vac))
    return out


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos");p.add_argument("spectra")
    p.add_argument("--reheat",type=float,nargs="+",default=[.02,.05,.1])
    p.add_argument("--evolved-snapshots")
    p.add_argument("--shell-test",action="store_true")
    a=p.parse_args()
    rows=json.loads(Path(a.spectra).read_text())["rows"]
    r=[next(row for row in rows if row['TR']==tr and abs(row['TF']-.9*tr)<1e-8) for tr in a.reheat]
    result=dict(forced_parent_angular_test=forced_parent_test(),rows=selftest(StandardModel(a.eos),r),scalar_thresholds_over_T=[.025,.05,.1,1.])
    if a.shell_test:
        result["shell_projection_tests"]=[shell_test(StandardModel(a.eos),row) for row in r]
    if a.evolved_snapshots:
        snapshots=json.loads(Path(a.evolved_snapshots).read_text())["snapshots"]
        result["evolved_fixed_snapshots"]=evolved_selftest(StandardModel(a.eos),snapshots)
    print(json.dumps(result,indent=2,allow_nan=False))
