"""Prescribed expanding-history test scalar, with separate decay and mesh losses."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.special import roots_legendre,expit
from tagged_scalar_nonthermal_bath_reference import Bath
from soft_scalar_elastic_graph_reference import graph
from three_species_decay_onsager_reference import mp,m
from sm_dark_entropy_exchange_reference import StandardModel,MPL,C
from fd_number_heat_state_reference import state
from scipy.optimize import brentq

def matched(row,eos):
    # Same number and mean energy as each collisionless injection snapshot.
    T=row["TF"]
    Y=np.array(row["stable_yields"])
    n=Y.sum()*eos.state(T)[2]
    e=(Y*np.array(row["redshifted_means"][:2])).sum()/Y.sum()
    t=brentq(lambda v:(state(v,n)["rho"]-C*v**4)/n-e,.1*T,10*T)
    eta=state(t,n)["eta"]
    b=Bath(row,eos)
    q=b.q
    from scipy.integrate import cumulative_trapezoid
    b.interp_integral=-cumulative_trapezoid(expit(eta-q[::-1]/t),q[::-1],initial=0.)[::-1]
    return b,lambda E:expit(eta-E/t)

def differential(eos,row,edges,a,orders,shape,return_model=False,factorized=False,MB_trial=False):
    T=row["TF"]
    p=edges/a
    H=np.sqrt(eos.state(T)[0]/3)/MPL
    if shape=="FD":
        bath,fn=matched(row,eos)
    else:
        bath=Bath(row,eos,smooth=shape=="spline")
        fn=bath
    trial_T=T if MB_trial else None
    A=graph(T,p,n=orders,tail=34.,pole=.82,tagged_bath=fn,tagged_trial_T=trial_T)
    n=len(edges)-1
    from factorized_prompt_return_graph_reference import shell_weights
    V=shell_weights(p,trial_T)
    K=A[:n,:n].copy()
    np.fill_diagonal(K,0.)
    scatter=K.T/(V*H)
    out=A[:n,-1]/(V*H)
    z,w=roots_legendre(32)
    # Uniform in p^3 corresponds to a constant occupation within a cell.
    x=(p[:-1,None]**3*(1-z)+p[1:,None]**3*(1+z))/2
    ys=np.cbrt(x)
    ww=w[None,:]*np.ones((n,1))/2
    if MB_trial:
        # The pure-killing diagnostic must use the same weighted trial cell.
        ww=ww*np.exp(-np.hypot(ys,mp)/T)
        ww/=ww.sum(axis=1)[:,None]
    rates=np.sum(np.array([[bath.rate(u) for u in r] for r in ys])*ww,axis=1)/H
    assert min(np.min(scatter),np.min(out),np.min(rates))>=0.
    R=np.zeros((n,n+1))
    C=np.zeros(n)
    if return_model:
        if factorized:
            from factorized_prompt_return_graph_reference import graph as return_graph
            R=return_graph(p,T=T,bath=bath,occupation=fn,trial_T=trial_T)/V[:,None]/H
            C=R.sum(axis=1)
        else:
            from prompt_resonant_scalar_redistribution_reference import matrix
            R,C=matrix(bath,p,occupation=fn)
            R/=H; C/=H
    return scatter,rates,out,R,C

def prompt_source(eos,T,edges,a):
    from empty_lepton_pair_spectrum_reference import Source
    from three_species_decay_onsager_reference import M,P,Es
    E,W=Source((8,22,16),dark_mass=M).particles(T)
    rs,_,ss,_=eos.state(T)
    H=np.sqrt(rs/3)/MPL
    # Exact angle probabilities, not interpolation in outgoing p.
    p1=np.sqrt(np.maximum(0.,(E-M)*(E+M)))
    ww=2*p1*P/M
    high=(E*Es+p1*P)/M
    lo=(Es*Es+(p1*mp/M)**2)/high
    c=np.hypot(edges/a,mp)
    # Direct channel overlaps, not subtraction of cumulative totals that
    # can round to the same value (or to either side) in the extreme tail.
    a1=np.maximum(lo[:,None],c)
    b1=np.minimum(high[:,None],np.r_[c[1:],np.inf])
    f=np.maximum(0.,b1-a1)/np.maximum(ww[:,None],1e-290)
    w=(W@f)/(H*ss)
    assert min(w)>=0 and abs(w.sum()*H*ss/W.sum()-1)<1e-9
    return w

def calc(eos,rows,nbins=17,support=.45,orders=(5,4,3,12),step=1,shape="spline",return_model=False,inject=False,reflect=False,matrix_check=False,factorized=False,outer_diag=0.,new_mesh=False,MB_trial=False,return_only=False):
    assert not MB_trial or factorized
    assert not return_only or return_model
    TR=rows[0]["TR"]
    n=nbins
    er=eos.state(TR)[2]
    chosen=sorted(rows,key=lambda r:r["TF"],reverse=True)
    chosen=[r for i,r in enumerate(chosen) if i%step==0 or i==len(chosen)-1]
    z=np.r_[0.,[np.log(er/eos.state(r["TF"])[2])/3 for r in chosen]]
    edges=np.r_[0.,np.geomspace(.5*mp,TR*support,n)]
    soft=.005*chosen[-1]["TF"]*np.exp(z[-1])
    if new_mesh:
        # Add knots around typical momentum and keep an actual exterior.
        knots=np.r_[.01,.02,.05,.1,.2,.4,.65,1.,1.5,2.,3.,4.,6.,8.,12.,16.,24.,32.,48.,64.,96.,128.,192.,256.]
        edges=np.unique(np.r_[edges,knots[knots<support]*TR])
    if inject:edges=np.unique(np.r_[edges,soft])
    if outer_diag:
        assert 0<outer_diag<support
        edges=np.unique(np.r_[edges,outer_diag*TR])
    n=len(edges)-1
    outer=np.flatnonzero(edges[:-1]>=outer_diag*TR) if outer_diag else np.empty(0,dtype=int)
    inner=np.flatnonzero(edges[1:]<=outer_diag*TR) if outer_diag else np.empty(0,dtype=int)
    extra=3 if outer_diag else 0
    seq=[(np.zeros((n,n)),np.zeros(n),np.zeros(n),np.zeros((n,n+1)),np.zeros(n))]
    for a,r in zip(np.exp(z[1:]),chosen):
        seq.append(differential(eos,r,edges,a,orders,shape,return_model,factorized,MB_trial))
    out=dict(TR=TR,final_visible_temperature=chosen[-1]["TF"],n_comoving_bins=n,
             a_end=np.exp(z[-1]),support_q_over_TR=support,orders=orders,
             snapshot_count=len(z),bath_shape=shape, factorized_return=factorized,interior_boundary_q_over_TR=outer_diag, within_cell_trial="MB" if MB_trial else "constant")
    src_function=None
    if inject:
        zs=np.unique(np.r_[z,np.linspace(0,z[-1],45)])
        temps=[TR if x==0. else chosen[-1]["TF"] if x==z[-1] else
               brentq(lambda T:np.log(er/eos.state(T)[2])/3-x,chosen[-1]["TF"],TR)
               for x in zs]
        source=np.array([prompt_source(eos,T,edges,a) for T,a in zip(temps,np.exp(zs))])
        norms=np.trapezoid(source,zs,axis=0)
        total=norms.sum()
        assert total>0
        out["source_yield_total"]=float(total)
        out["source_launch_outside_mesh_fraction"]=float(norms[-1]/total)
        out["source_count_relative_to_integrated_spectrum"]=(float(total/chosen[-1]["stable_yields"][1]))
        assert abs(out["source_count_relative_to_integrated_spectrum"]-1)<.015,out["source_count_relative_to_integrated_spectrum"]
        src_function=interp1d(zs,source/total,axis=0,assume_sorted=True)
    losses=[item[1] for item in seq]
    assert max(np.max(v) for v in losses)>1.
    def solve(mode="both",kind="linear",rtol=3e-9,start=0.,q0=0.,source_on=False,matrix=False):
        As,KS,RS,Back,Event=map(tuple,zip(*seq))
        cs=interp1d(z,np.stack(As),axis=0,kind=kind,assume_sorted=True)
        kd=interp1d(z,np.stack(KS),axis=0,kind=kind,assume_sorted=True)
        kr=interp1d(z,np.stack(RS),axis=0,kind=kind,assume_sorted=True)
        back=interp1d(z,np.stack(Back),axis=0,kind=kind,assume_sorted=True)
        ev=interp1d(z,np.stack(Event),axis=0,kind=kind,assume_sorted=True)
        def M(t):
            S=cs(t) if mode in ("both","scatter","return","reflect") else np.zeros((n,n))
            D=kd(t) if mode in ("both","decay") else np.zeros(n)
            R=kr(t) if mode in ("both","scatter","return","reflect") else np.zeros(n)
            if mode in ("return","reflect"):
                mapped=back(t)
                add=mapped[:,:n].T.copy()
                np.fill_diagonal(add,0.)
                S=S+add
                R=R+mapped[:,-1] if mode=="return" else np.zeros(n)
            B=np.zeros((n+3+extra,n+3+extra))
            B[:n,:n]=S-np.diag(S.sum(axis=0)+D+R)
            B[n,:n]=D
            B[n+1,:n]=R
            if mode=="return": B[n+2,:n]=ev(t)
            if extra:
                B[n+3,inner]=S[np.ix_(outer,inner)].sum(axis=0)
                B[n+4,outer]=S[np.ix_(inner,outer)].sum(axis=0)
                B[n+5,outer]=D[outer]+R[outer]
            assert np.max(abs(B[:n+2,:].sum(axis=0)))<4e-7*np.maximum(1,np.max(B))
            return B
        N=np.zeros(n+3+extra)
        i=int(np.searchsorted(edges,q0,side="right")-1)
        if not source_on:N[i]=1.
        def rhs(t,y):
            v=M(t)@y
            if source_on:
                s=src_function(t)
                v[:n]+=s[:n]
                v[n+1]+=s[-1]
            return v
        if matrix:
            from scipy.linalg import expm
            y=N
            for lo,hi in zip(z[:-1],z[1:]):
                h=(hi-lo)/28
                for dd in range(28):
                    y=expm(h*M(lo+(dd+.5)*h))@y
            assert np.min(y)>-1e-8 and abs(y[:n+2].sum()-1)<5e-7
            return dict(surviving_in_mesh=float(y[:n].sum()),inverse_decay_removed=float(y[n]),
                        outer_mesh_exit=float(y[n+1]))
        Q=solve_ivp(rhs,(start,z[-1]),N,jac=lambda t,y:M(t),rtol=rtol,
                    atol=rtol*3e-4,method="BDF",max_step=.021)
        assert Q.success and np.min(Q.y)>-2e-7
        y=Q.y[:,-1]
        maxabs=float(abs(y[:n+2].sum()-1) if source_on else np.max(abs(Q.y[:n+2].sum(axis=0)-1)))
        assert maxabs<(2e-5 if source_on else 5e-8),maxabs
        result=dict(surviving_in_mesh=float(y[:n].sum()),inverse_decay_removed=float(y[n]),
                    first_comoving_cell_retention=float(y[0]),
                    fraction_of_initial_or_injected_below_final_p_over_T_005=float(y[:n][edges[1:]<=soft].sum()),
                    outer_mesh_exit=float(y[n+1]),explicit_inverse_return_events=float(y[n+2]),number_residual=maxabs,
                    nfev=Q.nfev,start_log_a=start,initial_q_over_TR=q0/TR)
        if extra:
            result.update(currently_in_explicit_exterior=float(y[outer].sum()),
                          cumulative_crossings_into_explicit_exterior=float(y[n+3]),
                          cumulative_returns_into_interior=float(y[n+4]),
                          removed_after_being_in_explicit_exterior=float(y[n+5]))
            launch=float(norms[outer].sum()/total) if source_on else 0.
            residual=abs(y[outer].sum()-y[n+3]+y[n+4]+y[n+5]-launch)
            assert residual<3e-5,residual
            result["separate_exterior_ledger_residual"]=float(residual)
            if mode=="return":
                from scipy.linalg import solve as linear_solve
                B=M(z[-1])
                Oo=-B[np.ix_(outer,outer)].T
                inner_rate=B[np.ix_(inner,outer)].sum(axis=0)
                exit_rate=B[n+1,outer]
                if np.all(np.diag(Oo)>0.):
                    rec=linear_solve(Oo,inner_rate)
                    esc=linear_solve(Oo,exit_rate)
                    residence=linear_solve(Oo,np.ones(len(outer)))
                    assert np.max(abs(rec+esc-1.))<2e-7
                    influx=B[np.ix_(outer,inner)]@y[inner]
                    if influx.sum()>0.:
                        w=influx/influx.sum()
                        result.update(frozen_final_snapshot_return_probability=float(w@rec),
                            frozen_final_snapshot_distant_loss_probability=float(w@esc),
                            frozen_final_snapshot_mean_exterior_residence_H_time=float(w@residence))
        if source_on:
            ef,wf=roots_legendre(30)
            lo_p=edges[:-1,None]/np.exp(z[-1])
            hi_p=edges[1:,None]/np.exp(z[-1])
            pf=(lo_p*(1-ef)+hi_p*(1+ef))/2
            ef=np.hypot(pf,mp)
            wm=wf*pf*pf
            if MB_trial:wm*=np.exp(-ef/chosen[-1]["TF"])
            avg=np.sum(wm*ef,axis=1)/wm.sum(axis=1)
            rf=np.array(chosen[-1]["stable_yields"])@np.array(chosen[-1]["redshifted_means"][:2])
            result["in_mesh_scalar_energy_relative_to_prescribed_fermions"]=float(total*(y[:n]@avg)/rf)
            f=total*er*Q.y[:n]/(np.diff(edges**3)/(6*np.pi**2))[:,None]
            i,j=np.unravel_index(np.argmax(f),f.shape)
            result["max_time_and_cell_average_occupation"]=float(f[i,j])
            result["max_terminal_cell_average_occupation"]=float(f[:,-1].max())
            result["maximum_occupation_q_over_TR_bin"]=list(edges[i:i+2]/TR)
        return result
    if not return_only:
        out["rest"]=solve()
        if matrix_check:
            out["independent_exponential_control"]=solve(matrix=True)
            for k,j in out["independent_exponential_control"].items():
                assert abs(j-out["rest"][k])<1e-4,(k,j,out["rest"][k])
        for fraction in (.25,.5,.75):
            out[f"delayed_rest_{fraction}"]=solve(start=z[-1]*fraction)
        out["rest_decay_only"]=solve("decay")
        out["rest_scatter_only"]=solve("scatter")
        out["rest_no_reaction"]=solve("none")
        out["soft_warm_start"]=solve(q0=.005*TR)
        out["tighter_solver"]=solve(rtol=3e-10)
    if return_model:
        out["rest_prompt_return"]=solve("return")
        if not return_only:out["warm_prompt_return"]=solve("return",q0=.005*TR)
        assert out["rest_prompt_return"]["inverse_decay_removed"]==0.
        if reflect:
            out["rest_prompt_reflect"]=solve("reflect")
            assert abs(out["rest_prompt_reflect"]["surviving_in_mesh"]-1)<1e-8
    if not return_only: assert abs(out["tighter_solver"]["inverse_decay_removed"]-out["rest"]["inverse_decay_removed"])<1e-6
    if not return_only: assert out["rest_no_reaction"]["surviving_in_mesh"]==1.
    if not return_only: assert out["rest_scatter_only"]["inverse_decay_removed"]==0.
    if not return_only: assert out["rest_decay_only"]["outer_mesh_exit"]==0.
    if inject:
        if not return_only:
            out["injected_no_reaction"]=solve("none",source_on=True)
            out["injected_elastic_only"]=solve("scatter",source_on=True)
        if return_model:out["injected_prompt_return"]=solve("return",source_on=True)
        if return_model and reflect:out["injected_prompt_reflect"]=solve("reflect",source_on=True)
    return out

def gen(eos):
    from empty_lepton_pair_spectrum_reference import spectrum
    return dict(rows=[spectrum(eos,t,float(x),nout=14) for t in (.02,.05,.1)
                      for x in np.geomspace(.999*t,max(.01,.55*t),12)])

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos");p.add_argument("snapshots")
    p.add_argument("--TR",type=float,default=.02);p.add_argument("--support",type=float,default=.45)
    p.add_argument("--bins",type=int,default=17);p.add_argument("--step",type=int,default=1)
    p.add_argument("--orders",type=int,nargs=4,default=[5,4,3,12])
    p.add_argument("--shape",choices=["spline","piecewise","FD"],default="spline")
    p.add_argument("--create",action="store_true")
    p.add_argument("--prompt-return",action="store_true")
    p.add_argument("--inject",action="store_true")
    p.add_argument("--reflect",action="store_true")
    p.add_argument("--matrix-check",action="store_true")
    p.add_argument("--factorized",action="store_true")
    p.add_argument("--outer-diagnostic",type=float,default=0.)
    p.add_argument("--new-mesh",action="store_true")
    p.add_argument("--MB-trial",action="store_true")
    p.add_argument("--return-only",action="store_true")
    a=p.parse_args()
    eos=StandardModel(a.eos)
    if a.create:
        Path(a.snapshots).write_text(json.dumps(gen(eos),allow_nan=False))
    else:
        data=json.loads(Path(a.snapshots).read_text())
        rows=[r for r in data["rows"] if r["TR"]==a.TR]
        print(json.dumps(calc(eos,rows,a.bins,a.support,tuple(a.orders),a.step,a.shape,a.prompt_return,a.inject,a.reflect,a.matrix_check,a.factorized,a.outer_diagnostic,a.new_mesh,a.MB_trial,a.return_only),indent=2,allow_nan=False))
