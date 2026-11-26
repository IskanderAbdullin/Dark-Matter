"""Positive below-pole phi-chi2 equilibrium graph in a fixed fermion bath."""
import argparse,json
import numpy as np
from scipy.special import roots_legendre,expit
from scipy.linalg import eigh,solve
from scalar_rest_yukawa_elastic_ir_reference import cross as rest_cross
from three_species_decay_onsager_reference import M,m,mp,y,susceptibility
from decay_scalar_ir_relaxation_reference import log_rate
from sm_dark_entropy_exchange_reference import StandardModel,MPL

def spin(s,t):
    uu=2*(m*m+mp*mp)-s-t
    a=1/(s-M*M)
    b=1/(uu-M*M)
    S=M*(a+b)
    d=.5*(a*(s+m*m-mp*mp)+b*(uu+m*m-mp*mp))
    U=s*a*a+uu*b*b+2*(m*m-mp*mp)*a*b
    return y**4*(2*S*S*(2*m*m-t/2)+8*S*m*d+4*d*d+t*U)

def check_spin():
    err=[]
    c=np.array([-.99,-.8,-.2,0.,.7,.98])
    for E in (.021,.1,.7,10.,30.):
        s=m*m+mp*mp+2*mp*E
        pc2=(mp*mp*(E*E-m*m))/s
        t=-2*pc2*(1+c)
        ref=rest_cross(1.,E,c)*16*np.pi*s
        err.append(float(np.max(abs(spin(s,t)/ref-1))))
    assert max(err)<5e-13
    return max(err)

def nodes(a,b,n,*,log=False):
    x,w=roots_legendre(n)
    if len(a)==0:return np.empty((0,n)),np.empty((0,n))
    x=x[None,:]
    r=(a[:,None]*(1-x)+b[:,None]*(1+x))/2
    q=(b-a)[:,None]*w[None,:]/2
    if log:
        r=np.exp(r)
        q=q*r
    return r,q

def graph(T,edges,*,n=(5,4,2,12),pole=.82,tail=28.,reverse_order=False,eta=0.,eta_F=0.,tagged_bath=None,tagged_trial_T=None):
    assert eta<mp/T
    nU,nL,na,nf=n
    edges=np.asarray(edges)
    e_edges=np.hypot(edges,mp)
    bsz=len(edges)-1
    # Last category is an explicit thermal reservoir. Its start is the
    # last physical edge; its maximum is limited by total lab energy.
    result=np.zeros((bsz+1,bsz+1))
    u0=m+mp
    umax=np.sqrt(pole)*M
    ex=umax-u0
    cuts=np.unique(np.r_[1e-12, np.geomspace(1e-10,ex,14),ex])
    xx,vv=nodes(np.log(cuts[:-1]),np.log(cuts[1:]),nU,log=True)
    az=np.cos((np.arange(nf)+.5)*2*np.pi/nf)
    for gap,dU in zip(xx.ravel(),vv.ravel()):
        u=u0+gap
        s=u*u
        pc2=gap*(u+u0)*(u-m+mp)*(u+m-mp)/(4*s)
        pc=np.sqrt(pc2)
        ec=np.hypot(mp,pc)
        E_cap=u+T*tail
        le=(u*(ec*ec+edges*edges)/(e_edges*ec+edges*pc))
        he=u*(e_edges*ec+edges*pc)/(mp*mp)
        sub=np.unique(np.r_[u,E_cap,np.clip(le,u,E_cap),np.clip(he,u,E_cap),
                     np.clip(u+T*np.array([.01,.05,.25,.6,1,2,4,8,14,20]),u,E_cap)])
        labs,ws=nodes(sub[:-1],sub[1:],nL)
        for L,dL in zip(labs.ravel(),ws.ravel()):
            q=np.sqrt((L-u)*(L+u))
            A=L*ec/u
            B=q*pc/u
            if B<=0.:continue
            Emin=(ec*ec+(q*mp/u)**2)/(A+B)
            # Integrate both physical endpoints of each energy shell.
            l=np.clip((e_edges-Emin)/B-1,-1.,1.)
            l=np.r_[l,1.]
            labels=np.flatnonzero(np.diff(l)>0)
            if not len(labels):continue
            li,lf=(na+1,na) if reverse_order else (na,na+1)
            def subangle(num):
                z,w=roots_legendre(num)
                e0=Emin+B*(l[labels]+1)
                D=B*(l[labels+1]-l[labels])
                diff=np.log1p(D/e0)
                inc=e0[:,None]*np.expm1(diff[:,None]*(z+1)/2)
                ee=(e0[:,None]+inc).ravel()
                aa=np.clip((l[labels,None]+inc/B).ravel(),-1.,1.)
                if np.any(ee<mp*(1-1e-9)):raise AssertionError((u,L,ee.min(),mp))
                b=(diff[:,None]*w).ravel()*ee/(2*B)
                return aa,b,np.maximum(ee,mp),np.sqrt(np.maximum(0.,1-aa*aa))
            ci,wi,ei,si=subangle(li)
            cf,wf,ef,sf=subangle(lf)
            cc=ci[:,None]*cf+si[:,None]*sf*az[:,None,None]
            tv=-2*pc2*(1-cc)
            v=spin(s,tv).mean(axis=0)
            assert np.min(v)>=-1e-40
            # Evaluate the forward flux only. In/out interchange at each
            # event is a nontrivial test of the FD/BE factors and energies.
            if tagged_bath is None:
                Fi=expit(-(L-ei)/T+eta_F)/np.expm1(ei/T-eta)
                Ff=expit((L-ef)/T-eta_F)/(-np.expm1(eta-ef/T))
            else:
                # Tagged test scalar: no stimulated factor. A nonthermal
                # bath need not generate reciprocal scalar bin rates.
                Fi=tagged_bath(L-ei)
                if tagged_trial_T is not None:Fi=Fi*np.exp(-ei/tagged_trial_T)
                Ff=1-tagged_bath(L-ef)
            mat=(v*(wi*Fi)[:,None]*(wf*Ff)).reshape(len(labels),li,len(labels),lf).sum(axis=(1,3))
            result[np.ix_(labels,labels)]+=dU*dL*q*pc2/u*mat
    return result*4/(256*np.pi**5)

def metric(T,edges,eta=0.,eta_F=0.):
    lo,hi=edges[:-1].copy(),edges[1:].copy()
    # Use square root near zero to resolve the 86-eV mass scale.
    r,w=nodes(np.sqrt(lo),np.sqrt(hi),30)
    P=r*r; w*=2*r
    E=np.hypot(P,mp)
    f=1/np.expm1(E/T-eta)
    d=w*P**2*f*(1+f)/(2*np.pi**2)
    supp=d.sum(axis=1)
    rate=np.exp(np.array([log_rate(p,T,eta_F) for p in P.ravel()])).reshape(P.shape)
    return supp,(d*rate).sum(axis=1)

def independent_walk(P,out,metric,H,n=40000):
    rng=np.random.default_rng(711863)
    active=np.zeros(n,dtype=int)
    t=np.zeros(n)
    jump=np.zeros(n,dtype=int)
    cdf=np.cumsum(P,axis=1)
    for _ in range(400):
        ind=np.flatnonzero(active<len(out))
        if not ind.size:break
        old=active[ind].copy()
        t[ind]+=rng.exponential(size=len(ind))*metric[old]*H/out[old]
        xi=rng.random(len(ind))
        active[ind]=(xi[:,None]>cdf[old]).sum(axis=1)
        jump[ind]+=1
    assert not np.any(active<len(out))
    return t.mean(),t.std(ddof=1)/np.sqrt(n),jump.mean()

def trial(T,eos,bins=11,n=(5,4,2,12),pole=.82,cap=4.,eta=0.,first=.5,limit=.95,tail=28.,eta_F=0.):
    from scipy.optimize import brentq
    H=np.sqrt((eos.state(T)[0]+sum(susceptibility(T,j,g,b)[1] for j,g,b in ((M,4,False),(m,4,False),(mp,1,True))))/3)/MPL
    cross=brentq(lambda lp:log_rate(T*np.exp(lp),T,eta_F)-np.log(H),-17,0)*1.
    q=np.exp(cross)*T
    pmax=min(cap*q,limit*T)
    assert 0<first*mp<q<pmax
    ed=np.unique(np.r_[0.,q,np.geomspace(first*mp,pmax,bins)])
    G=graph(T,ed,n=n,pole=pole,eta=eta,tail=tail,eta_F=eta_F)
    R=graph(T,ed,n=n,pole=pole,reverse_order=True,eta=eta,tail=tail,eta_F=eta_F)
    Hg=np.linalg.norm(G-G.T)/np.linalg.norm(G)
    Qg=np.linalg.norm(G-R.T)/np.linalg.norm(G)
    asym_row=np.sum(abs(G-G.T),axis=1)/np.maximum(G.sum(axis=1),1e-290)
    assert max(Hg,Qg,float(asym_row.max()))<.02,(Hg,Qg,asym_row.tolist())
    G=(G+G.T)/2
    m0,k=metric(T,ed,eta,eta_F)
    L=np.diag(G.sum(axis=1)[:-1])-G[:-1,:-1]
    j=np.diag(G[:-1,:-1].sum(axis=1))-G[:-1,:-1]
    assert np.linalg.norm(j@np.ones(len(m0)))/np.linalg.norm(j)<5e-15
    B=np.diag(k)+L
    ev=eigh(B,np.diag(m0),eigvals_only=True)
    assert ev[0]>0
    tau=solve(B,m0,assume_a="pos")*H
    # Auxiliary first passage: absorbing physical bins at or above q*.
    # This is independent of the distant thermal boundary used for B.
    zq=np.searchsorted(ed,q)
    L0=L[:zq,:zq]
    L1=L0+np.diag(k[:zq])
    y0=-L[:zq,zq:].sum(axis=1)+G[:zq,-1]
    t_hit=solve(L0,m0[:zq],assume_a="pos")*H
    # Check an alternative, reflecting outer boundary with true decay
    # killing, and quantify the artificial-reservoir sensitivity.
    B_reflect=j+np.diag(k)
    t0_reflect=solve(B_reflect,m0,assume_a="pos")[0]*H
    ev_reflect=eigh(B_reflect,np.diag(m0),eigvals_only=True)
    outgoing=np.diag(L0).copy()
    n_jumps=solve(L0,outgoing,assume_a="pos")
    P_jump=G[:zq,:zq].copy()
    np.fill_diagonal(P_jump,0.)
    P_jump/=outgoing[:,None]
    t_MC,se_MC,j_MC=independent_walk(P_jump,outgoing,m0[:zq],H)
    assert abs(t_MC-t_hit[0]) < 5*se_MC,(t_MC,se_MC,t_hit[0])
    remain=np.ones(zq)
    jump_cdf={}
    for j_step in range(1,33):
        remain=P_jump@remain
        if j_step in (1,2,4,8,16,32):
            jump_cdf[str(j_step)]=float(1-remain[0])
    # Remove all but neighboring bin transitions; a conservative diagnostic
    # of the conductance available without large projected momentum jumps.
    adjacent=np.triu(G[:zq+1,:zq+1],1)-np.triu(G[:zq+1,:zq+1],2)
    adj=adjacent+adjacent.T
    A_adj=(np.diag(adj.sum(axis=1))-adj)[:zq,:zq]
    assert np.min(np.diag(adjacent,1))>0.
    neighbor_time=solve(A_adj,m0[:zq],assume_a="pos")[0]*H
    success=solve(L1,y0,assume_a="pos")
    assert np.all((success>=0)&(success<=1+1e-7))
    # Spectral integral of a positive perturbation in the full operator.
    mask=1/np.sqrt(m0)
    sp,V=eigh(B*mask[:,None]*mask)
    t_spec=((V/sp)@V.T@(1/mask))*mask*H
    assert np.max(abs(t_spec/tau-1))<3e-7
    # The closed graph dissipates perturbation free energy monotonically.
    h=np.linspace(-.17,.27,len(m0))**3
    dE=h@j@h
    balance=.5*np.sum(G[:-1,:-1]*(h[:,None]-h)**2)
    assert balance>=0 and abs(dE/balance-1)<3e-12
    z=solve(B,k,assume_a="pos")
    esc=solve(B,G[:-1,-1],assume_a="pos")
    assert np.max(abs(z+esc-1))<1e-7
    return dict(T_GeV=T,scalar_reference_eta=eta,common_fermion_eta=eta_F,below_pole_s_fraction=pole,thermal_pair_excess_cut_over_T=tail,maximum_support_p_over_T=limit,bins=len(m0),base_bins=bins,orders=n,
                independently_checked_spin_relative_error=check_spin(),
                first_cell_elastic_first_passage_to_decay_fast_H_time=t_hit[0],
                first_cell_probability_reaching_decay_fast_before_decay=success[0],
                first_cell_direct_handoff_over_H=y0[0]/m0[0]/H,
                first_cell_log10_center_decay_rate_over_H=(log_rate(.5*first*mp,T,eta_F)-np.log(H))/np.log(10.),
                bin_jumps_before_first_passage=n_jumps[0],
                first_passage_probability_by_bin_jump_count=jump_cdf,
                monte_carlo_mean_H_time=t_MC,monte_carlo_standard_error_H_time=se_MC,
                monte_carlo_average_bin_jumps=j_MC,
                neighboring_bin_only_first_passage_H_time=neighbor_time,
                maximum_row_relative_asymmetry=asym_row.max(),
                first_row_relative_asymmetry=asym_row[0],
                first_cell_reflecting_outer_mean_hubble_time=t0_reflect,
                slowest_mode_reflecting_outer_over_H=ev_reflect[0]/H,
                reciprocity_residual=Hg,swapped_quadrature_residual=Qg,
                soft_cell_limits_over_T=(ed[:2]/T).tolist(),upper_p_over_T=pmax/T,
                direct_decay_boundary_p_over_T=q/T,slowest_mode_over_H=ev[0]/H,
                first_cell_conditional_mean_hubble_time=tau[0],
                first_cell_eventual_decay_loss_probability=z[0],
                first_cell_reservoir_loss_probability=esc[0],
                first_cell_decay_only_over_H=k[0]/m0[0]/H,
                first_cell_diagonal_elastic_over_H=L[0,0]/m0[0]/H)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos");p.add_argument("--T",type=float,default=.1); p.add_argument("--bins",type=int,default=11);
    p.add_argument("--orders",nargs=4,type=int,default=[5,4,2,12]);p.add_argument("--cap",type=float,default=4.);p.add_argument("--pole",type=float,default=.82)
    p.add_argument("--eta-scalar",type=float,default=0.)
    p.add_argument("--first",type=float,default=.5)
    p.add_argument("--limit",type=float,default=.95)
    p.add_argument("--tail",type=float,default=28.)
    p.add_argument("--eta-fermion",type=float,default=0.)
    a=p.parse_args()
    print(json.dumps(trial(a.T,StandardModel(a.eos),a.bins,tuple(a.orders),a.pole,a.cap,a.eta_scalar,a.first,a.limit,a.tail,a.eta_fermion),indent=2,allow_nan=False))
