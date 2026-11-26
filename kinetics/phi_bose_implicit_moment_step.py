"""Small on-shell Yukawa χ1↔χ2+φ implicit entropy-dual moment time pilot.

NOT the four-species physical FRW simulation. Every event uses a chosen
positive numerical rate weight (not the physical Γrest), and one implicit
Euler moment step solves M(ηnext)-M(ηold)=dt*Q(ηnext). Enforces Bose
inverse fugacity positivity by line search, no particle/entropy clipping.
Charge/energy preservation and nondecreasing physical quantum entropy
follow on the admissible finite Galerkin manifold from the same hat
source, provided Newton residual is controlled.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.linalg import solve
from kinetics.phi_bose_moment_galerkin import (
    Mode, hats, eta_at, occupation_from_eta, physical_moments,
    yukawa_kinematics, yukawa_moment_source,
)


def pack(eta):
    return np.concatenate(eta)


def split(vec, modes):
    n = np.cumsum([0]+[len(m.nodes) for m in modes])
    return [vec[n[i]:n[i+1]] for i in range(len(modes))]


def source_with_jacobian(modes, eta, events):
    """All on-shell stoichiometric event source plus exact η Jacobian."""
    n = sum(len(m.nodes) for m in modes)
    offsets=np.cumsum([0]+[len(m.nodes) for m in modes])
    C=np.zeros(n); JC=np.zeros((n,n)); H=0.
    for p,c,weight in events:
        if not(math.isfinite(weight) and weight>=0):
            raise ValueError('physical positive event rate weight')
        energy = yukawa_kinematics(p,c,*[mo.mass for mo in modes])
        hats_by=[]
        for mo,E in zip(modes,energy):
            h=np.zeros(len(mo.nodes))
            j,a,b=hats(mo,E)
            h[j]=a;h[j+1]=b
            hats_by.append(h)
        fug=[float(np.dot(x,h)) for x,h in zip(eta,hats_by)]
        f1=float(occupation_from_eta(fug[0],True))
        f2=float(occupation_from_eta(fug[1],True))
        b=float(occupation_from_eta(fug[2],False))
        forward=f1*(1-f2)*(1+b)
        reverse=(1-f1)*f2*b
        J=forward-reverse
        affinity=-fug[0]+fug[1]+fug[2]
        H+=weight*J*affinity
        derivatives=(-f1*(1-f1)*(1+b-f2),
                     f2*(1-f2)*(f1+b),
                     -b*(1+b)*(f1-f2))
        for i,sign in enumerate((-1,1,1)):
            ri=slice(offsets[i],offsets[i+1])
            C[ri]+=sign*weight*J*hats_by[i]
            for j,dJ in enumerate(derivatives):
                rj=slice(offsets[j],offsets[j+1])
                JC[ri,rj]+=sign*weight*dJ*np.outer(hats_by[i],hats_by[j])
    return C,JC,H


def ledger(modes, eta, order=32):
    ints=[physical_moments(m,a,order=order) for m,a in zip(modes,eta)]
    n=[i['N'] for i in ints]
    u=[i['U'] for i in ints]
    return {'N1_plus_N2':n[0]+n[1],
            'N1_plus_Nphi':n[0]+n[2],
            'U_total':sum(u),
            'S_total':sum(i['S'] for i in ints),
            'N_phi':n[2],
            'min_phi_eta':float(min(eta[2])),
            'M':pack([i['M'] for i in ints]),
            'susceptibilities':[i['dM_deta'] for i in ints]}


def implicit_step(modes, eta_old, events, dt, *, quad_order=24,
                  tol=5e-13, max_newton=35):
    """Damped Newton full implicit Euler in integrated thermodynamic moments.

    This is a *local* finite-graph and finite-quadrature guarantee; no
    independent global time error estimate or real cosmological Γ_el.
    """
    if len(modes)!=3 or not dt>0 or not math.isfinite(dt):
        raise ValueError('three species and finite positive dt required')
    old=ledger(modes,eta_old,quad_order)
    initialM=old['M']
    start=pack([np.asarray(a,dtype=float) for a in eta_old])
    lengths=[len(m.nodes) for m in modes]
    phi_sl=slice(lengths[0]+lengths[1],sum(lengths))
    if np.any(start[phi_sl]<=0):
        raise ValueError('Bose finite inverse fugacity must be strictly positive for time step')
    state=start.copy()
    lastnorm=math.inf
    for it in range(max_newton):
        pieces=split(state,modes)
        ints=[physical_moments(m,a,order=quad_order) for m,a in zip(modes,pieces)]
        M=pack([x['M'] for x in ints])
        C,JC,H=source_with_jacobian(modes,pieces,events)
        residual=M-initialM-dt*C
        scale=max(float(np.max(abs(initialM))),1e-12)
        norm=float(np.max(abs(residual)))
        if norm<tol*scale:
            after=ledger(modes,pieces,quad_order)
            errors={k:after[k]-old[k] for k in
                    ('N1_plus_N2','N1_plus_Nphi','U_total','S_total','N_phi')}
            entropy_instant=H
            # Concavity of entropy on this moment manifold implies for
            # backward Euler ΔS >= eta_new·ΔM = dt*Σ W J ln(F/R).
            lower_bound=float(np.dot(state,M-initialM))
            if errors['S_total'] < lower_bound - 1e-10:
                raise ArithmeticError('discrete entropy inequality violated')
            return {'eta':pieces,'iterations':it,'residual_abs':norm,
                    'increments':errors,'entropy_instant':entropy_instant,
                    'entropy_concavity_lower_bound':lower_bound,
                    'max_spectral_phi_eta':max(pieces[2]),
                    'min_phi_eta':min(pieces[2]),
                    'not_cosmological':True}
        jac=-dt*JC
        off=0
        for block in ints:
            n=block['dM_deta'].shape[0]
            jac[off:off+n,off:off+n]+=block['dM_deta']
            off+=n
        try:
            direction=solve(jac,-residual,assume_a='gen',check_finite=True)
        except Exception as ex:
            raise ArithmeticError('moment susceptibility failed linear solve') from ex
        accepted=False
        alpha=1.
        for _ in range(45):
            candidate=state+alpha*direction
            if np.all(np.isfinite(candidate)) and np.all(candidate[phi_sl]>0):
                cp=split(candidate,modes)
                try:
                    new_m=pack([physical_moments(m,a,order=quad_order)['M']
                                for m,a in zip(modes,cp)])
                    new_C,_,_=source_with_jacobian(modes,cp,events)
                    newnorm=float(np.max(abs(new_m-initialM-dt*new_C)))
                    if newnorm < norm*(1-1e-4*alpha):
                        state=candidate
                        accepted=True
                        break
                except (ArithmeticError,ValueError,FloatingPointError):
                    pass
            alpha*=.5
        if not accepted:
            raise ArithmeticError('Newton cannot decrease residual within physical Bose domain')
        lastnorm=norm
    raise ArithmeticError('Newton maximum iteration reached; no false success')


def demonstration(bins=16,dt=1.,*,mphi=86e-9,weight=1e-4):
    m1,m2=.1,.02
    modes=[Mode(m,np.linspace(m,4.,bins+1),fermion=i<2)
           for i,m in enumerate((m1,m2,mphi))]
    mus=(.0,.01,-.01)
    eta=[(mo.nodes-mu)/1. for mo,mu in zip(modes,mus)]
    eta[0]=eta[0]+.12*np.sin(3*modes[0].nodes)
    events=[(p,c,weight) for p in (.0,.05,.15,.7) for c in (-.8,-.2,.3,.9)]
    res=implicit_step(modes,eta,events,dt,quad_order=24)
    return {k:v for k,v in res.items() if k!='eta'}


def evolve_moment_graph(modes, eta_initial, events, *, end_time=2., steps=8,
                        quad_order=24):
    """Sequential fully implicit finite-graph steps with audited ledger.

    No operator splitting, other SM elastic coupling or cosmological H;
    this remains first-order time with prescribed positive event weights.
    """
    if not (isinstance(steps,int) and steps>=1 and end_time>0):
        raise ValueError('positive steps and time')
    state=[np.array(x,copy=True) for x in eta_initial]
    first=ledger(modes,state,quad_order)
    total_iterations=0
    min_entropy=math.inf
    max_residual=0.
    for _ in range(steps):
        step=implicit_step(modes,state,events,end_time/steps,
                           quad_order=quad_order)
        state=step['eta']
        total_iterations+=step['iterations']
        min_entropy=min(min_entropy,step['increments']['S_total'])
        max_residual=max(max_residual,step['residual_abs'])
    last=ledger(modes,state,quad_order)
    return {'state':state,'steps':steps,'end_time':end_time,
            'delta_Nphi':last['N_phi']-first['N_phi'],
            'min_step_entropy_gain':min_entropy,
            'total_entropy_gain':last['S_total']-first['S_total'],
            'charge_N1N2_residual':last['N1_plus_N2']-first['N1_plus_N2'],
            'charge_N1Nphi_residual':last['N1_plus_Nphi']-first['N1_plus_Nphi'],
            'total_energy_residual':last['U_total']-first['U_total'],
            'total_newton_iterations':total_iterations,
            'max_nonlinear_residual_abs':max_residual,
            'min_final_phi_eta':float(min(state[2])),
            'only_synthetic_one_to_two_graph':True}


def time_grid_audit(*, bins=16, mphi=86e-9, end_time=2.,
                    time_steps=(2,4,8,16,32,64)):
    modes=[Mode(m,np.linspace(m,4.,bins+1),fermion=i<2)
           for i,m in enumerate((.1,.02,mphi))]
    mus=(.0,.01,-.01)
    initial=[mo.nodes-mu for mo,mu in zip(modes,mus)]
    initial[0]=initial[0]+.12*np.sin(3*modes[0].nodes)
    events=[(p,c,1e-4) for p in (0.,.05,.15,.7)
           for c in (-.8,-.2,.3,.9)]
    rows=[]
    for steps in time_steps:
        z=evolve_moment_graph(modes,initial,events,end_time=end_time,steps=steps)
        rows.append({k:v for k,v in z.items() if k!='state'})
    selferr=[]
    for x,y in zip(rows[:-1],rows[1:]):
        selferr.append({'steps':x['steps'],
                        'Nphi_self_diff':abs(x['delta_Nphi']-y['delta_Nphi'])})
    return {'mass_GeV_synthetic':mphi,'mode_bins':bins,
            'method':'fully implicit Euler on physical Bose/Fermi moment space',
            'rows':rows,'self_differences':selferr,'physical_dm_relic':False}


def physically_weighted_yukawa_events(*, m1=.1,m2=.02,mphi=86e-9,
                                       y=.001,pmax=2.,nr=16,nangle=12,g1=4):
    """Actual tree-level vacuum Yukawa Γ_rest with moving mother time dilation.

    Physical on-shell measure for χ1↔χ2+real φ weak number, in common
    natural GeV units; still vacuum amplitudes, prescribed f, no thermal
    masses / SM plasma / FRW / actual parameter constraints.
    """
    from numpy.polynomial.legendre import leggauss
    if not(m1>m2+mphi>=0 and y>0 and pmax>0 and nr>=2 and
           nangle>=2 and g1>0):raise ValueError('physical open Yukawa phase space')
    lam=(m1*m1-(m2+mphi)**2)*(m1*m1-(m2-mphi)**2)
    width=y*y*((m1+m2)**2-mphi*mphi)*math.sqrt(lam)/(16*math.pi*m1**3)
    xp,wp=leggauss(nr)
    ps=(xp+1)*pmax/2
    pw=wp*pmax/2
    cc,wc=leggauss(nangle)
    events=[]
    for p,weight_p in zip(ps,pw):
        gamma=width*m1/math.hypot(p,m1)
        for c,weight_c in zip(cc,wc):
            W=g1/(2*math.pi**2)*weight_p*p*p*gamma*weight_c/2
            events.append((float(p),float(c),float(W)))
    return events,width


def physical_weak_rate_control(*,mphi=86e-9,bins=16,nr=12,nangle=12,
                               pmax=2.,y=.001):
    """Matched full parent p/c weak scalar rate from dual mode source.

    Independently integrate the same on-shell F/R with a second direct
    quadrature loop (not using the hat-moment source deposition).
    """
    from scipy.special import expit
    from numpy.polynomial.legendre import leggauss
    m1,m2=.1,.02
    upper=math.hypot(m1,pmax)+.2
    modes=[Mode(m,np.linspace(m,upper,bins+1),degeneracy=(4 if i<2 else 1),
                fermion=i<2) for i,m in enumerate((m1,m2,mphi))]
    # Closed common-T with a slightly colder parent; negative real-phi net.
    eta=[(mo.nodes/T) for mo,T in zip(modes,(.9,1.,1.))]
    event_list,width=physically_weighted_yukawa_events(
        m1=m1,m2=m2,mphi=mphi,y=y,pmax=pmax,nr=nr,nangle=nangle,g1=4)
    projected=yukawa_moment_source(modes,eta,event_list)
    px,pw=leggauss(nr); cx,cw=leggauss(nangle)
    independent=0.
    for x,w in zip(px,pw):
        p=pmax*(x+1)/2
        wp=w*pmax/2
        e1=math.hypot(p,m1)
        a=expit(-eta_at(modes[0],eta[0],e1))
        gamma=width*m1/e1
        for c,wa in zip(cx,cw):
            _,e2,eph=yukawa_kinematics(p,c,m1,m2,mphi)
            d=expit(-eta_at(modes[1],eta[1],e2))
            b=1/math.expm1(eta_at(modes[2],eta[2],eph))
            F=a*(1-d)*(1+b)
            R=(1-a)*d*b
            independent+=4/(2*math.pi**2)*wp*p*p*gamma*wa/2*(F-R)
    return {'mphi':mphi,'vacuum_gamma_rest':width,
            'Rphi_from_moment_projection':projected['rate_phi'],
            'Rphi_independent_parent_pc_integral':independent,
            'relative_discrepancy':projected['rate_phi']/independent-1,
            'entropy_source':projected['entropy_from_event'],
            'full_physical_SM_FRW':False}


def physical_weight_time_audit(*,mphi=86e-9,pmax=2.,
                               bins=(8,16,32),times=(1e6,1e7,1e8),
                               nr=12,nangle=12,y=.001):
    """Conditional scalar Yukawa microphysics time, no finite SM or FRW.

    Uses actual Γrest mass+coupling and positive on-shell parent p/c
    measure but the benchmark masses, temperatures and Yukawa are only
    illustrative; not a viable model parameter fit or relic prediction.
    """
    m1,m2=.1,.02
    upper=math.hypot(m1,pmax)+.2
    events,width=physically_weighted_yukawa_events(m1=m1,m2=m2,
        mphi=mphi,y=y,pmax=pmax,nr=nr,nangle=nangle,g1=4)
    rows=[]
    for n in bins:
        modes=[Mode(m,np.linspace(m,upper,n+1),
                    degeneracy=4 if i<2 else 1,fermion=i<2)
               for i,m in enumerate((m1,m2,mphi))]
        eta=[mo.nodes/T for mo,T in zip(modes,(.9,1.,1.))]
        initial=ledger(modes,eta,order=24)
        C,JC,H=source_with_jacobian(modes,eta,events)
        initial_weak=float(C[2*(n+1):].sum())
        for dt in times:
            step=implicit_step(modes,eta,events,dt,quad_order=24)
            rows.append({'bins':n,'time_GeV_inverse_synthetic':dt,
                 'Nphi_increment':step['increments']['N_phi'],
                 'initial_weak_Rphi':initial_weak,
                 'initial_weak_entropy':H,
                 'S_increment':step['increments']['S_total'],
                 'entropy_lower_bound':step['entropy_concavity_lower_bound'],
                 'charge_errors':[step['increments']['N1_plus_N2'],
                                  step['increments']['N1_plus_Nphi']],
                 'energy_error':step['increments']['U_total'],
                 'min_phi_eta':step['min_phi_eta'],
                 'newton_iterations':step['iterations'],
                 'residual_abs':step['residual_abs']})
    return {'schema_version':1,'stage':'physical-yukawa-vacuum-width but synthetic masses/thermal populations',
       'mphi_GeV':mphi,'m1_GeV':m1,'m2_GeV':m2,'y_assigned':y,
       'vacuum_gamma_rest':width,'mother_pmax_GeV':pmax,
       'gauss_nr':nr,'gauss_nangle':nangle,'rows':rows,
       'not_real_FRW_or_physical_86eV_relic':True}


if __name__=='__main__':
    import argparse,json
    parser=argparse.ArgumentParser()
    parser.add_argument('--audit',action='store_true')
    args=parser.parse_args()
    if args.audit:
        data={'schema_version':1,'date_utc':'2026-11-22',
              'implicit_first_order_time_self_convergence':time_grid_audit(),
              'vacuum_width_physical_yukawa_conditional_time':physical_weight_time_audit(),
              'parent_weak_rate_radial_reference':[physical_weak_rate_control(mphi=86e-9,bins=16,nr=n,nangle=12)
                                               for n in (8,12,16,24,32)],
              'scope':'physical p-space moment implicit finite Yukawa 1to2 only, no SM/FRW/cosmological Ω/BEC/CLASS'}
    else:
        data={str(dt):demonstration(dt=dt) for dt in (.25,.5,1.,2.)}
    print(json.dumps(data,indent=2,default=float))
