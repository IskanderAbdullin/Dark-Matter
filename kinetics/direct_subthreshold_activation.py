"""Direct late cubic activation without a second nearly degenerate species.

V = M^2 A^2/2 + mS^2/2*(S + g0*x*A^2/(2*mS^2))^2
    + V_chi,eta; x=chi/v.
The square supplies the stabilizer lambda_A(x)=3*g0^2*x^2/mS^2,
zero at x=0. Renormalized trigger potential prescribed; compare the
two-loop MS threshold, including the one-loop S tadpole displacement,
at several matching scales. No claim of UV completion.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.integrate import trapezoid,solve_ivp
from scipy.special import zeta
from subthreshold_roll_scan import params,roll
from soft_mixing_activation import variance_table,bath_integrals
from subthreshold_transport_tradeoff import rayleigh_grid
from subthreshold_recombination import halo,zero_rate
from direct_activation_consistency import rho_tot,CGS_RHO_TO_EV4,T0,radiation_T4
from late_activation_overlap import H,Om,Or,critical_density
from resonant_mediator_inverse import M,C,HBAR,SEC_GYR
from quartic_consistency import lambda_from_sigma
from velocity_only_exposure import RHO0

VREF=np.sqrt(2)*44.1e5
SCALE=M*M*(VREF/C)**2
RHO=critical_density(VREF)

def sunset_threshold(g2,Q):
    """Near mS=2M limit of the Martin MS-bar vacuum sunset.
    I/M^2 at (1,1,4), Q/M arbitrary; Kallen discriminant zero.
    Relative error O(det/M^2 * log(det/M^2)) for our tiny detuning.
    """
    la=-2*np.log(Q/M);ls=np.log(4)+la
    I=M*M*(-4*la*ls+la*la+4*la+8*ls-15)
    tad=M*M*(la-1.)
    fig=3*g2*tad*tad/(8*(4*M*M))
    shift=-g2*tad*tad/(8*(4*M*M))
    sun=-g2*I/4
    return dict(Q_over_M=Q/M,I_eV2=float(I),
                figure8_eV4=float(fig/(16*np.pi**2)**2),
                relaxed_S_tadpole_eV4=float(shift/(16*np.pi**2)**2),
                sunset_eV4=float(sun/(16*np.pi**2)**2),
                total_eV4=float((fig+shift+sun)/(16*np.pi**2)**2))

def solve(row,p,scheme,g2=0.,seed=1e-8,n=14000,rtol=1e-8,medium_sign=1.):
    # Medium is a deliberately local Gaussian static-S diagnostic, not
    # exact real-time nonlocal response. Its extrapolation to early
    # relativistic epochs is used only for sensitivity.
    V,Tr=variance_table(row);eta=bath_integrals()
    mu=p['mu_eV'];X=p['X'];r_on=p['z_on']+1
    td0=p['Td_on_eV']/r_on;m_eta=p['m_eta_at_vev_eV'];cw=p['eta_CW_over_X']
    def fun(t,y):
        r=np.exp(-t);x,dx=y;hh=HBAR*H(r-1)
        dh=-.5*(3*Om*r**3+4*Or*r**4)/(Om*r**3+Or*r**4+1-Om-Or)
        therm=12*(r/r_on)**2*eta(m_eta*abs(x)/(td0*r))
        med=medium_sign*g2*V(r)**2/(2*4*M*M*X)
        xx=x*x;lx=np.log(max(xx,1e-300))
        cwterm=4*cw*(xx*lx if scheme=='origin-vev' else 1+xx*(lx-1))
        return [dx,-(3+dh)*dx-(mu/hh)**2*x*(xx-1+therm+med+cwterm)]
    t0=-np.log(1e14)
    ts=np.r_[np.linspace(t0,-np.log(101.),1000),np.linspace(-np.log(101.),0.,n)[1:]]
    sol=solve_ivp(fun,(t0,0.),[seed,0.],t_eval=ts,rtol=rtol,atol=rtol*1e-10,method='DOP853',max_step=.045,first_step=.01)
    assert sol.success,sol.message
    r=np.exp(-ts);x,dx=sol.y
    v0=row['final']['peak_bin_q']*Tr/M*C
    return dict(t=ts,r=r,x=x,dx=dx,nfev=sol.nfev,v0=v0,
                max_early_abs_x=float(max(abs(x[r>=1001]))),
                max_abs_x=float(max(abs(x))),today_x=float(x[-1]))

def diagnose(row,p,rr,aa,QQ,formation_z=2.,efold=1.):
    t,r,x,dx=[rr[k] for k in ['t','r','x','dx']]
    clock=1/H(r-1);h=r<=formation_z+1;late=r<=101
    A4=float(trapezoid(x[h]**4*clock[h],t[h]))
    if A4<=1e8:return dict(active=False)
    cosmic=float(trapezoid(r[late]**4*x[late]**4*clock[late],t[late]))
    sig0=efold/(RHO*VREF*QQ*A4);g2=lambda_from_sigma(sig0)*aa*SCALE
    tau=sig0*RHO0*rr['v0']*cosmic
    loop=[sunset_threshold(v,M) for v in g2]
    loopM=np.array([2*abs(x['total_eV4'])/p['X'] for x in loop])
    loopall=np.array([max(2*abs(sunset_threshold(v,q*M)['total_eV4'])/p['X'] for q in [.5,1.,2.]) for v in g2])
    # Compare tree trigger energy and a simple upper bound on eta energy
    # with the *actual* background during the roll, not nominal z_on.
    # Vacuum offset chosen at x=1; include specified CW potential shapes
    # separately in a later calculation if needed.
    rho=rho_tot(r-1);X=p['X'];mu=p['mu_eV']
    tree=X*(.5*(HBAR*H(r-1)/mu)**2*dx*dx+.25*(x*x-1.)**2)
    td=p['Td_on_eV']*r/(p['z_on']+1);nn=p['Neta']*zeta(3,1)*td**3/np.pi**2
    eta_upper=p['Neta']*np.pi**2/30*td**4+nn*p['m_eta_at_vev_eV']*abs(x)
    emax=float(max((tree+eta_upper)[late]/rho[late]))
    hhalf=np.flatnonzero(abs(x)**4>=.5)
    # The direct pole has no H conversion. Evaluate LL scattering at
    # cluster speed, s-channel dominant; no velocity-independent limit assumed.
    cs=sig0*(aa/(aa+(1000e5/VREF)**2))**2
    rows=[]
    for i,a in enumerate(aa):
        rows.append(dict(a=float(a),Q=float(QQ[i]),sigma0=float(sig0[i]),g2_eV2=float(g2[i]),
           lambda_stabilizer=float(3*g2[i]/(4*M*M-a*SCALE)),LL_tau=float(tau[i]),
           cosmic_proxy=float(5.756e7*tau[i]),loop_over_X_at_M=float(loopM[i]),
           max_loop_over_X_scales=float(loopall[i]),LL_sigma_m_1000kms=float(cs[i]),
           max_tree_plus_eta_fraction_zle100=emax))
    passr=[x for x in rows if x['cosmic_proxy']<=.01 and x['LL_sigma_m_1000kms']<=.2 and x['max_tree_plus_eta_fraction_zle100']<=.04]
    return dict(active=True,formation_z=formation_z,active_time_Gyr=A4/SEC_GYR,
       z_half=float(r[hhalf[0]]-1) if len(hhalf) else None,rows=rows,
       best=min(passr,key=lambda x:x['loop_over_X_at_M']) if passr else None,
       best_robust=min(passr,key=lambda x:x['max_loop_over_X_scales']) if passr else None,
       feasible_at_M=sum(x['loop_over_X_at_M']<=1. for x in passr),
       feasible_across_scales=sum(x['max_loop_over_X_scales']<=1. for x in passr),
       max_tree_plus_eta_fraction_zle100=emax)

def direct_recombination(row,p,rr,b):
    a=b['a'];g=np.sqrt(b['g2_eV2']);ms=np.sqrt(4*M*M-a*SCALE)
    z=zero_rate(1.,g,3*g*g/ms**2,M,M+1.,ms)
    # amplitude proportional to x^3 at fixed masses; cold zero-momentum
    # cosmic gain includes final A occupancy sampled at recoil q.
    r,t,x=rr['r'],rr['t'],rr['x'];l=r<=101
    r,t,x=r[l],t[l],x[l]
    V,Tr=variance_table(row);mom=z['recoil_momentum_eV']
    edges=np.asarray(row['frozen_cell_spectrum']['edges']);f=np.maximum(row['frozen_cell_spectrum']['f_phi'],0.)
    j=np.searchsorted(edges,mom/(Tr*r))-1;valid=(j>=0)&(j<len(f));ff=np.zeros_like(r);ff[valid]=f[j[valid]]
    n=RHO0*CGS_RHO_TO_EV4/M*r**3
    rate=z['coefficient']*x**6*n*n/(HBAR*H(r-1))*(1+ff)
    fS=2*np.pi**2*rate*n/mom**3
    return dict(zero_L=z,cosmic_S_fraction=float(trapezoid(rate,t)),max_cosmic_line_fS=float(max(fS)),
                max_final_A_f=float(max(ff)),halo_zero_S_fraction_10Gyr=float(z['coefficient']*(RHO*CGS_RHO_TO_EV4/M)**2*10*SEC_GYR/HBAR),
                scope='Cold zero-incoming-momentum initial source; no inverse, S stimulation, or hot tail.')

def scan(row,aa,QQ,output,group):
    cases=[]
    for z in [8,12,16,20,25,30]:
        for dn in [.1,.2]:
            for xf,muH in [(.5,10),(.7,10),(.7,20),(.9,10)]:
                cases.append((z,dn,xf,muH))
    out=[]
    for i,(z,dn,xf,muH) in enumerate(cases):
        if i%4!=group:continue
        p=params(z,dn,xf,1.,muH)
        p['delta']=0.
        for scheme in ['origin-vev','vev-curvature']:
            rr=solve(row,p,scheme,n=8000,rtol=3e-7)
            q=diagnose(row,p,rr,aa,QQ)
            rec=dict(params=p,scheme=scheme,roll={k:v for k,v in rr.items() if k not in ['t','r','x','dx']},screen=q)
            out.append(rec)
            Path(output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
            print('DIRECT',z,dn,xf,muH,scheme,q.get('z_half'),q.get('feasible_at_M'),q.get('feasible_across_scales'),q.get('best'),flush=True)
    return out

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--tradeoff',required=True)
    ap.add_argument('--output',required=True);ap.add_argument('--group',type=int,required=True);a=ap.parse_args()
    row=json.loads(Path(a.input).read_text())[0];q=json.loads(Path(a.tradeoff).read_text())
    scan(row,np.array(q['aa']),np.array(q['Q96']),a.output,a.group)
