"""Conditional number-conserving Hartree quasiparticle gas + zero mode.

Real phi^4/4! scalar; assumes a saturated effective Bose chemical potential
equal to its quasiparticle gap, NOT a conserved exact U(1) or a derived
coherent-condensate dispersion. Vacuum-subtracted Hartree functional:
U = integral E f + m_eff N0 - Delta^2/(2 lambda).
Stationarity in Delta=m_eff^2-m_vac^2 yields the gap equation.
P = kinetic pressure + Delta^2/(2 lambda).
Then U-3P=m_vac^2*(N0/m_eff+integral f/E). All in Tref units.
The equilibrium limits are conditional, not time-dependent nucleation.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from scipy.optimize import brentq
from scipy.interpolate import PchipInterpolator
from scipy.integrate import cumulative_trapezoid
from quartic_consistency import shell_nodes,lambda_from_sigma,gap,variance

def thermal(T,m,order=96):
    if T<=0:return dict(n=0.,kin=0.,p=0.,V=0.,D=0.,s=0.)
    z=m/T
    u,w=roots_legendre(order)
    v=16*(u+1)/2;wv=8*w
    scale=np.sqrt(1+z);y=scale*v
    E=np.hypot(y,z);K=y*y/(E+z)
    f=np.exp(-K)/(-np.expm1(-K))
    weight=wv*scale*y*y*f/(2*np.pi**2)
    n=T**3*weight.sum()
    kin=T**4*(weight@K)
    p=T**4*(weight@(y*y/(3*E)))
    V=T**2*(weight@(1/E))
    D=T**3*(weight@(K/E))
    return dict(n=float(n),kin=float(kin),p=float(p),V=float(V),D=float(D),s=float((kin+p)/T))

def saturation(N,T,m0,lam,order=96):
    if lam==0:
        m=m0;delta=0.;th=thermal(T,m,order)
    else:
        upper=(lam*N/2)**(2/3)
        if m0>0:upper=min(upper,lam*N/(2*m0))
        def eval_gap(u):
            m=np.sqrt(m0*m0+u*upper)
            if m==0:return -1e300
            th=thermal(T,m,order)
            return u-(lam/2)*(N-th['D'])/(m*upper)
        if eval_gap(0)>0:raise ValueError('No saturated gap with specified number and temperature')
        root=brentq(eval_gap,0.,1+1e-12,xtol=1e-14,rtol=1e-12)
        delta=upper*root;m=np.sqrt(m0*m0+delta)
        th=thermal(T,m,order)
    B=delta*delta/(2*lam) if lam else 0.
    N0=N-th['n'];V=N0/m+th['V'] if m>0 else 0.
    U=m*N+th['kin']-B
    P=th['p']+B
    return dict(T=T,m_eff=m,delta=delta,N0=N0,excited_fraction=th['n']/N,U=U,P=P,
                entropy=th['s'],V=V,mean_field_B=B,
                trace_residual=float((U-3*P-m0*m0*V)/max(U,1e-300)))

def initial_U(edges,f,m0,lam,order=8):
    q,w=shell_nodes(edges,f,order)
    m=gap(q,w,m0,lam)
    V=variance(q,w,m)
    B=lam*V*V/8
    return dict(N=float(w.sum()),U=float(w@np.hypot(q,m)-B),m_eff=float(m),mean_field_B=float(B),
                U_vac=float(w@np.hypot(q,m0)),lambda_f_peak=float(lam*max(f)))

def match(N,U,m0,lam,order=96):
    ground=saturation(N,0.,m0,lam,order)
    if U<ground['U']*(1-1e-10):raise ValueError('Energy below saturated Hartree ground state')
    T0=(30*U/np.pi**2)**.25
    def residual(T):return saturation(N,T,m0,lam,order)['U']-U
    hi=T0*1.5
    try:
        while residual(hi)<0:hi*=2
        T=brentq(residual,0,hi,xtol=1e-13,rtol=1e-12)
    except ValueError:
        return dict(saturated=False,reason='Saturated branch unavailable in root bracket')
    st=saturation(N,T,m0,lam,order)
    if st['N0']<0:return dict(saturated=False,reason='Thermal capacity exceeds N',trial_excited_fraction=st['excited_fraction'])
    st.update(saturated=True,ground_energy_fraction=ground['U']/U,energy_match_relative=st['U']/U-1,
        condensate_fraction=st['N0']/N,entropy_per_particle=st['entropy']/N)
    # First law at fixed particle count and bare mass.
    eps=1e-4
    lo=saturation(N,T*(1-eps),m0,lam,order);hi=saturation(N,T*(1+eps),m0,lam,order)
    st['first_law_relative_error']=float((hi['U']-lo['U'])/(T*(hi['entropy']-lo['entropy']))-1)
    assert abs(st['trace_residual'])<2e-9 and abs(st['energy_match_relative'])<1e-8
    assert abs(st['first_law_relative_error'])<1e-5
    return st

def adiabatic(row,N,S,lam,ngrid=160,order=96):
    Tgamma0=2.7255*8.617333262e-5
    theta=row['final']['theta'];tr=(3.909/106.75)**(1/3)
    Tref0=Tgamma0*tr/theta
    a=np.geomspace(Tgamma0/30000,1.,ngrid)
    m0=86*a/Tref0
    states=[]
    Tguess=(45*S/(2*np.pi**2))**(1/3)
    for ai,m in zip(a,m0):
        def f(T):return saturation(N,T,m,lam,order)['entropy']-S
        hi=max(Tguess*2,1e-30)
        while f(hi)<0:hi*=2
        T=brentq(f,0,hi,xtol=max(hi*1e-13,1e-30),rtol=1e-12)
        st=saturation(N,T,m,lam,order)
        assert st['N0']>=0
        states.append(st);Tguess=T
    U=np.array([v['U'] for v in states]);P=np.array([v['P'] for v in states])
    # Derivatives along the comoving adiabatic saturated branch.
    ln=np.log(a)
    du=PchipInterpolator(ln,np.log(U)).derivative()(ln)
    dp=PchipInterpolator(ln,np.log(P)).derivative()(ln)
    cont=du-(U-3*P)/U
    c2=P*(4-dp)/(3*(U+P))
    assert np.max(abs(cont))<.01,np.max(abs(cont))
    assert np.min(c2)>0 and np.max(c2)<1/3+.01
    h=.674;H0=100*h;Ob=.02237/h**2;Od=row['final']['omega_phi_h2_m86']/h**2
    Or=2.4728e-5/h**2*(1+.22710731766*3.044);Ol=1-Ob-Od-Or
    H=H0*np.sqrt(Or/a**4+Ob/a**3+Od*U/a**4/U[-1]+Ol)
    ae=Or/(Ob+Od)
    af=np.geomspace(a[0],1.,6000);cf=np.exp(PchipInterpolator(ln,np.log(c2))(np.log(af)))
    HF=np.exp(PchipInterpolator(ln,np.log(H))(np.log(af)))
    r=cumulative_trapezoid(299792.458*np.sqrt(cf)/(af*af*HF),af,initial=0.)
    eq=float(np.interp(np.log(ae),np.log(af),r))
    sample=[]
    for z in (1e8,1e7,1e6,1e5,1e4,3400,0):
        i=np.argmin(abs(ln-np.log(1/(1+z))))
        sample.append(dict(z=float(1/a[i]-1),m0=float(m0[i]),w=float(P[i]/U[i]),cs2=float(c2[i]),**states[i]))
    return dict(ngrid=ngrid,order=order,maximum_continuity_residual=float(max(abs(cont))),
        sound_horizon_equality_Mpc=eq,sound_horizon_today_Mpc=float(r[-1]),
        today_excited_fraction=states[-1]['excited_fraction'],samples=sample,
        scope='Idealized self-consistent Hartree saturated number-conserving single fluid; no true condensate excitation dispersion, finite kinetics, or chemical evolution.')

def main():
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--output',required=True)
    p.add_argument('--sigmas',nargs='+',type=float,default=[1e-7,.001,.1,1.25])
    p.add_argument('--adiabatic',action='store_true');a=p.parse_args()
    row=json.loads(Path(a.input).read_text())[0];e=np.array(row['frozen_cell_spectrum']['edges'])
    hist=row.get('scalar_history_psd',[dict(x=row['final']['x'],f_phi=row['frozen_cell_spectrum']['f_phi'])])
    results=[]
    for sigma in a.sigmas:
        lam=float(lambda_from_sigma(sigma));states=[]
        for target in [.175,.2,.5,2.,80.]:
            i=int(np.argmin([abs(s['x']-target) for s in hist]))
            if any(s['x']==hist[i]['x'] for s in states):continue
            h=hist[i];f=np.array(h['f_phi']);m0=1e-12*h['x']
            ini=initial_U(e,f,m0,lam);r=match(ini['N'],ini['U'],m0,lam)
            states.append(dict(x=h['x'],initial=ini,equilibrium=r))
        out=dict(sigma=sigma,lambda_quartic=lam,states=states)
        if a.adiabatic and states[-1]['equilibrium']['saturated']:
            N=states[-1]['initial']['N'];S=states[-1]['equilibrium']['entropy']
            ad=[adiabatic(row,N,S,lam,n) for n in [96,192]]
            out['adiabatic']=ad
            out['sound_grid_relative_drift']=abs(ad[0]['sound_horizon_equality_Mpc']/ad[1]['sound_horizon_equality_Mpc']-1)
        results.append(out)
        Path(a.output).write_text(json.dumps(results,indent=2,allow_nan=False)+'\n')
        print('SIGMA',sigma,'states',json.dumps([dict(x=s['x'],gap=s['initial']['m_eff'],equilibrium=s['equilibrium']) for s in states]),flush=True)
        if 'adiabatic' in out:print('SOUND',sigma,out['adiabatic'][-1]['sound_horizon_equality_Mpc'],out['sound_grid_relative_drift'],flush=True)
if __name__=='__main__':main()
