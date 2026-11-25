"""Leading-density real forward-scattering screen for the derivative portal.

Massless external A, tree s/u intermediate off-diagonal fermion, spin averaged:
 M=-8 D (p.k)^2/[F^2(D^2-4(p.k)^2)], D=Mi^2-Mj^2.
Principal value at the decay/inverse-decay pole. For a dilute gas the
on-shell mass-squared diagnostic is -integral f M d^3p/[(2pi)^3 2E].
No self-consistent thermal spectral function or imaginary-part evolution.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.special import expit,roots_legendre
from derivative_portal_consistency import M,R,K,rate_mapping
from soft_mixing_activation import scales

def amp_primitive(y):
    y=np.asarray(y);out=np.empty_like(y,dtype=float)
    a=y<.02;b=~a
    t=y[a];out[a]=t**3*(1/3+t*t*(1/5+t*t*(1/7+t*t*(1/9+t*t/11))))
    v=y[b]
    # Exact PV primitive of y^2/(1-y^2); log singular only at 1.
    out[b]=.5*(np.log1p(v)-np.log(np.maximum(abs(1-v),1e-300)))-v
    return out

def angle_amplitude(k,p,mi,mj,F):
    k,p=np.broadcast_arrays(np.asarray(k,float),np.asarray(p,float))
    E=np.hypot(mi,p);D=mi*mi-mj*mj;d=abs(D)
    yl=2*k*mi*mi/(d*(E+p));yh=2*k*(E+p)/d
    wid=yh-yl
    out=np.empty_like(k);small=wid<1e-8*np.maximum(1.,yh)
    mid=(yl[small]+yh[small])/2
    out[small]=mid*mid/(1-mid*mid)
    out[~small]=(amp_primitive(yh[~small])-amp_primitive(yl[~small]))/wid[~small]
    return -2*D/F**2*out

def fd_pi(k,T,mi,mj,F,g=4,order_rtol=2e-8):
    D=abs(mi*mi-mj*mj);s=D/(2*k)
    root=abs(s*s-mi*mi)/(2*s*T)
    def ff(q):
        E=np.hypot(mi,T*q)
        return q*q*expit(-E/T)*float(angle_amplitude(k,T*q,mi,mj,F))/E
    points=[root] if 0<root<60 else None
    ans,err=quad(ff,0,60,points=points,epsabs=1e-40,epsrel=order_rtol,limit=200)
    pi=-g*T**3/(4*np.pi**2)*ans
    N=g/(2*np.pi**2)*quad(lambda q:q*q*expit(-np.hypot(mi/T,q)),0,60,epsabs=1e-40,epsrel=1e-9)[0]
    return pi,N,abs(g*T**3/(4*np.pi**2)*err)

def fermion_pi(e,f,Tref,k,mi,mj,F,g,order=32):
    # Piecewise-constant actual heavy distribution. Split the shell with
    # the log endpoint singularity, and refine order separately.
    D=abs(mi*mi-mj*mj);t=D/(2*k)
    root=abs(t*t-mi*mi)/(2*t*Tref)
    e=np.asarray(e);f=np.maximum(f,0.)
    ed=np.unique(np.r_[e,[root] if e[0]<root<e[-1] else []])
    j=np.clip(np.searchsorted(e,(ed[1:]+ed[:-1])/2)-1,0,len(f)-1)
    u,w=roots_legendre(order);l,h=ed[:-1],ed[1:]
    q=(l[:,None]+h[:,None])/2+(h-l)[:,None]*u/2
    a=angle_amplitude(k,q*Tref,mi,mj,F)
    return float(-g*Tref**3/(4*np.pi**2)*np.sum(f[j]*(h-l)/2*((q*q*a/np.hypot(mi,q*Tref))@w)))

def scalar_pi(e,f,Tref,p,mi,mj,F,order=24):
    E=np.hypot(mi,p);D=abs(mi*mi-mj*mj)
    roots=[D/(2*(E+p)*Tref),D*(E+p)/(2*mi*mi*Tref)]
    e=np.asarray(e);f=np.maximum(f,0.)
    ed=np.unique(np.r_[e,[r for r in roots if e[0]<r<e[-1]]])
    j=np.clip(np.searchsorted(e,(ed[1:]+ed[:-1])/2)-1,0,len(f)-1)
    u,w=roots_legendre(order);l,h=ed[:-1],ed[1:]
    q=(l[:,None]+h[:,None])/2+(h-l)[:,None]*u/2
    a=angle_amplitude(q*Tref,p,mi,mj,F)
    # A massless, 1 dof, f in Tref units.
    return float(-Tref*Tref/(4*np.pi**2)*np.sum(f[j]*(h-l)/2*((q*a)@w)))

def forward_trace(k,p,mi,mj,F,axial):
    # Independent spin trace and s/u diagrams, away from poles.
    I=np.eye(2);Z=np.zeros((2,2))
    pauli=[np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1])]
    g=[np.block([[I,Z],[Z,-I]])]+[np.block([[Z,s],[-s,Z]]) for s in pauli]
    g5=1j*g[0]@g[1]@g[2]@g[3];I4=np.eye(4)
    slash=lambda p:sum(v*x for v,x in zip(np.asarray(p)*[1,-1,-1,-1],g))
    a=slash(k)@(g5 if axial else I4)/F;sp=slash(p);ans=0.
    for sg in [1,-1]:
        q=np.asarray(p)+sg*np.asarray(k);den=q[0]**2-np.dot(q[1:],q[1:])-mj*mj
        ans-=np.trace((sp+mi*I4)@a@(slash(q)+mj*I4)@a).real/(2*den)
    return ans

def check_trace():
    out=[]
    for mi,mj in [(1.,.9),(.9,1.)]:
        for k in [.001,.03,.3]:
            p=[np.hypot(mi,.4),0,0,.4];kv=[k,k,0,0]
            pk=k*p[0];D=mi*mi-mj*mj
            expected=-8*D*pk*pk/(D*D-4*pk*pk)
            a=[forward_trace(kv,p,mi,mj,1.,x) for x in [False,True]]
            out.append(float(max(abs(np.array(a)/expected-1))))
    return max(out)

def screen(row,order=24):
    M1=row['parent_GeV'];M2=row['fermion_ratio']*M1
    F=2*rate_mapping(M1,row['fermion_ratio'],row['kappa'],row['kappa2'],row['kappa_sm1'],row['g_bath'])['f_over_C_GeV']
    hist=row['history'];p0=row['final']['peak_bin_q'];g=row['g_fermion'];e=row['frozen_cell_spectrum']['edges']
    states=[];D=M1*M1-M2*M2
    split=scales()['Delta']*1e-18 # eV^2->GeV^2
    for h,s in zip(hist,row['scalar_history_psd']):
        x=h['x'];Tr=M1/x;T=Tr*h['theta'];k=p0*Tr
        a,Na,ea=fd_pi(k,T,M1,M2,F,g);b,Nb,eb=fd_pi(k,T,M2,M1,F,g)
        # Shape surrogate only: rescale equilibrium by actual number density.
        r1=max(h['mother_number'],0.)/(Na*h['theta']**3)
        r2=max(h['fermion_number'],0.)/(Nb*h['theta']**3)
        pv=a*r1+b*r2
        actual=None
        if 'f_parent' in s:
            e1=row['frozen_cell_spectrum']['parent_edges'];e2=row['frozen_cell_spectrum']['chi2_edges']
            p1=fermion_pi(e1,s['f_parent'],Tr,k,M1,M2,F,g,order)
            p2=fermion_pi(e2,s['f_chi2'],Tr,k,M2,M1,F,g,order)
            p1f=fermion_pi(e1,s['f_parent'],Tr,k,M1,M2,F,g,2*order)
            p2f=fermion_pi(e2,s['f_chi2'],Tr,k,M2,M1,F,g,2*order)
            actual=dict(Pi1_GeV2=p1f,Pi2_GeV2=p2f,Pi_total_GeV2=p1f+p2f,
              Pi_total_over_k2=(p1f+p2f)/k**2,
              Pi_total_over_B_split=(p1f+p2f)/split,
              order_absolute_error_over_k2=(abs(p1f-p1)+abs(p2f-p2))/k**2)
        ferm=[]
        for pt in [1.,3.,6.]:
            p=pt*T
            pa=scalar_pi(e,s['f_phi'],Tr,p,M1,M2,F,order)
            pb=scalar_pi(e,s['f_phi'],Tr,p,M2,M1,F,order)
            qa=scalar_pi(e,s['f_phi'],Tr,p,M1,M2,F,2*order)
            qb=scalar_pi(e,s['f_phi'],Tr,p,M2,M1,F,2*order)
            ferm.append(dict(p_over_T=pt,Pi1_GeV2=qa,Pi2_GeV2=qb,
                splitting_change_over_vacuum=(qa-qb)/D,
                absolute_quadrature_difference_over_gap=max(abs(pa-qa),abs(pb-qb))/D))
        states.append(dict(x=x,T_GeV=T,peak_k_GeV=k,Pi_A_eq_1_GeV2=a,Pi_A_eq_2_GeV2=b,
          Pi_A_number_rescaled_GeV2=pv,Pi_A_number_rescaled_over_k2=pv/k**2,
          Pi_A_eq_over_k2=(a+b)/k**2,Pi_A_number_rescaled_over_soft_B_split=pv/split,
          parent_number_over_FD=r1,daughter_number_over_FD=r2,
          eq_quad_error_GeV2=ea+eb,actual_fermion_PSD=actual,fermion_from_actual_scalar_PSD=ferm))
    # Locate sign/level crossings of the number-rescaled shape surrogate.
    crossings=[]
    for a,b in zip(states[:-1],states[1:]):
        d1=(a['actual_fermion_PSD']['Pi_total_GeV2'] if a['actual_fermion_PSD'] else a['Pi_A_number_rescaled_GeV2'])-split
        d2=(b['actual_fermion_PSD']['Pi_total_GeV2'] if b['actual_fermion_PSD'] else b['Pi_A_number_rescaled_GeV2'])-split
        if d1*d2<0:
            # Only secant diagnostic across sparse x snapshots, not a bound
            # on conversion (real thermal shapes and damping can differ).
            l1,l2=np.log(a['x']),np.log(b['x']);slope=(d2-d1)/(l2-l1)
            lx=l1-d1/slope;xx=np.exp(lx);E=p0*M1/xx
            HH=1.66*np.sqrt(row['g_bath'])*(M1/xx)**2/1.220890e19
            off=scales()['delta']*1e-8*1e-18
            crossings.append(dict(x_secant=xx,linearized_crossing_slope_GeV2=slope,
               LZ_exponent_secant=np.pi*off**2/(E*HH*abs(slope))))
    return dict(trace_relative_error=check_trace(),F_GeV=F,scalar_pv_order=order,states=states,
       level_crossing_secant_diagnostics=crossings,
       scope='Vacuum massless derivative tree forward amplitude and leading-density real PV shifts. Actual recorded fermion PSD used when present; otherwise number-rescaled FD surrogate for scalar shift. Actual scalar PSD for fermion shift. No thermal widths, dressed denominators, finite-T spin kinetics, coherence, or nonlinear corrections. Secant LZ is not a conversion bound.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--output',required=True);p.add_argument('--order',type=int,default=24);a=p.parse_args()
    r=screen(json.loads(Path(a.input).read_text())[0],a.order)
    Path(a.output).write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
    print('trace',r['trace_relative_error'],'crossings',r['level_crossing_secant_diagnostics'])
    for v in r['states']:
        print(v['x'],v['actual_fermion_PSD'],v['Pi_A_number_rescaled_over_k2'],v['Pi_A_number_rescaled_over_soft_B_split'],max(abs(f['splitting_change_over_vacuum']) for f in v['fermion_from_actual_scalar_PSD']),max(f['absolute_quadrature_difference_over_gap'] for f in v['fermion_from_actual_scalar_PSD']),flush=True)
