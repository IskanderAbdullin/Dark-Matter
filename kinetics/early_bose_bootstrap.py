"""Conditional first-exit bounds for early heavy-state Bose seeding.

Before a specified portal-off epoch, B^4 acts on the B flavor only; allow
arbitrary diagonal thermal dispersion, coherent mixing and diagonal A source.
Afterwards assume the vacuum eigenbasis and all six reversible quartic
channels. No direct B source except mixing/scattering. We bound both
integrated heavy number and maximum mode occupancy; the latter is essential
because a tiny number fraction alone does not exclude Bose stimulation.

Bounds use phase-space rearrangement at fixed number and occupancy cap,
not a frozen A momentum shape. They require a positive quasiparticle
kinetic description with diagonal populations for collisions, the stated
mass floor, number/source caps and homogeneous prescribed chi history.
They do not cover off-diagonal collision kernels near thermal crossings,
condensates/coherent B fields, or arbitrary primordial B.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid,trapezoid,quad
from mixing_crossing_bound import early_bound,solve_roll
from soft_mixing_activation import scales,frac_B
from late_activation_overlap import H,H0,Or
from direct_activation_consistency import CGS_RHO_TO_EV4
from resonant_mediator_inverse import M,HBAR
from velocity_only_exposure import RHO0

PAIRS=[(0,0),(0,1),(1,1)]
N0=RHO0*CGS_RHO_TO_EV4/M

def rearrange(n,F,m):
    """Upper I=int f dE and V=int d^3p f/((2pi)^3 E).

    Fill a sphere at f=F, because the integrand per particle decreases
    with momentum. I and V bounds use the smaller NR/UR bound, avoiding
    loss of precision in the exact V integral at p<<m.
    """
    n=np.maximum(n,0.);F=max(float(F),1e-300)
    p=np.cbrt(6*np.pi**2*n/F)
    I=F*p*p/(np.hypot(p,m)+m)
    V=np.minimum(n/m,F*p*p/(4*np.pi**2))
    return I,V

def ci(y,t):
    return cumulative_trapezoid(y,t,initial=0.)

def pre_portal(row,p,seed,cap,FA,rcut,n,m):
    eb=early_bound(row,p,seed,rcut=rcut)
    maxchi=eb['assumed_seed_envelope']
    hc=.5*H0*np.sqrt(Or)*HBAR
    # Analytic upper amplitude at every r from infinity; worst p=0.
    # n_A <= 2*N0*r^3 and |H_AB| <= delta chi/m.
    A0=abs(p['delta'])*maxchi/(2*hc*m)
    t=np.linspace(-np.log(rcut)-35.,-np.log(rcut),n);r=np.exp(-t)
    amp=A0/r**2;y=2*amp**2
    I,V=rearrange(y*N0*r**3,cap,m)
    # Exact Bose quartic cancellation. Per fixed external mode, angular
    # interval <=2p, incoming distinct-pair factor <=2, E>=m.
    # HH->HH is the only diagonal B collision before portal-off.
    G=p['lambda_B']**2*(1+2*cap)*I*I/(64*np.pi**3*m*hc*r*r)
    integrated=float(trapezoid(G,t))
    # Missing r=inf tail: I <= F^(2/3)*(6pi^2*n)^(1/3),
    # n_B <= const/r, so G <= const*r^(-8/3).
    pstar=np.cbrt(6*np.pi**2*y[0]*N0*r[0]**3/cap)
    Gfirst=p['lambda_B']**2*(1+2*cap)*(cap*pstar)**2/(64*np.pi**3*m*hc*r[0]**2)
    integrated+=3*Gfirst/8
    # If F'<=2 a sqrt(F)+G, a comparison solution is
    # F <= (int a + sqrt(int G))^2, even if G depends on the cap.
    fmax=(np.sqrt(FA)*amp[-1]+np.sqrt(integrated))**2
    return dict(maxchi=maxchi,heavy_fraction_bound=float(y[-1]),mode_occupation_bound=float(fmax),
       pre_collisional_mode_gain=integrated,pre_cap_closed=bool(fmax<cap),missing_tail_bound=float(3*Gfirst/8),
       scope='B-flavor diagonal B^4 collisions only; arbitrary diagonal detuning; A number and occupation bounded; no primordial B.')

def gain_mode(s,lam,u,I,V,F,m,h):
    # Incoming unordered species pair i,j. Sum over possible outgoing
    # partner k for fixed external species s. Identical final particle
    # multiplicities already absorbed by choosing conservative coefficient.
    out=np.zeros_like(h)
    for i,j in PAIRS:
        for k in [0,1]:
            nh=i+j+s+k
            amp=lam**2*u**(4-nh)*(1-u)**nh
            out+=amp*I[i]*I[j]*(1+F[s]+F[k])
    return out/(64*np.pi**3*m*h)

def post_portal(row,p,sol,seed,cap,FA,rcut,n,m,pre):
    t=np.linspace(-np.log(rcut),-np.log(1001.),n);r=np.exp(-t)
    x,dx=sol.sol(t);u=frac_B(x,p['delta'],p['Delta'])
    g=abs(p['delta']*p['Delta']*dx/(p['Delta']**2+4*p['delta']**2*x*x))
    h=.5*HBAR*H(r-1)
    nr=N0*r**3;IA,VA=rearrange(2*nr,FA,m)
    # Project B flavor into the vacuum heavy state with triangle inequality.
    entry=(np.sqrt(pre['heavy_fraction_bound'])+np.sqrt(2*u[0]))**2
    fentry=(np.sqrt(pre['mode_occupation_bound'])+np.sqrt(FA*u[0]))**2
    # Number source AA->AH and AA->HH. Ignore threshold, reverse flux,
    # and use conservative symmetry coefficient (distinct final maximum).
    lam=p['lambda_B']
    amp1=lam**2*u**3*(1-u);amp2=lam**2*u*u*(1-u)**2
    amp3=lam**2*u*(1-u)**3
    S=(amp1*(1+FA+cap)+2*amp2*(1+2*cap))*VA*VA/(32*np.pi*nr*h)
    # AH->HH catalytic increase. V_H <= n_H/m for the global bound.
    C=amp3*VA*(1+2*cap)/(32*np.pi*m*h)
    intS=ci(S,t);intC=ci(C,t);intg=ci(g,t)
    y=np.exp(np.minimum(intC,700))*(np.sqrt(entry)+np.sqrt(2)*intg+np.sqrt(intS))**2
    # Even if the number envelope grows, number conservation caps it at 2.
    y=np.minimum(y,2.)
    IH,VH=rearrange(y*nr,cap,m)
    F=[FA,cap];I=[IA,IH];V=[VA,VH]
    GH=gain_mode(1,lam,u,I,V,F,m,h)
    GA=gain_mode(0,lam,u,I,V,F,m,h)
    gainH=float(trapezoid(GH,t));gainA=float(trapezoid(GA,t))
    fH=(np.sqrt(fentry)+np.sqrt(FA)*intg[-1]+np.sqrt(gainH))**2
    # Initial A cap assumed <= FA/2 before the vacuum projection.
    # Bound the projected A amplitude and all positive collision gains.
    aentry=(np.sqrt(FA/2)+np.sqrt(pre['mode_occupation_bound']*u[0]))**2
    fA=(np.sqrt(aentry)+np.sqrt(cap)*intg[-1]+np.sqrt(gainA))**2
    # Bound L1 change of A occupations relative to Nref due to all
    # diagonal 2->2 channels. Sum ordered incoming pair/outgoing pairs.
    # Each event changes at most four A occupations; omit reverse by
    # including both orientations. This is conservative, not net relaxation.
    Gtot=np.zeros_like(r)
    for i,j in PAIRS:
        for k,l in PAIRS:
            nh=i+j+k+l
            a=lam**2*u**(4-nh)*(1-u)**nh
            Gtot+=a*V[i]*V[j]*(1+F[k]+F[l])/(32*np.pi*nr*h)
    # A-B coherent number transfer contributes additional L1 change.
    l1_coll=4*float(trapezoid(Gtot,t))
    l1_mix=2*float(trapezoid(g*np.sqrt(2*y),t))
    return dict(entry_fraction_bound=float(entry),entry_occupation_bound=float(fentry),
       max_heavy_fraction_bound=float(max(y)),end_heavy_fraction_bound=float(y[-1]),
       integrated_direct_number_source=float(intS[-1]),integrated_catalytic_exponent=float(intC[-1]),
       integrated_rotation=float(intg[-1]),heavy_mode_gain_upper=gainH,light_mode_gain_upper=gainA,
       max_heavy_occupation_bound=float(fH),max_light_occupation_bound=float(fA),
       heavy_cap_closed=bool(fH<cap),light_cap_closed=bool(fA<FA),
       collision_L1_over_Nref_upper=l1_coll,coherent_L1_over_Nref_upper=l1_mix,
       min_r=1001.,max_r=rcut,
       max_B_Hartree_mass2_over_m2=float(max(lam*VH/(2*m*m))),
       max_B_Hartree_mass2_over_vacuum_split=float(max(lam*VH/(2*p['Delta']))))

def run(row,p,scheme,seed,n,cap,sol=None):
    if sol is None:sol,_=solve_roll(row,p,scheme,seed)
    fmax=float(max(row['frozen_cell_spectrum']['f_phi']))
    FA=2.1*fmax
    m=.99*M
    rcut=1e14
    pre=pre_portal(row,p,seed,cap,FA,rcut,n,m)
    post=post_portal(row,p,sol,seed,cap,FA,rcut,n,m,pre)
    return dict(seed=seed,cap_H=cap,cap_A=FA,mass_floor_eV=m,n=n,params=p,scheme=scheme,
        pre=pre,post=post,first_exit_caps_closed=bool(pre['pre_cap_closed'] and post['heavy_cap_closed'] and post['light_cap_closed']),
        assumptions='Pre-portal A cap externally assumed, B diagonal B^4 collisions; post-portal vacuum mixing, all diagonal quartic channels, equal mass angular bound, total number <=2 Nref, A entry <= cap_A/2; H lower bound 0.5 LCDM. No direct primordial B, coherent condensate, offdiagonal collision term, arbitrary thermal mass corrections or 2->4.')

def main():
    a=argparse.ArgumentParser();a.add_argument('--input',required=True);a.add_argument('--output',required=True)
    a.add_argument('--sigma',nargs='+',type=float,default=[.1,.3,1.25]);a.add_argument('--seeds',nargs='+',type=float,default=[1e-8,1e-7,1e-6,1e-5]);a.add_argument('--n',type=int,default=8000)
    args=a.parse_args();row=json.loads(Path(args.input).read_text())[0];out=[]
    for s in args.sigma:
        z=5. if s==.1 else 8.;mu=100. if s==.1 else (30. if s==.3 else 60.)
        p0=scales(s,z_on=z,dneff=.1,muH=mu)
        p=scales(s,z_on=z,dneff=.1,muH=mu,X=1.05*p0['X_loop_floor'])
        for seed in args.seeds:
            sol,_=solve_roll(row,p,'vev-curvature',seed)
            for cap in [1e-5,1e-4,1e-3,.01,.1,1.]:
                result=run(row,p,'vev-curvature',seed,args.n,cap,sol)
                out.append(result)
                print(s,seed,cap,result['first_exit_caps_closed'],result['pre']['mode_occupation_bound'],result['post']['max_heavy_occupation_bound'],result['post']['end_heavy_fraction_bound'],result['post']['collision_L1_over_Nref_upper'],flush=True)
                if result['first_exit_caps_closed']:break
            Path(args.output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
if __name__=='__main__':main()
