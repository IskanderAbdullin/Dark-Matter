"""Exact homogeneous chi + collisionless eta energy budget for direct switch.
Use the same prescribed LCDM H as the field equation and check continuity.
The two-loop A/S term is already absorbed into renormalized chi parameters;
its magnitude is a separate matching sensitivity, not counted twice.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from scipy.integrate import cumulative_trapezoid
from direct_subthreshold_activation import solve,diagnose
from subthreshold_roll_scan import params
from late_activation_overlap import H
from resonant_mediator_inverse import HBAR
from direct_activation_consistency import rho_tot

def energy(rr,p,scheme,nq=32):
    r,t,x,dx=[rr[k] for k in ['r','t','x','dx']]
    l=r<=101
    r,t,x,dx=r[l],t[l],x[l],dx[l]
    td=p['Td_on_eV']*r/(p['z_on']+1)
    mass=p['m_eta_at_vev_eV']*abs(x);j=mass/td
    u,w=roots_legendre(nq);e=np.array([0.,1e-5,1e-4,.001,.01,.1,1.,5.,15.,50.])
    q=((e[:-1,None]+e[1:,None])/2+np.diff(e)[:,None]*u/2).ravel()
    ww=(np.diff(e)[:,None]*w/2).ravel()
    f=q*q*ww/np.expm1(q)
    E=np.hypot(q[None,:],j[:,None])
    norm=p['Neta']*td**4/(2*np.pi**2)
    rhoeta=norm*(E@f);peta=norm*((1/E)@(f*q*q))/3
    X=p['X'];xx=x*x;log=np.log(np.maximum(xx,1e-300))
    cw=p['eta_CW_over_X']
    V=X*(.25*(xx-1)**2+cw*(xx*xx*log+(-.5*xx*xx+.5 if scheme=='origin-vev' else -1.5*xx*xx+2*xx-.5)))
    K=.5*X*(HBAR*H(r-1)*dx/p['mu_eV'])**2
    rho=K+V+rhoeta;pr=K-V+peta
    # comoving radiation-scaled continuity check for the stipulated
    # collisionless eta distribution and the specified CW potential.
    R=rho/r**4;target=R[0]+cumulative_trapezoid((rho-3*pr)/r**4,t,initial=0.)
    i=np.argmax(rho/rho_tot(r-1))
    return dict(max_fraction=float(rho[i]/rho_tot(r[i]-1)),redshift_of_max=float(r[i]-1),
        max_abs_comoving_continuity_residual=float(max(abs(R-target))/max(abs(R))),
        today_chi_fraction=float((K[-1]+V[-1])/rho_tot(0.)),today_eta_fraction=float(rhoeta[-1]/rho_tot(0.)),
        min_potential_over_X=float(min(V)/X),initial_eta_radiation_fraction=float(rhoeta[0]/rho_tot(r[0]-1)),
        quadrature_order=nq,scope='Exact homogeneous chi+eta energy at prescribed H, no A/S medium energy or metric perturbations.')

def run(row,table,group,output):
    # Focus on near-margin cases. Vev-curvature matching was found more
    # favorable; test nearby X and slower field response explicitly.
    cases=[(z,dn,xf,mu) for z,dn in [(30,.1),(30,.2),(35,.1),(35,.2),(40,.1),(40,.2)]
           for xf,mu in [(.65,7),(.75,7),(.7,9),(.8,9)]]
    aa=np.array(table['aa']);QQ=np.array(table['Q96']);out=[]
    for i,(z,dn,xf,mu) in enumerate(cases):
        if i%4!=group:continue
        p=params(z,dn,xf,1.,mu);p['delta']=0.
        rr=solve(row,p,'vev-curvature',n=12000,rtol=1e-7)
        q=diagnose(row,p,rr,aa,QQ)
        en=energy(rr,p,'vev-curvature')
        for b in q.get('rows',[]):
            b['exact_chi_eta_fraction']=en['max_fraction']
        candidates=[b for b in q.get('rows',[]) if b['cosmic_proxy']<=.01 and b['LL_sigma_m_1000kms']<=.2 and en['max_fraction']<=.04]
        best=min(candidates,key=lambda b:b['max_loop_over_X_scales']) if candidates else None
        out.append(dict(params=p,scheme='vev-curvature',roll={k:v for k,v in rr.items() if k not in ['r','t','x','dx']},screen=q,energy=en,best_exact_energy=best))
        Path(output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
        print('BG',z,dn,xf,mu,en['max_fraction'],en['max_abs_comoving_continuity_residual'],
              None if best is None else (best['a'],best['max_loop_over_X_scales'],best['cosmic_proxy']),flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--tradeoff',required=True);ap.add_argument('--output',required=True);ap.add_argument('--group',type=int,required=True);a=ap.parse_args()
    run(json.loads(Path(a.input).read_text())[0],json.loads(Path(a.tradeoff).read_text()),a.group,a.output)
