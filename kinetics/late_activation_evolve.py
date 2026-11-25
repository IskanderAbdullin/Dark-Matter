"""Late contact switch-on, homogeneous NR kinetic pilot.
Comoving momentum fixed, cosmology enters through integrated rho*sigma*v.
Conservative graph; excludes condensate, gravitational clustering and 4->2.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.integrate import BDF,quad
from elastic_implicit_trial import ImplicitElastic,merge_grid
from late_activation_overlap import H

def run(row,z,sigma,n,qmin,wall=100.):
    h=row['final'];p0=h['peak_bin_q']
    e0=np.asarray(row['frozen_cell_spectrum']['edges']);f0=np.asarray(row['frozen_cell_spectrum']['f_phi'])
    e,f=merge_grid(e0,f0,n,qmin=qmin,qmax=e0[-1])
    e=e/p0
    # Exact NR angular geometry to relative O((p/m)^2). Dummy m only
    # defines a numerically stable limit of the existing relativistic graph.
    m=1e10;model=ImplicitElastic(e,m);model.rate*=(2*m)**3
    v=model.vol;N=float(v@f);E=float((v*model.K)@f)
    k2=model.p**2
    M0=4*np.pi**2*N
    norm=np.maximum(f,1.);ff=np.maximum(f,1e-38)
    scale=1/(32*np.pi**2*M0)/1e8
    v0=p0*2.7255*8.617333262e-5*(3.909/row['g_bath'])**(1/3)/h['theta']/86*2.99792458e10
    rho0=.12*1.8783416e-29
    tau=sigma*quad(lambda zz:rho0*v0*(1+zz)**3/H(zz),0,z,epsrel=1e-10)[0]
    tmax=tau*1e8
    def rhs(t,y):return model.rates_jac(y*norm,False)*scale/norm
    def jac(t,y):return model.rates_jac(y*norm,True)*scale*norm[None,:]/norm[:,None]
    weights=[v*k2*k2]
    for c in [2.,6.,20.]:
        l,u=e[:-1],e[1:]
        weights.append(v*np.maximum(u**3-np.maximum(l,c)**3,0.)/(u**3-l**3))
    W=np.array(weights);M=W@f
    initial=(W@(rhs(0.,ff/norm)*norm))/M*1e8
    start=time.monotonic();sol=BDF(rhs,0.,ff/norm,tmax,jac=jac,rtol=1e-5,atol=1e-12,max_step=tmax/20)
    hist=[];steps=0;msg='finished';last=start
    while sol.status=='running':
        tprev=sol.t;yold=sol.y.copy();note=sol.step();steps+=1
        g=sol.y*norm
        if v@np.maximum(-g,0.)/N>1e-9:
            sol.t=tprev;sol.y=yold;msg='negative accepted number';break
        if sol.status=='failed':msg=str(note);break
        if steps%10==0 or sol.status=='finished':
            hist.append(dict(t=float(sol.t),moment_changes=(W@g/M-1).tolist(),N_error=float(v@g/N-1),
              kinetic_E_error=float((v*model.K)@g/E-1),low_two_fraction=float(v[:2]@g[:2]/N),fmax=float(g.max())))
        if time.monotonic()-last>10:
            print('progress',sol.t,tmax,steps,flush=True);last=time.monotonic()
        if time.monotonic()-start>wall:msg='wall budget';break
    g=sol.y*norm
    # Dilute NR mean-field chemical potential / smallest resolved kinetic E.
    lam=float(np.sqrt(128*np.pi*86**3*sigma*1.7826619216278976e-33/(1.973269804e-5)**2))
    # g_NR = lambda/(8 m^2); n in eV^3 and p_peak at activation.
    pphys=p0*2.7255*8.617333262e-5*(3.909/row['g_bath'])**(1/3)/h['theta']*(1+z)
    mf_over_peak_kin=lam*N*pphys/(4*86.)
    return dict(success=sol.status=='finished' and msg=='finished',message=msg,z=z,sigma=sigma,n=len(f),qmin=qmin,
        requested_tau=tau,achieved_tau=float(sol.t/1e8),steps=steps,initial_rates_over_n_sigma_vpeak=initial.tolist(),
        final_moment_changes=(W@g/M-1).tolist(),N_error=float(v@g/N-1),kinetic_E_error=float((v*model.K)@g/E-1),
        mean_field_over_peak_kinetic=mf_over_peak_kin,mean_field_over_minimum_kinetic=float(mf_over_peak_kin/k2[0]),
        final_low_two_fraction=float(v[:2]@g[:2]/N),final_fmax=float(g.max()),wall_seconds=time.monotonic()-start,
        history=hist,edges_over_peak=e.tolist(),initial_f=f.tolist(),final_f=g.tolist(),scope='homogeneous NR elastic pilot with fixed cross section after activation, no zero mode, no full physical viability')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--output',required=True)
    p.add_argument('--z',type=float,default=3.);p.add_argument('--sigma',type=float,default=.1)
    p.add_argument('--n',type=int,default=64);p.add_argument('--qmin',type=float,default=1e-5)
    p.add_argument('--wall',type=float,default=90.);a=p.parse_args()
    r=run(json.loads(Path(a.input).read_text())[0],a.z,a.sigma,a.n,a.qmin,a.wall)
    Path(a.output).write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
    print('FINAL',json.dumps({k:v for k,v in r.items() if k not in ('history','edges_over_peak','initial_f','final_f')}),flush=True)
