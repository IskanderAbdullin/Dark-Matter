"""Independent BEST v1.3.5 check on smooth and burst PSDs.

External packages and checked-out solver live under ignored .runtime, no global
installation. Unit matrix element: multiply all rates by lambda**2 afterwards.
This is a fixed-geometry comparison, not an interacting relic calculation.
"""
import argparse,json,time,sys,os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
DEPS=ROOT/'.runtime'/'best-deps'
if DEPS.exists():
    sys.path.insert(0,str(DEPS))
    dll=DEPS/'Library'/'bin'
    if dll.exists():
        os.environ.setdefault('I_MPI_ROOT',str(DEPS/'Library'))
        _dll_handle=os.add_dll_directory(str(dll)) if os.name=='nt' else None
sys.path.insert(0,str(ROOT/'.runtime'/'best'))
sys.path.insert(0,str(ROOT/'kinetics'))
import numpy as np
from scipy.special import expit
from besthep import BEST,grid_weights
import gvar
from quartic_consistency import gap,lambda_from_sigma,shell_nodes

def load_candidate(path,sigma,x):
    row=json.loads(Path(path).read_text())[0]
    if 'scalar_history_psd' in row:
        i=int(np.argmin([abs(v['x']-x) for v in row['scalar_history_psd']]))
        h=row['history'][i];f=np.asarray(row['scalar_history_psd'][i]['f_phi'])
    else:
        h=row['final'];f=np.asarray(row['frozen_cell_spectrum']['f_phi'])
    e=np.asarray(row['frozen_cell_spectrum']['edges'])
    q,w=shell_nodes(e,f,8)
    lam=float(lambda_from_sigma(sigma));m=gap(q,w,row['scalar_ratio']*h['x'],lam)
    def f0(q):
        j=np.searchsorted(e,q,side='right')-1
        if j<0 or j>=len(f):return 0.
        return max(0.,float(f[j]))
    factor=lam**2*row['parent_GeV']*h['x']/(1.66*np.sqrt(row['g_bath'])*row['parent_GeV']**2/1.220890e19)*np.sqrt(np.pi**2/30*row['g_bath']/h['total_rho_Tref4'])
    return f0,m,factor,dict(x=h['x'],reference_number=float(w.sum()),reference_q2=float(w@(q*q)),reference_energy=float(w@np.hypot(q,m)),lambda_quartic=lam,sigma=sigma)

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--input');p.add_argument('--x',type=float,default=80.)
    p.add_argument('--sigma',type=float,default=1e-7)
    p.add_argument('--n',type=int,default=80);p.add_argument('--neval',type=int,default=20000)
    p.add_argument('--nitn',type=int,default=3);p.add_argument('--seed',type=int,default=192)
    p.add_argument('--qmin',type=float);p.add_argument('--qmax',type=float)
    p.add_argument('--equilibrium',action='store_true');p.add_argument('--step',type=float,default=0.)
    p.add_argument('--method',default='exprb_full');p.add_argument('--output',required=True)
    a=p.parse_args()
    if a.input:
        f,m,factor,meta=load_candidate(a.input,a.sigma,a.x)
        qmin=a.qmin or 1e-5;qmax=a.qmax or 50.
    else:
        m=1.;factor=1.;meta={}
        f=(lambda q:1/np.expm1(np.hypot(q,m)/2)) if a.equilibrium else (lambda q:expit((3-q)/2))
        qmin=a.qmin or .02;qmax=a.qmax or 35.
    s=BEST(q_min=qmin,q_max=qmax,n_grid=a.n)
    gvar.ranseed(a.seed+s.world_rank)
    s.initialize_species('phi',f,stat='boson',mass=m)
    # Unit coupling avoids underflow in rate-density variance. Restore units in moments.
    s.add_process('elastic',['phi','phi'],['phi','phi'],lambda p:np.ones(p.shape[2]),neval=a.neval,nitn=a.nitn)
    q=s.r_grids['phi'];f0=s.distributions_1d['phi'].copy();w=grid_weights(q)*q*q/(2*np.pi**2)
    E=np.hypot(q,m)
    t=time.monotonic();rates=s.compute_rates(['elastic'])['phi']
    if s.world_rank==0:
        N=float(w@f0);en=float(w@(E*f0));Q2=float(w@(q*q*f0))
        n0=float(s.zero_mode.get('phi',0.))
        out=dict(input=a.input,equilibrium=a.equilibrium,n=a.n,neval=a.neval,nitn=a.nitn,seed=a.seed,
            qmin=qmin,qmax=qmax,mass=m,world_size=s.world_size,unit_rate_multiplier=factor,meta=meta,
            initial_number=N+n0,initial_grid_number=N,initial_zero_mode=n0,initial_energy=en+m*n0,initial_q2=Q2,
            q2_derivative_per_unit=float(w@(q*q*rates)),relative_q2_derivative=float(factor*(w@(q*q*rates))/Q2),
            number_source=float(w@rates),energy_source=float(w@(E*rates)),
            number_net_over_absolute=float((w@rates)/max(w@abs(rates),1e-300)),
            energy_net_over_absolute=float((w@(E*rates))/max(w@(E*abs(rates)),1e-300)),
            rate_wall_seconds=time.monotonic()-t,quality={k:v for k,v in s.error_stats.items() if k!='rel_errs'},
            q=q.tolist(),f=f0.tolist(),rates=rates.tolist())
        if a.input:
            out['remesh_number_relative_error']=(N+n0)/meta['reference_number']-1
            out['remesh_energy_relative_error']=(en+m*n0)/meta['reference_energy']-1
            out['remesh_q2_relative_error']=Q2/meta['reference_q2']-1
        Path(a.output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
        print('BENCHMARK',json.dumps({k:v for k,v in out.items() if k not in ('q','f','rates')}),flush=True)
    if a.step>0:
        # a.step is in Hubble times for input, unit time otherwise.
        t=time.monotonic()
        s.evolve_step(dt=a.step*factor,method=a.method,adapt_dt=False)
        if s.world_rank==0:
            f1=s.distributions_1d['phi'];n01=float(s.zero_mode.get('phi',0.))
            out.update(step_requested=a.step,method=a.method,step_actual=s.current_time/factor,
                number_drift=float((w@f1+n01)/(N+n0)-1),
                energy_drift=float((w@(E*f1)+m*n01)/(en+m*n0)-1),
                q2_drift=float(w@(q*q*f1)/Q2-1),final_zero_mode=n01,
                final_min_f=float(f1.min()),final_f=f1.tolist(),step_wall_seconds=time.monotonic()-t)
            Path(a.output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
            print('STEP',json.dumps({k:out[k] for k in ('step_actual','number_drift','energy_drift','q2_drift','final_zero_mode','step_wall_seconds')}),flush=True)
if __name__=='__main__':main()
