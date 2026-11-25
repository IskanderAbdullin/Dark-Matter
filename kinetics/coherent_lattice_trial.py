"""Small independent classical real-field lattice diagnostic.

Dimensionless equation y''=Laplacian(y)-alpha*y-beta*y^3.
Periodic cubic nearest-neighbour spatial operator and kick-drift leapfrog.
Tests the conditional coherent-field Floquet prediction and nonlinear
backreaction in flat spacetime. No quantum vacuum UV noise, expansion,
decays or identification of the full dark-sector PSD with this field.
"""
import argparse,json,time,sys
from pathlib import Path
import numpy as np
DEPS=Path(__file__).resolve().parents[1]/'.runtime'/'elastic-deps'
if DEPS.exists():sys.path.insert(0,str(DEPS))
from numba import njit

@njit(cache=True)
def accel(y,out,idx2,alpha,beta):
    n=y.shape[0]
    for i in range(n):
        ip=(i+1)%n;im=(i-1)%n
        for j in range(n):
            jp=(j+1)%n;jm=(j-1)%n
            for k in range(n):
                kp=(k+1)%n;km=(k-1)%n
                z=y[i,j,k]
                out[i,j,k]=(y[ip,j,k]+y[im,j,k]+y[i,jp,k]+y[i,jm,k]+y[i,j,kp]+y[i,j,km]-6*z)*idx2-alpha*z-beta*z*z*z

@njit(cache=True)
def advance(y,v,a,idx2,alpha,beta,dt,nsteps):
    for _ in range(nsteps):
        accel(y,a,idx2,alpha,beta)
        v+=dt*a
        y+=dt*v

def lattice_k2(n,L):
    k=2*np.pi*np.fft.fftfreq(n,L/n)
    v=4*np.sin(k*L/(2*n))**2/(L/n)**2
    return v[:,None,None]+v[None,:,None]+v[None,None,:]

def stats(y,v,a,idx2,alpha,beta,dt,k2,mask):
    accel(y,a,idx2,alpha,beta)
    vel=v+dt*a/2
    kin=.5*np.mean(vel*vel)
    grad=.5*idx2*sum(np.mean((np.roll(y,1,axis=ax)-y)**2) for ax in range(3))
    pot=np.mean(alpha*y*y/2+beta*y**4/4)
    mean=float(y.mean());vm=float(vel.mean())
    zero=.5*vm*vm+alpha*mean**2/2+beta*mean**4/4
    # Discrete energy partition; potential cross terms assigned to nonzero component.
    var=float(np.mean((y-mean)**2));vvar=float(np.mean((vel-vm)**2))
    fk=np.fft.fftn(y-mean,norm='ortho')
    vk=np.fft.fftn(vel-vm,norm='ortho')
    # Use positive fixed frequency to avoid a zero crossing in the mode diagnostic.
    bandpower=float(np.sum((abs(vk[mask])**2+(k2[mask]+1)*abs(fk[mask])**2)/(k2[mask]+1))/y.size)
    # UV-shell energy proxy for detecting lattice-end contamination.
    uv=k2>(.75*np.sqrt(k2.max()))**2
    fenergy=abs(vk)**2+(k2+1)*abs(fk)**2
    return dict(U=float(kin+grad+pot),zero_energy=float(zero),zero_fraction=float(zero/(kin+grad+pot)),
                field_variance=var,velocity_variance=vvar,bandpower=bandpower,
                uv_quadratic_energy_fraction=float(fenergy[uv].sum()/max(fenergy.sum(),1e-300)),
                gradient=float(grad),kinetic=float(kin),potential=float(pot))

def run(n,L,dt,end,seed,noise,alpha=0.,sample=2.,mode='coherent',init=None):
    beta=1-alpha;dx=L/n;k2=lattice_k2(n,L)
    assert dt*np.sqrt(k2.max()+alpha+3*beta)<1.5
    rng=np.random.default_rng(seed);s=beta/2
    lo=3-3*s;hi=lo+3*s*s/(2*np.sqrt(1-s+s*s)+2-s)
    mask=(k2>=lo)&(k2<=hi)
    if mode=='coherent':
        assert mask.sum()>0
        f=np.fft.fftn(rng.standard_normal((n,n,n)),norm='ortho')
        f*= (k2>.7)&(k2<3.)
        perturb=np.fft.ifftn(f,norm='ortho').real
        perturb*=noise/np.std(perturb)
        y=1+perturb;vel=np.zeros_like(y)
    elif mode=='gas':
        y,vel=init(k2,rng)
    else:raise ValueError(mode)
    a=np.empty_like(y);accel(y,a,1/dx**2,alpha,beta)
    v=vel-dt*a/2
    hist=[];Nstep=int(round(end/dt));every=max(1,int(round(sample/dt)))
    t0=time.monotonic()
    for i in range(0,Nstep+1,every):
        h=stats(y,v,a,1/dx**2,alpha,beta,dt,k2,mask);h['t']=i*dt
        hist.append(h)
        if i>=Nstep:break
        advance(y,v,a,1/dx**2,alpha,beta,dt,min(every,Nstep-i))
    if hist[-1]['t']<Nstep*dt:
        h=stats(y,v,a,1/dx**2,alpha,beta,dt,k2,mask);h['t']=Nstep*dt;hist.append(h)
    U=hist[0]['U']
    ans=dict(n=n,L=L,dt=dt,seed=seed,noise=noise,alpha=alpha,mode=mode,band_modes=int(mask.sum()),
             maximum_energy_drift=max(abs(h['U']/U-1) for h in hist),
             final=hist[-1],min_zero_fraction=min(h['zero_fraction'] for h in hist),
             wall_seconds=time.monotonic()-t0,history=hist)
    if mode=='coherent':
        # Fit amplitude envelope in an intermediate linear interval.
        t=np.array([h['t'] for h in hist]);p=np.array([h['bandpower'] for h in hist])
        good=(t>=60)&(t<=180)&(p<1e-4)
        if good.sum()>3:ans['fitted_mu']=float(np.polyfit(t[good],.5*np.log(p[good]),1)[0])
    return ans

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=40);p.add_argument('--L',type=float,default=40.)
    p.add_argument('--dt',type=float,default=.04);p.add_argument('--end',type=float,default=450.);p.add_argument('--seed',type=int,default=42)
    p.add_argument('--noise',type=float,default=1e-6);p.add_argument('--output',required=True);a=p.parse_args()
    r=run(a.n,a.L,a.dt,a.end,a.seed,a.noise)
    Path(a.output).write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
    print('FINAL',json.dumps({k:v for k,v in r.items() if k!='history'}),flush=True)
