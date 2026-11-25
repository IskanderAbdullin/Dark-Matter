"""Late conversion by nearly degenerate bosonic parent decay.

A -> B + massless eta, same inertial mass in NR kinematics, isotropic
fixed recoil. B has large constant self interaction, not included here.
First solve the restricted spontaneous + B-stimulated decay system
(f_eta=0, no inverse reaction). This is NOT a physical lower bound when
finite recoil, eta trapping, or B redistribution are important.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.integrate import BDF,quad,trapezoid
from late_activation_overlap import H,H0,Or
from resonant_mediator_inverse import M,C,HBAR,SEC_GYR
from direct_activation_consistency import T0
from quartic_consistency import shell_nodes

def mesh(row,n=112,cut=20.,max_b=0.):
    p0=row['final']['peak_bin_q'];e=np.array(row['frozen_cell_spectrum']['edges'])/p0
    f=np.maximum(row['frozen_cell_spectrum']['f_phi'],0.)
    ee=np.unique(np.r_[0.,e[(e>0)&(e<cut)],cut])
    j=np.searchsorted(e,(ee[:-1]+ee[1:])/2)-1
    ff=np.where((j>=0)&(j<len(f)),f[np.clip(j,0,len(f)-1)],0.)
    N=ff*np.diff(ee**3)/(6*np.pi**2)
    ix=np.unique(np.rint(np.linspace(0,len(ff),min(n,len(ff))+1)).astype(int))
    a=ee[ix];vol=np.diff(a**3)/(6*np.pi**2)
    g=np.array([N[l:h].sum() for l,h in zip(ix[:-1],ix[1:])])/vol
    fullN=float(shell_nodes(e,f,12)[1].sum())
    # Extend daughter support beyond recoil and parent cut. No parent in tail.
    tail=np.geomspace(cut,cut+max(1.,1.01*max_b),max(2,int(n/8))+1)[1:]
    b=np.r_[a,tail]
    return a,b,g,dict(retained_number_fraction=float(vol@g/fullN),fullN=fullN)

def vgrid(e):
    l,h=e[:-1],e[1:]
    return np.diff(e**3)/(6*np.pi**2),np.sqrt(.6*(h**5-l**5)/(h**3-l**3))

def transition(q,b,edges):
    # P(q_B|q_A)=q_B/(2 q_A b) within [|q_A-b|,q_A+b].
    big=np.maximum(q,b)[:,None];small=np.minimum(q,b)[:,None]
    cs=((edges[None,:]-big)*(edges[None,:]+big)-small*small)/(2*q[:,None]*b)
    return np.diff(np.clip(cs,-1.,1.),axis=1)/2

def decay(row,recoil,lifetime=10.,n=112,cut=20.,rend=10001.,rtol=2e-5,wall=130.):
    start=time.monotonic()
    Tref=T0*(3.909/row['g_bath'])**(1/3)/row['final']['theta'];p0=Tref*row['final']['peak_bin_q']
    rmax=min(1e8,.03*M/(cut*p0));rmin=rend
    assert rmax>rmin
    eA,eB,f0,meta=mesh(row,n,cut,recoil/(p0*rmin))
    vA,qA=vgrid(eA);vB,qB=vgrid(eB);N0=float(vA@f0)
    na=len(qA);nb=len(qB)
    # Both densities scaled in occupancy; B initially zero.
    normA=np.maximum(f0,1.)
    normB=np.maximum(np.r_[f0,np.zeros(nb-na)],1.)
    scale=np.r_[normA,normB]
    data=np.r_[f0,np.zeros(nb)]
    rate=1/(lifetime*SEC_GYR)
    last_t=None;W=None;g=None;maxerr=0.;bad_geometry=[];end=np.log(rmax/rmin)
    def geo(t):
        nonlocal last_t,W,g,maxerr
        if t!=last_t:
            r=rmax*np.exp(-np.clip(t,0.,end))
            W=transition(qA,recoil/(p0*r),eB)
            err=float(max(abs(W.sum(axis=1)-1)))
            maxerr=max(maxerr,err)
            if err>1e-8 and len(bad_geometry)<10:bad_geometry.append([float(t),float(r),err])
            g=rate/H(r-1);last_t=t
    def rhs(t,y):
        geo(t);f=y*scale;A=f[:na];B=f[na:]
        return g*np.r_[-A*(W@(1+B)),(W.T@(vA*A))/vB*(1+B)]/scale
    def jac(t,y):
        geo(t);f=y*scale;A=f[:na];B=f[na:]
        J=np.zeros((na+nb,na+nb))
        J[:na,:na]=np.diag(-(W@(1+B)))
        J[:na,na:]=-A[:,None]*W
        J[na:,:na]=(1+B)[:,None]*W.T*vA[None,:]/vB[:,None]
        J[na:,na:]=np.diag((W.T@(vA*A))/vB)
        return g*J*scale[None,:]/scale[:,None]
    y0=data/scale
    sol=BDF(rhs,0.,y0,end,jac=jac,rtol=rtol,atol=1e-11,max_step=.03,first_step=1e-5)
    hist=[];first={};steps=0;msg='finished';last=start
    def stat(t,y):
        f=y*scale;A=f[:na];B=f[na:]
        return dict(r=float(rmax*np.exp(-t)),B_fraction=float(vB@B/N0),
           number_error=float((vA@A+vB@B)/N0-1),B_max_f=float(B.max()),
           negative_fraction=float((vA@np.maximum(-A,0)+vB@np.maximum(-B,0))/N0))
    hist.append(stat(0.,sol.y))
    while sol.status=='running':
        tp=sol.t;yp=sol.y.copy();note=sol.step();steps+=1
        s=stat(sol.t,sol.y)
        if s['negative_fraction']>1e-8:
            sol.t=tp;sol.y=yp;msg='negative accepted number';break
        if sol.status=='failed':msg=str(note);break
        for lim in [.001,.01,.1,.5,.9]:
            if s['B_fraction']>=lim and str(lim) not in first:first[str(lim)]=s
        if steps%5==0 or sol.status=='finished':hist.append(s)
        if time.monotonic()-last>12:
            print('progress',json.dumps(s),flush=True);last=time.monotonic()
        if time.monotonic()-start>wall:msg='wall budget';break
    if hist[-1]['r']!=float(rmax*np.exp(-sol.t)):hist.append(stat(sol.t,sol.y))
    return dict(success=sol.status=='finished' and msg=='finished',message=msg,nA=na,nB=nb,cut=cut,rmax=rmax,rend=rend,
        recoil_eV=recoil,kick_kms=recoil/M*C/1e5,lifetime_Gyr=lifetime,peak_momentum_today_eV=p0,
        meta=meta,max_transition_row_error=maxerr,bad_geometry=bad_geometry,first_crossings=first,final=hist[-1],history=hist,
        wall_seconds=time.monotonic()-start,steps=steps,edges_A=eA.tolist(),edges_B=eB.tolist(),
        initial_A=f0.tolist(),final_A=(sol.y[:na]*scale[:na]).tolist(),final_B=(sol.y[na:]*scale[na:]).tolist(),
        scope='Frozen recoil NR decay network; f_eta set zero, no inverse, no B elastic or mean field, no hot parent modes. Not a self-consistent exclusion.')

def collinear(row,lifetime=10.):
    p0=row['final']['peak_bin_q'];e=np.array(row['frozen_cell_spectrum']['edges']);f=np.maximum(row['frozen_cell_spectrum']['f_phi'],0.)
    w=f*np.diff(e**3)/(6*np.pi**2);N=w.sum()
    out=[]
    for r in [1e6,3e5,1e5,3e4,1e4,1001.]:
        t=quad(lambda l:1/H(np.exp(l)-1),np.log(r),np.log(1e11),epsrel=1e-9)[0]
        g=t/(lifetime*SEC_GYR)
        # fB=(f0+1)/(1+f0 exp(-(1+f0)g))-1, stable even for large f0.
        z=-np.minimum((f+1)*g,745.)
        b=f*(-np.expm1(z))/(1+f*np.exp(z))
        frac=float((b*np.diff(e**3)/(6*np.pi**2)).sum()/N)
        out.append(dict(r=r,t_seconds=t,converted_fraction=frac,spontaneous_fraction=-np.expm1(-g)))
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--output',required=True)
    p.add_argument('--recoil',nargs='+',type=float,default=[1e-5,1e-4,1e-3,.01]);p.add_argument('--lifetime',type=float,default=10.)
    p.add_argument('--n',type=int,default=112);p.add_argument('--rend',type=float,default=10001.);p.add_argument('--wall',type=float,default=100.);a=p.parse_args()
    row=json.loads(Path(a.input).read_text())[0];out=dict(collinear_limit=collinear(row,a.lifetime),cases=[])
    for k in a.recoil:
        r=decay(row,k,a.lifetime,a.n,rend=a.rend,wall=a.wall);out['cases'].append(r)
        Path(a.output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
        print('FINAL',json.dumps({k:v for k,v in r.items() if k not in ['history','edges_A','edges_B','initial_A','final_A','final_B']}),flush=True)
