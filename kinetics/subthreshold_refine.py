"""Refine small-mixing subthreshold candidates with actual time integration.

Conditional homogeneous roll, fixed near-critical thermal halo since z=2.
Transport normalized to one local Rayleigh e-fold, not complete halo
thermalization. Bose conversion is a frozen initial source, not evolved
abundance. No claim about primordial heavy/mediator population.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.integrate import trapezoid
from scipy.interpolate import PchipInterpolator
from subthreshold_roll_scan import params,roll
from subthreshold_transport_tradeoff import bose_source
from subthreshold_recombination import cosmic_zero,halo
from soft_mixing_activation import frac_B
from late_activation_overlap import H,critical_density
from resonant_mediator_inverse import M,C,SEC_GYR
from quartic_consistency import lambda_from_sigma
from velocity_only_exposure import RHO0

def get_kernels(d,a,n):
    return bose_source(d,.5,a,n,n,return_kernels=True)

def refine(row,p0,b,scheme,table,formation_z=2.,n=14000,rtol=1e-8):
    p=params(p0['z_on'],p0['dneff'],p0['Xfrac'],p0['loop_margin'],
             p0['mu_over_H_on'],b['u'])
    rr=roll(row,p,scheme,rtol=rtol,n=n)
    r,t,x=rr['r'],rr['t'],rr['x']
    vT=np.sqrt(2)*44.1e5;rho=critical_density(vT);scale=M*M*(vT/C)**2
    D,delta,u,a=p['Delta'],p['delta'],p['u'],b['a']
    shift=lambda x:-2*delta**2*x*x/(np.hypot(D,2*delta*x)+D)
    aval=a+4*(shift(x)-shift(1))/scale
    assert min(aval)>0
    mix=frac_B(x,delta,D)
    gap=np.hypot(D,2*delta*x)/scale
    shape=(mix/u)**4*(a/aval)**2
    qfun=PchipInterpolator(np.log(table['aa']),np.log(table['Q96']))
    Q=np.exp(qfun(np.log(np.clip(aval,min(table['aa']),max(table['aa'])))))
    h=r<=formation_z+1;late=r<=101
    clock=1/H(r-1)
    integ=lambda y,mask:float(trapezoid(y[mask]*clock[mask],t[mask]))
    norm=rho*vT*integ(shape*Q,h)
    sig=1/norm;g2=float(lambda_from_sigma(sig))*a*scale/u**2
    S2=4*(M*M+shift(1))-a*scale
    lam=3*g2/S2
    tau=sig*RHO0*rr['v0']*integ(r**4*shape,late)
    w1=(mix/u)**3*(1-mix)/u*(a/aval)**2
    w2=(mix/u)**2*((1-mix)/u)**2*(a/aval)**2
    # d,a vary only at order u times the field excursion. Compute full
    # kernels at three actual x^2 nodes; compare n=96,144,192 at full x.
    z2=np.unique(np.r_[min(x[h]**2),1.,max(x[h]**2)])
    kvals=[]
    for zz in z2:
        dp=np.hypot(D,2*delta*np.sqrt(zz))/scale
        ap=a+4*(shift(np.sqrt(zz))-shift(1))/scale
        kvals.append(get_kernels(dp,ap,144))
    kvals=np.asarray(kvals)
    k1=np.interp(x*x,z2,kvals[:,0]);k2=np.interp(x*x,z2,kvals[:,1])
    source1=sig*rho*vT*integ(w1*k1,h)
    source2=sig*rho*vT*integ(w2*k2,h)
    conv=[]
    for order in (96,144,192):
        k=get_kernels(np.hypot(D,2*delta)/scale,a,order)
        conv.append(dict(order=order,kernels=k,constant_kernel_source= sig*rho*vT*(k[0]*integ(w1,h)+k[1]*integ(w2,h))))
    # Wilsonian estimate, not regulator-independent bound: require the
    # hard stabilizer tadpole and cubic threshold <= physical splitting.
    cubic=g2/(16*np.pi**2)
    cutoff=np.sqrt(max(0.,(D-cubic)*32*np.pi**2/lam))
    # S propagator one-loop slope below heavy-heavy threshold.
    from scipy.integrate import quad
    B2=M*M+D+shift(1)*(-1)
    slope=g2/(32*np.pi**2)*quad(lambda y:y*(1-y)/(B2-y*(1-y)*S2),0,1,epsabs=1e-11)[0]
    slopeL=g2*u*u/(32*np.pi**2)*quad(lambda y:y*(1-y)/(M*M+shift(1)-y*(1-y)*S2),0,1,epsabs=1e-11)[0]
    ret=dict(params=p,scheme=scheme,candidate=b,formation_z=formation_z,rtol=rtol,n=n,
       sigma0_for_one_efold=sig,g_eV=float(np.sqrt(g2)),lambda0=lam,D=D,delta=delta,u=u,a=a,
       S_mass_eV=float(np.sqrt(S2)),pole_detuning_eV2=float(a*scale),
       gap_over_T=float(np.hypot(D,2*delta)/scale),min_a=float(min(aval)),max_a=float(max(aval)),
       halo_active_Gyr=integ(shape,h)/SEC_GYR,halo_Q_time_Gyr=integ(shape*Q,h)/SEC_GYR,
       LL_tau=tau,cosmic_proxy=5.756e7*tau,
       cold_v2_over_detuning=float(max((12*rr['v0']*r[late]/vT)**2/aval[late])),
       LH_initial_source=source1,HH_initial_source=source2,initial_H_source=source1+source2,
       conversion_convergence=conv,kernel_x2=z2.tolist(),kernel_values=kvals.tolist(),
       split_loop_ratio=float((cubic+lam*M*M/(32*np.pi**2))/D),
       pole_loop_ratio=float(g2/(32*np.pi**2*a*scale)),
       quartic_cutoff_eV=float(cutoff),S_bubble_derivative_HH=float(slope),S_bubble_derivative_LL=float(slopeL),
       roll={k:v for k,v in rr.items() if k not in ['t','r','x','dx']})
    return ret,rr

def select(scans):
    criteria=[(8,.2,.7,1.1,20,'vev-curvature','factor2_margin'),
      (10,.1,.7,1.1,20,'vev-curvature','factor2_margin'),
      (12,.1,.7,1.1,20,'vev-curvature','factor10_margin'),
      (10,.2,.5,1.1,20,'vev-curvature','factor2_margin'),
      (8,.3,.7,1.1,20,'origin-vev','factor2_margin')]
    for z,dn,xf,margin,muH,scheme,lab in criteria:
        v=[s for s in scans if all(s['params'][k]==val for k,val in zip(
            ['z_on','dneff','Xfrac','loop_margin','mu_over_H_on'],[z,dn,xf,margin,muH])) and s['scheme']==scheme]
        if not v or lab not in v[0]['screen']['best']:continue
        yield v[0]['params'],v[0]['screen']['best'][lab],scheme

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--input',required=True);ap.add_argument('--input2',required=True)
    ap.add_argument('--tradeoff',required=True);ap.add_argument('--scan',nargs='+',required=True)
    ap.add_argument('--output',required=True);ap.add_argument('--index',type=int,default=-1)
    a=ap.parse_args();rows=[json.loads(Path(p).read_text())[0] for p in [a.input,a.input2]]
    table=json.loads(Path(a.tradeoff).read_text())
    scans=sum([json.loads(Path(p).read_text()) for p in a.scan],[])
    out=[]
    for i,(p,b,scheme) in enumerate(select(scans)):
        if a.index>=0 and i!=a.index:continue
        case=[]
        for j,row in enumerate(rows):
            ret,rr=refine(row,p,b,scheme,table,n=14000 if j==0 else 22000,rtol=1e-8 if j==0 else 3e-9)
            if j==0:
                ret['cosmic_recombination']=cosmic_zero(row,ret,rr)
                ret['halo_recombination_at_full']=halo({k:ret[k] for k in ['u','D','delta','a','g_eV']},powers=(16,18),replicates=2)
            case.append(ret)
            Path(a.output).write_text(json.dumps(out+[dict(index=i,runs=case)],indent=2,allow_nan=False)+'\n')
            print('REFINE',i,j,{k:ret[k] for k in ['sigma0_for_one_efold','g_eV','initial_H_source','LL_tau','cosmic_proxy','split_loop_ratio','pole_loop_ratio','quartic_cutoff_eV']},flush=True)
        out.append(dict(index=i,runs=case))
        Path(a.output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
