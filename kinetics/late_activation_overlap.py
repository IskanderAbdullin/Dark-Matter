"""Conditional late activation overlap: frozen homogeneous PSD vs local halo.
No mediator model, backreaction, structure formation or observational bound.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import zeta

H0=.674*100*1e5/3.085677581e24
Om=.315;Or=9.2e-5
def H(z):return H0*np.sqrt(Om*(1+z)**3+Or*(1+z)**4+1-Om-Or)
def elapsed(z):return quad(lambda x:1/((1+x)*H(x)),0,z,epsrel=1e-10)[0]
def critical_density(v):
    m=86*1.7826619216278976e-33
    return m*zeta(1.5,1)*(m*v/(2*np.sqrt(np.pi)*1.054571817e-27))**3

def screen(data,halo_rate):
    rates=np.array(data['orders'][-1]['relative_rate_over_n_sigma_vpeak'])
    rho=data['rho0_g_cm3'];v=data['peak_velocity_today_cm_s']
    cases=[('illustrative dilute',1e-24,1e7,.851076865),('critical original F583-4 temperature',critical_density(np.sqrt(2)*44.1e5),np.sqrt(2)*44.1e5,halo_rate)]
    def exposure(z):return quad(lambda x:rho*v*(1+x)**3/H(x),0,z,epsrel=1e-10)[0]*rates
    def upper(z,cap):return cap/np.max(abs(exposure(z)))
    out=[]
    for name,d,vh,c in cases:
        def lower(z):return 1/(c*d*vh*elapsed(z))
        rows=[]
        for z in (.1,.5,1.,2.,3.,5.,10.,30.,100.):
            rows.append(dict(z=z,elapsed_Gyr=elapsed(z)/(365.25*86400*1e9),sigma_min_one_local_efold=lower(z),
                        sigma_max_frozen_log_change_001=upper(z,.01),sigma_max_frozen_log_change_01=upper(z,.1),
                        moment_exposures_per_cm2_g=exposure(z).tolist()))
        cross={}
        for cap in (.01,.1):
            a,b=1e-4,100.
            vroot=brentq(lambda z:np.log(upper(z,cap)/lower(z)),a,b) if upper(a,cap)>lower(a) and upper(b,cap)<lower(b) else None
            cross[str(cap)]=dict(z_cross=vroot,sigma_at_cross=lower(vroot) if vroot else None)
        out.append(dict(case=name,rho_g_cm3=d,v_T_cm_s=vh,rayleigh_over_n_sigma_vT=c,rows=rows,crossings=cross))
    return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--nr',required=True);p.add_argument('--rayleigh',required=True);p.add_argument('--output',required=True);a=p.parse_args()
    r=screen(json.loads(Path(a.nr).read_text()),json.loads(Path(a.rayleigh).read_text())[-1]['rayleigh_over_n_sigma_vT'])
    Path(a.output).write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
    for d in r:
        print(d['case'],d['rho_g_cm3'],d['crossings'])
        for row in d['rows']:print(row['z'],row['sigma_min_one_local_efold'],row['sigma_max_frozen_log_change_001'])
