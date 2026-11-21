"""High-occupation *six-dimensional* phase volume of an exactly critical Bose halo.
Mathematical local Poisson/Vlasov model, not a physical collisionless halo-formation demonstration.
"""
import json
import numpy as np
from scipy.special import zeta
from scipy.integrate import solve_ivp,quad
nuc=float(zeta(1.5))
z_half=zeta(.5);z_mhalf=zeta(-.5);z_m3half=zeta(-1.5)

def density_of_y(y):
    y=max(0,y)
    # Li3/2(exp(-y)) series valid and accurate at the small y used by F>=10
    return (nuc-2*np.sqrt(np.pi*y)-z_half*y+
             z_mhalf*y*y/2-z_m3half*y**3/6)/nuc

def phase_volume(F):
    L=np.log1p(1/F)
    eps=1e-7
    def rhs(r,w):return[w[1]/r**2,r*r*density_of_y(w[0])]
    def stop(r,w):return w[0]-L
    stop.terminal=True
    stop.direction=1
    sol=solve_ivp(rhs,[eps,4],[eps**2/6,eps**3/3],
       rtol=1e-10,atol=1e-13,max_step=.006,dense_output=True,events=stop)
    assert sol.success and len(sol.t_events[0])==1
    rt=float(sol.t_events[0][0])
    fn=lambda xi:xi*xi*max(0,L-sol.sol(xi)[0])**1.5
    vol=(16*np.pi**2/3)*np.sqrt(2)**3*quad(fn,0,rt,epsabs=1e-13,epsrel=1e-11)[0]
    leading=4*np.sqrt(3)*np.pi**3*L**3
    return dict(F=F,L=L,xi_max=rt,
        volume_over_R3_mV03=vol,leading_6D_ball=leading,ratio=vol/leading)

def run():
    data=[phase_volume(F) for F in (10,30,100,300,1000,3000,10000,100000)]
    vals=data[-4:]
    slope=float(np.polyfit(np.log([x["F"] for x in vals]),
              np.log([x["volume_over_R3_mV03"] for x in vals]),1)[0])
    B=128/(35*np.sqrt(np.pi)*nuc)
    assert abs(data[-1]["ratio"]-1)<.005
    assert abs(slope+3)<.03
    assert abs(data[-1]["ratio"]-(1+B*np.sqrt(data[-1]["L"])))<1e-4
    return dict(nu_critical=nuc,leading_coefficient_4root3pi3=4*np.sqrt(3)*np.pi**3,
       analytic_first_correction_B=B,loglog_slope=slope,
       rows=data,checks_passed=3,
       caution="Fine-grained closed-patch Casimir analysis, not empirical halo phase-space data or proof of condensate production.")
if __name__=="__main__":print(json.dumps(run(),indent=2))
