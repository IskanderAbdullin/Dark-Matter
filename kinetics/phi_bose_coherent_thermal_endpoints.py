"""Conditional scalar ideal-gas endpoints at fixed energy and (optionally) action count.
eV units. Does not calculate equilibration rates or a physical BEC.
"""
import math,json
from phi_bose_narrow_sweep import thermal
from scipy.integrate import quad
from scipy.optimize import brentq
m=86.; n0=.12*1.05375e4/m*(1.973269804e-5)**3

def action_correction(A):
    return 4/math.pi*quad(lambda t:math.cos(t)**2*(A*A*(1+math.sin(t)**2))/(math.sqrt(1+A*A*(1+math.sin(t)**2))+1),0,math.pi/2,epsabs=1e-15,epsrel=1e-11)[0]

def action(A,l):
    return 6*m**3/l*A*A*(1+action_correction(A))

def getA(a,l):
    return math.exp(brentq(lambda t:math.log(action(math.exp(t),l)*a**3/n0),-30,25))

def moments(T,critical=False):
    x=m/T; scale=math.sqrt(1+x)
    def fun(z,i):
        y=z*scale; e=math.hypot(y,x); ke=y*y/(e+x); d=ke if critical else e
        f=math.exp(-d)/(-math.expm1(-d))
        return scale*[y*y*f,y*y*e*f,y**4*f,y*y*ke*f,y**4/(3*e)*f][i]
    ns,r,ps,u,pr=[quad(lambda z:fun(z,i),0,80,epsabs=1e-10,epsrel=1e-9)[0]/(2*math.pi**2) for i in range(5)]
    return dict(n=ns*T**3,rho=r*T**4,pr=pr*T**4,
                prms=T*math.sqrt(ps/ns) if ns>0 else 0,u=u*T**4)

def solve(T,l):
    a,H=thermal(T); A=getA(a,l); J=action(A,l); rho=6*m**4/l*(A*A+A**4)
    tc=math.exp(brentq(lambda x:math.log(max(moments(math.exp(x))['rho'],1e-300)/rho),math.log(m/100),math.log(1e10)))
    chem=moments(tc); n=chem['n']; u=6*m**4/l*A*A*(A*A-action_correction(A))
    tk=math.exp(brentq(lambda x:math.log(moments(math.exp(x),True)['u']/u),math.log(m/1e9),math.log(1e10)))
    kin=moments(tk,True); qcut=1.14255e-6; p=qcut/a; Ecut=math.hypot(m,p)
    assert abs(chem['rho']/rho-1)<1e-9 and abs(kin['u']/u-1)<1e-9
    return dict(T=T,A=A,delta=u/(m*J),chem_T_over_m=tc/m,chem_Omega=.12*n/J,
                chem_qrms=a*chem['prms'],chem_gas_n_over_J=n/J,
                chem_s_over_J=(rho+chem['pr'])/(tc*J),kin_T_over_m=tk/m,
                kin_nexc_over_J=kin['n']/J,kin_rms_total=a*kin['prms']*math.sqrt(kin['n']/J),
                hot_number_fraction_bound_if_J_conserved=min(1,u/(J*p*p/(Ecut+m))),
                qcut_phys_over_m=p/m,thermal_mass_over_m_sq=l*(tc/m)**2/24)

if __name__=='__main__':
    for l,Ti in [(4.698884e-13,413349.430520685),(3.82587e-8,5692.465733402297)]:
        for T in [Ti,2e7,1e5,1e3]: print(l,json.dumps(solve(T,l)),flush=True)
    r=moments(m*1e4)
    assert abs(r['rho']/(m*1e4)**4/(math.pi**2/30)-1)<1e-7
    from scipy.special import zeta
    T=m*1e-6; n=moments(T,True)['n']
    assert abs(n/((m*T/(2*math.pi))**1.5*zeta(1.5))-1)<1e-5
    assert abs(action_correction(1e-6)/1e-12-5/8)<1e-10
