"""Adiabatic real scalar background and linear self-resonance, not a relic solver."""
import json
import numpy as np
from scipy.special import ellipk, ellipj
from scipy.integrate import quad, solve_ivp, simpson
from scipy.optimize import brentq

M=86e-9  # GeV
RHO0=.12*1.05375e-5*(1.973269804e-14)**3
TG0=2.34865e-13
GS0=3.91
MPL=1.22089e19
LAMBDAS=[4.698883527e-13,3.825874e-8]

def integrals(A):
    def d(t): return 1+A*A*(1+np.sin(t)**2)
    def integ(f): return quad(f,0,np.pi/2,epsabs=1e-12,epsrel=1e-12)[0]
    I0=integ(lambda t:1/np.sqrt(d(t)))
    Ic=integ(lambda t:np.cos(t)**2*np.sqrt(d(t)))
    D0=integ(lambda t:-A*A*(1+np.sin(t)**2)*d(t)**(-1.5))/I0
    Dc=integ(lambda t:np.cos(t)**2*A*A*(1+np.sin(t)**2)/np.sqrt(d(t)))/Ic
    return I0,Ic,D0,Dc

def action(A):
    return (4/np.pi)*A*A*integrals(A)[1]

def amp(j):
    return np.exp(brentq(lambda l:np.log(action(np.exp(l)))-np.log(j),-50,50,xtol=1e-12))

def background(a,lam):
    escale=6*M**4/lam
    enow=RHO0/escale
    A0=np.sqrt(2*enow/(1+np.sqrt(1+4*enow)))
    j=action(A0)/a**3
    A=amp(j)
    I0,Ic,D0,Dc=integrals(A)
    e=A*A*(1+A*A)
    I4=quad(lambda t:np.sin(t)**4/np.sqrt(1+A*A*(1+np.sin(t)**2)),0,np.pi/2,epsabs=1e-12,epsrel=1e-12)[0]
    return dict(A=A,j=j,rho=escale*e,w=A**4*I4/(I0*e),cs2=-D0/(2+Dc),omega=M*np.pi/(2*I0),phi=A*np.sqrt(12/lam)*M)

def band(A):
    o=np.sqrt(1+2*A*A)
    p=(A/o)**2
    lo=3*(1+A*A)
    width=o**2*3*p*p/(2*np.sqrt(1-p+p*p)+2-p)
    return p,o,lo,width

def floquet(A,k2,rtol=2e-12):
    p,o,lo,width=band(A)
    L=2*ellipk(p)
    h=(k2+1+6*A*A)/o**2
    def f(z,v):
        sn=ellipj(z,p)[0]
        b=h-6*p*sn*sn
        return [v[2],v[3],-b*v[0],-b*v[1]]
    s=solve_ivp(f,[0,L],[1,0,0,1],method='DOP853',rtol=rtol,atol=rtol*.01)
    assert s.success
    mat=s.y[:,-1].reshape(2,2)
    tr=np.trace(mat)
    det=np.linalg.det(mat)
    d=abs(tr)/(2*np.sqrt(det))-1
    mu=np.arccosh(1+max(0,d))/L*o
    return dict(mu_over_m=mu,trace=float(tr),det=float(det),disc_excess=float(d))

def independent(A):
    I0,Ic,D0,Dc=integrals(A)
    p,o,l,w=band(A)
    k2=l+w/2
    def rhs(t,v):
        x,dx,a,b,da,db,ip=v
        B=1+k2+6*x*x
        return dx,-x-2*x**3,da,db,-B*a,-B*b,dx*dx-x*x-x**4
    s=solve_ivp(rhs,[0,2*I0],[A,0,1,0,0,1,0],rtol=2e-13,atol=1e-14,method='DOP853')
    assert s.success
    x,dx,a,b,da,db,ip=s.y[:,-1]
    det=a*db-b*da
    d=abs(a+db)/(2*np.sqrt(det))-1
    mu=np.arccosh(1+max(0,d))/(2*I0)
    return dict(w=ip/(2*I0*A*A*(1+A*A)),mu=mu,closure=x/A+1,det=det)

def cosmo(T,gs,gr,tnu,lam):
    a=TG0/T*(GS0/gs)**(1/3)
    b=background(a,lam)
    rho_r=np.pi**2/30*gr*T**4
    H=np.sqrt(8*np.pi/3*(rho_r+b['rho']))/MPL
    p,o,lo,width=band(b['A'])
    k2=lo+width/2
    f=floquet(b['A'],k2)
    return dict(T_MeV=T*1e3,gs=gs,gr=gr,a=a,**b,H=H,H_over_omega=H/b['omega'],neff_equiv=b['rho']/(7*np.pi**2/120*(T*tnu)**4),rho_over_rad=b['rho']/rho_r,klo_over_m=np.sqrt(lo),khi_over_m=np.sqrt(lo+width),qmid_eV=a*M*1e9*np.sqrt(k2),mu_over_H=f['mu_over_m']*M/H,**f)

def mode_growth(T0,lam,gs=10.75,gr=10.75):
    a0=TG0/T0*(GS0/gs)**(1/3)
    b0=background(a0,lam)
    _,_,lo,width=band(b0['A'])
    k20=lo+width/2
    def band_diff(lna):
        b=background(a0*np.exp(lna),lam)
        _,_,l,w=band(b['A'])
        return k20*np.exp(-2*lna)-l
    end=brentq(band_diff,0,10,xtol=1e-12)
    def fn(lna):
        a=a0*np.exp(lna)
        T=T0*np.exp(-lna)
        b=background(a,lam)
        H=np.sqrt(8*np.pi/3*(np.pi**2/30*gr*T**4+b['rho']))/MPL
        k2=k20*np.exp(-2*lna)
        _,_,l,w=band(b['A'])
        if not l<k2<l+w:return 0.
        return floquet(b['A'],k2)['mu_over_m']*M/H
    vals=[]
    for n in [65,129]:
        grid=np.linspace(0,end,n)
        vals.append(float(simpson([fn(t) for t in grid],x=grid)))
    return dict(T0_MeV=T0*1e3,ln_a_exit=end,T_exit_MeV=T0*np.exp(-end)*1e3,linear_exponent=vals[-1],rel_grid=abs(vals[-1]-vals[0])/vals[-1],q_eV=a0*M*1e9*np.sqrt(k20))

def main():
    checks=[]
    for A in [.05,.1,.3,1,3,10,100]:
        p,o,l,w=band(A)
        res=[floquet(A,l+w*x) for x in [-.2,0,.5,1,1.2]]
        I0,Ic,_,_=integrals(A)
        assert abs(I0-ellipk(p)/o)<2e-12
        alt=independent(A)
        assert abs(alt['w']-(2*Ic/((1+A*A)*I0)-1))<1e-11
        assert abs(alt['mu']/res[2]['mu_over_m']-1)<1e-5
        checks.append(dict(A=A,lo=l,width=w,values=res,independent=alt))
    cases=[]
    for lam in LAMBDAS:
        rows=[cosmo(T*1e-3,gs,gr,tn,lam) for T,gs,gr,tn in [(20,10.75,10.75,1),(5,10.75,10.75,1),(1,10.75,10.75,1),(.1,3.91,3.36,(4/11)**(1/3)),(.05,3.91,3.36,(4/11)**(1/3))]]
        cases.append(dict(lam=lam,rows=rows,mode_growth=mode_growth(.02,lam)))
    print(json.dumps(dict(m_GeV=M,checks=checks,cases=cases),indent=2))

if __name__=='__main__':main()
