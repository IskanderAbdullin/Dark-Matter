"""Signed two-temperature FD/BE heat on one exact below-pole event measure."""
import json
import numpy as np
from scipy.special import expit

from continuum_dark_normal_capacity_reference import fugacity, moments, rule
from soft_scalar_elastic_graph_reference import spin
from three_species_decay_onsager_reference import m,mp,M


def bdist(z):
    # log(f) and log(1+f); z>0 for the normal BE branch.
    l=-np.log(-np.expm1(-z))
    return l-z,l


def fdist(z):
    return -np.logaddexp(0.,z),-np.logaddexp(0.,-z)


def nodes(T,mass,n,tail):
    x,w=rule(n)
    e=np.sqrt(np.unique(np.clip([0.,.001,.015,.1,.5,1.5,4.,9.,18.,36.,tail],0.,tail)))
    y=(e[:-1,None]*(1-x)+e[1:,None]*(1+x))/2
    k=(T*y*y).ravel()
    E=k+mass
    p=np.sqrt(k*(k+2*mass))
    W=(w*(e[1:]-e[:-1])[:,None]*T*y).ravel()*p*E/(2*np.pi**2)
    return k,E,p,W


def heat(Tp,Tc,ep,ec,orders=(8,8,10,4),cut=.82,tail=60.,stats=False):
    """Return Q_phi (>0 is scalar heating); includes FD/BE quantum factors."""
    n,im,om,na=orders
    phi=nodes(Tp,mp,n,tail)
    K,E,p,Wc=nodes(Tc,m,n,tail)
    le,be=fdist(K/Tc-ec)
    cm,wm=rule(im)
    co,wo=rule(om)
    sin=np.sqrt(np.maximum(0.,1-co*co))
    cs=np.cos((np.arange(na)+.5)*2*np.pi/na)
    temp_difference=1/Tp-1/Tc
    res=np.zeros(3)
    for k,F,u,Wp in zip(*phi):
        lpi,bpi=bdist(k/Tp-ep)
        d=k*E[:,None]+mp*K[:,None]-u*p[:,None]*cm
        s=(m+mp)**2+2*d
        ok=s<cut*M*M
        q2=d*(2*m*mp+d)/s
        q=np.sqrt(np.maximum(0.,q2))
        root=np.sqrt(s)
        L=F+E[:,None]
        P=np.sqrt(np.maximum(0.,L*L-s))
        rest=(mp*(m+mp)+d)/root
        center=L*rest/root
        B=P*q/root
        lo=(rest*rest+(P*mp/root)**2)/(center+B)
        # Avoid a subtraction of lab rest energies in the outgoing BE tail.
        outk=((lo-mp)[:,:,None]+B[:,:,None]*(co+1))
        outk=np.maximum(outk,np.finfo(float).tiny)
        delta=outk-k
        ecprime=(K[:,None,None]-delta)/Tc-ec
        lpf,bpf=bdist(outk/Tp-ep)
        ll,bb=fdist(ecprime)
        lf=lpi+le[:,None,None]+bpf+bb
        lr=lpf+ll+bpi+be[:,None,None]
        log_ratio=delta*temp_difference
        # Check detailed balance pointwise where both logs are representable.
        assert np.max(abs((lf-lr)-log_ratio))<6e-10
        logmax=np.maximum(lf,lr)
        net=np.sign(log_ratio)*np.exp(logmax)*(-np.expm1(-abs(log_ratio)))
        ci=np.clip(np.divide(F-center,B,out=np.zeros_like(B),where=B>0),-1.,1.)
        si=np.sqrt(np.maximum(0.,1-ci*ci))
        amplitude=0.
        s0=np.where(ok,s,(m+mp)**2)
        for az in cs:
            cos=ci[:,:,None]*co+si[:,:,None]*sin*az
            tr=-2*q2[:,:,None]*(1-cos)
            amplitude+=spin(s0[:,:,None],np.where(ok[:,:,None],tr,0.))/na
        coeff=Wp*4*Wc[:,None]*(wm/2)*q*root/(F*E[:,None])/(32*np.pi*s)
        raw=coeff[:,:,None]*amplitude*wo*ok[:,:,None]
        res[0]+=np.sum(raw*net*delta)/2
        if stats:
            res[1]+=np.sum(raw*np.exp(lf)*delta*delta)/(2*Tp*Tp)
            res[2]+=np.sum(raw*np.exp(lf))
    if stats: return tuple(float(x) for x in res)
    return float(res[0])


def test():
    small=(8,8,10,4)
    fine=(13,12,16,6)
    cases=[]
    for Tp,Tc in [(1e-4,1e-4),(7e-5,1.4e-4),(5e-5,1e-5),(1.8e-5,4e-8)]:
        n=np.array([2e-17,5e-15])
        ep=fugacity(Tp,mp,n[0],True,20)
        ec=fugacity(Tc,m,n[1],False,20)
        a=heat(Tp,Tc,ep,ec,small)
        b=heat(Tp,Tc,ep,ec,fine)
        assert (a==b==0.) if Tp==Tc else np.sign(a)==np.sign(Tc-Tp)==np.sign(b)
        if Tp!=Tc: assert abs(a/b-1)<.005
        cases.append(dict(scalar_T=Tp,fermion_T=Tc,scalar_heat_GeV5=b,quadrature_change=float(abs(a/b-1.)) if b else 0.))
    T=1e-4
    n=np.array([2e-17,5e-15])
    ep=fugacity(T,mp,n[0],True,20)
    ec=fugacity(T,m,n[1],False,20)
    g=heat(T,T,ep,ec,fine,stats=True)[1]
    def f(p):
        return heat(p,T,fugacity(p,mp,n[0],True,20),ec,fine)
    h=T*2e-5
    err=abs((f(T-h)-f(T+h))/(2*h*g)-1.)
    from late_scalar_fermion_heat_rate_reference import rates
    ref=rates(T,n,ep,ec,14,12,16,6)["equilibrium_energy_conductance_GeV4"]
    assert max(err,abs(g/ref-1.))<1e-5
    return dict(cases=cases,finite_difference_error=err,previous_equilibrium_error=abs(g/ref-1.))


if __name__=="__main__":
    print(json.dumps(test(),indent=2,allow_nan=False))
