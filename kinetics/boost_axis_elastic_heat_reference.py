"""Independent elastic heat integral in invariant mass and boost-axis angles."""
import json
import numpy as np
from scipy.special import roots_legendre,expit
from timelike_chemical_heat_eigenmodes_reference import leptons
from timelike_portal_onsager_coverage_reference import V,mchi,gprod

def integrate(Ts,Td,eta,m,*,n=20,cut=29.,linear=False):
    assert Ts>0 and Td>0 and (not linear or Ts==Td)
    # Dimensionless momenta in powers of Ts; restore Ts^9 for Q.
    M,u0=mchi/Ts,m/Ts
    eps=Td/Ts
    scale=max(1.,eps)
    mv2=(V/Ts)**2
    G=gprod/V**2
    r,w=roots_legendre(n)
    z=(1+r)*np.sqrt(cut*scale)/2
    dz=w*np.sqrt(cut*scale)*z
    # Orient both CM polar axes along the total lab boost. Then dcos(out)
    # dphi is rotation invariant, and t depends only on their relative angle.
    c=r
    si=np.sqrt(1-c*c)
    az=(np.arange(2*n)+.5)*np.pi/n
    cosg=c[:,None,None]*c[None,:,None]+si[:,None,None]*si[None,:,None]*np.cos(az)
    polar_w=w[:,None]*w[None,:]
    result=0.
    for u,du in zip(M+u0+z*z,dz):
        s=u*u
        pcm2=((u-M-u0)*(u+M+u0)*(u-M+u0)*(u+M-u0))/(4*s)
        pcm=np.sqrt(pcm2)
        t=-2*pcm2*(1-cosg)
        S=s-M*M-u0*u0
        spin=2*G*G*(S*S+(S+t)**2+2*(M*M+u0*u0)*t)/(1-t/mv2)**2
        avg=spin.mean(axis=2)*polar_w
        L=u+z*z
        P=np.sqrt((L-u)*(L+u))
        E=(L[:,None]*(s+M*M-u0*u0)/(2*u)+P[:,None]*pcm*c)/u
        W=L[:,None]-E
        fa=expit(eta-E/eps)
        fb=expit(-W)
        fwd=(fa*fb)[:,:,None]*((1-fa)*(1-fb))[:,None,:]
        d=E[:,None,:]-E[:,:,None]
        if linear:
            signed=fwd*d*d
        else:
            signed=(fwd+fwd.transpose(0,2,1))*d*np.tanh(d*(1/eps-1)/2)
        jac=P*pcm2/u
        result+=du*np.sum((dz*jac)*np.sum(signed*avg,axis=(1,2)))
    # Both species count particle and antiparticle spin, 4 x 4. The term
    # for a time-reversed event is integrated over the same domain.
    return result*Ts**9*16/(512*np.pi**5)

def test():
    xs=[]
    for T,r,eta in ((.03,.71,1.1),(.15,.85,.25),(.6882487043,1.,0.),(.6882487043,1.05,.13)):
        a=np.array([integrate(T,T*r,eta,m,n=18,linear=r==1.) for m in leptons])
        b=np.array([integrate(T,T*r,eta,m,n=25,linear=r==1.) for m in leptons])
        ref=abs(a.sum()/b.sum()-1)
        bc=np.array([integrate(T,T*r,eta,m,n=25,cut=36.,linear=r==1.) for m in leptons])
        cut=abs(bc.sum()/b.sum()-1)
        assert max(ref,cut)<2e-5,(T,ref,cut)
        if r==1.:
            j=3e-5
            up=sum(integrate(T,T*(1-j),eta,m,n=25) for m in leptons)
            down=sum(integrate(T,T*(1+j),eta,m,n=25) for m in leptons)
            linear=abs((up-down)/(2*j*b.sum())-1)
            assert linear<1e-6
            assert sum(integrate(T,T,eta,m) for m in leptons)==0.
        else:
            linear=None
        assert (b.sum()>0) if r<=1. else (b.sum()<0)
        xs.append(dict(T_SM=T,ratio=r,eta=eta,quadrature_shift=ref,
                       cutoff_shift=cut,thermal_linear_response_error=linear))
    return dict(rows=xs,scope="Charged leptons only; finite-vector elastic vacuum heat.")

if __name__=="__main__":
    print(json.dumps(test(),indent=2,allow_nan=False))
