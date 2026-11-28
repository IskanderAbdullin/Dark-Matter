"""Conservative three-species number/heat response of the decay/inverse pair."""
import argparse
import json
import numpy as np
from scipy.special import expit, roots_legendre, kve
from scipy.linalg import block_diag
from sm_dark_entropy_exchange_reference import StandardModel,MPL

M,m,mp,y=.1,.02,86e-9,.001
P=np.sqrt((M-m-mp)*(M+m+mp)*(M-m+mp)*(M+m-mp))/(2*M)
Es=(M*M-m*m+mp*mp)/(2*M)
width=y*y*((M+m)**2-mp*mp)*P/(8*np.pi*M*M)

def events(T,ne=80,na=40):
    z,w=roots_legendre(ne)
    u=(z+1)*np.sqrt(65.)/2
    e1=M+T*u*u
    p=u*np.sqrt(T*(e1+M))
    ephi_hi=(e1*Es+p*P)/M
    # Rational exact lower endpoint avoids subtracting two almost equal
    # boosts for phi close to massless and chi1 highly relativistic.
    # E_- E_+ = E_*^2 + (p m_phi / M)^2.
    ephi_lo=(Es*Es+(p*mp/M)**2)/ephi_hi
    x,v=roots_legendre(na)
    d=np.log(ephi_hi/ephi_lo)
    ep=ephi_lo[:,None]*np.exp(d[:,None]*(x+1)/2)
    # Lab-angle factor: d cos(theta*)/2 = M dEphi/(2 p p*).
    a_w=v[None,:]*d[:,None]*ep*M/(4*P*p[:,None])
    E1=np.broadcast_to(e1[:,None],ep.shape)
    E=np.stack((E1,E1-ep,ep),axis=-1)
    # g1=4 (two independent Dirac charge sectors); parent time dilation.
    moment=4/(2*np.pi**2)*M*width*p*(w*T*np.sqrt(65.)*u)
    weights=moment[:,None]*a_w
    assert np.max(abs(a_w.sum(axis=1)-1))<5e-12
    return E.reshape(-1,3), weights.ravel()

def flow(E,W,T,eta,Ts):
    x=E/np.asarray(Ts)-np.asarray(eta)
    assert np.all(x[:,2]>0.)
    f=expit(-x[:,:2])
    G=f[:,0]*(1-f[:,1])/(-np.expm1(-x[:,2])) # BE 1+f_phi from thermal affinity
    A=x[:,1]+x[:,2]-x[:,0]
    net=np.sign(A)*(-np.expm1(-abs(A)))*np.where(A>=0,G,G*np.exp(-A))
    v=vectors(E,T)
    return (W*net)@v, np.sum(W*net*A)

def vectors(E,T):
    v=np.stack((np.ones(len(E)),E[:,0]/T,np.ones(len(E)),E[:,1]/T,
                np.ones(len(E)),E[:,2]/T),axis=1)
    v[:,:2]*=-1
    return v

def susceptibility(T,mi,g,BE=False):
    z,w=roots_legendre(180)
    q=(z+1)*np.sqrt(65.)/2
    e=np.hypot(q*q,mi/T)
    measure=g/(2*np.pi**2)*T**3*w*np.sqrt(65.)*q**5
    if BE:
        f=1/np.expm1(e)
        n=f*(1+f)
    else:
        f=expit(-e)
        n=f*expit(e)
    return np.array([[measure@n,(measure*n)@e],[(measure*n)@e,(measure*n)@(e*e)]]),T*(measure*f)@e

def main(table):
    eos=StandardModel(table)
    rows=[]
    null=np.array([[1,0,1,0,0,0],[1,0,0,0,1,0],[0,1,0,1,0,1]]).T
    for T in (.010,.014,.020,.04,.1,.3,.6882487043,.95):
        E,W=events(T)
        v=vectors(E,T)
        G=expit(-E[:,0]/T)*expit(E[:,1]/T)/(-np.expm1(-E[:,2]/T))
        F=(v.T*(W*G))@v
        eq=float(np.linalg.norm(F@null)/np.linalg.norm(F))
        species=[susceptibility(T,mi,g,BE) for mi,g,BE in ((M,4,False),(m,4,False),(mp,1,True))]
        X=block_diag(*(q[0] for q in species))
        lhs=np.linalg.cholesky(np.linalg.inv(X))
        eigen=np.linalg.eigvalsh(lhs.T@F@lhs)
        assert eigen[3]>0 and max(abs(eigen[:3]))/eigen[3]<3e-11
        rh=eos.state(T)[0]+sum(q[1] for q in species)
        H=np.sqrt(rh/3)/MPL
        # Independently perturb both temperature and chemical forces; keep
        # the BE chemical potential below its upper limit in both directions.
        d=np.array([.4,.7,-.6,.2,-.03,-.55])*min(3e-6,mp/T)
        f0,p0=flow(E,W,T,d[::2],T*np.exp(d[1::2]))
        f1,p1=flow(E,W,T,-d[::2],T*np.exp(-d[1::2]))
        err=np.linalg.norm((f0-f1)/2+F@d)/np.linalg.norm(F@d)
        Er,Wr=events(T,108,56)
        Vr=vectors(Er,T)
        Pr=expit(-Er[:,0]/T)*expit(Er[:,1]/T)/(-np.expm1(-Er[:,2]/T))
        Fr=(Vr.T*(Wr*Pr))@Vr
        ref=np.linalg.norm(F-Fr)/np.linalg.norm(Fr)
        mb=4*M*M*T/(2*np.pi**2)*width*kve(1,M/T)*np.exp(-M/T)
        vac=expit(-E[:,0]/T)@W
        classical=np.exp(-E[:,0]/T)@W
        assert max(eq,err,ref,abs(classical/mb-1))<8e-7,(T,eq,err,ref,classical/mb-1)
        assert p0>=0 and p1>=0
        # Put the exact scalar heat column in the three-dimensional
        # energy-only, fixed-chemical-affinity subspace as a separate diagnostic.
        rows.append(dict(T_GeV=T,nonzero_modes_per_H=(eigen[3:]/H).tolist(),
                        boosted_vacuum_FD_over_MB=vac/mb,
                        total_equilibrium_forward_over_vacuum_FD=float((W@G)/vac),
                        invariance_residual=eq,two_transport_coordinates_error=float(err),
                        phase_space_refinement_error=float(ref)))
    return dict(width_GeV=width,rows=rows,scope="Local FD/BE decay-inverse response near each common equilibrium; no portal, new kinetic closure, condensate, or cosmological trajectory.")

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos")
    print(json.dumps(main(p.parse_args().eos),indent=2,allow_nan=False))
