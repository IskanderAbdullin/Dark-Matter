"""Finite-width vacuum e/mu/tau with FD blocking and spin correlations."""
import json
import numpy as np
from scipy.special import expit, roots_legendre
from scipy.integrate import quad_vec
from timelike_portal_onsager_coverage_reference import width, V, mchi, factor
from timelike_chemical_heat_eigenmodes_reference import leptons, rates

def integrator(Td,Ts,eta,ny=36,nx=44,*,equilibrium_spatial_moment=False):
    if equilibrium_spatial_moment and (Td!=Ts or eta!=0.):
        raise ValueError("Equilibrium spatial response requires detailed balance.")
    hot=max(Td,Ts)
    z,w=roots_legendre(ny)
    y=(z+1)*np.sqrt(55.)/2
    dw=w*np.sqrt(55.)*hot*y
    x,wx=roots_legendre(nx)
    wx/=2
    p2=(3*x*x-1)/2
    half=(x[None,:]/2)
    def vals(E,_):
        s=E*E
        vi2=1-4*mchi*mchi/s
        if vi2<=0.: return np.zeros(12)
        qq=E+hot*y*y
        p=np.sqrt(hot*y*y*(qq+E))
        e=qq[:,None]/2
        boost=p[:,None]*half
        a=expit((e+np.sqrt(vi2)*boost)/Td-eta)*expit((e-np.sqrt(vi2)*boost)/Td-eta)
        a0,a2=a@wx,a@(wx*p2)
        ld=2*eta-qq/Td
        ls=-qq/Ts
        delta=2*eta+qq*((Td-Ts)/(Td*Ts))
        q_w=dw*p/E
        b=factor*(s+2*mchi*mchi)*E**3*np.sqrt(vi2)
        out=np.zeros((4,3))
        weights=np.array([np.exp(ls)*q_w,
                          np.sign(delta)*np.exp(np.maximum(ls,ld))*(-np.expm1(-abs(delta)))*q_w])
        # The opt-in third column uses |p_pair|^2 instead of E_lab^2.
        # Build it from the invariant constraint in this independent quadrature.
        last=(qq-E)*(qq+E) if equilibrium_spatial_moment else qq*qq
        qpowers=np.stack((np.ones_like(qq),qq,last),axis=1)
        for m in leptons:
            vf2=1-4*m*m/s
            if vf2<=0.: continue
            b_m=b*(1+2*m*m/s)*np.sqrt(vf2)
            bl=expit((e+np.sqrt(vf2)*boost)/Ts)*expit((e-np.sqrt(vf2)*boost)/Ts)
            c0,c2=bl@wx,bl@(wx*p2)
            f2=(2/3)*vi2*vf2/(3-vi2-vf2+vi2*vf2/3)
            iso=a0*c0
            cor=f2*a2*c2
            out[0]+=b_m*weights[0]@qpowers  # equilibrium vacuum MB normalization
            out[1]+=b_m*(weights[0]*(iso+cor))@qpowers
            out[2]+=b_m*(weights[0]*cor)@qpowers
            out[3]+=b_m*(weights[1]*(iso+cor))@qpowers
        return out.ravel()
    ab=(2*mchi,max(60*hot,V+30*hot))
    points=[*ab,V,*[2*m for m in leptons],*[V+j*width for j in (-300,-100,-30,-10,-3,-1,1,3,10,30,100,300)],*[hot*j for j in (.1,.2,.5,1,2,3,5,7,10,15,22,30,50)]]
    bw=V*width
    # Log-E keeps low-mass thresholds well separated at low T.
    boundaries=np.log(sorted(set(x for x in points if ab[0]<=x<=ab[1])))
    def f(log_E):
        E=np.exp(log_E)
        return 2*E*E/((E*E-V*V)**2+bw*bw)*vals(E,hot)
    return sum((quad_vec(f,l,r,epsrel=4e-9,epsabs=1e-180)[0]
                for l,r in zip(boundaries[:-1],boundaries[1:])),np.zeros(12)).reshape(4,3)

def test():
    results=[]
    for T in (.06,.1,.3,.45,.6882487043,1.,2.):
        v=integrator(T,T,0.)
        ref=integrator(T,T,0.,ny=52,nx=58)
        err=float(np.max(abs(v[:2]/ref[:2]-1)))
        ermb=float(np.max(abs(v[0]/rates(T)-1)))
        assert err<2e-6 and ermb<2e-7,(T,err,ermb)
        assert np.all(v[3]==0.)
        dt=2e-6
        uu=integrator(T*(1+dt),T,.3*dt)[3,:2]
        ud=integrator(T*(1-dt),T,-.3*dt)[3,:2]
        jac=np.array([[v[1,0],v[1,1]],[v[1,1],v[1,2]]])
        got=(uu-ud)/(2*dt)
        ex=jac@np.array([.6,1/T])
        ed=float(np.max(abs(got/ex-1)))
        assert ed<5e-7,(T,ed)
        assert np.linalg.det(jac)>0
        results.append(dict(T_GeV=T,fd_over_mb=(v[1]/v[0]).tolist(),
                            spin_correlation_fraction=(v[2]/v[1]).tolist(),
                            q0_and_angle_refinement_max_relative_error=err,
                            vacuum_mb_independent_bessel_max_relative_error=ermb,
                            two_force_response_max_relative_error=ed))
    return dict(rows=results,scope="Direct e/mu/tau channels, spin-correlated 2-to-2 vacuum matrix element and exact FD occupancies; no medium self-energy, real-mediator double counting, hadrons, spacelike exchange, or evolved FD trajectory.")
if __name__=="__main__":
    print(json.dumps(test(),indent=2,allow_nan=False))
