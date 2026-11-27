"""Bare virtual longitudinal Ward checks with two unequal-mass legs."""
import json

import numpy as np

from early_pole_hard_yukawa_photoproduction_reference import G,I,dot,slash


def inverse(p,m):
    return slash(p)-m*I


def numerator(p,m):
    return slash(p)+m*I


def propagator(p,m):
    return numerator(p,m)/(dot(p,p)-m*m)


def trace(a):
    return np.trace(a)


def carry(P,b,M,m):
    Pspin,bspin=numerator(P,M),numerator(b,m)
    born=trace(bspin@Pspin)
    def probe(k):
        K=slash(k)
        vertex=trace(bspin@K@propagator(b-k,m)@propagator(P-k,M)@K@Pspin)/born
        def leg(p,m):
            s=numerator(p,m)
            norm=trace(s@G[0])
            algebra=[]
            def f(dp,sign):
                X=p+np.array([dp,0,0,0])
                L=X-sign*k
                D=inverse(X,m)
                S=propagator(L,m)
                direct=K@S@K
                expand=D@S@D-D-sign*K
                algebra.append(float(np.max(abs(direct-expand))))
                return trace(s@direct)/norm
            def even(dp):
                return (f(dp,1)+f(dp,-1))/2
            residue=[(even(h)-even(-h))/(2*h) for h in (3.e-5,1.e-5,3.e-6)]
            return even(0.),residue,max(algebra)
        pm,pd=(leg(P,M),leg(b,m))
        ward=max(abs(vertex+(pm[1][j]+pd[1][j])/2) for j in (1,2))
        def extrapolate(z):
            h1,h2=1.e-5,3.e-6
            return (h1*h1*z[2]-h2*h2*z[1])/(h1*h1-h2*h2)
        refined=abs(vertex+(extrapolate(pm[1])+extrapolate(pd[1]))/2)
        if refined>2e-8:
            raise ArithmeticError('Extrapolated pole residue did not cancel the common vertex.')
        if max(abs(vertex-1),abs(pm[0]),abs(pd[0]),pm[2],pd[2])>2e-10 or ward>2e-6:
            raise ArithmeticError('Longitudinal virtual and residue Ward algebra failed.')
        return dict(loop_momentum=k.tolist(),normalized_longitudinal_vertex=vertex.real,
                    symmetrized_onshell_shift=[pm[0].real,pd[0].real],
                    parent_residue_derivative=[z.real for z in pm[1]],
                    daughter_residue_derivative=[z.real for z in pd[1]],
                    maximum_matrix_ward_defect=max(pm[2],pd[2]),
                    common_vertex_and_residue_defect=float(ward),
                    extrapolated_common_defect=float(refined))
    def cut_check(momentum,momentum_mass,partner,partner_mass,is_parent,sign,v):
        L=np.r_[sign*np.hypot(np.linalg.norm(v),momentum_mass),v]
        k=momentum-L
        K=slash(k)
        N=numerator(L,momentum_mass)
        U=propagator(partner-k,partner_mass)
        T=K@U@N@K if is_parent else K@N@U@K
        cut=trace(bspin@T@Pspin)/born
        leg=trace(numerator(momentum,momentum_mass)@K@N@K)
        D=inverse(momentum,momentum_mass)
        residual=float(np.max(abs(K@N@K-D@N@D)))
        if max(abs(cut),abs(leg),residual)>2e-10:
            raise ArithmeticError('One bare longitudinal charged cut did not vanish.')
        return dict(cut_leg='parent' if is_parent else 'daughter',energy_sign=sign,
                    internal_charge_spatial_momentum=v.tolist(),
                    projected_vertex_cut_abs=float(abs(cut)),projected_self_energy_cut_abs=float(abs(leg)),
                    unprojected_longitudinal_cut_identity_abs=residual)
    loops=[probe(np.asarray(k)) for k in ((.18,.21,-.16,.06),(.54,-.18,.34,.27),
                                          (-.09,.07,.29,-.35),(.09,.2,-.1,-.75))]
    cuts=[cut_check(p,a,q,b,is_parent,sg,np.asarray(v))
          for p,a,q,b,is_parent in ((P,M,b,m,True),(b,m,P,M,False))
          for sg in (-1,1) for v in ((.14,-.08,.22),(-.43,.19,.08),(.074,.66,-.16))]
    return dict(parent=P.tolist(),daughter=b.tolist(),vacuum_longitudinal_loops=loops,charged_on_shell_cuts=cuts)


def run():
    M,m,ms=1.,.06,.0174
    e=(M*M+m*m-ms*ms)/(2*M)
    v=np.sqrt((e-m)*(e+m))*np.array([.6,.0,.8])
    P=np.array([M,0.,0.,0.])
    b=np.r_[e,v]
    eta=.54
    def boost(a):
        out=a.copy()
        out[0]=np.cosh(eta)*a[0]+np.sinh(eta)*a[3]
        out[3]=np.sinh(eta)*a[0]+np.cosh(eta)*a[3]
        return out
    return dict(scope='Ward algebra with a common smooth longitudinal loop regulator and bare on-shell charged cuts; not a real-time double-pole or screened medium calculation.',
                vertex_and_residue=[carry(P,b,M,m),carry(boost(P),boost(b),M,m)])


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
