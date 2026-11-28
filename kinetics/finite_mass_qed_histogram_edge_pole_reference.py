"""Local coherent spacelike subtraction at an on-shell histogram step."""
import gzip
import json
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import roots_legendre

from early_pole_hard_yukawa_photoproduction_reference import I,G,dot,slash,spinors
from finite_mass_qed_nonnull_photon_ward_poles_reference import setup


def analytic_moments(a,b,eps,n):
    m0=(np.arctan(b/eps)-np.arctan(a/eps))/eps
    m1=.5*np.log((b*b+eps*eps)/(a*a+eps*eps))
    z=[m0,m1]
    for j in range(2,n+1):z.append((b**(j-1)-a**(j-1))/(j-1)-eps**2*z[j-2])
    return np.array(z)


def run_case(grid,d):
    a,b=-.35,.44
    k,w,c=.00065,.00013,-.5
    M,m,q,K,T,Phi=setup(d,k,w,c)
    E0_name,e_idx=('fermion',32) if grid=='fine' else ('parent',48)
    E0=np.hypot(M if grid=='third' else m,d[f'endpoint_{E0_name}_momentum_edges'][e_idx])
    if E0_name=='fermion':E0+=T[0]
    h=2*dot(q,K)
    K0=dot(K,K)
    def f(phi,E):
        p,_=Phi(E,phi)
        return (K0+2*dot(p,K))/h
    ph=brentq(lambda x:f(x,E0),0,np.pi,xtol=2e-15)
    V=np.linalg.norm(T[1:])
    axis=T[1:]/V
    ax=np.array([axis[2],0.,-axis[0]])
    Q0=(M*M+dot(T,T)-m*m)/2
    def row(x):
        e=E0 if x==0. else brentq(lambda y:f(ph,y)-x,E0-.002,E0+.002,xtol=3e-16)
        p,bv=Phi(e,ph)
        par=(e*T[0]-Q0)/V
        transverse=np.sqrt((e-M)*(e+M)-par**2)
        dp=axis[2]*T[0]/V+ax[2]*np.cos(ph)*(e-par*T[0]/V)/transverse
        jac=h/(2*(w-k*dp))
        if jac<=0:raise ArithmeticError('Nonmonotone pole chart')
        U=spinors(p,M);Vbar=spinors(bv,m).conj().T@G[0]
        A=slash(p+K)+M*I;B=slash(bv-K)+m*I
        vL=np.ravel([Vbar@A@slash(z)@U for z in np.eye(4)[1:3]])
        vR=np.ravel([Vbar@slash(z)@B@U for z in np.eye(4)[1:3]])
        C=vL@vR.conj()
        val=jac*np.r_[np.sum(abs(vL)**2),np.sum(abs(vR)**2),C.real,C.imag]
        i=[int(np.searchsorted(d[f'endpoint_{nm}_momentum_edges'],np.linalg.norm(v[1:]))-1)
           for nm,v in zip(('parent','fermion'),(p,bv),strict=True)]
        return val,i,e
    scale=row(0.)[0][0]
    def raw(x):
        return row(x)[0]/scale
    wgt=[]
    indices=[]
    for l,r in ((a,0.),(0.,b)):
        tests=[row(x)[1] for x in np.linspace(l+1e-7,r-1e-7,19)]
        if any(z!=tests[0] for z in tests):
            raise ArithmeticError('Another occupation edge enters diagnosed side')
        indices.append(tests[0])
        fp,fb=[d[f'endpoint_{name}_occupations'][j] for j,name in zip(tests[0],('parent','fermion'),strict=True)]
        wgt.append(fp*(1-fb))
    if min(wgt)<=0 or wgt[0]==wgt[1] or abs(f(ph,E0))>3e-10:raise ArithmeticError('Missing empirical step')
    # Normalize to mean weight, which makes the diagnostic coefficients portable.
    wgt=np.asarray(wgt)/np.mean(wgt)
    sample=np.linspace(a,b,40)
    r0=raw(0.)
    def fit(n):
        x=(a+b+(b-a)*np.cos(np.pi*(np.arange(n+4)+.5)/(n+4)))/2
        z=np.array([raw(u) for u in x])
        p=[Polynomial.fit(x,z[:,i],n).convert() for i in range(4)]
        p=[v+(r0[i]-v(0.)) for i,v in enumerate(p)]
        defect=max(abs((np.array([v(sample) for v in p]).T-np.array([raw(x) for x in sample])).ravel()))
        if defect>3e-8:raise ArithmeticError(('Smooth side numerator accuracy',defect))
        return p,float(defect)
    p0,e0=fit(6);p1,e1=fit(10)
    xg,wg=roots_legendre(64)

    def scalar_quad(f,l,r):
        return (r-l)/2*(wg@f((r+l)/2+xg*(r-l)/2))

    def finite(p):
        f,g,c,i=p
        deriv=f.deriv()(0.)+2*c(0.)
        logcoef=(wgt[1]-wgt[0])*deriv
        tails=-f(0.)*(wgt[0]/(-a)+wgt[1]/b)
        B=deriv*(wgt[1]*np.log(b)-wgt[0]*np.log(-a))
        im=np.pi*sum(wgt)*i(0.)
        fp=Polynomial(f.coef[2:])
        cp=Polynomial(c.coef[1:])
        R=sum(W*scalar_quad(lambda x:fp(x)+2*cp(x)+g(x)/(1-x)**2+2*c(x)/(1-x),l,r)
              for W,(l,r) in zip(wgt,((a,0.),(0.,b)),strict=True))
        return logcoef,tails+B+R+im
    climit,F=finite(p1)
    eps_list=(.015,.004,.001,2.5e-4,1e-5,4e-7,2e-9)
    def smooth_integral(p,e,edge=0.,remove_area=True):
        f,g,c,i=p
        total=0.
        for wt,(l,r) in zip(wgt,((a,edge),(edge,b)),strict=True):
            near=f+2*c*Polynomial([0.,1.])+2*i*e
            v=analytic_moments(l,r,e,len(near.coef))
            if remove_area:
                L=abs(l) if r==0. else abs(r)
                v[0]=-np.arctan(e/L)/e
                # Only remove f(0) from the leading Lorentzian; retain finite eps*Im terms.
                total+=wt*(near.coef@v[:len(near.coef)]+np.pi*i(0.))
            else:
                total+=wt*(near.coef@v[:len(near.coef)])
            total+=wt*scalar_quad(lambda x:(g(x)+2*(c(x)*(1-x)-i(x)*e))/((1-x)**2+e*e),l,r)
        return total
    def direct(e):
        def numer(x):
            f,g,c,i=raw(x)
            return f/(x*x+e*e)+g/((1-x)**2+e*e)+2*(complex(c,i)/((x+1j*e)*(1-x-1j*e))).real
        return (sum(W*quad(numer,l,r,epsabs=2e-9,epsrel=6e-12,limit=220)[0]
                    for W,(l,r) in zip(wgt,((a,0.),(0.,b)),strict=True))-np.pi/e)
    finite_rows=[]
    for e in eps_list:
        X=smooth_integral(p1,e)
        finite_rows.append(dict(width=e,area_subtracted=X,area_and_log_subtracted=X-climit*np.log(1/e),
                                difference_vs_finite_limit=X-climit*np.log(1/e)-F))
    checkmax=max(abs(direct(e)-smooth_integral(p1,e)) for e in (.012,.004))
    if (max(abs(np.array(finite(p0))-np.array([climit,F])))>2e-7 or checkmax>5e-7
            or abs(climit)<1e-3 or abs(finite_rows[-1]['difference_vs_finite_limit'])>2e-5):
        raise ArithmeticError(('Unresolved occupations at an edge',grid,climit,F,checkmax,finite_rows))
    moved=[]
    epsilon=2.5e-7
    for lam in (-4.,-1.,0.,1.,4.):
        x=epsilon*smooth_integral(p1,epsilon,edge=epsilon*lam,remove_area=False)
        y=wgt[0]*(np.pi/2+np.arctan(lam))+wgt[1]*(np.pi/2-np.arctan(lam))
        if abs(x-y)>2e-5:raise ArithmeticError(('Step/root displacement law',x,y))
        moved.append(dict(edge_displacement_in_widths=lam,scaled_measured=x,predicted_leading=y))
    def bridge_eps(e,ramp,cubic=False):
        f,g,c,i=p1
        X=Polynomial([0.,1.])
        z=(f+2*c*X+2*i*e).coef
        odd=Polynomial(np.where(np.arange(len(z))%2==1,2*z,0.))
        def reg(x):return (g(x)+2*(c(x)*(1-x)-i(x)*e))/((1-x)**2+e*e)
        def integrand(y):
            x=e*y
            z=x/ramp
            h=(1.5*z-.5*z**3) if cubic else z
            jump=(h-1.)*(wgt[1]-wgt[0])/2
            return jump*(odd(x)/(x*x+e*e)+reg(x)-reg(-x))*e
        r=ramp/e
        edges=np.unique(np.r_[0.,r,np.geomspace(.04,max(r,.04),90)])
        edges=edges[(edges>=0)&(edges<=r)]
        return sum(scalar_quad(integrand,l,b) for l,b in zip(edges[:-1],edges[1:],strict=True))
    shape=[]
    eps=8e-8
    for ratio in (.1,.5,1.,4.,16.,64.):
        v=bridge_eps(eps,ratio*eps)
        h=climit*(1.-np.arctan(ratio)/ratio-.5*np.log1p(ratio*ratio))
        if abs(v-h)>1e-7:
            raise ArithmeticError(('Linear smoothing limit failed',ratio,v,h))
        cub=bridge_eps(eps,ratio*eps,True)
        cubref=climit*(1.5*(ratio-np.arctan(ratio))/ratio
                      -.5*(ratio**3/3-ratio+np.arctan(ratio))/ratio**3-.5*np.log1p(ratio*ratio))
        if abs(cub-cubref)>1e-7:
            raise ArithmeticError(('Cubic smoothing limit failed',ratio,cub,cubref))
        shape.append(dict(smoothing_width_in_regulator_widths=ratio,measured_change=v,
                          cubic_measured_change=cub,cubic_exact_joint_limit=cubref,
                          cubic_absolute_defect=abs(cub-cubref),
                          exact_joint_limit=h,absolute_defect=abs(v-h)))
    sigma=.002
    vals=[]
    for e in (.003,.0003,3e-5,3e-6,3e-7):
        v=smooth_integral(p1,e)+bridge_eps(e,sigma)
        vals.append(dict(regulator=e,after_area_at_fixed_smoothing=v))
    final_cubic=smooth_integral(p1,3e-7)+bridge_eps(3e-7,sigma,True)
    cubic_minus_linear_offset=final_cubic-vals[-1]['after_area_at_fixed_smoothing']
    r=sigma/3e-7
    finite_profile_theory=.5*climit*((r-np.arctan(r))/r-(r**3/3-r+np.arctan(r))/r**3)
    if abs(cubic_minus_linear_offset-finite_profile_theory)>3e-6:
        raise ArithmeticError(('Smoothing profile difference',grid,cubic_minus_linear_offset,finite_profile_theory))
    finite_sigma_prediction=F+climit*(1+np.log(1/sigma))
    if abs(vals[-1]['after_area_at_fixed_smoothing']-finite_sigma_prediction)>4e-4:
        raise ArithmeticError(('Fixed smoothing saturation not seen',grid,vals,finite_sigma_prediction))
    return dict(grid=grid,k=k,w=w,cosine=c,actual_histogram=E0_name,edge_index=e_idx,edge_parent_energy=E0,
                exact_azimuth=ph,pole_chart_separation=h,root_x_residual=abs(f(ph,E0)),occupations_by_side=indices,
                normalized_forward_weights=wgt.tolist(),local_coordinate_range=[a,b],
                smooth_polynomial_error=max(e0,e1),log_coefficient_after_half_residue=climit,finite_after_two_subtractions=F,
                diagonal_log_coefficient=(wgt[1]-wgt[0])*p1[0].deriv()(0.),
                interference_log_coefficient=2*(wgt[1]-wgt[0])*p1[2](0.),
                independent_direct_width_error=checkmax,ladder=finite_rows,edge_displacement_fixture=moved,
                linear_smoothed_edge_joint_limit=shape,fixed_smooth_halfwidth=sigma,
                fixed_smoothing_area_subtracted=vals,small_smoothing_finite_limit_prediction=finite_sigma_prediction,
                cubic_minus_linear_at_fixed_smoothing=cubic_minus_linear_offset,sharp_limit_offset_difference=climit/3,finite_ratio_profile_offset_prediction=finite_profile_theory)


if __name__=='__main__':
    p=Path(__file__).parent/'results'
    a=[]
    for g in ('fine','third'):
        with gzip.open(p/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz','rt') as file:d=json.load(file)
        a.append(run_case(g,d))
    print(json.dumps(dict(scope='Diagnostic equal-width local pole edge; no matched finite-width plasma kernel.',data=a),
                     indent=2,allow_nan=False))
