"""Resolved bare self-energy controls with fixed external positive-energy spinors."""
import gzip
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import roots_legendre

from early_pole_decay_mode_clock_reference import intervals,Q
from early_pole_hard_yukawa_photoproduction_reference import G,I,MET,dot,slash

D=np.longdouble


def charged(m,E,delta,sign):
    m,E,delta=map(D,(m,E,delta))
    t=E+delta
    p=np.sqrt((E-m)*(E+m))
    r=np.sqrt((t-m)*(t+m))
    R=p*r
    if sign==1:
        # (E t - m² - R)(E t - m² + R) = m² (E-t)².
        wide=E*t-m*m+R
        narrow=m*m*delta*delta/wide
        near,far=-2*narrow,-2*wide
        sep=delta
    else:
        near=2*(E*t+m*m+R)
        far=2*m*m*(E+t)**2/(E*t+m*m+R)
        sep=-t-E
    J1=np.log(near/far)/(2*R)
    J2=2/(far*near)
    F=(2+2*m*m*J1)/E
    deriv=2*sep*(J1+2*m*m*J2)/E
    # Reduce the apparent double inverse analytically before subtraction.
    if sign==1:
        regular=(t+E)/((r+p)*E)+r*sep*J1/E
    else:
        regular=r*(1/sep+sep*J1)/E
    def offshell(h):
        s=-2*sep*h+h*h
        return (2+(2*m*m-s)*np.log((near+s)/(far+s))/(2*R))/E
    return F,deriv,regular,offshell,near,far


def independent_angular(m,E,delta,sign):
    t=E+delta
    p=np.sqrt((E-m)*(E+m))
    r=np.sqrt((t-m)*(t+m))
    # Different from the discriminant factorization: positive cut near
    # collinearity is derived independently using the rapidity difference.
    b=2*m*m*np.sinh((np.arcsinh(p/m)-np.arcsinh(r/m))/2)**2 if sign==1 else -E*t-p*r-m*m
    R=p*r
    width=abs(b/R)
    cuts=[0.,*[j*width for j in (1,4,16,64,256,1024) if j*width<2],2.]
    def kernel(u):
        gap=-2*(b+R*u)
        return 2*(sign*t-E)*(gap+2*m*m)/(E*gap**2)
    chunks=[quad(kernel,0,cuts[1],epsabs=3e-11,epsrel=1e-8)[0]]
    for a,b1 in zip(cuts[1:-1],cuts[2:],strict=True):
        chunks.append(quad(lambda x: np.exp(x)*kernel(np.exp(x)),np.log(a),np.log(b1),
                           epsabs=3e-11,epsrel=1e-8)[0])
    angular=fsum(chunks)
    P=np.array([E,0.,0.,p])
    trace=[]
    for u in (-.6,.2,.97):
        L=np.array([sign*t,r*np.sqrt(1-u*u),0.,r*u])
        N=slash(L)+m*I
        actual=np.trace((slash(P)+m*I)@sum(g*z@N@z for g,z in zip(MET,G,strict=True))).real/(4*E)
        target=(-2*dot(P,L)+4*m*m)/E
        trace.append(abs(actual-target)/max(1e-10,abs(target)))
    return angular,max(trace)


def photon_finite_remainder(m,E):
    beta=np.sqrt((E-m)*(E+m))/E
    if beta<1e-5:
        return 4*(1+beta**2/3+beta**4/5)/(E*E)
    return 2*(np.log1p(beta)-np.log(m*m/(E*E)/(1+beta)))/(E*E*beta)


def photon_checks(m,E):
    p=np.sqrt((E-m)*(E+m))
    beta=p/E
    L=np.log(m*m/E/(E+p));R=np.log1p(beta)
    flex=photon_finite_remainder(m,E)
    results=[]
    for k in (1e-6,5e-4,.01,.8,4.):
        def eval_direct(y,h):
            x=k*E*np.exp(y)
            return sum((2*(m*m+sg*x)/E-2*h)/
                       (-2*sg*x+2*(E-sg*k)*h+h*h) for sg in (-1,1))*np.exp(y)/beta
        h=1e-10*k*m*m/E**2
        proj=quad(lambda y:eval_direct(y,0.),L,R,epsabs=1e-10)[0]
        slope=quad(lambda y:eval_direct(y,1j*h).imag/h,L,R,epsabs=1e-10)[0]
        error=max(abs(proj/(-4/E)-1),abs(slope/(-4/k**2+flex)-1))
        if error>2e-9:
            raise ArithmeticError('Paired null photon projection or derivative did not match.')
        results.append(dict(k=k,angular_onshell=proj,angular_energy_derivative=slope,relative_defect=error))
    # Rayleigh-Jeans singular derivative is omitted here by explicit
    # subtraction, not assigned a counterterm in a real-time rate.
    thermal_moment=np.pi**2/6
    numeric=quad(lambda k:k*np.exp(-k)/(-np.expm1(-k)),0,np.inf,epsabs=1e-10)[0]
    if abs(numeric-thermal_moment)>2e-10:
        raise ArithmeticError('The finite remainder failed its independent Bose moment.')
    return dict(angular_remainder_after_minus4_over_k2=flex,
                bare_onshell_at_T1=-(2/E)*thermal_moment,
                bare_regular_derivative_at_T1=flex*thermal_moment/2,controls=results)


def run_radial(m,E,edges,occup,sign,tuned):
    p=np.sqrt((E-m)*(E+m))
    A=p/E if sign==1 else 0.
    rho=3e-4
    grid=np.asarray(edges)-E
    if tuned:
        idx=np.argmin(abs(grid))
        if not 0<idx<len(grid)-1 or abs(grid[idx])>1e-14:
            raise ArithmeticError('Expected declared archived step in the external node.')
        grid[idx]=0.
    else:
        idx=None
    jump=float(occup[idx]-occup[idx-1]) if idx else 0.
    predicted=A*jump
    def func(z, subtract):
        F,der,reg,*_=charged(m,E,z,sign)
        r=np.sqrt((E+z-m)*(E+z+m))
        return np.array([float(r*F/2),float(reg if subtract else r*der/2)])
    def accumulate(n,eps):
        nodes,weights=roots_legendre(n)
        logs=[]
        cuts=[]
        for a,b,f in zip(grid[:-1],grid[1:],occup,strict=True):
            if not f:
                continue
            if eps is None and sign==1:
                logs.append(A*f*np.log((abs(b) or rho)/(abs(a) or rho)))
            intervals=[(a,b)]
            if eps is not None and sign==1:
                intervals=[(a,min(b,-eps)),(max(a,eps),b)]
            for a,b in intervals:
                if b<=a:
                    continue
                geo=[z for r in (1e-12,1e-10,1e-8,1e-6,1e-4,.002,.02,.2)
                     for z in (-r,r) if a<z<b]
                bounds=sorted({a,b,*geo,*([0.] if a<0<b else [])})
                for lo,hi in zip(bounds[:-1],bounds[1:],strict=True):
                    phi=np.pi*(nodes+1)/2
                    zs=lo+(hi-lo)*np.sin(phi/2)**2
                    jac=np.pi*(hi-lo)*np.sin(phi)/4
                    vals=np.array([func(z,eps is None) for z in zs])
                    cuts.append(f*((jac*weights)@vals))
        result=np.sum(cuts,axis=0)
        result[1]+=fsum(logs)
        return result,fsum(logs)
    a,lead=accumulate(26,None)
    b,_=accumulate(44,None)
    resolved=float(max(abs(a-b)))
    if resolved>3e-6:
        raise ArithmeticError(f'Regularized archive quadrature did not resolve: {a,b}.')
    ladder=[]
    for eps in (1e-6,1e-8,1e-10):
        v,_=accumulate(36,eps)
        ladder.append(dict(guard=eps,raw_guarded_derivative=v[1],
                           subtracted_guarded_derivative=v[1]-predicted*np.log(rho/eps)))
    defect=abs(ladder[-1]['subtracted_guarded_derivative']-b[1])
    slope=(ladder[-1]['raw_guarded_derivative']-ladder[-2]['raw_guarded_derivative'])/np.log(100.)
    if abs(slope-predicted)>2e-7:
        raise ArithmeticError('The independently cut radial log slope did not match the archive jump.')
    if defect>2e-5:
        raise ArithmeticError(f'Exact odd extraction did not match guarded archive: {ladder,b}.')
    return dict(sign=sign,onshell_projection=b[0],fixed_spin_energy_derivative=b[1],
                analytic_odd_component=lead,radial_step_jump=jump,expected_step_log_coefficient=predicted,
                guard_ladder=ladder,measured_step_log_coefficient=float(slope),quadrature_defect=resolved,guard_replay_defect=float(defect))


def run():
    rows=[]
    rest=max(abs(photon_finite_remainder(m,m)*m*m/4-1) for m in (.05,.003))
    if rest>1e-14:
        raise ArithmeticError('The null-photon angular remainder lost its rest limit.')
    for grid in ('fine','third'):
        path=Path(__file__).parent/f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz'
        with gzip.open(path) as f:d=json.load(f)
        M=d['x_end'];m=d['mass_ratio']*M
        ms=d['checkpoints'][-1]['scalar_pole_over_T']
        Es,lo,hi,*_=intervals(M,ms,m,Q)
        central=(lo+hi)/2
        for name,ma,shift in (('parent',M,0.),('fermion',m,Es)):
            E=central-shift
            edges=np.hypot(d[f'endpoint_{name}_momentum_edges'],ma)
            occ=np.asarray(d[f'endpoint_{name}_occupations'])
            good=np.flatnonzero((edges[1:-1]+shift>lo+4e-4)&(edges[1:-1]+shift<hi-4e-4))
            step=edges[1+good[np.argmax(abs(np.diff(occ)[good]))]] if len(good) else None
            independent=[]
            for sign in (-1,1):
                for delta in (-1e-3,-1e-5,1e-5,1e-3):
                    F,C,R,off,ne,fa=charged(ma,E,delta,sign)
                    eta=min(abs(ne),abs(fa))*D('1e-6')/E
                    via_complex=off(1j*eta).imag/eta
                    direct,trace=independent_angular(ma,E,delta,sign)
                    independent.append(float(max(abs(direct/C-1),abs(via_complex/C-1),trace)))
            if max(independent)>2e-8:
                raise ArithmeticError('Charged angular derivative did not match its independent reconstruction.')
            rows.append(dict(grid=grid,species=name,energy=E,
                             maximum_independent_angular_and_trace_defect=max(independent),
                             central_charged=[run_radial(ma,E,edges,occ,s,False) for s in (-1,1)],
                             tuned_energy=float(step) if step is not None else None,
                             tuned_charged=run_radial(ma,step,edges,occ,1,True) if step is not None else None,
                             bare_null_photon=photon_checks(ma,E)))
    return dict(scope='Bare Feynman-gauge on-shell self-energy projection and fixed-spin energy derivative; no thermal sign, medium pole normalization, width or matched rate assigned.',photon_rest_limit_defect=rest,rows=rows)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
