"""Charge-even hard shell projection of the real on-shell log-discrepancy term."""
import json

import mpmath as mp
import numpy as np
from scipy.special import expit, roots_legendre


def gauss(n):
    t,w=roots_legendre(n)
    return (t+1)/2,w/2


def integrated_kernel(edges,x,n=24):
    """Symmetric K and dK/dln(x); K includes x² but not gauge prefactor."""
    nc=len(edges)-1
    t,w=gauss(n)
    A,B=t[:,None]**2,t[None,:]
    W=(2*t*w)[:,None]*w[None,:]
    K=np.empty((nc,nc))
    D=np.empty_like(K)
    v=np.diff(edges)
    def integral(p,r,gap,measure):
        ep,er=np.hypot(p,x),np.hypot(r,x)
        base=(x*x*measure)*(p/ep)*(r/er)*(np.log(p+r)-np.log(gap))
        return np.sum(base),np.sum(base*(2-(x/ep)**2-(x/er)**2))
    for i in range(nc):
        l,u=edges[i:i+2]
        wi=v[i]
        p=l+wi*A
        length=wi*A
        for j in range(i+1):
            wj=v[j]
            if i==j:
                gap=length*B**4
                vals=np.asarray(integral(p,p-gap,gap,2*wi*W*length*4*B**3))
            elif i==j+1:
                # Rectangle sharing the boundary l: scale distances to [0,1],
                # split along normalized equal distances, then set small=A*B.
                left=A*B
                vals=np.asarray(integral(l+wi*A,l-wj*left,A*(wi+wj*B),W*A*wi*wj))
                vals+=integral(l+wi*left,l-wj*A,A*(wi*B+wj),W*A*wi*wj)
            else:
                r=edges[j]+wj*B
                vals=integral(p,r,p-r,W*wi*wj)
            K[i,j]=K[j,i]=vals[0]
            D[i,j]=D[j,i]=vals[1]
    return K,D


def matrices(edges,a,masses,g,n=24,logs=True):
    """H maps cell occupations to derivative of a⁴ rho; T maps to its scale derivative."""
    nc=len(edges)-1
    t,w=gauss(n)
    q=edges[:-1,None]+np.diff(edges)[:,None]*t
    weights=np.diff(edges)[:,None]*w
    x=masses*a
    XX=x[:,None,None]
    e=np.hypot(q,XX)
    j=(weights*q*q/e).sum(axis=-1)
    dj=-(weights*q*q/e*(XX/e)**2).sum(axis=-1)
    H=np.zeros((len(x)*nc,len(x)*nc))
    T=np.zeros_like(H)
    gpu=g*g/np.pi**4
    for i in range(1,len(x)):
        f=slice(i*nc,(i+1)*nc)
        p=slice(0,nc)
        H[f,f]=gpu*np.outer(j[i],j[i])
        T[f,f]=gpu*(np.outer(j[i],dj[i])+np.outer(dj[i],j[i]))
        H[p,f]=gpu*np.outer(j[0],j[i]); H[f,p]=H[p,f].T
        T[p,f]=gpu*(np.outer(j[0],dj[i])+np.outer(dj[0],j[i])); T[f,p]=T[p,f].T
        if logs and x[i]>0:
            K,D=integrated_kernel(edges,x[i],n)
            H[f,f]-=gpu*K/2
            T[f,f]-=gpu*D/2
    return H,T


def shape(z):
    return np.array([1/np.expm1(z/.95+.4),expit(-z/.9-.3),expit(-z/.45-1.)])


def continuum_log(a,m,g,n=152):
    t,w=gauss(n)
    p=.25+7.75*t[:,None]
    L=p-.25
    gap=L*t**4
    r=p-gap
    W=(7.75*w[:,None])*(L*4*t**3*w)
    # Two identical halves; log kernel symmetric in p,r.
    pole=0.
    for i in (1,2):
        x=a*m[i]
        pole-=g*g*x*x/(2*np.pi**4)*np.sum(
                  W*(p/np.hypot(p,x))*shape(p)[i]*(r/np.hypot(r,x))*shape(r)[i]*
                  (np.log(p+r)-np.log(gap)))
    return pole


def angular_check():
    # Contract the two signed-energy spin traces before the angular average.
    # Check the integral independently at high precision, including unequal
    # but near-coincident momenta; equal momenta have an integrable radial log.
    out=[]
    with mp.workdps(65):
        for masses in ('.003','.05','.5'):
            for p,k in (('.25','7.9'),('.6','.7'),('.25','.250002'),('1.','1.00000001')):
                p,k,m=map(mp.mpf,(p,k,masses))
                e1,e2=mp.sqrt(p*p+m*m),mp.sqrt(k*k+m*m)
                m2=m*m
                def angular(j):
                    scal=e1*e2-p*k*(1-2*j)
                    return sum((-8*z+16*m2)/(2*(m2-z)) for z in (scal,-scal))/8
                scale=((e1*e2-p*k)-m2)/(2*p*k)
                # Resolve the endpoint layer explicitly in the direct quadrature.
                breaks=[mp.mpf(0)]+sorted({z for z in (scale,mp.sqrt(scale),mp.mpf('.1')) if 0<z<1})+[mp.mpf(1)]
                ref=mp.quad(angular,breaks)
                got=1-m2/(2*p*k)*mp.log((p+k)/abs(p-k))
                error=abs(ref-got)
                if error>mp.mpf('2e-45'):
                    raise ValueError(f'Angular branch sum mismatch {error}.')
                out.append(dict(m=float(m),p=float(p),k=float(k),absolute_error=float(error)))
    return out


def run():
    a,g=1.,.005
    angular=angular_check()
    output=[]
    for masses in (np.array([0.,.05,.003]),np.array([0.,.5,.25])):
        exact=continuum_log(a,masses,g)
        rows=[]
        for nc in (8,16,32,64):
            edges=np.linspace(.25,8.,nc+1)
            H,T=matrices(edges,a,masses,g)
            H0,T0=matrices(edges,a,masses,g,logs=False)
            vol=(edges[1:]**3-edges[:-1]**3)/3
            u,w=gauss(24)
            q=edges[:-1,None]+np.diff(edges)[:,None]*u
            f=((shape(q)*q*q*w).sum(axis=2)*np.diff(edges)/vol).ravel()
            N=(np.array([2.,4.,4.])[:,None]*vol/(2*np.pi**2)).ravel()
            energy=lambda F,mat: (F@mat@F)/2
            logged=energy(f,H-H0)
            trace=energy(f,T)
            explicit_log_trace=energy(f,T-T0)
            direction=.35*np.sin(np.arange(f.size)+.4)*f
            chi=1e-3
            lhs=(energy(f+chi*direction,H)-energy(f-chi*direction,H))/(2*chi)
            rhs=(H@f)@direction
            n0=f*N
            delta=direction*chi*N
            energychange=energy((n0+delta)/N,H)-energy(n0/N,H)
            midpoint=((H@((n0+delta/2)/N))/N)@delta
            if max(abs(H-H.T).ravel())>0 or abs(lhs-rhs)>1e-17 or abs(energychange-midpoint)>1e-19:
                raise ValueError('Shell energy derivative/reciprocity check failed.')
            row=dict(shells=nc,interaction_energy=float(energy(f,H)),
                     log_interaction_energy=float(logged),
                     relative_log_error=float(abs(logged/exact-1)),
                     fixed_comoving_trace=float(trace),
                     fixed_comoving_log_trace=float(explicit_log_trace),
                     energy_derivative_abs_error=float(abs(lhs-rhs)),
                     midpoint_identity_abs_error=float(abs(energychange-midpoint)),
                     max_relative_fermion_shift=float(max(abs(((H-H0)@f)[nc:]/(H0@f)[nc:]))))
            if nc==32:
                Hfine,Tfine=matrices(edges,a,masses,g,n=40)
                coeff_scale=max(np.max(abs(H-H0)),np.max(abs(T-T0)))
                qerr=max(np.max(abs(H-Hfine)),np.max(abs(T-Tfine)))/coeff_scale
                if qerr>5e-8:
                    raise ValueError('Same-shell diagonal quadrature precision not certified.')
                row['scaled_matrix_quadrature_defect']=float(qerr)
                hh=2e-4
                minus=energy(f,matrices(edges,a*np.exp(-hh),masses,g)[0])
                plus=energy(f,matrices(edges,a*np.exp(hh),masses,g)[0])
                err=abs((plus-minus)/(2*hh)-trace)
                if err>1e-13:
                    raise ValueError('Fixed-comoving trace derivative failed.')
                row['trace_finite_difference_abs_error']=float(err)
            rows.append(row)
        output.append(dict(masses=masses.tolist(),continuous_log_reference=float(exact),rows=rows))
        if rows[-1]['relative_log_error']>=.005 or not all(3 < r1['relative_log_error']/r2['relative_log_error'] < 4.2
                for r1,r2 in zip(rows,rows[1:])):
            raise ValueError('Expected shell-shape refinement approaching second order.')
    return dict(angular_kernel_checks=angular,scope='Isotropic charge-even hard shapes; reciprocal on-shell exchange candidate, not dynamic screening or collision matching.',
                cases=output)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
