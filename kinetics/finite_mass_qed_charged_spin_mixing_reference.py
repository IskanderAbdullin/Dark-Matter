"""Frozen free-line charged cut and pole spinor projector pilot."""
import gzip,json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from scipy.integrate import quad
from early_pole_decay_mode_clock_reference import Q,intervals
from early_pole_hard_yukawa_photoproduction_reference import G,I,MET,slash,dot
from finite_mass_qed_dispersive_cut_bridge_reference import cuts

D=np.longdouble
BETA=G[0]
ALPHA=np.array([G[0]@x for x in G[1:]])


def angular(m,p,k,sg):
    m,p,k=map(D,(m,p,k))
    E,t=np.hypot(p,m),np.hypot(k,m)
    R=p*k
    d=(k-p)*(k+p)/(E+t)
    if sg>0:
        big=E*t-m*m+R
        small=m*m*d*d/big
        near,far=-2*small,-2*big
    else:
        big=E*t+m*m+R
        small=m*m*(E+t)**2/big
        near,far=2*big,2*small
    A=(near+far)/2
    B=2*R
    z=B/A
    if abs(z)<.005:
        J=2/A*(1+z*z/3+z**4/5+z**6/7)
        L=-2/A*(z/3+z**3/5+z**5/7+z**7/9)
    else:
        J=np.log(near/far)/B
        L=(2-A*J)/B
    return np.array([-2*sg*t*J,2*k*L,4*m*J],dtype=D)


def direct_angular(m,p,k,sg):
    E=np.hypot(p,m)
    t=sg*np.hypot(m,k)
    P=np.array([E,0.,0.,p])
    def value(u,j):
        L=np.array([t,k*np.sqrt((1-u)*(1+u)),0.,k*u])
        K=P-L
        U=slash(L)+m*I
        N=BETA@sum(z*R@U@R for z,R in zip(MET,G,strict=True))
        operator=(I,ALPHA[2],BETA)[j]
        return (np.trace(N@operator).real/4)/dot(K,K)
    return np.array([quad(lambda u:value(u,j),-1,1,epsabs=3e-11,epsrel=3e-10,limit=120)[0]
                     for j in range(3)])

def integrate(m,p,e,f,sg,n):
    x,w=roots_legendre(n)
    y=np.pi*(x+1)/2
    vals=[]
    for a,b,fi in zip(e[:-1],e[1:],f,strict=True):
        if not fi:
            continue
        bins=sorted({a,b,*[z for z in (p,p-1e-3,p+1e-3,p-.01,p+.01) if a<z<b]})
        for l,r in zip(bins[:-1],bins[1:],strict=True):
            k=l+(r-l)*np.sin(y/2)**2
            W=fi*(r-l)*np.pi/4*w*np.sin(y)*k*k/(2*np.hypot(k,m))
            v=np.array([angular(m,p,q,sg) for q in k])
            vals.append(W@v)
    return np.sum(vals,axis=0).astype(float)


def core(pvec,m,e,occ):
    p=np.linalg.norm(pvec)
    E=np.hypot(p,m)
    A=np.einsum('i,ijk->jk',pvec/p,ALPHA)
    H=p*A+m*BETA
    plus=(I+H/E)/2
    comps=np.array([integrate(m,p,e,occ,s,64) for s in (-1,1)])
    alt=np.array([integrate(m,p,e,occ,s,43) for s in (-1,1)])
    err=float(max(np.ravel(abs(comps-alt))))
    if err>3e-8:
        raise ArithmeticError('Matrix component independent angular mesh failure.')
    check=[]
    ret=[]
    for z in comps:
        d=z[0]*I+z[1]*A+z[2]*BETA
        pole=(np.trace(plus@d)/2).real
        dP=((I-plus)@d@plus+plus@d@(I-plus))/(2*E)
        # Compare a derivative of independently diagonalized Hermitian
        # Hamiltonians at the frozen on-shell matrix, no scalar trace fitting.
        def diag(h):
            evals,U=np.linalg.eigh(H+h*d)
            V=U[:,2:]
            return V@V.conj().T,np.mean(evals[2:])
        d1,y=diag(1e-5)
        d2,z1=diag(-1e-5)
        num=(d1-d2)/2e-5
        eerr=abs((y-z1)/2e-5-pole)
        check.append(float(max(np.max(abs(num-dP)),eerr)))
        ret.append((pole,dP,d))
    c=cuts(p,m,e,occ,N=56)
    if max(abs(c[i]-ret[i][0]) for i in (0,1))>2e-8 or max(check)>2e-8:
        raise ArithmeticError('Matrix pole projection or independent eigensystem failed.')
    spin_trace=[]
    for k in (.12,.54,3.2):
        for s in (-1,1):
            x=direct_angular(m,p,k,s)
            y=angular(m,p,k,s)
            spin_trace.append(float(np.max(abs(x-y)/np.maximum(1.,abs(x)))))
    if max(spin_trace)>3e-8:
        raise ArithmeticError('Direct angular gamma projection did not match the signed matrix cut.')
    return plus,ret,dict(components_negative_positive=comps.tolist(),pole_projection=c.tolist(),
                         refinement_error=err,independent_diagonalization_error=max(check),direct_matrix_angular_defect=max(spin_trace))


def case(d):
    M=d['x_end'];m=M*d['mass_ratio'];ms=d['checkpoints'][-1]['scalar_pole_over_T']
    Es,lo,hi,*_=intervals(M,ms,m,Q)
    E=(lo+hi)/2
    z=(E*Es-(M*M+ms*ms-m*m)/2)/Q
    P=np.array([np.sqrt((E-M)*(E+M)-z*z),0,z])
    B=P-np.array([0.,0.,Q])
    a,A,ca=core(P,M,d['endpoint_parent_momentum_edges'],d['endpoint_parent_occupations'])
    b,Bs,cb=core(B,m,d['endpoint_fermion_momentum_edges'],d['endpoint_fermion_occupations'])
    def trace(X,Y):
        return np.trace(X@BETA@Y@BETA).real
    base=trace(b,a)
    direct=np.array([[trace(b,u[1])/base,trace(v[1],a)/base] for u,v in zip(A,Bs,strict=True)])
    ea,eb=np.hypot(np.linalg.norm(P),M),np.hypot(np.linalg.norm(B),m)
    ha, hb=(2*a-I)*ea,(2*b-I)*eb
    dh1=A[0][2]+A[1][2]
    dh2=Bs[0][2]+Bs[1][2]
    def split(h):
        def pin(H):
            _,U=np.linalg.eigh(H)
            return U[:,2:]@U[:,2:].conj().T
        return trace(pin(hb+h*dh2),pin(ha+h*dh1))/base
    rate=float(np.sum(direct))
    half=(split(2e-6)-split(-2e-6))/4e-6
    eps=.005**2/(4*np.pi**2)
    actual=split(eps)-1
    if abs(rate-half)>5e-7 or abs(actual-eps*rate)>7e-12:
        raise ArithmeticError(f'Held trace: {rate}, {half}, {actual/eps}.')
    return dict(parent=ca,daughter=cb,Born_normalized_spin_trace=base,
                charged_spin_trace_derivatives_by_branch_and_leg=direct.tolist(),
                sum_before_gauge=float(np.sum(direct)),
                fraction_at_g005=rate*eps,exact_frozen_matrix_at_g005=actual,
                combined_centered_trace_derivative_defect=float(abs(rate-half)))



def run():
    outputs=[]
    for grid in ('fine','third'):
        with gzip.open(Path(__file__).parent/f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz') as f:
            d=json.load(f)
        outputs.append(dict(grid=grid,**case(d)))
    return dict(scope='Frozen charge-symmetric bare loop Hamiltonian spin-projector correction at held Born momenta. No derivative residue, phase shift or complete thermal rate.',rows=outputs)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
