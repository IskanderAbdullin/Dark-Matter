"""Finite, kinematically on-shell chi1 <-> chi2+phi resonance graph.

For energy lattice e1=.1+i*h, e2=.052+j*h and ephi=.048+k*h,
the reaction i=j+k is exactly energy preserving; every retained event
also obeys physical boosted two-body angle bounds. We prove exactly
three global left invariants, FD/BE fixed point and a bounded implicit
per-event extent. The positive lattice weights and rates are ILLUSTRATIVE,
not an IR-complete physical collision operator or a relic prediction.
"""
import json
import math
import numpy as np

M1, M2, MPH = .1, .02, 86e-9


def graph(n=8, h=.01):
    if n < 2 or h <= 0: raise ValueError('nontrivial positive lattice')
    e1=.1+h*np.arange(n+1)
    e2=.052+h*np.arange(n+1)
    ep=.048+h*np.arange(n+1)
    lam=(M1*M1-(M2+MPH)**2)*(M1*M1-(M2-MPH)**2)
    qstar=math.sqrt(lam)/(2*M1)
    ephstar=math.hypot(qstar,MPH)
    events=[]
    for i in range(n+1):
        p=math.sqrt(max(0.,e1[i]**2-M1*M1))
        for k in range(i+1):
            j=i-k
            if p == 0:
                if abs(ep[k]-ephstar)>1e-10: raise ArithmeticError('at-rest shell mismatch')
                c=0.
            else:
                a=e1[i]/M1*ephstar
                b=p/M1*qstar
                c=(a-ep[k])/b
                if abs(c)>1+1e-10: raise ArithmeticError('kinematically forbidden event')
            if abs(e1[i]-e2[j]-ep[k])>1e-12:
                raise ArithmeticError('nonresonant lattice')
            events.append((i,j,k,float(c)))
    return (e1,e2,ep),events


def incidence(n,events):
    S=np.zeros((3*(n+1),len(events)))
    for col,(i,j,k,c) in enumerate(events):
        S[i,col]=-1.
        S[n+1+j,col]=1.
        S[2*(n+1)+k,col]=1.
    return S


def invariant_audit(n=8):
    energies, events=graph(n)
    S=incidence(n,events)
    r=np.linalg.matrix_rank(S,tol=1e-10)
    size=3*(n+1)
    ones=np.ones(n+1);zeros=np.zeros(n+1)
    charges=np.vstack((np.r_[ones,ones,zeros],np.r_[ones,zeros,ones],np.concatenate(energies)))
    err=float(np.max(abs(charges@S)))
    return {'n':n,'n_modes':size,'n_events':len(events),'stoichiometric_rank':int(r),
            'left_nullity':int(size-r),'expected_physical_invariants':3,
            'invariant_residual':err,'cos_min':float(min(e[3] for e in events)),
            'cos_max':float(max(e[3] for e in events))}


def equilibrium(n=8, T=.7, mu2=-.01, mup=-.01):
    energies,events=graph(n)
    mu1=mu2+mup
    e1,e2,ep=energies
    f1=1/(1+np.exp((e1-mu1)/T))
    f2=1/(1+np.exp((e2-mu2)/T))
    bp=1/np.expm1((ep-mup)/T)
    F=[];R=[]
    for i,j,k,_ in events:
        F.append(f1[i]*(1-f2[j])*(1+bp[k]))
        R.append((1-f1[i])*f2[j]*bp[k])
    return {'max_forward_reverse_difference':float(np.max(abs(np.array(F)-np.array(R)))),
            'max_relative_difference':float(np.max(abs(np.array(F)-np.array(R))/np.maximum(F,1e-30)))}


def entropy(f1,f2,b):
    def ferm(f):
        return -f*math.log(f)-(1-f)*math.log1p(-f)
    def bos(z):
        return (1+z)*math.log1p(z)-z*math.log(z) if z>0 else 0.
    return ferm(f1)+ferm(f2)+bos(b)


def implicit_event(f1,f2,b,dt=1., rate=1.,capacities=(1.,1.,1.)):
    """Bracket at least one positive Pauli/Bose implicit reaction extent."""
    w1,w2,wp=capacities
    if not(0<f1<1 and 0<f2<1 and b>0 and dt>0 and rate>0 and min(capacities)>0):
        raise ValueError('interior Pauli/Bose state, positive capacities and rate')
    lo=max(w1*(f1-1),-w2*f2,-wp*b)
    hi=min(w1*f1,w2*(1-f2))
    def state(x):return f1-x/w1,f2+x/w2,b+x/wp
    def flux(x):
        a,d,z=state(x)
        return a*(1-d)*(1+z)-(1-a)*d*z
    def residual(x): return x-dt*rate*flux(x)
    if residual(lo)>1e-12 or residual(hi)<-1e-12:
        raise ArithmeticError('physical extent does not bracket a root')
    for _ in range(90):
        mid=(lo+hi)/2
        if residual(mid)<=0:lo=mid
        else:hi=mid
    extent=(lo+hi)/2
    new=state(extent)
    def sf(f): return -f*math.log(f)-(1-f)*math.log1p(-f)
    def sb(z): return (1+z)*math.log1p(z)-z*math.log(z) if z>0 else 0.
    initial=w1*sf(f1)+w2*sf(f2)+wp*sb(b)
    final=w1*sf(new[0])+w2*sf(new[1])+wp*sb(new[2])
    F=new[0]*(1-new[1])*(1+new[2])
    R=(1-new[0])*new[1]*new[2]
    H=(F-R)*math.log(F/R)
    return {'extent':extent,'old':[f1,f2,b],'new':list(new),
            'entropy_gain':final-initial,'implicit_residual':residual(extent),
            'entropy_lower_bound':dt*rate*H,
            'invariant_charge1_error':(-extent+extent),
            'invariant_charge2_error':(-extent+extent)}


def network_run(n=8, time=2., nsteps=4, rate=.02, reverse=False):
    """Synthetic-rate Lie composition on one connected exact-energy event lattice."""
    if nsteps<1 or time<=0 or rate<=0:raise ValueError('positive time and rate')
    (e1,e2,ep),events=graph(n)
    f1=1/(1+np.exp(e1/.9))
    f2=1/(1+np.exp(e2/.7))
    b=1/np.expm1(ep/.7)
    start=(f1.copy(),f2.copy(),b.copy())
    def inv(a,d,z):
        return np.array([sum(a)+sum(d),sum(a)+sum(z),
                         np.dot(a,e1)+np.dot(d,e2)+np.dot(z,ep)])
    prior=inv(*start)
    def total_entropy(a,d,z):
        return float(sum(-a*np.log(a)-(1-a)*np.log1p(-a))+
                     sum(-d*np.log(d)-(1-d)*np.log1p(-d))+
                     sum((1+z)*np.log1p(z)-z*np.log(z)))
    S0=total_entropy(*start)
    for _ in range(nsteps):
        for i,j,k,_ in reversed(events) if reverse else events:
            q=implicit_event(float(f1[i]),float(f2[j]),float(b[k]),time/nsteps,rate)
            f1[i],f2[j],b[k]=q['new']
    delta=inv(f1,f2,b)-prior
    return {'n':n,'time':time,'nsteps':nsteps,'reverse':reverse,
            'invariants_residual':delta.tolist(),
            'entropy_gain':total_entropy(f1,f2,b)-S0,
            'delta_Nphi':float(sum(b-start[2])),
            'bounds':{'f1_min':float(min(f1)),'f1_max':float(max(f1)),
                      'f2_min':float(min(f2)),'f2_max':float(max(f2)),
                      'phi_min':float(min(b))},
            'class_ready':False}


def run():
    return {'schema_version':1,
            'graph':[invariant_audit(n) for n in (2,4,8,16)],
            'thermal':[equilibrium(n) for n in (2,8,16)],
            'implicit_event':implicit_event(.5,.2,.7,dt=3),
            'network_time':[network_run(8,nsteps=i) for i in (1,2,4,8)],
            'limitations':'Kinematically allowed E1=E2+Ephi lattice, not full phase-space weighted continuous physics, no IR rest node, no finite SM/FRW/BEC/CLASS.',
            'class_ready':False}


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
