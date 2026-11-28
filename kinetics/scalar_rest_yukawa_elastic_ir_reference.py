"""Below-pole chi2-phi elastic damping of a scalar rest-mode perturbation."""
import argparse
import json
import numpy as np
from scipy.special import expit, roots_legendre
from three_species_decay_onsager_reference import M,m,mp,y,susceptibility
from sm_dark_entropy_exchange_reference import StandardModel,MPL

def dot(a,b):
    return a[...,0]*b[...,0]-np.sum(a[...,1:]*b[...,1:],axis=-1)

def cross(T,E,xc,*,exact=True):
    p=np.sqrt((E-m)*(E+m))
    s=m*m+mp*mp+2*mp*E
    u=np.sqrt(s)
    o=mp*p/u
    z=(mp*mp+mp*E)/u
    ec=u-z
    if exact:
        sn=np.sqrt(np.maximum(0.,1-xc*xc))
        a=np.array([ec,0.,0.,o])
        b=np.array([np.full_like(xc,ec),-o*sn,np.zeros_like(xc),-o*xc]).T
        kp=np.array([np.full_like(xc,z),o*sn,np.zeros_like(xc),o*xc]).T
        sc=np.array([u,0.,0.,0.])
        uuvec=a-kp
        d1=s-M*M
        d2=dot(uuvec,uuvec)-M*M
        F=M*(1/d1+1/d2)
        Q=sc/d1+uuvec/d2[:,None]
        b_a=dot(b,a)
        Tr=2*F*F*(b_a+m*m)+4*F*m*dot(a+b,Q)+2*(2*dot(b,Q)*dot(a,Q)+(m*m-b_a)*dot(Q,Q))
        return y**4*Tr/(16*np.pi*s)
    return np.full_like(xc,y**4*m*m/(np.pi*(M-m)**2*(m+mp)**2))

def trace_control():
    I=np.eye(2)
    O=np.zeros((2,2))
    S=(np.array([[0.,1.],[1.,0.]]),np.array([[0.,-1j],[1j,0.]]),np.diag([1.,-1.]))
    g=[np.block([[I,O],[O,-I]])]+[np.block([[O,t],[-t,O]]) for t in S]
    unit=np.eye(4)
    def slash(q):
        return g[0]*q[0]-sum(q[k]*g[k] for k in range(1,4))
    errors=[]
    for E in (.021,.1,.7,10.,30.):
        s=m*m+mp*mp+2*mp*E
        u=np.sqrt(s)
        p=mp*np.sqrt(E*E-m*m)/u
        z=mp*(mp+E)/u
        c=(s+m*m-mp*mp)/(2*u)
        p1=np.array([c,0,0,p])
        total=np.array([u,0,0,0])
        for v in (-.97,-.5,0.,.7,.99):
            q=np.array([z,p*np.sqrt(1-v*v),0.,p*v])
            p2=total-q
            b=p1-q
            A=(slash(total)+M*unit)/(s-M*M)+(slash(b)+M*unit)/(dot(b,b)-M*M)
            t=np.trace((slash(p2)+m*unit)@A@(slash(p1)+m*unit)@g[0]@A.conj().T@g[0]).real/2
            x=y**4*t/(16*np.pi*s)
            y_ref=cross(1.,E,np.array([v]))[0]
            errors.append(abs(x/y_ref-1))
    assert max(errors)<1e-12
    return max(errors)

def calc(T,hard,*,n=52,cut=36.,exact=True):
    v,w=roots_legendre(n)
    k=(np.hypot(hard,mp)-mp)/(2*mp)
    e_threshold=k*mp+np.sqrt((k*mp)**2+m*m+k*(m*m+mp*mp))
    # Put the opening of the hard-to-reach outgoing gate exactly on
    # a phase-space panel boundary to independently converge its tail.
    bounds=np.unique(np.clip([0.,np.sqrt(max(0.,(e_threshold-m)/T)),np.sqrt(cut)],0.,np.sqrt(cut)))
    xx=(bounds[:-1,None]+(bounds[1:]-bounds[:-1])[:,None]*(v+1)/2).ravel()
    xd=((bounds[1:]-bounds[:-1])[:,None]*w).ravel()
    u=m+T*xx*xx
    dp=T*xx*xd
    total=survivor=vac=0.
    for E,dE in zip(u,dp):
        p=np.sqrt((E-m)*(E+m))
        s=m*m+mp*mp+2*mp*E
        assert (s-m*m)/(M*M-m*m)<.006  # neither elastic channel approaches the pole
        d=mp*p*p/s
        full=mp+2*d
        cut_lab=np.hypot(hard,mp)
        def angle(low,bose):
            if low>=full:return 0.
            k=np.log(full/low)
            ef=low*np.exp(k*(v+1)/2)
            out=(ef-mp)/d-1
            dw=.25*k*ef*w/d
            sf=cross(T,E,np.clip(out,-1,1),exact=exact)
            a=sf*dw
            if bose:
                a*=expit((E+mp-ef)/T)/(-np.expm1(-ef/T))
            return a.sum()
        density=4*dE*p*p/(2*np.pi**2)*expit(-E/T)
        total+=density*angle(mp,True)
        survivor+=density*angle(cut_lab,True)
        vac+=density*angle(mp,False)
    eq=-np.expm1(-mp/T)
    return np.array([total*eq,survivor*eq,vac])

def main(file):
    from decay_scalar_ir_relaxation_reference import main as ir
    source=ir(file)
    eos=StandardModel(file)
    out=[]
    for row in source["rows"]:
        T=row["T_GeV"]
        p=row["first_crossing_p_GeV"]
        x=calc(T,p)
        z=calc(T,p,n=76,cut=42)
        u=calc(T,p,n=76,cut=42,exact=False)
        err=np.max(abs(x[[0,2]]/z[[0,2]]-1))
        hard_err=abs(x[1]/z[1]-1) if z[1]>0 and x[1]>0 else None
        assert err<1e-4 and (hard_err is None or hard_err<1e-3),(err,hard_err)
        H=np.sqrt((eos.state(T)[0]+sum(susceptibility(T,mm,g,be)[1] for mm,g,be in ((M,4,False),(m,4,False),(mp,1,True))))/3)/MPL
        rows=dict(T_GeV=T,local_rest_damping_over_H=z[0]/H,tail_handoff_damping_over_H=z[1]/H,
                      vacuum_loss_clock_over_H=z[2]/H,thermal_angle_and_tail_error=err,handoff_relative_refinement=hard_err,
                      low_invariant_contact_control=abs(u[0]/z[0]-1))
        out.append(rows)
    return dict(rows=out,independent_gamma_trace_max_relative_error=trace_control(),scope="Below-pole elastic phi-chi2 only, local diagonal scalar damping with other scalar modes held at equilibrium; finite incoming-energy cutoff and no full transport/eigenmode, RIS or BEC inference.")

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("eos")
    print(json.dumps(main(p.parse_args().eos),indent=2,allow_nan=False))
