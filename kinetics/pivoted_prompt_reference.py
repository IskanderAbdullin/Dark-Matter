"""Parent-local affine centering for the prompt collision and analytic JVP."""
import numpy as np
from joint_resonant_centered_stiff_reference import Centered


def pivoted(model,z,qref,R0):
    k=model.k
    if not len(k.y):
        return np.zeros(model.d),lambda v:np.zeros(model.d)
    f,g=model.arrays(qref+z)
    f=k.c[0].f(f,k.ph,k.I)
    g=k.c[1].f(g,k.ch,k.J)
    b=k.pre.astype(np.longdouble)*(1+f)*(1-g)
    active=np.r_[0,np.flatnonzero(np.diff(k.y))+1]
    counts=np.diff(np.r_[active,len(k.y)])
    def unfold(x):return np.repeat(x,counts)
    def accumulate(x):return unfold(np.add.reduceat(x,active))
    s,t=model.arrays(np.asarray(z,dtype=np.longdouble))
    uu=model.us.astype(np.longdouble); vv=model.uc.astype(np.longdouble)
    x=s[0,k.I]+s[1,k.I]*uu+t[0,k.J]+t[1,k.J]*vv
    pivot=unfold(np.minimum.reduceat(x,active))
    shift=x-pivot
    e=np.exp(-shift)
    md=np.expm1(-shift)
    Rstar=np.asarray(R0,dtype=np.longdouble)*np.exp(-pivot)
    near=unfold(np.maximum.reduceat(shift,active))<.75
    B=accumulate(b)
    A=accumulate(b*e)
    J=accumulate(b*md)
    D=B+Rstar*A
    omega=np.divide(k.fac*Rstar,D*model.H,out=np.zeros_like(D),where=D>0)
    v=np.where(near,J-B*md,A-B*e)
    result=model.proj(np.asarray(omega*b*v,dtype=float))
    def derivative(d):
        p,q=model.arrays(np.asarray(d,dtype=np.longdouble))
        h=p[0,k.I]+p[1,k.I]*uu
        j=q[0,k.J]+q[1,k.J]*vv
        db=b*(-f*h+g*j)
        de=-e*(h+j)
        dB=accumulate(db)
        dA=accumulate(db*e+b*de)
        dJ=accumulate(db*md+b*de)
        # Holding the pivot fixed yields the exact JVP, since the physical
        # rate is algebraically independent of this common parent pivot.
        dV=np.where(near,dJ-dB*md-B*de,dA-dB*e-B*de)
        dD=dB+Rstar*dA
        v1=omega*(v*db+b*dV-b*v*np.divide(dD,D,out=np.zeros_like(D),where=D>0))
        return model.proj(np.asarray(v1,dtype=float))
    return result,derivative


class Pivoted(Centered):
    def centered(self,z,qref,R0):
        return pivoted(self,z,qref,R0)


def test(eos,row):
    from joint_resonant_kinetic_linear_modes_reference import configuration,metric
    T,H,k,N,t,Q,sh=configuration(row,eos,8,6.5)
    q=np.concatenate([x.T.ravel() for x in sh])
    model=Pivoted(k,H)
    ref,R,_=model.reference(q,T)
    direction=np.cos(np.arange(len(q))*.37)
    result=[]
    for amp in (.07,1e-7):
        z=amp*np.sin(np.arange(len(q))*.31)
        a,ja=Centered.centered(model,z,ref,R)
        b,jb=model.centered(z,ref,R)
        flux=np.linalg.norm(a-b)/np.linalg.norm(a)
        jvp=np.linalg.norm(ja(direction)-jb(direction))/np.linalg.norm(ja(direction))
        e=2e-6
        fd=(model.centered(z+e*direction,ref,R)[0]-model.centered(z-e*direction,ref,R)[0])/(2*e)
        jfd=np.linalg.norm(fd-jb(direction))/np.linalg.norm(fd)
        assert max(flux,jvp,jfd)<5e-8
        _,_,inv=metric(k,sh)
        normalized=inv/np.linalg.norm(inv,axis=0)
        defect=float(np.linalg.norm(normalized.T@b)/np.linalg.norm(b))
        jdefect=float(np.linalg.norm(normalized.T@jb(direction))/np.linalg.norm(jb(direction)))
        assert max(defect,jdefect)<3e-13
        assert np.linalg.norm(model.centered(np.zeros_like(z),ref,R)[0])==0.
        result.append(dict(amplitude=amp,flux_relative_error=float(flux),jvp_relative_error=float(jvp),
                           finite_difference_defect=float(jfd),flux_invariant_defect=defect,
                           JVP_invariant_defect=jdefect))
    return result


if __name__=="__main__":
    import argparse,json
    from pathlib import Path
    from sm_dark_entropy_exchange_reference import StandardModel
    p=argparse.ArgumentParser();p.add_argument("eos");p.add_argument("spectrum")
    a=p.parse_args()
    row=next(x for x in json.loads(Path(a.spectrum).read_text())["rows"] if x["TR"]==.05 and abs(x["TF"]-.045)<1e-8)
    print(json.dumps(dict(rows=test(StandardModel(a.eos),row)),indent=2,allow_nan=False))
