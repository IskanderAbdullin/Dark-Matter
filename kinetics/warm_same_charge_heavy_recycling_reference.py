"""Spin-antisymmetrized same-charge heavy production on passive warm archives."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from archived_heavy_recycling_gate_reference import angular_integral
from archived_scalar_pair_reaction_reference import cropped, points, snapshot_functions
from late_radiation_background_reference import TabulatedRadiation, MPL
from warm_ballistic_heavy_exposure_reference import redshift
from three_species_decay_onsager_reference import m, M, mp, y


def interference_integral(s):
    """Integral of Re(A_t A_u*) (spin average, without y**4)."""
    t0=s/2-M*M-m*m
    center=mp*mp+t0
    span=2*np.sqrt(np.maximum(0.,(s/4-m*m)*(s/4-M*M)))
    zmin=mp*mp+(M*M-m*m)**2/(t0+span)
    dint=np.divide(np.log1p(2*span/zmin),span*center,
                   out=np.asarray(2/center**2).copy(),where=span>0.)
    d0=s/2-m*m
    d1=s/2-M*M
    C0=s*s/8-d0*d1+M*m*s+m*m*d1+M*M*d0+M*M*m*m
    return (C0+center**2/2)*dint-1.


def exact_angular(s):
    ic=interference_integral(np.maximum(s,4*M*M))
    I=angular_integral(s)
    same_both=I-ic
    result=np.where(s>=4*M*M,same_both,0.)
    assert np.all(result>=0.) and np.all(result<=2*I)
    return result


def independent_check():
    I2=np.eye(2); Z=np.zeros((2,2)); one=np.eye(4)
    pauli=(np.array([[0,1],[1,0]]), np.array([[0,-1j],[1j,0]]), np.diag([1,-1]))
    gm=[np.block([[I2,Z],[Z,-I2]])]+[np.block([[Z,P],[-P,Z]]) for P in pauli]
    A=lambda E,v,mass:gm[0]*E-sum(g*x for g,x in zip(gm[1:],v))+mass*one
    x,w=roots_legendre(180)
    err=[]
    for s in np.geomspace(4*M*M*(1+1e-7),.9,28):
        E=np.sqrt(s)/2; p=np.sqrt(s/4-m*m); k=np.sqrt(s/4-M*M)
        B0=s/2-M*M-m*m+mp*mp; span=2*p*k
        zlo,zhi=np.log([B0-span,B0])
        z=np.exp((zhi+zlo)/2+(zhi-zlo)/2*x)
        angle=(B0-z)/span
        a,b=[A(E,[0.,0.,sign*p],m) for sign in (1.,-1.)]
        trace=[]
        for c in angle:
            v=k*np.array([np.sqrt(np.maximum(0.,1-c*c)),0.,c])
            d,e=A(E,v,M),A(E,-v,M)
            trace.append(np.trace(d@a@e@b).real/4)
        C=(zhi-zlo)/span*np.sum(w*np.array(trace)/(2*B0-z))
        I=angular_integral(s)
        err.append(abs((I-C)/exact_angular(s)-1.))
    assert max(err)<1e-9
    return max(err)


def folded(edges,occ,nx,nz):
    p,e,w=points(edges,m,occ,nx)
    E,F=e[:,None],e[None,:]
    P=p[:,None]*p[None,:]
    cap=np.clip((E*F+m*m-2*M*M)/P,-1.,1.)
    x,z=roots_legendre(nz)
    mu=((cap-1)/2)[:,:,None]+((cap+1)/2)[:,:,None]*x
    s=2*m*m+2*(E[:,:,None]*F[:,:,None]-P[:,:,None]*mu)
    k=np.sqrt(np.maximum(0.,s/4-M*M))
    rel=y**4*2*k/np.sqrt(s)/(64*np.pi*(E*F)[:,:,None])
    # For equal occupations per charge, the two identical initial-channel
    # factors of 1/2 sum to one. exact_angular includes 1/2 for identical
    # final states: 0.5 integral(A_t^2 + A_u^2 - 2 Re A_t A_u*) .
    weight=(w[:,None]*w[None,:]*(cap+1)/4)[:,:,None]*z
    rel*=weight
    return float(4*np.sum(rel*angular_integral(s))),float(4*np.sum(rel*exact_angular(s)))



def read(p):
    return json.loads(Path(p).read_bytes())


def integrate(eos, paths, positive_exposure, crop_to=16.):
    defect = independent_check()
    s = [read(p)["snapshots"][-1] for p in paths]
    s.append(cropped(s[-1], crop_to))
    ref = read(positive_exposure)["cases"]
    answers = []
    for state, row in zip(s, ref, strict=True):
        p = np.asarray(state["lab_momentum_edges_GeV"])
        _, _, occ, _ = snapshot_functions(state)
        R = TabulatedRadiation(eos, state["T_visible_GeV"])
        den = row["normalization_at_anchor_GeV3"]
        assert R.state(1.)[0] == row["anchor_T_GeV"]
        e = np.log(row["b_panels"])

        def sampled(b, order):
            H = np.sqrt(R.state(b)[1]/3)/MPL
            return 2*b**3/(den*H)*np.asarray(folded(p/b, redshift(occ, m, b), *order))

        controls = []
        for order, nt in (((7, 21), 8), ((9, 27), 10)):
            x, w = roots_legendre(nt)
            a = (e[:-1, None]+e[1:, None])/2
            d = np.diff(e)[:, None]/2
            vals = (d*w).ravel() @ np.asarray([
                sampled(b, order) for b in np.exp(a+d*x).ravel()])
            against = next(k for k in row["controls"] if k["log_b_order"] == nt)
            phi, opp = against["integrated_phi_and_opposite_stable"]
            assert abs(vals[0]/opp-1.) < 2e-13
            total = float(phi+sum(vals))
            controls.append(dict(order=order, time=nt, integrated_phi_opp_same=[phi, *vals],
                                 total=total, same_vs_opposite=vals[1]/vals[0],
                                 exact_over_optimistic=total/against["optimistic_including_same_charge"]))
        assert abs(controls[-2]["total"]/controls[-1]["total"]-1.) < .002
        answers.append(dict(snapshot=row["snapshot"], controls=controls))
    return dict(scope="Pauli-free unpolarized forward folding on the unchanged passive "
                      "redshifting shapes; final heavy Pauli blocking and all later redistribution are omitted.",
                max_spin_trace_relative_error=defect,
                sha256=[hashlib.sha256(Path(p).read_bytes()).hexdigest()
                        for p in [eos, *paths, positive_exposure]],
                cases=answers)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos")
    p.add_argument("positive_exposure")
    p.add_argument("paths", nargs=2)
    a = p.parse_args()
    print(json.dumps(integrate(a.eos, a.paths, a.positive_exposure),
                     indent=2, allow_nan=False))
