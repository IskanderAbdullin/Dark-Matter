"""Construct a positive two-temperature template matching mean q and q²."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.special import expit,zeta
from archived_collisionless_streaming_reference import Background
from three_species_decay_onsager_reference import mp


def moments(q,f):
    w=q*q*f/np.trapezoid(q*q*f,q)
    return [float(np.trapezoid(w*q**k,q)) for k in (1,2,4)]


def fluid(q,f,mob):
    w=q*q*f/np.trapezoid(q*q*f,q)
    e=np.hypot(mob,q)
    r=float(np.trapezoid(w*e,q))
    p=float(np.trapezoid(w*q*q/(3*e),q))
    z=float(np.trapezoid(w*q**4/(3*e**3),q))
    return dict(scaled_rho=r,pressure_to_rho=p/r,adiabatic_sound_speed_squared=(5*p-z)/(3*(r+p)))


def run(eos,folder):
    bg=Background(eos)
    q,f=np.loadtxt(folder/"direct-sig07.dat").T
    p1,p2,p4=moments(q,f)
    c1=7*np.pi**4/(180*zeta(3))
    c2=15*zeta(5)/zeta(3)
    mu=p1/c1
    sd=np.sqrt(p2/c2-mu**2)
    assert 0.<sd<mu
    T1,T2=mu-sd,mu+sd
    Tq=np.sqrt(p2/c2)
    x=np.geomspace(1e-7,45,5000)
    fdl=lambda T:expit(-x/T)/(1.5*zeta(3)*T**3)
    fx=(fdl(T1)+fdl(T2))/2
    single=fdl(Tq)
    a,b,c=moments(x,fx)
    assert max(abs(a/p1-1),abs(b/p2-1))<2e-8
    assert abs(moments(x,single)[1]/p2-1)<2e-8
    file=folder/"fd-p12-mix.dat"
    np.savetxt(file,np.c_[x,fx],fmt="%.16e")
    ini=(folder/"std-amp3.ini").read_text()
    ini=ini.replace("direct-sig07.dat","fd-p12-mix.dat").replace("std-amp3_","std-p12_")
    (folder/"std-p12.ini").write_text(ini)
    out=[]
    fd_err=[]
    for z in (1e-7,2e-6,1e-5):
        x0=z*mp/bg.T0
        dt=1e-5
        for Q,F in ((q,f),(x,single),(x,fx)):
            N=Q*Q*F/np.trapezoid(Q*Q*F,Q)
            def direct(h):
                E=np.hypot(Q,x0*np.exp(h))
                return np.exp(-4*h)*np.asarray([np.trapezoid(N*E,Q),np.trapezoid(N*Q*Q/(3*E),Q)])
            dv=direct(dt)-direct(-dt)
            fd_err.append(abs((dv[1]/dv[0])/fluid(Q,F,x0)["adiabatic_sound_speed_squared"]-1))
    assert max(fd_err)<3e-8
    for s in (1e-8,1e-7,1e-6,2e-6,5e-6,1e-5,bg.aeq,1/21):
        z=s*mp/bg.T0
        v=[fluid(Q,F,z) for Q,F in ((q,f),(x,single),(x,fx))]
        out.append(dict(scale=s,controls=v))
    return dict(scope="Thermal-shaped controls are synthetic scalar PSDs at the same amplified abundance, not FD scalar "
                      "microphysics. Number, mean q and mean q² are matched by the positive two-temperature mixture.",
                source_moments_q1_q2_q4=[p1,p2,p4],reference_single_moments=moments(x,single),
                two_temperature_controls_in_photon_units=[T1,T2],
                matched_mixture_moments=[a,b,c],psd_sha256=hashlib.sha256(file.read_bytes()).hexdigest(),
                max_sound_speed_defect_via_log_a_finite_difference=float(max(fd_err)),column_order=["direct","q2_matched_single_FD_shape","q1_q2_matched_double_FD_shape"],rows=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path);p.add_argument("native",type=Path)
    a=p.parse_args();print(json.dumps(run(a.eos,a.native),indent=2,allow_nan=False))
