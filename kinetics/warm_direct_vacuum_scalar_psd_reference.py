"""Construct positive scalar comoving shell counts from the two vacuum direct sources."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from archived_collisionless_streaming_reference import Background
from empty_lepton_pair_spectrum_reference import Source
from kloe_pion_band_portal_reference import POINT_MUON,table
from late_radiation_background_reference import TabulatedRadiation
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import Es,M,P,mp
from timelike_portal_onsager_coverage_reference import V,width,mpi,factor
from warm_light_occupation_stage_audit_reference import read


def pion_parents(T,data,n,nq,na):
    grid,cs,_=data
    a,w=roots_legendre(n)
    s=(grid[:,None]+.005*a).ravel()
    R=np.sqrt(s)
    wd=(np.repeat(cs,n)*np.tile(.005*w,len(grid)))
    v=(R-2*M)*(R+2*M)/s
    vp=(R-2*mpi)*(R+2*mpi)/s
    z,d=roots_legendre(nq)
    t=(z+1)*np.sqrt(75.)/2
    Q=R[:,None]+T*t*t
    wg=T*np.sqrt(75.)*t*d
    p=np.sqrt((Q-R[:,None])*(Q+R[:,None]))
    Z=(2*factor*(s+2*M*M)*s*np.sqrt(v)*wd*s/POINT_MUON/((s-V*V)**2+(V*width)**2))
    W=Z[:,None]*np.exp(-Q/T)*p*wg
    x,w=roots_legendre(na)
    h=p*np.sqrt(v)[:,None]/2
    Ep=Q[:,:,None]/2+h[:,:,None]*x
    u=p[:,:,None]*np.sqrt(vp)[:,None,None]*x/(2*T)
    z=Q[:,:,None]/(2*T)
    B=1./((-np.expm1(-z+u))*(-np.expm1(-z-u)))
    P2=(3*x*x-1)/2
    A=B@(w/2)
    B2=(-2*v/(3-v))[:,None]*(B@(w*P2/2))
    N=(A[:,:,None]+B2[:,:,None]*P2)*W[:,:,None]*(w/2)
    assert np.all(N>=0.)
    return Ep.ravel(),N.ravel()


def scalar_cumulative(E,W,cut):
    p=np.sqrt(np.maximum((E-M)*(E+M),0))/M
    upper=E*Es/M+P*p
    lower=(Es*Es+(mp*p)**2)/upper
    # Exact uniform scalar lab energy interval per individual parent.
    span=upper-lower
    assert np.all(span>0.)
    ci=np.hypot(mp,cut)
    return np.array([W@np.clip((c-lower)/span,0.,1.) for c in ci])


def run(eos,kloe,fd,pion):
    bg=Background(eos)
    data=table(kloe)
    xr=read(fd); xp=read(pion)
    match=np.array([xr["rows"][-1]["controls"][-1]["FD_yield"],xp["controls"][-1]["BE"]])
    R=TabulatedRadiation(eos,.05)
    edge=np.unique(np.r_[0.,np.geomspace(.001,18.,75),.1,np.inf])
    rows=[]
    for nt,orders,porders in ((16,(9,28,24,30),(5,24,18)),(26,(14,46,32,48),(8,36,26))):
        source=Source(orders[:3],tail=75.,dark_mass=M)
        x,w=roots_legendre(nt)
        a,b=np.log([.012,.05]);h=(b-a)/2
        c=np.zeros((2,len(edge)))
        for T,dln in zip(np.exp((b+a)/2+h*x),h*w,strict=True):
            A=bg.scale(T)
            _,rho=R.state(A/bg.scale(.05))
            ent=2*np.pi**2/45*T**3*np.exp(bg.gs(np.log(T)))
            dt=(3+bg.gs.derivative()(np.log(T)))/3
            fold=dln*dt/(ent*np.sqrt(rho/3)/MPL)
            for sp,(E,N) in enumerate((source.particles(T,orders[-1]),pion_parents(T,data,*porders))):
                pcut=edge*bg.T0/A
                C=scalar_cumulative(E,N,pcut)
                assert np.all(np.diff(C)>=-1e-14*C[-1])
                c[sp]+=fold*C
        j=int(np.flatnonzero(edge==.1)[0])
        assert np.max(abs(c[:,j]/match-1))<.018
        hist=np.maximum(np.diff(c,axis=1),0.)
        np.testing.assert_allclose(hist.sum(axis=1),c[:,-1],rtol=1e-12)
        rows.append(dict(time_order=nt,lepton_orders=orders,pion_orders=porders,
                         conditional_total_scalar_yield=c[:,-1].tolist(),
                         cold_yield_check_over_previous=(c[:,j]/match).tolist(),
                         shell_yields=hist.tolist(),normalized_shell_counts=(hist/c[:,-1,None]).tolist()))
    u,v=[np.asarray(r["normalized_shell_counts"]) for r in rows]
    error=np.max(abs(np.cumsum(u,axis=1)-np.cumsum(v,axis=1)))
    assert error<.008
    return dict(scope="12–50 MeV direct vacuum lepton FD and KLOE measured-window pion BE, unpolarized open decay. "
                      "Momenta redshift with the reference entropy clock. This is an independently produced "
                      "positive finite-shell input spectrum, not a relic prediction or hierarchy output.",
                KLOE_sha256=hashlib.sha256(kloe.read_bytes()).hexdigest(),species_source_order=["lepton","charged pion KLOE window"],
                conserved_momentum_edges_q_over_T0=[float(t) if np.isfinite(t) else None for t in edge],
                max_normalized_CDF_change=float(error),controls=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","kloe","fd","pion"):
        p.add_argument(k,type=Path)
    a=p.parse_args();print(json.dumps(run(a.eos,a.kloe,a.fd,a.pion),indent=2,allow_nan=False))
