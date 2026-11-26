"""Exact visible FD source for initially empty vacuum dark pairs."""
import argparse,json
import numpy as np
from positive_fd_portal_quadrature_reference import PositivePortal,rule
from timelike_chemical_heat_eigenmodes_reference import leptons
from sm_dark_entropy_exchange_reference import StandardModel,MPL
from fd_number_heat_state_reference import evaluate,state
from scipy.optimize import brentq
from scipy.special import expit,kve
from timelike_portal_onsager_coverage_reference import integrated,factor,V
from sm_dark_entropy_exchange_reference import C
from three_species_decay_onsager_reference import M,m,mp,P,Es

class Source(PositivePortal):
    def event_moments(self,T,visible_fd=True,*,include_boltzmann=True):
        """Return angular coefficients; optionally factor out exp(-q/T)."""
        E,q,p,W,v=self.geometry(T)
        x,w=rule(self.orders[2])
        b0=np.zeros_like(q)
        b2=b0.copy()
        for ml in leptons:
            sel=E>2*ml
            vl=1-4*ml*ml/E[sel]**2
            root=p[sel,:,None]*np.sqrt(vl)[:,None,None]*x/(2*T)
            a=q[sel,:,None]/(2*T)
            bb=1/(1+np.exp(-a+root)+np.exp(-a-root)+np.exp(-2*a)) if visible_fd else np.ones_like(root)
            massw=(1+2*ml*ml/E[sel]**2)*np.sqrt(vl)
            vv=v[sel]
            c2=(2/3)*vv*vl/(3-vv-vl+vv*vl/3)
            b0[sel]+=massw[:,None]*(bb@(w/2))
            b2[sel]+=(massw*c2)[:,None]*(bb@(w*(3*x*x-1)/4))
        phase=W*np.exp(-q/T) if include_boltzmann else W
        return (q,p,v,phase*b0,phase*b2)

    def particles(self,T,n=24):
        q,p,v,N,A=self.event_moments(T)
        x,w=rule(n)
        ef=q[:,:,None]/2+p[:,:,None]*np.sqrt(v)[:,None,None]*x/2
        g=(N[:,:,None]+A[:,:,None]*(3*x*x-1)/2)*w
        assert g.min()>=-1e-200
        assert abs(g.sum()/(2*N.sum())-1)<1e-12
        assert abs(np.sum(ef*g)/np.sum(q*N)-1)<1e-12
        return ef.ravel(),g.ravel()

def reduced_moments(T,src):
    q,p,v,N,A=src.event_moments(T)
    # Pair count N. E and E² of one final member, integrated and multiplied by 2.
    return np.array([N.sum(),np.sum(q*N),
        np.sum((q*q/2+p*p*v[:,None]/6)*N+(p*p*v[:,None]/15)*A)])

def spectrum(eos,TR,TF,order=(9,24,16),nout=18):
    srcs=[Source(order,dark_mass=j) for j in (m,M)]
    x,w=rule(nout)
    lo,hi=np.log([TF,TR])
    xs=lo+(x+1)*(hi-lo)/2
    dT=w*(hi-lo)/2
    sF=eos.state(TF)[2]
    dn=np.zeros(4)
    rho=np.zeros(4)
    pp=np.zeros(4)
    before=after=0.
    edges=np.r_[0.,np.geomspace(.005,90.,97)]
    hist=np.zeros((3,len(edges)-1))
    inscale=np.zeros(3)
    pcuts=np.array([.005,.01,.025,.05,.1])
    fractional_soft=np.zeros(len(pcuts))
    ds=rule(16)
    def accum(k,mi,ei,g,red):
        p=np.sqrt(np.maximum(0.,(ei-mi)*(ei+mi)))*red
        value=np.hypot(p,mi)
        dn[k]+=g.sum()
        rho[k]+=np.sum(g*value)
        pp[k]+=np.sum(g*p*p)
        if k<3:
            hist[k]+=np.histogram(p/TF,bins=edges,weights=g)[0]
            inscale[k]+=g[p<edges[1]*TF].sum()
    for v,d in zip(xs,dT):
        T=np.exp(v)
        rs,_,ss,cs=eos.state(T)
        H=np.sqrt(rs/3)/MPL
        red=(sF/ss)**(1/3)
        a=cs*d/(3*H*ss**2)
        for src in srcs:
            E,W=src.particles(T)
            W*=a
            if src.dark_mass==m:
                accum(0,m,E,W,red)
                continue
            # Compress the produced-parent energy to its preserved measure
            # before integrating two-body decay. Each parent gives exactly
            # one stable heavy fermion and one scalar.
            b=np.r_[M,np.geomspace(M*(1+1e-7),max(90*T,M*1.2),203),np.inf]
            i=np.searchsorted(b,E,side="right")-1
            W0=np.bincount(i,weights=W,minlength=len(b)-1)
            W1=np.bincount(i,weights=W*E,minlength=len(b)-1)
            ok=W0>0.
            E=W1[ok]/W0[ok]
            W=W0[ok]
            delta=np.sqrt(np.maximum(0.,(E-M)*(E+M)))*P/M
            center=E*Es/M
            low=(Es*Es+(mp/M)**2*(E-M)*(E+M))/(center+delta)
            angular_soft=np.clip((np.hypot(pcuts[None,:]*TF/red,mp)-low[:,None])/(2*np.maximum(delta[:,None],1e-100)),0.,1.)
            fractional_soft+=W@angular_soft
            Ephi=E[:,None]*Es/M+delta[:,None]*ds[0]
            E2=E[:,None]-Ephi
            gw=W[:,None]*ds[1]/2
            before+=np.sum(W*E)
            after+=np.sum(gw*(E2+Ephi))
            accum(1,m,E2,gw,red)
            accum(2,mp,Ephi,gw,red)
            accum(3,M,E,W,red)
    assert abs(before/after-1)<1e-13 and abs(dn[1]/dn[2]-1)<1e-13
    nt=dn[:2].sum()*sF
    em=rho[:2].sum()/dn[:2].sum()
    def therm(t):
        d=state(t,nt)
        return (d['rho']-C*t**4)/nt
    target=brentq(lambda t:therm(t)-em,.1*TF,10*TF,xtol=2e-13)
    fitting=state(target,nt)
    zx,zw=rule(160)
    z=(zx+1)*4
    u=z*z
    f=expit(fitting['eta']-np.hypot(u,m/target))
    w=8*zw*z*u*u*f
    p2=np.sum(w*u*u)/w.sum()*target*target
    return dict(TR=TR,TF=TF,stable_yields=dn[:2].tolist(),
        stable_energy_per_particle=em,FD_moment_matched_T_over_SM=target/TF,
        FD_moment_matched_eta=fitting['eta'],
        stable_actual_p2_over_moment_matched_FD=pp[:2].sum()/dn[:2].sum()/p2,
        stable_yield_to_neutral_thermal=(dn[:2].sum()*sF/evaluate(TF,0.)['n']),
        redshifted_means=(rho[:3]/dn[:3]).tolist(),
        redshifted_p_rms_over_T=(np.sqrt(pp[:3]/dn[:3])/TF).tolist(),
        scalar_fraction_below_005_histogram=inscale[2]/dn[2],
        scalar_exact_decay_angle_soft_fractions=(fractional_soft/dn[2]).tolist(),
        scalar_soft_threshold_p_over_T=pcuts.tolist(),
        boltzmann_tail_fraction=(1-hist.sum(axis=1)/dn[:3]).tolist(),
        scalar_differential_bins=edges.tolist(),scalar_normalized_counts=(hist[2]/dn[2]).tolist(),stable_normalized_counts=(hist[:2].sum(axis=0)/dn[:2].sum()).tolist())

def normalization_check(T,source):
    mD=source.dark_mass
    def fun(E,t):
        s=E*E
        if E<=2*mD:return np.zeros(3)
        K1=kve(1,E/t)*np.exp(-E/t)
        K0=kve(0,E/t)*np.exp(-E/t)
        dynamic=0.
        for ml in leptons:
            if E>2*ml:
                dynamic+=(1+2*ml*ml/s)*np.sqrt(1-4*ml*ml/s)
        a=factor*(s+2*mD*mD)*E**3*np.sqrt(1-4*mD*mD/s)*dynamic
        return a*np.array([t*K1,2*t*t*K1+t*E*K0,
               (6*t**3+t*s)*K1+3*t*t*E*K0])
    r=integrated(T,2*mD,max(60*T,V+30*T),kernel=fun)
    q,p,v,n,b=source.event_moments(T,visible_fd=False)
    return max(abs(np.array([n.sum(),np.sum(q*n),np.sum(q*q*n)])/r-1))

def test(file):
    eos=StandardModel(file)
    p=Source()
    rows=[]
    norms=[]
    for T in (.02,.05,.1,.3,.5):
        a=reduced_moments(T,p)
        b=-p.flow(T,T,-35.)[:2]
        assert np.max(abs(a[:2]/b-1))<2e-11
        for mass in (m,M):
            src=Source((9,24,16),dark_mass=mass)
            z=Source((12,36,26),dark_mass=mass)
            v=reduced_moments(T,src)
            r=reduced_moments(T,z)
            ef,w=src.particles(T)
            er=normalization_check(T,src)
            assert er<3e-6,er
            norms.append(er)
            assert np.max(abs(v/r-1))<3e-6
            assert abs((w@(ef**2))/v[2]-1)<2e-12
    for T in (.02,.03,.05,.075,.1,.15,.3,.5):
        for f in (.9,.65):
            rows.append(spectrum(eos,T,max(.01,T*f)))
    fine=spectrum(eos,.1,.065,(12,32,20),24)
    mid=[r for r in rows if r['TR']==.1 and abs(r['TF']-.065)<1e-12][0]
    assert abs(fine['stable_yield_to_neutral_thermal']/mid['stable_yield_to_neutral_thermal']-1)<2e-5
    return dict(rows=rows,fine_control=fine,max_independent_MB_normalization_error=max(norms),
        scope="Injection-only vacuum e/mu/tau pairs, prompt chi1 decay, collisionless redshift; no later scattering, reverse decay, depletion, or actual abundance claim.")

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument('eos')
    print(json.dumps(test(p.parse_args().eos),indent=2,allow_nan=False))
