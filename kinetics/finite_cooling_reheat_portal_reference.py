"""Conditional reheating-cooling extension of the vacuum portal count source."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from scipy.special import roots_legendre

from conditional_relic_normalization_reference import conversion
from empty_lepton_pair_spectrum_reference import Source,reduced_moments
from late_radiation_background_reference import TabulatedRadiation
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import M,m,mp,P,Es

TF=.001


def instantaneous(T,C):
    d,h=[reduced_moments(T,j) for j in C]
    nD,nH=2*d[0],2*h[0]
    pP2=h[2]-M*M*nH
    phi=Es**2/M**2*h[2]+P*P/(3*M*M)*pP2-mp*mp*nH
    stable=(M-Es)**2/M**2*h[2]+P*P/(3*M*M)*pP2-m*m*nH
    return np.array([nH,nD+nH,phi,d[2]-m*m*nD+stable])

def control():
    # Explicitly integrate both cross-section daughter angular distributions
    # and a second isotropic heavy-decay angle as an independent check.
    r,w=roots_legendre(32)
    S=[Source((9,23,18),dark_mass=mm) for mm in (m,M)]
    errs=[]
    for T in (.006,.012,.027):
        q,d=S[0].particles(T,32)
        h,z=S[1].particles(T,32)
        beta=np.sqrt(np.maximum((h-M)*(h+M),0.))
        def decay_sq(dmass,estar):
            lab=(h[:,None]*estar+beta[:,None]*P*r)/M
            return (z@((lab*lab-dmass*dmass)@(w/2)))
        ref=np.array([z.sum(),d.sum()+z.sum(),decay_sq(mp,Es),
                      d@((q-m)*(q+m))+decay_sq(m,M-Es)])
        m0=instantaneous(T,S)
        errs.append(np.max(abs(ref/m0-1.)))
    assert max(errs)<2e-12
    r,w=roots_legendre(24)
    L=np.log(3)
    residual=[]
    for n in (4.,8.,12.,14.):
        a=n-12.
        I=(8/3)*L/2*(w@np.exp((r+1)*L/2*a))
        J=(8/3)*(np.expm1(a*L)/a if a!=0 else L)
        residual.append(abs(I/J-1.))
    assert max(residual)<2e-13
    return dict(independent_angular_relative_defect=float(max(errs)),
                matter_cooling_power_law_source_quadrature_defect=float(max(residual)))


class SourceFold:
    def __init__(self,table,order,ng,nt):
        self.bg=TabulatedRadiation(table,.01)
        self.grad=self.bg.sg.derivative()
        self.C=[Source(order,dark_mass=mm) for mm in (m,M)]
        self.edges=np.linspace(np.log(TF),np.log(.051),ng)
        self.nodes,self.weights=roots_legendre(nt)
        y=np.array([instantaneous(np.exp(t),self.C) for t in self.edges])
        assert y.min()>0
        self.logR=PchipInterpolator(self.edges,np.log(y),axis=0)
        points=np.linspace(np.log(.0028),np.log(.048),27)
        err=[]
        for l in points:
            r=instantaneous(np.exp(l),self.C)
            err.append(np.max(abs(np.exp(self.logR(l))/r-1)))
        self.interpolation_error=float(max(err))

    def s(self,u):
        return 2*np.pi**2/45*np.exp(3*u+self.bg.sg(u))

    def H(self,u):
        return np.sqrt(np.pi**2/90*np.exp(4*u+self.bg.rg(u)))/MPL

    def fold(self,lo,hi,extra,fn):
        if hi==lo:
            return np.zeros(4)
        e=np.unique(np.clip(np.r_[lo,hi,extra],lo,hi))
        l,h=e[:-1,None],e[1:,None]
        a=((h+l)+(h-l)*self.nodes)/2
        b=(h-l)*self.weights/2
        return np.sum(b[:,:,None]*fn(a),axis=(0,1))

    def post(self,tr):
        r=np.log(tr)
        l=np.log(TF)
        def k(u):
            r2=(self.s(l)/self.s(u))**(2/3)
            a=(3+self.grad(u))/(3*self.H(u)*self.s(u))
            r=np.exp(self.logR(u))
            r[:,:,2:]*=r2[:,:,None]
            return r*a[:,:,None]
        return self.fold(l,r,np.r_[self.edges,self.bg.sg.x],k)

    def before(self,tr,ratio,ln_a=False,alpha=.375,beta=1.5):
        r=np.log(tr)
        sf=(self.s(np.log(TF))/self.s(r))**(2/3)
        def f(ell):
            # ell=ln(T/TR); H/H_R=exp(beta*ell/alpha), a/a_R=exp(-ell/alpha).
            r0=np.exp(self.logR(r+ell))
            r0[:,:,2:]*=sf*np.exp(-2*ell[:,:,None]/alpha)
            return r0*np.exp(-(3+beta)*ell[:,:,None]/alpha)/(self.H(r)*self.s(r))
        edge=self.edges-r
        p=np.log(ratio)
        if ln_a:
            return self.fold(-p/alpha,0.,-edge/alpha,lambda l: f(-alpha*l))
        return self.fold(0.,p,edge,f)/alpha

def energy_envelope(F,data,temperature,rescale_p2=1.):
    if np.max(data[:2])==0.:
        return [0.,0.]
    u=np.log(temperature)
    mass=np.array([mp,m])
    q2= data[2:]/data[:2]*rescale_p2
    # Jensen with concave sqrt: number and one p^2 moment bound energy
    # at this epoch without choosing a distribution.
    d=F.s(u)*data[:2]
    rho=3*MPL*MPL*F.H(u)**2
    return [float(d@mass/rho),float(d@np.sqrt(mass*mass+q2)/rho)]


def run(table):
    from archived_collisionless_streaming_reference import Background
    bg=Background(table)
    anchor=float(bg.scale(TF))
    masses=np.array([mp,m])
    def upper(rms):
        q=np.asarray(rms)*anchor
        a=np.array([bg.streams(q,mass,anchor,bg.aeq,n) for mass in masses
                    for n in (6,10)]).reshape(2,2,2)
        selected=np.array([a[j,:,j] for j in range(2)])
        assert np.max(abs(selected[:,0]/selected[:,1]-1.))<1e-8
        return selected[:,1].tolist()
    F=conversion(table)["Omega_h2_per_GeV_times_yield"]
    variants=[]
    for order,n,nt in (((7,18,14),130,4),((10,26,20),230,6)):
        K=SourceFold(table,order,n,nt)
        assert K.interpolation_error<0.004
        def val(tr,rat):
            a=K.post(tr)
            b=K.before(tr,rat)
            return a,b
        result=[]
        for f in (1.,1.5,2.,3.,5.):
            def root(tr):
                a,b=val(tr,f)
                return F*m*(a+b)[1]-.120
            TR=brentq(root,.0028,.0096,xtol=1e-12)
            a,b=val(TR,f)
            c=a+b
            if f>1.:
                v=K.before(TR,f,ln_a=True)
                assert np.max(abs(v/b-1.))<2e-12
            result.append(dict(Tmax_over_TR=f,T_R_GeV=TR,Tmax_GeV=f*TR,
                pre_fraction_of_scalar_and_stable_number=(b[:2]/c[:2]).tolist(),
                scalar_prompt_fixed_count_Omega_h2=F*mp*c[0],
                pre_dark_to_visible_at_TR_energy_interval=energy_envelope(
                    K,b,TR,(K.s(np.log(TR))/K.s(np.log(TF)))**(2/3)),
                combined_dark_to_visible_at_1MeV_energy_interval=energy_envelope(K,c,TF),
                pre_versus_post_rms_momentum_over_1MeV=[
                    None if f==1. else (np.sqrt(b[2:]/b[:2])/TF).tolist(),
                    (np.sqrt(a[2:]/a[:2])/TF).tolist()],
                combined_rms_momentum_over_1MeV=(np.sqrt(c[2:]/c[:2])/TF).tolist(),
                post_1MeV_mean_streaming_distance_upper_to_equality_Mpc=upper(np.sqrt(c[2:]/c[:2]))))
        fixed=[]
        for Tmax in (.010,.020,.040):
            root=brentq(lambda T: F*m*(K.post(T)+K.before(T,Tmax/T))[1]-.120,
                       .0028,min(.0096,Tmax),xtol=1e-12)
            a,b=val(root,Tmax/root)
            c=a+b
            fixed.append(dict(fixed_Tmax_GeV=Tmax,T_R_GeV=root,
                 prompt_scalar_fixed_count_omega_h2=F*mp*c[0],
                 pre_fraction_stable=float(b[1]/c[1])))
        variants.append(dict(production_quadrature=order,rate_log_grid_nodes=n,fixed_Tmax_cases=fixed,
                 max_independent_rate_interpolation_relative_error=K.interpolation_error,rows=result))
    d=[abs(a["T_R_GeV"]/b["T_R_GeV"]-1) for a,b in zip(variants[0]["rows"],variants[1]["rows"])]
    assert max(d)<3e-5
    # Additional deliberately prescribed clocks, at fixed beta=3/2,
    # vary only the visible cooling law, not the microscopic source.
    scan=[]
    for alpha in (.25,.375,.50):
        rat=3.
        r=brentq(lambda T: F*m*(K.post(T)+K.before(T,rat,alpha=alpha))[1]-.120,
                 .0028,.0096,xtol=1e-12)
        y=K.before(r,rat,alpha=alpha)
        j=K.before(r,rat,alpha=alpha,ln_a=True)
        assert np.max(abs(y/j-1))<2e-12
        z=K.post(r)+y
        scan.append(dict(T_power_exponent=alpha,Tmax_over_TR=rat,T_R_GeV=r,
                          stable_pre_fraction=float(y[1]/z[1]),scalar_prompt_omega_h2=F*mp*z[0]))
    return dict(controls=variants,variant_cooling_clock_at_TR_ratio3=scan,
        analytical_cross_checks=control(),max_relative_root_change=max(d),target_stable_Omega_h2=.120,
        streaming_bound_scope="In the existing reference radiation+matter background, from the 1 MeV snapshot to equality only. Jensen bounds mean speed using the collisionless redshifted RMS momentum. Not a transfer function or a claimed actual mean.",
        eos_sha256=hashlib.sha256(Path(table).read_bytes()).hexdigest(),
        scope="Prescribed cooling branch of matter-dominated reheating with T~a^-3/8, H~a^-3/2, glued continuously in H to entropy-conserving radiation at T_R. No heating branch, inflaton dynamics, nonthermal SM, backreaction, inverse, scalar pairs or deformation of the vacuum matrix elements. Momentum is collisionless after the moment of instantaneous parent decay. Fixed-count calibration, not full reheating/BBN or scalar abundance fit.")


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("eos");a=p.parse_args()
    print(json.dumps(run(a.eos),indent=2,allow_nan=False))
