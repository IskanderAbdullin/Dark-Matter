"""Positive on-shell lab-frame 2↔2 quantum Yukawa event integration.

Explicit Dirac chi2 (spin2), anti-chi2 (spin2) -> two IDENTICAL real phi,
with intermediate off-diagonal Dirac chi1 in t/u channels. All four
laboratory event energies are exactly on shell and energy conserving.

Physical event weight:
  dW = [gchi*ganti/(8π^4)] p^2 dp q^2 dq dcos(angle_in)
       * vMoller * [dσ_forward/dcos(theta_CM)] dcos(theta_CM)
       * d(azimuth_CM)/(2π).
dσ/dcos includes outgoing identical phi 1/2!, initial spin average 1/4.
This is the first real-mass pair continuum quantum WEAK collision source,
NOT a joint 1↔2 +2↔2 finite-time scheme or physical FRW/SM freezeout.
"""
import json
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import expit
from kinetics.phi_bose_pair_annihilation_mb import (
    M1 as MV,M2 as M,MPHI as MP,Y
)

def differential_pair_sigma(s,cos_theta,m_in=M,m_virtual=MV,y=Y):
    """Vectorized positive spin-average identical-real-phi dσ/dcos CM."""
    M=m_in;MV=m_virtual;Y=y
    if min(m_in,m_virtual)<=0 or y<0:raise ValueError('Positive Dirac masses and coupling required')
    cos_theta=np.atleast_1d(np.asarray(cos_theta,float))
    if s<=4*M*M:raise ValueError("CM pair threshold")
    E=math.sqrt(s)/2
    pin=math.sqrt(s/4-M*M)
    pout=math.sqrt(s/4-MP*MP)
    ab=E*E+pin*pin
    qz=[pin-pout*cos_theta,pin+pout*cos_theta]
    qsq=[-pin*pin-pout*pout+2*pin*pout*cos_theta,
         -pin*pin-pout*pout-2*pin*pout*cos_theta]
    accum=np.zeros_like(cos_theta)
    for i in (0,1):
        for j in (0,1):
            aq=pin*qz[i]
            ar=pin*qz[j]
            qr=qsq[i] if i==j else pout*pout-pin*pin
            four=-2*aq*ar-ab*qr
            two=2*M*MV*(aq+ar)+MV*MV*ab-M*M*qr
            numerator=four+two-M*M*MV*MV
            accum+=numerator/((qsq[i]-MV*MV)*(qsq[j]-MV*MV))
    dsigma=Y**4*accum*pout/(64*math.pi*s*pin)
    if np.min(dsigma)<-1e-18*np.max(np.abs(dsigma)):
        raise ArithmeticError("Negative squared tree amplitude")
    return np.maximum(dsigma,0.)

class PhysicalPairEvents:
    """One finite 5D positive physical event table for two Dirac charge sectors.

    A psi average is used to project CM scattering orientation against the
    moving LAB plasma frame. Outputs 4 physical LAB E and one shared
    nonnegative W per event; no arbitrary cell deposits or hidden sources.
    """
    def __init__(self,P=20.,nrad=16,nrelative=12,ncm=48,naz=8,
                 m_in=M,m_virtual=MV,y=Y):
        M=m_in; MV=m_virtual; Y=y
        self.m_in=m_in;self.m_virtual=m_virtual;self.y=y
        if m_in<=0 or m_virtual<=0 or y<0:raise ValueError('Invalid Yukawa pair masses or y')
        if P<=0 or min(nrad,nrelative,ncm,naz)<1:
            raise ValueError("Positive cutoff and Gaussian orders required")
        self.P=P
        self.orders=(nrad,nrelative,ncm,naz)
        x,wx=leggauss(nrad)
        p=P*(1+x)/2
        wp=P*wx/2
        ci,wci=leggauss(nrelative)
        u,wu=leggauss(ncm)
        u=14*u
        c=np.tanh(u)
        wc=14*wu*(1-c*c)
        sint=np.sqrt(1-c*c)
        cospsi=np.cos(2*math.pi*np.arange(naz)/naz)
        out=[[],[],[],[],[]]
        pref=4/(8*math.pi**4)
        for pa,wpa in zip(p,wp):
            ea=math.hypot(pa,M)
            for pb,wpb in zip(p,wp):
                eb=math.hypot(pb,M)
                etot=ea+eb
                for cij,wcij in zip(ci,wci):
                    s=2*M*M+2*(ea*eb-pa*pb*cij)
                    pcm=math.sqrt(max(0.,s/4-M*M))
                    ecm=math.sqrt(s)/2
                    kout=math.sqrt(s/4-MP*MP)
                    vm=math.sqrt(max(0.,s*(s-4*M*M)))/(2*ea*eb)
                    ds=differential_pair_sigma(s,c,m_in=M,m_virtual=MV,y=Y)
                    w=pref*pa*pa*wpa*pb*pb*wpb*wcij*vm*wc*ds
                    gamma=etot/math.sqrt(s)
                    bz=(ea/gamma-ecm)/pcm
                    bt=math.sqrt(max(0.,1-1/(gamma*gamma)-bz*bz))
                    e3=gamma*(ecm+kout*(bz*c[:,None]
                         +bt*sint[:,None]*cospsi[None,:]))
                    e4=etot-e3
                    if np.min(e3)<MP*(1-1e-8) or np.min(e4)<MP*(1-1e-8):
                        raise ArithmeticError("Off-shell scalar event")
                    out[0].append(np.full(e3.size,ea))
                    out[1].append(np.full(e3.size,eb))
                    out[2].append(e3.ravel())
                    out[3].append(e4.ravel())
                    out[4].append(np.repeat(w/naz,naz))
        self.e_chi,self.e_anti,self.e_phi3,self.e_phi4,self.W=[
            np.concatenate(x) for x in out]
        if np.min(self.W)<=0:raise ArithmeticError("Nonpositive physical pair event W")
        self.event_count=len(self.W)
        self.max_kinematic_energy_defect=float(np.max(np.abs(
            self.e_phi3+self.e_phi4-self.e_chi-self.e_anti)))

    def MB_forward(self,T):
        """Independent LAB p,q,relative angle+CM angle integral, no quantum."""
        if T<=0:raise ValueError("T>0")
        return float(np.dot(self.W,np.exp(-(self.e_chi+self.e_anti)/T)))

    def quantum_weak(self,Tchi=1.,Tanti=1.,Tphi=1.,
                     muchi=0.,muanti=0.,muphi=0.,bose_curvature=0.):
        """Microscopic FD+BE shared-event source and H for arbitrary shown eta."""
        if min(Tchi,Tanti,Tphi)<=0:raise ValueError("T>0")
        a=(self.e_chi-muchi)/Tchi
        b=(self.e_anti-muanti)/Tanti
        c=(self.e_phi3-muphi)/Tphi+bose_curvature*(1-np.cos(self.e_phi3))
        d=(self.e_phi4-muphi)/Tphi+bose_curvature*(1-np.cos(self.e_phi4))
        if min(float(np.min(c)),float(np.min(d)))<=0:
            raise ValueError("Bosonic eta must be positive")
        f=expit(-a)
        g=expit(-b)
        h=1/np.expm1(c)
        k=1/np.expm1(d)
        F=f*g*(1+h)*(1+k)
        R=(1-f)*(1-g)*h*k
        J=F-R
        affinity=c+d-a-b
        weightJ=self.W*J
        rate=float(np.sum(weightJ))
        q1=-float(np.dot(weightJ,self.e_chi))
        q2=-float(np.dot(weightJ,self.e_anti))
        qp=float(np.dot(weightJ,self.e_phi3+self.e_phi4))
        H=float(np.dot(weightJ,affinity))
        return dict(R_events_GeV4=rate,dN_chi=-rate,dN_anti=-rate,
                    dN_phi=2*rate,dQ_chi=q1,dQ_anti=q2,dQ_phi=qp,
                    total_energy_defect=q1+q2+qp,
                    conserved_U1_defect=-rate-(-rate),
                    quantum_H_GeV4=H,
                    min_micro_entropy_JA=float(np.min(J*affinity)),
                    event_count=self.event_count,
                    min_Bose_eta=min(float(np.min(c)),float(np.min(d))),
                    max_energy_shell_defect=self.max_kinematic_energy_defect,
                    class_ready=False)

def audit():
    data=[]
    for nrad,nrel,ncm,naz in ((12,12,32,4),(16,12,48,8),
                              (20,16,64,8),(24,24,64,6)):
        m=PhysicalPairEvents(20.,nrad,nrel,ncm,naz)
        data.append(dict(grid=[nrad,nrel,ncm,naz],
                         MB_forward_T1=m.MB_forward(1.),
                         common_T_equilibrium=m.quantum_weak(muchi=.05,muanti=-.05),
                         chemical_mu_phi_negative=m.quantum_weak(
                            muchi=.05,muanti=-.05,muphi=-.2),
                         curved_nonthermal=m.quantum_weak(
                            Tchi=.9,Tanti=1.1,Tphi=1.,bose_curvature=.01)))
    return dict(schema_version=1,scope="on-shell positive finite-cutoff 5D pair quantum weak events and rate N/E/H, no finite-time/SM/FRW/zero-mode/86eV relic",
                rows=data,class_ready=False)

if __name__=="__main__":
    print(json.dumps(audit(),indent=2,allow_nan=False))
