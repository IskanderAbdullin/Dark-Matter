"""Conditional instantaneous low-reheat MB upper-source reference, GeV units."""
from math import pi,sqrt,exp
from scipy.integrate import quad
from scipy.special import k1,kve
from scipy.optimize import brentq
import json
MPL=1.2209e19
M=10.
GD=.1
EPS=1e-5
ALPHA=1/137.035999084
COUP2=GD**2*EPS**2*4*pi*ALPHA
GAMMA=.005305164626087453
GSTAR=10.75
S0=2891.2
RHOC=1.05375e-5
MPHI=86e-9
M1=.1
M2=.02
LEPTONS=(.00051099895,.1056583755,1.77686)


def gamma(T,mchi,ml=0.,cut=100.):
    """MB forward l+l- -> chi+chi- reaction count/volume/time.
    Spin 2 for each incoming charge, sigma spin-averaged. Vacuum BW.
    The integral is truncated at u=threshold+cut; sufficient at low T.
    """
    u0=2*max(mchi,ml)/T
    def f(x):
        u=u0+x
        s=(T*u)**2
        a=ml**2/s;b=mchi**2/s
        pref=(1+2*a)*(1+2*b)*sqrt(max(0.,1-4*a))*sqrt(max(0.,1-4*b))
        # lambda(s,ml^2,ml^2) / sqrt(s) times sigma; incoming flux cancels.
        return u**6*pref*kve(1,u)*exp(-x)/((s-M*M)**2+(M*GAMMA)**2)
    return COUP2*T**8/(48*pi**5)*exp(-u0)*quad(f,0,cut,epsabs=1e-30,epsrel=1e-8)[0]


def yield_pair(TR,mchi,ml=0.,cut=100.,gstar=GSTAR):
    """Per-charge Y from injection only. g*S=g*rho fixed; T_R to 0."""
    def f(z):
        T=TR*exp(-z)
        s=2*pi*pi/45*gstar*T**3
        h=1.66*sqrt(gstar)*T*T/MPL
        return gamma(T,mchi,ml,cut)/(s*h)
    return quad(f,0,5,epsabs=1e-50,epsrel=1e-7)[0]


def main():
    rows=[]
    for TR in (.005,.008,.01,.015,.02,.025,.03,.05,.1):
        ys=[sum(yield_pair(TR,m,ml) for ml in LEPTONS) for m in (M1,M2)]
        es=[yield_pair(TR,m,0.) for m in (M1,M2)]
        # Each chi1 or antichi1 decay yields one phi, one stable charged chi2.
        # All yields include particle and antiparticle (factor 2).
        Nheavy=2*sum(ys); Nphi_direct=2*ys[0]
        rows.append({'TR_MeV':TR*1000,'Y1_per_charge_MB_all_leptons':ys[0],'Y2_per_charge_MB_all_leptons':ys[1],
           'Y1_per_charge_e_massless_MB':es[0],'Y2_per_charge_e_massless_MB':es[1],
           'omega_heavy_injection_only':M2*Nheavy*S0/RHOC,
           'omega_phi_decay_only':MPHI*Nphi_direct*S0/RHOC,
           'omega_phi_if_all_injected_heavy_convert_1to1':MPHI*(Nheavy+Nphi_direct)*S0/RHOC})
    r=MPHI/M2
    depletions=[]
    for F in (.5,.9,.99):
        q=(1-F)/F*r
        depletions.append({'phi_fraction_min':F,'heavy_over_phi_number_max':q,'heavy_survival_fraction_max':2*q/(1+q)})
    return {'rows':rows,'depletions':depletions,'gstar':GSTAR,'assumptions':'instantaneous reheat, empty dark sector; MB forward vacuum leptons, no backreaction, no hadrons; constant gstar illustrative above 20 MeV; no kinetic solution'}

if __name__=='__main__':
    print(json.dumps(main(),indent=2))
