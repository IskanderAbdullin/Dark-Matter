"""One-loop fermion CW matching, zero-momentum potential, not full physical amplitude."""
import math,json
M=.1;m=.02;y=.001;mu_phi=86e-9;delta=M-m;MPL=1.2209e19
GEV_INV_CM=1.973269804e-14;GEV_G=1.78266192e-24;AU_CM=1.495978707e13;KPC_CM=3.0856775814913673e21
T0=2.34865e-13;G0=3.91;RHOC=1.05375e-5
def coefficients(scale):
    def fp(z):return 4*z**3*(math.log(z*z/(scale*scale))-1)
    def fpp(z):return 4*z*z*(3*math.log(z*z/(scale*scale))-1)
    def fppp(z):return 8*z*(3*math.log(z*z/(scale*scale))+2)
    a=fp(M)-fp(m);b=fpp(M)+fpp(m);c=fppp(M)-fppp(m)
    c1=a/delta;c2=-a/delta**3+b/(2*delta**2);c3=2*a/delta**5-b/delta**4+c/(6*delta**3)
    return [-2*y*y*c1/(16*math.pi**2),-24*y**4*c2/(16*math.pi**2),-720*y**6*c3/(16*math.pi**2)]

def numeric(scale):
    """Optional independent 70-digit differentiation (requires mpmath)."""
    import mpmath as mp
    mp.mp.dps=70
    A=mp.mpf(str(M));B=mp.mpf(str(m));Y=mp.mpf(str(y));S=mp.mpf(str(scale))
    def V(z):
        d=mp.sqrt((A-B)**2+4*Y**2*z*z);u=(A+B+d)/2;v=(A+B-d)/2
        return -sum(w**4*(mp.log(w*w/(S*S))-mp.mpf('1.5')) for w in (u,v))/(16*mp.pi**2)
    return [float(mp.diff(V,mp.mpf('0'),n)) for n in (2,4,6)]

def halo(lam):
    sigma=lam*lam/(128*math.pi*mu_phi**2)
    sigma_per_mass=sigma*GEV_INV_CM**2/(mu_phi*GEV_G)
    tf=math.pi*math.sqrt(lam/(32*math.pi))*MPL/mu_phi**2*GEV_INV_CM/AU_CM if lam>0 else None
    return dict(lambda_total=lam,sigma_per_mass_cm2_g=sigma_per_mass,TF_radius_AU=tf)

def cold_diagnostic(T,g,lam,c6):
    n0=.12*RHOC/mu_phi*GEV_INV_CM**3
    a=T0/T*(G0/g)**(1/3);n=n0/a**3
    r=lam*n/(6*mu_phi**3)
    H=1.66*math.sqrt(g)*T*T/MPL
    amp=lam*lam/mu_phi**2-c6
    rate=math.sqrt(3)/(6144*math.pi)*n**3/mu_phi**4*amp**2
    return dict(T_GeV=T,g_assumed=g,a=a,n_GeV3=n,phiamp_GeV=math.sqrt(2*n/mu_phi),ratio_quartic_quadratic=r,rate_over_H_incoherent_dilute_diagnostic=rate/H)

def main(high_precision=False):
    rows=[]
    for s in (.01,.02,.05,.1,.2,1.):
        c=coefficients(s)
        if high_precision:
            n=numeric(s)
            assert max(abs(a/b-1) for a,b in zip(c,n))<1e-11
        rows.append(dict(scale_GeV=s,delta_m2_GeV2=c[0],delta_lambda=c[1],c6_GeV_minus2=c[2],abs_delta_m2_over_mphi2=abs(c[0])/mu_phi**2,delta_m_curvature_eV=math.sqrt(abs(c[0]))*1e9))
    assert max(abs(r['c6_GeV_minus2']/rows[0]['c6_GeV_minus2']-1) for r in rows)<1e-10
    dm2_dlog=(rows[4]['delta_m2_GeV2']-rows[3]['delta_m2_GeV2'])/math.log(2)
    dlambda_dlog=(rows[4]['delta_lambda']-rows[3]['delta_lambda'])/math.log(2)
    expect_m=y*y*(M*M+m*m+M*m)/math.pi**2
    expect_l=6*y**4/math.pi**2
    assert abs(dm2_dlog/expect_m-1)<1e-12 and abs(dlambda_dlog/expect_l-1)<1e-12
    lam_loop=abs(coefficients(.1)[1]);c6=coefficients(.1)[2]
    lam_kpc=(KPC_CM/GEV_INV_CM*mu_phi**2/(math.pi*MPL))**2*32*math.pi
    sigmax=1.25;lam_sidm=math.sqrt(sigmax*128*math.pi*mu_phi**3*GEV_G/GEV_INV_CM**2)
    n0=.12*RHOC/mu_phi*GEV_INV_CM**3
    tcrit=T0*(6*mu_phi**3/lam_loop/n0*(G0/10.75))**(1/3)
    return dict(rows=rows,dm2_dlog=dm2_dlog,dlambda_dlog=dlambda_dlog,loop_size_halo=halo(lam_loop),lambda_for_1kpc_TF=lam_kpc,lambda_for_sigma_1p25=lam_sidm,halo_sigma_1p25=halo(lam_sidm),cold_Tcrit_GeV=tcrit,cold_diagnostic=[cold_diagnostic(T,10.75,lam_loop,c6) for T in (.02,.005,.001)],high_precision_checked=high_precision)

if __name__=='__main__':
    import sys
    print(json.dumps(main('--high-precision' in sys.argv),indent=2))
