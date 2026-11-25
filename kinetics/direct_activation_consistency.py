"""Direct quartic activation: density, thermal and radiative consistency.

Potential V= m^2 A^2/2 + chi^2 A^4/(24 Lambda^2)
             + lambda_chi (chi^2-v^2)^2/4.
lambda_late=v^2/Lambda^2. Gaussian NR medium, standard two-loop
double bubble at MS scale m, and explicit illustrative tuning diagnostics.
No universal EFT no-go: counterterms/symmetries may change threshold.
"""
import argparse,json,math
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from quartic_consistency import lambda_from_sigma
from late_activation_overlap import H,critical_density
from resonant_mediator_inverse import HBAR,HBARC,M

RHO0_CGS=.12*1.8783416e-29
# 1 g in eV divided by the volume (1 cm in eV^-1)^3.
CGS_RHO_TO_EV4=HBARC**3/1.7826619216278976e-33
MPL_REDUCED=2.435e27
T0=2.7255*8.617333262e-5
I=M*M/(16*np.pi**2)
KPC=3.085677581e21/HBARC

def rho_tot(z):return 3*(H(z)*HBAR)**2*MPL_REDUCED**2
def radiation_T4(z,dneff):return (7/4)*(4/11)**(4/3)*dneff*(T0*(1+z))**4

def sphere(k):
    # Boundary chi(R)=v gives an upper envelope to the interior solution.
    if k<1e-5:return 1.,0.
    if k>1e4:return 3/(4*k),math.log10(2*k)-k/math.log(10)
    log_sinh=lambda x:x-math.log(2)+math.log(-math.expm1(-2*x))
    log_center=math.log(k)-log_sinh(k)
    # Integrate in distance from boundary in units of screening length.
    def integrand(s):
        x=1-s/k
        if x<=0:return 0.
        u=math.exp(log_sinh(k*x)-log_sinh(k)-math.log(x))
        return 3*x*x*u**4/k
    avg=quad(integrand,0,min(k,40),epsabs=1e-12,epsrel=1e-9)[0]
    return avg,log_center/math.log(10)

def direct(sigma,z=3.):
    lam=float(lambda_from_sigma(sigma))
    rho=RHO0_CGS*(1+z)**3*CGS_RHO_TO_EV4;V=rho/M**2
    # I is magnitude of renormalized tadpole at MS scale m.
    double_bubble=lam*I*I/8
    # X=mu^2 v^2, critical density mass squared = V^2/(4 Lambda^2).
    Xden=lam*V*V/4
    rows=[]
    rh=critical_density(np.sqrt(2)*44.1e5)*CGS_RHO_TO_EV4;Vh=rh/M**2
    for L in [M,10*M,1e6,8.6e13,1e15,1e16,1e17,1e18,1e19,MPL_REDUCED]:
        v=np.sqrt(lam)*L;mu=V/(2*L)
        m_in=np.sqrt(max(Vh*Vh-V*V,0.))/(2*L)
        k=m_in*.6*KPC;avg,center=sphere(k)
        rows.append(dict(cutoff_eV=L,vev_eV=v,mu_eV=mu,mu_over_H=mu/(HBAR*H(z)),
            halo_kappa=m_in*.6*KPC,halo_cross_section_volume_upper_ratio=avg,halo_log10_central_field_upper_ratio=center))
    return dict(sigma=sigma,z=z,lambda_late=lam,rho_DM_ev4=rho,rho_total_ev4=rho_tot(z),variance_NR_ev2=V,
        vacuum_two_loop_threshold_ev4=double_bubble,
        threshold_over_1percent_total=double_bubble/(.01*rho_tot(z)),
        density_trigger_tree_X_ev4=Xden,density_trigger_tree_vacuum_ev4=Xden/4,
        threshold_over_density_tree=double_bubble/(Xden/4),
        induced_A_mass_squared_ev2=lam*I/2,
        induced_lower_operator_medium_over_quartic_medium=2*I/V,
        thermal=[dict(delta_neff=d,T_D_on_eV=radiation_T4(z,d)**.25,
            X_quartic_loop_upper=np.pi**2/9*radiation_T4(z,d),
            threshold_over_thermal_X_over4=4*double_bubble/(np.pi**2/9*radiation_T4(z,d)))
            for d in [.03,.1,.2,.5]],
        halo_density_cgs=rh/CGS_RHO_TO_EV4,halo_radius_kpc=.6,
        density_trigger_unscreened_cutoff_min_eV=np.sqrt(max(Vh*Vh-V*V,0.))*.6*KPC/2,
        density_trigger_mu_over_H10_cutoff_max_eV=V/(20*HBAR*H(z)),cutoffs=rows,
        scope='Gaussian NR density trigger; static linear sphere upper envelope; MS two-loop naturalness diagnostic, not invariant observable or universal no-go.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args()
    r=[direct(s,z) for s in [.1,1.25] for z in [1.,3.,5.,10.]]
    Path(a.output).write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
    for v in r:
        print(v['sigma'],v['z'],{k:v[k] for k in ('vacuum_two_loop_threshold_ev4','threshold_over_1percent_total','threshold_over_density_tree','induced_lower_operator_medium_over_quartic_medium','density_trigger_unscreened_cutoff_min_eV','density_trigger_mu_over_H10_cutoff_max_eV')},flush=True)
