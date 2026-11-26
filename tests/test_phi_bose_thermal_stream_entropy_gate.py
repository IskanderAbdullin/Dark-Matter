"""Independent order of integration and cold entropy benchmark."""
from math import exp,sqrt
from scipy.integrate import quad
from kinetics.phi_bose_thermal_stream_entropy_gate import main,AEQ,C,H0,OR,OM,OL,weight,moment

def test_thermal_stream_entropy_gate():
    r=main()
    for k in ('scalar_BE','reference_FD'):
        d=r[k]
        fd=d['fd']
        q=d['T0_eV']/d['mass_eV']
        norm=moment(2,fd)
        def f(x):
            a=AEQ*exp(x)
            v=quad(lambda y:weight(y,fd)*(q*y)/sqrt((q*y)**2+a*a),0,60,epsabs=1e-11,epsrel=1e-8)[0]/norm
            return a*C/H0*v/sqrt(OR+OM*a+OL*a**4)
        val=quad(f,-35,0,epsabs=1e-9,epsrel=1e-8)[0]
        assert abs(val/d['free_stream_to_eq_Mpc']-1)<1e-7
        assert d['free_stream_to_eq_Mpc']<d['pure_RD_to_eq_Mpc']
    assert 0<r['cold_entropy_gate']['retained_fraction_max_from_initial_thermal_benchmark']<.00022
    assert 0<r['cold_entropy_gate']['saturating_cold_count_fraction']<1
