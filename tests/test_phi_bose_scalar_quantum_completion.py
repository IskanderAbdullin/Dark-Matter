import math
from kinetics.phi_bose_scalar_quantum_completion import coefficients,main,halo

def test_matching_and_scale():
    r=main()
    a=coefficients(.1)
    assert math.isclose(a[0],6.118843231538268e-10,rel_tol=1e-10)
    assert math.isclose(a[2],coefficients(.02)[2],rel_tol=1e-10)
    assert .46<r["loop_size_halo"]["TF_radius_AU"]<.48
    assert r["lambda_for_1kpc_TF"]>9e4
    assert r["cold_diagnostic"][0]["ratio_quartic_quadratic"]>20
