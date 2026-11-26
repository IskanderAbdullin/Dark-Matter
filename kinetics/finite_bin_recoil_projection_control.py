"""Finite-bin detailed-balance projection negative control; NOT physical chi-e 5D scattering."""
import math
from numpy.polynomial.legendre import leggauss
import numpy as np

def paired_integrals(n, width, T=1., Ei=2., Ej=1.):
    x,w=leggauss(n)
    a=Ei+width*x/2
    b=Ej+width*x/2
    A,B=np.meshgrid(a,b,indexing="ij")
    weights=np.outer(w,w)/4
    # Positive, symmetric microscopic event kernel; intentionally synthetic.
    kernel=np.exp(-0.17*(A+B))*(1+0.12*(A-B)**2)
    forward=np.sum(weights*kernel*np.exp((A-B)/(2*T)))
    reverse=np.sum(weights*kernel*np.exp(-(A-B)/(2*T)))
    return float(forward/reverse/math.exp((Ei-Ej)/T)-1)

def test_projection():
    expected={.8:.02277004610530553,.4:.00570907623234973,.2:.0014282495281294771,.1:.00035712279618960174,.05:.0000892844613802346}
    for width, ref in expected.items():
        vals=[paired_integrals(n,width) for n in (8,16,32,64)]
        assert max(abs(v-ref) for v in vals)<2e-12, (width,vals)
    # Error shrinks by ~4 when widths halve: an O(width^2) projection bias.
    assert 3.9<paired_integrals(32,.4)/paired_integrals(32,.2)<4.1
if __name__=="__main__":
    test_projection()
    for width in (.8,.4,.2,.1,.05):
        print(f"width={width:.2f} relative_balance_defect={paired_integrals(32,width):.12g}")
