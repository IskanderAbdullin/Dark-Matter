import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'kinetics'))
import numpy as np
from zero_mode_cutoff import exchange
from condensate_chemistry_screen import amplitude,rate
from hartree_condensate_limits import saturation,match

def test_cutoff_strip_exact_box_and_subdivision():
    A,L=4.7,2.3
    for k in (0.,1e-12,.1,.8,1.7,2.3,3.):
        v=exchange([0,L],[A],k)
        expect=A*A*(max(L-k,0)**2-max(L-2*k,0)**2)
        B=A/2*max(L-2*k,0)**2
        assert np.isclose(v['positive_cutoff_strip'],expect,atol=3e-14,rtol=1e-5)
        assert np.isclose(v['spontaneous_loss'],B,atol=1e-14)
        vv=exchange(np.linspace(0,L,31),np.full(30,A),k)
        assert np.isclose(vv['net'],v['net'],atol=1e-12,rtol=1e-12)

def test_threshold_amplitude_and_published_prefactor():
    m,lam,N=1.7,.03,.6
    p=np.array([[m,0,0,0]]*4+[[-2*m,np.sqrt(3)*m,0,0],[-2*m,-np.sqrt(3)*m,0,0]])
    assert np.isclose(amplitude(p,m,lam),-lam**2/m**2,rtol=1e-14)
    v,n,f=rate(N,m,lam,0.)
    assert np.isclose(v/2/N**4,np.sqrt(3)*lam**4/(12288*np.pi*m**8))
    assert n==v and f==0
    v,n,f=rate(N,m,lam,4.)
    assert np.isclose(n/v,1+2*f)

def test_saturated_hartree_thermodynamics():
    N,m,lam,T=.25,1e-5,1e-8,.13
    b=saturation(N,T,m,lam,128)
    assert b['N0']>0 and b['excited_fraction']<.01
    assert abs(b['delta']-lam*b['V']/2)<1e-15
    assert abs(b['trace_residual'])<1e-12
    c=match(N,b['U'],m,lam,128)
    assert c['saturated'] and abs(c['T']/T-1)<1e-8
    assert abs(c['first_law_relative_error'])<1e-5
