"""Independent explicit gamma-matrix check of unequal-mass Dirac-vector trace."""
import math
import numpy as np


def gamma_trace(M,m,s,c):
    I=np.eye(2,dtype=complex); z=np.zeros((2,2),complex)
    pauli=[np.array([[0,1],[1,0]],complex),
           np.array([[0,-1j],[1j,0]],complex),
           np.array([[1,0],[0,-1]],complex)]
    gam=[np.block([[I,z],[z,-I]])]
    gam += [np.block([[z,v],[-v,z]]) for v in pauli]
    eta=[1,-1,-1,-1]; eye=np.eye(4,dtype=complex)
    def slash(v):
        return sum(eta[i]*v[i]*gam[i] for i in range(4))
    p=math.sqrt((s-(M+m)**2)*(s-(M-m)**2))/(2*math.sqrt(s))
    E=(s+M*M-m*m)/(2*math.sqrt(s))
    w=(s+m*m-M*M)/(2*math.sqrt(s))
    si=math.sqrt(1-c*c)
    a=np.array([E,0,0,p]);b=np.array([w,0,0,-p])
    ap=np.array([E,p*si,0,p*c]);bp=np.array([w,-p*si,0,-p*c])
    ca=[[np.trace((slash(ap)+M*eye)@gam[u]@(slash(a)+M*eye)@gam[v])
         for v in range(4)] for u in range(4)]
    cb=[[np.trace((slash(bp)+m*eye)@gam[u]@(slash(b)+m*eye)@gam[v])
         for v in range(4)] for u in range(4)]
    exact=sum(eta[u]*eta[v]*ca[u][v]*cb[u][v]
              for u in range(4) for v in range(4))/4
    t=-2*p*p*(1-c)
    D=s-M*M-m*m
    predicted=2*(D*D+(-D-t)**2+2*(M*M+m*m)*t)
    return exact,predicted


def test_vector_spin_trace_for_unequal_masses():
    for M,m,s,c in [(.1,1.,2.,-.6),(.1,1.,5.,.8),
                     (1.1,.4,10.,.15),(.8,.2,3.,-.9)]:
        measured,expected=gamma_trace(M,m,s,c)
        assert abs(measured.imag)<1e-12
        assert math.isclose(measured.real,expected,rel_tol=5e-13)
