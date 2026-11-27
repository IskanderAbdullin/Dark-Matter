"""Compare bulk drag from timelike annihilation and spacelike elastic scattering."""
import json

import numpy as np

from charged_lepton_elastic_fd_bulk_drag_reference import response, leptons
from positive_fd_portal_quadrature_reference import PositivePortal
from finite_width_vector_pauli_rate_reference import integrator


def run():
    coarse,fine=PositivePortal(),PositivePortal((13,36,26),tail=65.)
    rows=[]
    for T in (.015,.05,.1,.3,.6882487043):
        C=coarse.moments(T,T,0.,common_momentum=True)
        Q,P,V=fine.moments(T,T,0.,common_momentum=True)
        if abs(Q-3*P-V)>3e-14*Q or not 0<P<Q/3:
            raise ValueError("Timelike lab-to-invariant momentum decomposition failed.")
        if np.max(abs(np.array([Q,P,V])/C-1))>8e-6:
            raise ValueError("Timelike momentum moments changed under quadrature refinement.")
        old=fine.moments(T,T,0.)[5]/T
        if abs(Q/old-1)>2e-14:
            raise ValueError("Existing equilibrium heat moment did not replay.")
        adaptive=integrator(T,T,0.)[1,2]/T
        adaptive_P=integrator(T,T,0.,equilibrium_spatial_moment=True)[1,2]/(3*T)
        if max(abs(Q/adaptive-1),abs(P/adaptive_P-1))>3e-6:
            raise ValueError("Independent adaptive heat reference did not replay.")
        qe,ts,pe,_=sum((response(T,0.,m,26,36.) for m in leptons),np.zeros(4))
        if pe<=qe/3 or abs(3*pe-qe-ts)>3e-14*(qe+ts):
            raise ValueError("Spacelike bulk response failed its invariant decomposition.")
        rows.append(dict(T=T,timelike_heat=Q,timelike_bulk=P,spacelike_heat=qe,spacelike_bulk=pe,
                         timelike_bulk_over_heat=P/Q,spacelike_bulk_over_heat=pe/qe,
                         combined_bulk_over_heat=(P+pe)/(Q+qe),
                         elastic_increment_relative_to_timelike_bulk=pe/P,
                         elastic_increment_relative_to_timelike_heat=qe/Q,
                         adaptive_heat_difference=abs(Q/adaptive-1),
                         adaptive_bulk_difference=abs(P/adaptive_P-1)))
    return dict(scope="Local charged-lepton vacuum heat and momentum at full equilibrium; no expansion or transfer model.",
                rows=rows)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
