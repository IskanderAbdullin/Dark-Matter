"""Independent physical-measure test of a ONE-eta-per-energy-cell Bose closure.

A cellwise constant eta is monotone in the cell particle number; however
its sole energy moment is locked to Ebar*N. A continuum on-shell decay
samples variable daughter energies inside cells, so energy and thermal
fixed points cannot both be exactly retained. Mathematical/numerical gate,
not a proposed physical dark-matter solver or a condensate model.
"""
import json
import math
import numpy as np
from numpy.polynomial.legendre import leggauss

M1, M2, MPH, YUKAWA = .1, .02, 86e-9, .001


def cells(n, upper, order=64):
    """Exact cell capacity W, Gauss physical thermal Bose N/U, constant eta."""
    if n < 4 or upper <= MPH or order < 12:
        raise ValueError("positive cells and physical quadrature required")
    e = MPH + (upper-MPH)*np.linspace(0., 1., n+1)**2
    p = np.sqrt(e*e-MPH*MPH)
    w = np.diff(p**3)/(6*math.pi**2)
    x, wx = leggauss(order)
    nt = np.zeros(n)
    ut = np.zeros(n)
    q = np.zeros(n)
    for i in range(n):
        pp = p[i]+(x+1)*(p[i+1]-p[i])/2
        weights = wx*(p[i+1]-p[i])/2*pp**2/(2*math.pi**2)
        ee = np.hypot(pp,MPH)
        occ = 1/np.expm1(ee)
        nt[i] = weights@occ
        ut[i] = weights@(ee*occ)
        q[i] = weights@ee
    eta = np.log1p(w/nt)
    ebar = q/w
    assert np.all(eta > 0)
    return e, w, nt, ut, eta, ebar


def mapped_events(pmax=2., nr=24, nangle=48):
    """Independent positive boosted log(Ephi) angular quadrature and physical rate."""
    if pmax <= 0 or nr < 2 or nangle < 2:
        raise ValueError("positive parent domain and quadrature required")
    x, wx = leggauss(nr)
    cx, cw = leggauss(nangle)
    lam = (M1*M1-(M2+MPH)**2)*(M1*M1-(M2-MPH)**2)
    qstar = math.sqrt(lam)/(2*M1)
    ephstar = math.hypot(qstar,MPH)
    width = YUKAWA**2*((M1+M2)**2-MPH**2)*math.sqrt(lam)/(16*math.pi*M1**3)
    result = []
    for z, wz in zip(x, wx):
        mom = pmax*(z+1)/2
        dp = wz*pmax/2
        en1 = math.hypot(mom,M1)
        a = en1/M1*ephstar
        b = mom/M1*qstar
        e_low, e_high = a-b, a+b
        if not (e_low > 0 and e_high > e_low):
            raise ArithmeticError("nonpositive physical daughter shell")
        L = math.log(e_high/e_low)
        eph = e_low*np.exp(L*(cx+1)/2)
        wg = cw*L/2*eph/b
        assert abs(float(wg.sum())/2-1.) < 1e-9
        en2 = en1-eph
        rate = 4/(2*math.pi**2)*dp*mom**2*width*M1/en1*wg/2
        result.extend((en1,float(e2),float(ep),float(w)) for e2,ep,w in zip(en2,eph,rate))
    return np.asarray(result)


def audit(n, pmax=2., nr=24, nangle=48, order=64):
    """Initially all three species at EXACT common T=1, mu=0.

    Fermion logits = physical E1,E2, Bose logit fixed per cell to reproduce
    its exact thermal NUMBER, not energy. A physical exact-thermal state has
    identically zero weak J; this projected state need not.
    """
    edges, W, N, U, eta, ebar = cells(n,math.hypot(pmax,M1)+.2,order)
    e1,e2,ep,w = mapped_events(pmax,nr,nangle).T
    ind = np.clip(np.searchsorted(edges,ep,side='right')-1,0,n-1)
    eta_p = eta[ind]
    f1 = 1/(1+np.exp(e1))
    f2 = 1/(1+np.exp(e2))
    b = 1/np.expm1(eta_p)
    gain = f1*(1-f2)*(1+b)
    loss = (1-f1)*f2*b
    j = gain-loss
    dN = w*j
    affinity = -e1+e2+eta_p
    # Exact thermal pointwise equilibrium control, independent of cell closure.
    be_true = 1/np.expm1(ep)
    jt = f1*(1-f2)*(1+be_true)-(1-f1)*f2*be_true
    # For every reaction parent E1=E2+Ephi. Using the physical ONE-DOF
    # scalar stored energy U=Ebar*N instead leaves this residual.
    defect = float(np.sum(dN*(ebar[ind]-ep)))
    abs_flux_energy = float(np.sum(abs(dN)*ep))
    return {
        'n_cells':n, 'parent_pmax_GeV':pmax,'nrad':nr,'log_angle_nodes':nangle,
        'boson_mass_GeV':MPH,'T_common_GeV':1.,
        'scalar_thermal_number':float(np.sum(N)),
        'scalar_thermal_U_relative_error':float(np.sum(ebar*N-U)/np.sum(U)),
        'spurious_equilibrium_scalar_rate':float(dN.sum()),
        'exact_thermal_scalar_rate':float(np.sum(w*jt)),
        'spurious_entropy_rate':float(np.dot(dN,affinity)),
        'energy_closure_residual':defect,
        'residual_over_absolute_energy_throughput':defect/abs_flux_energy,
        'max_event_affinity':float(np.max(abs(affinity))),
        'cell_source_sign_monotone':True,
        'class_ready':False,
    }


def run():
    return {'schema_version':1,
            'claim':'one local scalar eta cannot simultaneously match thermal number, continuum energy and exact equilibrium',
            'input':'tree chi1->chi2+phi; exact common-T thermal fermions; scalar exact cell N and one constant eta; mapped outgoing log angle',
            'results':[audit(n) for n in (8,16,32,64,128)],
            'parent_cutoff_control':[audit(32,pmax=x,nr=48) for x in (2.,8.,16.)],
            'interpretation':'Thermal common-T is a strict fixed point for the continuous source, but not this 1-eta cell projection. Positive entropy in discrete approximation does not fix lost continuum energy. No BEC, freezeout or Omega claim.'}


if __name__ == '__main__':
    print(json.dumps(run(), indent=2, allow_nan=False))
