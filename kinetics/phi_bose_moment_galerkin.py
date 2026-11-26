"""Physical p-space Bose/Fermi entropy-dual Galerkin moment pilot.

A mathematically conservative *instantaneous* basis for the 1<->2 Yukawa
chi1 <-> chi2 + real phi event; not a physical kinetic time solver.

Instead of assigning BE(E=mphi) to a whole first positive-volume cell,
interpolate the inverse fugacity eta(E) in piecewise-linear energy hats:
  b(E)=1/expm1(eta(E)), f(E)=1/(exp(eta(E))+1).
The state moment per hat is M_i = g/(2pi^2) int p^2 dp h_i(E) n(E).
Since sum h_i=1 and sum E_i h_i(E)=E, the total number and energy
are exactly sum_i M_i and sum_i E_i*M_i respectively. Thermal eta(E)
is exactly affine, so its physical integrated number does not acquire
a spurious low-mass (mphi << first cell) rest-node condensate.

The same hats deposit one on-shell Yukawa collision extent into
M_chi1, M_chi2, M_phi; the integrated entropy gradient is exactly eta_i
(fermion AND boson), so Sdot = W*(F-R)*ln(F/R) >= 0 with no clipping.

This is not a finite-step Pauli/Bose-bound-preserving solver: to evolve,
update the integrated moments and invert their nonlinear map to eta.
A true BEC delta(p) degree of freedom is separate and absent here.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import expit


@dataclass(frozen=True)
class Mode:
    mass: float
    nodes: np.ndarray
    degeneracy: float = 1.0
    fermion: bool = False

    def __post_init__(self):
        e = np.asarray(self.nodes, dtype=float)
        if (not math.isfinite(self.mass) or self.mass < 0 or
                not math.isfinite(self.degeneracy) or self.degeneracy <= 0 or
                e.ndim != 1 or len(e) < 3 or
                not np.all(np.isfinite(e)) or
                abs(e[0] - self.mass) > 1e-10 or
                not np.all(np.diff(e) > 0)):
            raise ValueError('physical, ordered E nodes beginning at rest mass required')
        object.__setattr__(self, 'nodes', e)


def hats(mode: Mode, energy):
    """Positive first-order hats: partition unity + reproduce physical E.

    A scalar energy is required; collisional arrays are passed one by one.
    """
    E = float(energy)
    z = mode.nodes
    if not math.isfinite(E) or E < z[0] - 2e-12 or E > z[-1] + 2e-12:
        raise ValueError('on-shell energy outside mode grid')
    E = min(max(E, z[0]), z[-1])
    j = min(int(np.searchsorted(z, E, side='right')) - 1, len(z) - 2)
    j = max(j, 0)
    t = (E - z[j])/(z[j+1] - z[j])
    return j, 1-t, t


def eta_at(mode: Mode, eta: np.ndarray, E: float):
    j, a, b = hats(mode, E)
    return float(a*eta[j] + b*eta[j+1])


def occupation_from_eta(eta, fermion=False):
    eta = np.asarray(eta)
    if fermion:
        return expit(-eta)
    # eta=0 at the single p=0 massless node is allowed; it has zero
    # continuum measure. No quadrature point lies exactly at p=0.
    with np.errstate(divide='ignore', invalid='raise', over='raise'):
        return 1/np.expm1(eta)


def physical_moments(mode: Mode, eta: np.ndarray, order=32):
    """Hat-integrated N_i, actual energy, quantum entropy, Jacobian.

    Integrate p^2 dp on each E-hat interval in physical p, avoiding the
    massive E-Jacobian threshold singularity. The symmetric negative
    matrix J_ij=dM_i/deta_j is an independent entropy-dual control.
    """
    z = mode.nodes
    lam = np.asarray(eta, dtype=float)
    if lam.shape != z.shape or not np.all(np.isfinite(lam)) or order < 6:
        raise ValueError('finite inverse fugacity on same mode grid')
    if not mode.fermion and (np.any(lam < 0) or
                             np.all(lam == 0)):
        raise ValueError('Bose inverse fugacity must be nonnegative and nontrivial')
    x,w = leggauss(order)
    M = np.zeros(len(z))
    J = np.zeros((len(z),len(z)))
    S = 0.
    pref = mode.degeneracy/(2*math.pi**2)
    for i in range(len(z)-1):
        plo = math.sqrt(max(0,z[i]**2-mode.mass**2))
        phi = math.sqrt(max(0,z[i+1]**2-mode.mass**2))
        p = plo + (x+1)*(phi-plo)/2
        E = np.hypot(p,mode.mass)
        a = (z[i+1]-E)/(z[i+1]-z[i])
        b = 1-a
        e = a*lam[i] + b*lam[i+1]
        if not mode.fermion and np.any(e <= 0):
            raise ValueError('Bose inverse fugacity nonpositive inside physical cell')
        f = occupation_from_eta(e, mode.fermion)
        dW = pref*(phi-plo)/2*w*p*p
        M[i] += np.dot(dW,a*f)
        M[i+1] += np.dot(dW,b*f)
        # Fermi/Bose physical entropy; stable even for eta << 1.
        if mode.fermion:
            s = np.logaddexp(0,-e) + e*f
            var = f*(1-f)
        else:
            s = -np.log(-np.expm1(-e)) + e*f
            var = f*(1+f)
        S += float(np.dot(dW,s))
        J[i,i] -= np.dot(dW,a*a*var)
        J[i,i+1] -= np.dot(dW,a*b*var)
        J[i+1,i] -= np.dot(dW,a*b*var)
        J[i+1,i+1] -= np.dot(dW,b*b*var)
    return {'M':M,'N':float(M.sum()),'U':float(np.dot(z,M)),
            'S':float(S),'dM_deta':J,
            'cond_number':float(np.linalg.cond(J))}


def yukawa_kinematics(p,c,m1,m2,mphi):
    if m1 <= m2+mphi or p < 0 or not -1 <= c <= 1:
        raise ValueError('on-shell open 1->2+phi required')
    lam = (m1*m1-(m2+mphi)**2)*(m1*m1-(m2-mphi)**2)
    e1 = math.hypot(p,m1)
    e2star = (m1*m1+m2*m2-mphi*mphi)/(2*m1)
    ephstar = (m1*m1-m2*m2+mphi*mphi)/(2*m1)
    pstar = math.sqrt(lam)/(2*m1)
    e2 = e1/m1*e2star + p/m1*pstar*c
    eph = e1/m1*ephstar - p/m1*pstar*c
    return e1,e2,eph


def yukawa_moment_source(modes, eta, events):
    """Accumulate synthetic positive-measure on-shell 1<->2 events.

    modes=(chi1 Fermi,chi2 Fermi,phi Bose); eta is inverse fugacity
    arrays of each species. Events=(mother momentum,CM cosangle,weight).
    Weights are nonnegative *chosen quadrature weights*, not a claim of
    cosmological Yukawa rate, no source time integrated here.
    """
    if len(modes)!=3 or len(eta)!=3 or (not modes[0].fermion or
            not modes[1].fermion or modes[2].fermion):
        raise ValueError('chi1,chi2 Fermi and phi Bose mode order required')
    dM = [np.zeros(len(mode.nodes)) for mode in modes]
    entropy_exact = 0.
    max_shell_defect = 0.
    max_detailed_balance = 0.
    rate = 0.
    for p,c,W in events:
        if not(math.isfinite(W) and W>=0):
            raise ValueError('positive event measure required')
        energies = yukawa_kinematics(p,c,*[mode.mass for mode in modes])
        e1,e2,eph = energies
        vals = [eta_at(mod, et, e) for mod,et,e in zip(modes,eta,energies)]
        a,d,b = (float(occupation_from_eta(vals[0],True)),
                 float(occupation_from_eta(vals[1],True)),
                 float(occupation_from_eta(vals[2],False)))
        F = a*(1-d)*(1+b)
        R = (1-a)*d*b
        flux = F-R
        affinity = -vals[0]+vals[1]+vals[2]
        if abs(flux*affinity) < 1e-25:
            # May be slightly negative only from roundoff, measured below.
            pass
        max_detailed_balance=max(max_detailed_balance,abs(flux))
        max_shell_defect=max(max_shell_defect,abs(e1-e2-eph))
        amount = W*flux
        rate += amount
        entropy_exact += amount*affinity
        for species, sign, en in ((0,-1,e1),(1,1,e2),(2,1,eph)):
            mod = modes[species]
            j,u,v = hats(mod,en)
            dM[species][j] += sign*amount*u
            dM[species][j+1] += sign*amount*v
    dN = [float(x.sum()) for x in dM]
    dU = [float(np.dot(mod.nodes,x)) for mod,x in zip(modes,dM)]
    ent = sum(float(np.dot(np.asarray(lam),x)) for lam,x in zip(eta,dM))
    return {'dM':dM,'dN':dN,'dU':dU,'dark_fermion_charge':dN[0]+dN[1],
            'parent_plus_phi_charge':dN[0]+dN[2],
            'total_energy_residual':sum(dU),
            'entropy_from_dual_gradient':ent,
            'entropy_from_event':entropy_exact,
            'max_shell_defect':max_shell_defect,'rate_phi':rate,
            'max_micro_flux_abs':max_detailed_balance,
            'not_a_physical_cosmological_rate':True}


def thermal_bose_cell_average_at_rest(mass, upper, bins, temperature=1.):
    """Independent physical FIRST cell average, diagnostic only.

    Replacing the rest-node occupation with this average does NOT make a
    microreversible node-based scheme: see audit() detailed-balance gate.
    """
    if not (mass > 0 and upper > mass and bins >= 4 and temperature > 0):
        raise ValueError('physical positive massive thermal test')
    from scipy.integrate import quad
    h = (upper-mass)/bins
    Ehi = mass+h/2
    pmax = math.sqrt(Ehi**2-mass**2)
    val = quad(lambda p:p*p/math.expm1(math.hypot(mass,p)/temperature),
               0,pmax,epsabs=2e-13,epsrel=2e-11,limit=250)[0]
    return 3*val/pmax**3


def audit(*, mass=86e-9, T=1., upper=12.1, bins=(16,64,256)):
    """Self-contained compact numerical artifact, no physical DM model."""
    from scipy.integrate import quad
    from scipy.special import zeta
    pmax = math.sqrt(upper*upper-mass*mass)
    reference = quad(lambda p:p*p/math.expm1(math.hypot(mass,p)/T),
                     0,pmax,epsabs=2e-13,epsrel=2e-11)[0]/(2*math.pi**2)
    rows=[]
    for n in bins:
        nodes = np.linspace(mass,upper,n+1)
        mode=Mode(mass,nodes)
        phi=physical_moments(mode,nodes/T,order=36)
        h = (upper-mass)/n
        # Explicitly reproduce flawed original positive-volume rest-node rule.
        left=np.maximum(nodes-h/2,mass)
        right=np.minimum(nodes+h/2,upper)
        w=(np.maximum(right*right-mass*mass,0)**1.5-
           np.maximum(left*left-mass*mass,0)**1.5)/(6*math.pi**2)
        bad=float(np.dot(w,1/np.expm1(nodes/T)))
        rows.append({'bins':n,'Galerkin_number':phi['N'],
                     'physical_reference_number':reference,
                     'Galerkin_relative_error':phi['N']/reference-1,
                     'bad_rest_node_number':bad,
                     'bad_rest_node_relative_factor':bad/reference,
                     'Galerkin_energy':phi['U'],
                     'Galerkin_entropy':phi['S'],
                     'first_dual_moment':float(phi['M'][0]),
                     'first_cell_thermal_average':thermal_bose_cell_average_at_rest(
                         mass,upper,n,T)})
    m1,m2=.1,.02
    modes=[Mode(m,np.linspace(m,upper,17),fermion=i<2)
           for i,m in enumerate((m1,m2,mass))]
    eq=[md.nodes/T for md in modes]
    events=[(.0,0.,1.),(.03,-.7,1.),(.03,.7,1.),(.5,.25,1.)]
    exact=yukawa_moment_source(modes,eq,events)
    compromised=[x.copy() for x in eq]
    bbar=thermal_bose_cell_average_at_rest(mass,upper,16,T)
    compromised[2][0]=math.log1p(1/bbar)
    faulty=yukawa_moment_source(modes,compromised,events)
    massless=Mode(0.,np.linspace(0,50,65))
    massless_n=physical_moments(massless,massless.nodes/T,order=48)['N']
    return {'schema_version':1,'scope':'synthetic equilibrium IR and instantaneous Yukawa on-shell moment basis; NOT cosmological 86eV dark matter or condensate',
            'mass_unit':'GeV if mass=86e-9, T and upper illustrative',
            'mphi':mass,'T':T,'Eupper':upper,'rows':rows,
            'm0_number_from_moment_basis':massless_n,
            'massless_zeta3_T3_pi2':float(zeta(3)/(math.pi**2)*T**3),
            'consistent_thermal_Rphi':exact['rate_phi'],
            'consistent_thermal_entropy':exact['entropy_from_event'],
            'cell_average_substituted_Rphi':faulty['rate_phi'],
            'cell_average_substituted_entropy':faulty['entropy_from_event'],
            'mathematical_conservation':{'sum_hats':'1','sum_Ehats':'physical E',
                'dS_dM':'eta_i','dM_deta':'symmetric negative definite'},
            'no_time_solver_or_zero_mode':True,'class_ready':False}


if __name__=='__main__':
    import json
    print(json.dumps(audit(),indent=2))
