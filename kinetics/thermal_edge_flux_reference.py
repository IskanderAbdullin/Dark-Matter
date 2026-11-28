"""Boundary flux diagnostic for isotropic F -> F + B with thermal occupations."""
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import expit
from scipy.integrate import quad, quad_vec
from massive_shell_moments import shell_weights


def grid(m, qfirst, qmax):
    # Interior cells are not used: only the means of the two edge cells.
    q = [0., qfirst, qmax/1.1, qmax]
    n, e, _, _ = shell_weights(q, m)
    return np.array([m, e[0]/n[0], e[-1]/n[-1], np.hypot(qmax, m)])


def integral(M, temp, qfirst, qmax, order, pmax, chemical_shift=0.):
    # All energies are in common units. The thermal tests use T=1.
    m, b = .08*M, .18*M
    u = (M*M+m*m-b*b)/(2*M)
    v = M-u
    k = np.sqrt((M-m-b)*(M+m+b)*(M-m+b)*(M+m-b))/(2*M)
    bos, fer = grid(m,qfirst,qmax), grid(b,qfirst,qmax)
    nodes, weights = leggauss(order)
    # Segment the parent momentum at kinematic crossings with every
    # energy cut (for either daughter). The angular integrand is then smooth
    # inside each interval, including for arbitrarily narrow edge regions.
    breaks = [0., pmax, .025*temp, .1*temp, .4*temp, temp, 2*temp, 4*temp,
              8*temp, 16*temp, 32*temp]
    for cuts, own, star in [(bos,m,u), (fer,b,v)]:
        for x in cuts:
            h = np.sqrt(max((x-own)*(x+own),0))
            for root in (M*(star*star+h*h)/(star*x+k*h),
                         M*(star*x+k*h)/(own*own)):
                if root > M:
                    breaks.append(np.sqrt((root-M)*(root+M)))
    bounds = np.unique(np.clip(breaks,0,pmax))
    lo, hi = bounds[:-1,None], bounds[1:,None]
    p = (lo+(hi-lo)/2*(nodes+1)).ravel()
    wp = ((hi-lo)/2*weights).ravel()
    E = np.hypot(p,M)
    emin = (m*m*E*E+k*k*M*M)/(M*(u*E+k*p))
    emax = (u*E+k*p)/M
    z = np.stack([emin,emax,*[np.full_like(E,w) for w in bos[1:]],*[E-w for w in fer[1:]]],axis=1)
    z = np.sort(np.clip(z,emin[:,None],emax[:,None]),axis=1)
    zl, zh = z[:,:-1,None],z[:,1:,None]
    # Log energy integration is stable even when the Bose factor varies rapidly.
    dt = np.log1p((zh-zl)/zl)
    phi = zl*np.exp(dt*(nodes+1)/2)
    ferm = E[:,None,None]-phi
    jac = phi*dt*weights/2
    norm = wp*p*p*M/E/(2*k*p/M)
    P = expit(-(E[:,None,None]-chemical_shift)/temp)
    F = expit(-ferm/temp)
    B = 1/np.expm1(phi/temp)
    gain = P*(1-F)*(1+B)
    loss = (1-P)*F*B
    measure = norm[:,None,None]*jac
    stats = np.array([np.sum(measure*gain),np.sum(measure*loss)])
    residual = np.sum(measure*abs(gain-loss))/sum(stats)
    balance = np.sum(measure*abs(gain-loss*np.exp(chemical_shift/temp)))/sum(stats)
    def fractions(x, bounds):
        low, high, overflow = x<bounds[1],(x>bounds[2])&(x<=bounds[3]),x>bounds[3]
        return [float(np.sum(measure*(gain+loss)*mask)/sum(stats))
                for mask in (low,high,overflow)]
    return dict(gain=stats[0], loss=stats[1], absolute_net_fraction=float(residual),
                detailed_balance_identity_residual=float(balance),
                phi_fractions_low_high_overflow=fractions(phi,bos),
                fermion_fractions_low_high_overflow=fractions(ferm,fer))


def independent_low_edge():
    # Direct angle integration on one difficult early-onset case, using
    # adaptive quadrature rather than the segmented log-energy rule.
    M, qfirst = .005, .00625
    m, b = .08*M, .18*M
    v = (M*M+m*m-b*b)/(2*M)
    k = np.sqrt(v*v-m*m)
    boundary = grid(m,qfirst,24.)[1]
    def outer(p):
        e = np.hypot(p,M)
        mid, span = v*e/M, k*p/M
        def fun(c):
            phi = mid+span*c
            return p*p*M/e*expit(-e)*expit(e-phi)/(-np.expm1(-phi))/2
        cut = np.clip((boundary-mid)/span,-1.,1.) if span else -1.
        return np.array([quad(fun,-1.,1.,epsabs=1e-11,epsrel=1e-9)[0],
                         quad(fun,-1.,cut,epsabs=1e-11,epsrel=1e-9)[0]])
    val, err = quad_vec(outer,0.,60.,epsabs=1e-10,epsrel=1e-8,
                        points=[.01,.05,.1,.5,1,2,4,8,16,32])
    return {"gain":float(val[0]), "low_fraction":float(val[1]/val[0]),
            "scipy_outer_error_estimate":float(err)}


def run():
    report={}
    for M in (1.,.1,.005):
        main=[]
        scan = [(.1,6.),(.025,12.),(.00625,24.)]
        if M == .005:
            scan += [(.0015625,24.),(.000390625,24.),(.00009765625,24.)]
        for a,b in scan:
            r=integral(M,1.,a,b,48,80.)
            t=integral(M,1.,a,b,96,100.)
            main.append(dict(first_q=a,last_q=b,**t,
                             max_relative_total_flux_drift=max(abs(r[s]/t[s]-1) for s in ("gain","loss")),
                             max_absolute_boundary_fraction_drift=max(abs(x-y) for key in
                                ("phi_fractions_low_high_overflow","fermion_fractions_low_high_overflow")
                                for x,y in zip(r[key],t[key]))))
        n=integral(M,1.,.025,12.,96,100.,.2)
        assert max(row["max_relative_total_flux_drift"] for row in main) < 1e-8
        assert max(row["max_absolute_boundary_fraction_drift"] for row in main) < 1e-8
        assert max(row["absolute_net_fraction"] for row in main) < 1e-13
        assert n["detailed_balance_identity_residual"] < 1e-13
        assert abs(n["gain"]/n["loss"]-np.exp(.2)) < 1e-12
        lows = [row["phi_fractions_low_high_overflow"][0] for row in main]
        assert all(x>y for x,y in zip(lows,lows[1:]))
        report[str(M)]={"equilibrium":main,"parent_chemical_shift":.2,"out_of_equilibrium":n}
    check = independent_low_edge()
    baseline = report["0.005"]["equilibrium"][2]
    assert abs(check["gain"]/baseline["gain"]-1) < 2e-8
    expected = baseline["phi_fractions_low_high_overflow"][0]
    assert abs(check["low_fraction"]-expected) < 2e-8
    report["independent_angle_check"] = check
    return report

if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=False))
