"""Independently factor the scalar-daughter pole residue on frozen archive steps."""
import argparse
import gzip
import json

import numpy as np
from scipy.special import roots_legendre

from early_pole_archive_scalar_below_pole_finite_q_reference import integrate
from early_pole_decay_mode_clock_reference import G, KAPPA
from early_pole_native_aperture_witness_reference import GRID, mode


def step(value, edges, density):
    k = np.searchsorted(edges, value, side='right') - 1
    return np.where((k >= 0) & (k < len(density)),
                    density[np.clip(k, 0, len(density) - 1)], 0.)


def residue(a, q, order, exposures=None):
    M, m = .05, .003
    ms = a['checkpoints'][-1]['scalar_pole_over_T']
    H = np.sqrt(np.pi**2*106.75/90)*(86./a['x_end'])/2.435e18
    vertex = (M + m)**2 - ms*ms
    pn = np.asarray(a['endpoint_parent_momentum_edges'])
    pf = np.asarray(a['endpoint_parent_occupations'])
    w0 = np.hypot(q, ms)
    a0 = (M*M - m*m - ms*ms)/2
    mr = np.sqrt((a0 - m*ms)*(a0 + m*ms))
    star_q = mr/M
    star_e = np.hypot(star_q, ms)
    residue_spin = (((M + m)**2 - ms*ms)/(M*M))**2
    P = np.asarray(a['endpoint_fermion_momentum_edges'])
    F = np.asarray(a['endpoint_fermion_occupations'])
    E = np.hypot(m, P)
    Q = np.asarray(a['scalar_momentum_edges'])
    B = np.maximum(0., np.asarray(a['endpoint_scalar_occupations']))
    cut = 1e-5
    kcut = cut*cut/(np.hypot(cut, ms) + ms)
    K = Q*Q/(np.hypot(Q, ms) + ms)
    tagged = float(step(q, Q, B)) if q else B[0]
    pmin, pmax = np.hypot(m, [.45, 5.])
    if q:
        a, b = (w0*a0 + np.array([-1., 1.])*q*mr)/(ms*ms)
        a, b = max(a, pmin), min(b, pmax)
        assert b > a
        cuts = np.unique(np.r_[a, E[(a < E) & (E < b)], b])
        if exposures is not None:
            parent_edges = np.hypot(pn,M)-w0
            cuts = np.unique(np.r_[cuts,parent_edges[(a<parent_edges)&(parent_edges<b)]])
        h = np.diff(cuts)[:, None]/2
        x, y = roots_legendre(order)
        es = ((cuts[:-1] + cuts[1:])[:, None]/2 + h*x).ravel()
        es_weight = (h*y).ravel()
    else:
        a = b = a0/ms
        assert pmin < a < pmax
        es, es_weight = np.array([a]), np.array([1.])
    out = np.zeros((2, 3))
    eliminated = np.zeros((2, 3))
    stationary_sources = np.zeros(3)
    exact_sources = np.zeros(3)
    geom = []
    fac = H*(16*np.pi*KAPPA*M)**2
    gamma_vac = H*KAPPA*M*(2*star_q/M)*(vertex/(M*M))
    source_prefactor = (G*KAPPA*M*vertex/(q*w0) if q
                        else G*KAPPA*M*vertex*2*mr/ms**3)
    max_thermal = 0.
    if exposures is not None:
        instantaneous=np.zeros((3,len(exposures)))
        integrated=np.zeros_like(instantaneous)
        old_est=np.zeros(3)
    for e, weight in zip(es, es_weight, strict=True):
        et = e + w0
        lab = np.sqrt((et - M)*(et + M))
        speed, boost = lab/et, et/M
        w_min = boost*(star_e - speed*star_q)
        klow = (boost*(star_q - speed*star_e))**2/(w_min + ms)
        khigh = klow + 2*lab*star_q/M
        lo = max(klow, et - ms - pmax)
        hi = min(khigh, et - ms - pmin)
        if hi <= lo:
            continue
        edge = np.r_[K, kcut, et - ms - E]
        edge = np.unique(np.r_[lo, edge[(lo < edge) & (edge < hi)], hi])
        k = (edge[1:] + edge[:-1])/2
        y = np.diff(edge)
        fin = float(step(e, E, F))
        fout, bout = step(et - ms - k, E, F), step(k, K, B)
        D, I = (1-fout)*(1+bout), fout*bout
        loss = fin*D
        gain = (1-fin)*(1+1/tagged)*I
        net = bout*(fout-fin) + bout*fout*(1-fin)/tagged - fin*(1-fout)
        pref = (M*M/(64*np.pi**3*q*w0*lab) if q
                else M*M/(32*np.pi**3*ms*ms))
        width_k = khigh-klow
        damping = float(y @ (1-fout+bout))/width_k
        fill = float(y @ I)/(width_k*damping)
        assert damping>0. and 0.<=fill<=1.
        assert abs(float(y @ (D+I))/width_k-damping)<1e-10*damping
        Pfin = float(step(lab,pn,pf))
        geom.append([damping,fill,Pfin,gamma_vac*damping,gamma_vac*damping*M*M/(H*et)])
        # Compare the separately eliminated parent source and the actual
        # archived parent using a numerically stable difference.
        def direct(P):
            return np.array([P*(1-fin)*(1+1/tagged),(1-P)*fin,
                             P-fin+P*(1-fin)/tagged])
        stationary_sources += source_prefactor*weight*direct(fill)
        if not q:
            exact_sources += source_prefactor*weight*direct(Pfin)
        if exposures is not None:
            labrate=gamma_vac*damping*M*M/(H*et)
            endpoint_fraction=-np.expm1(-labrate*exposures)
            avg_fraction=mean_relaxation(labrate*exposures)
            da=source_prefactor*weight*direct(Pfin)
            df=source_prefactor*weight*(fill-Pfin)*np.array([ (1-fin)*(1+1/tagged),-fin,1+(1-fin)/tagged])
            old_est+=da
            instantaneous+=da[:,None]+df[:,None]*endpoint_fraction
            integrated+=(da[:,None]+df[:,None]*avg_fraction)*exposures
        # Independent thermal local check: common-T chemical equilibrium
        # implies this stationary parent occupation at each lab energy.
        be = 1/np.expm1(ms+k)
        fd = 1/(1+np.exp(et-ms-k))
        therm = (fd*be)/(1-fd+be)
        max_thermal = max(max_thermal,float(np.max(abs(therm-1/(1+np.exp(et))))))
        for j, panel in enumerate((k < kcut, k >= kcut)):
            packet = (pref*weight*residue_spin)*(np.array([gain, loss, net])[:,panel] @ y[panel])
            out[j] += packet
            eliminated[j] += fac*np.pi/(gamma_vac*damping)*packet
    # Entrywise signs are the signed net after division by the tagged count.
    err = np.max(abs(eliminated.sum(axis=0)-stationary_sources)/np.maximum(abs(stationary_sources),1e-50))
    if err>7e-7:
        raise ValueError(f'Matched occupation-weighted elimination mismatch: {err:.3e}')
    assert max_thermal<1e-15
    if q:
        # Archived parent populations add their own lab-energy edges, not
        # present in the stationary elimination or pole scattering.
        pe=np.hypot(pn,M)-w0
        ec=np.unique(np.r_[a,b,E[(a<E)&(E<b)],pe[(a<pe)&(pe<b)]])
        em=ec[:-1]+np.diff(ec)/2
        ff,pp=step(em,E,F),step(np.sqrt((em+w0-M)*(em+w0+M)),pn,pf)
        vv=np.array([pp*(1-ff)*(1+1/tagged),(1-pp)*ff,
                     pp-ff+pp*(1-ff)/tagged])
        exact_sources=source_prefactor*(vv@np.diff(ec))
        _, moment, _ = mode(M,ms,m,q,pn,pf,pf,P,F,0.)
        native_check=source_prefactor*(moment[0]-moment[1]+moment[2]/tagged)
        assert abs(native_check-exact_sources[-1])<1e-8*abs(native_check)
    g=np.array([np.min(geom,axis=0),np.max(geom,axis=0)])
    min_et=a+w0
    min_lab=np.sqrt((min_et-M)*(min_et+M))
    damping_bound=1.+float(np.diff(K)@B)/(2*min_lab*star_q/M)
    relaxation_bound=gamma_vac*damping_bound*M*M/(H*min_et)
    assert g[1,4]<=relaxation_bound
    spectral=dict(occupied_damping=float(g[0,0]),
                  max_occupied_damping=float(g[1,0]),
                  parent_stationary_and_archived_range=g[:,1:3].tolist(),
                  dimensionless_width_in_s_over_M2_range=g[:,3].tolist(),
                  lab_relaxation_over_reference_H_range=g[:,4].tolist(),
                  lab_relaxation_over_reference_H_upper_bound=float(relaxation_bound),
                  real_pole_weight_model_over_H=eliminated.tolist(),
                  stationary_parent_direct_over_H=stationary_sources.tolist(),
                  archived_parent_direct_over_H=exact_sources.tolist(),
                  stationary_matching_error=float(err),maximum_thermal_equilibrium_error=max_thermal)
    if exposures is not None:
        np.testing.assert_allclose(old_est, exact_sources, rtol=3e-13, atol=1e-12)
        return dict(tagged_q=q,order=order,exposure=exposures.tolist(),
                    exact_archived_initial_gross_and_net_rate=exact_sources.tolist(),
                    stationary_direct_gross_and_net_rate=stationary_sources.tolist(),
                    instant_gross_and_net_rate=instantaneous.T.tolist(),
                    integrated_gross_and_net_gain=integrated.T.tolist(),
                    maximum_frozen_parent_relaxation=float(relaxation_bound))
    return out, [a, b], residue_spin, spectral


def mean_relaxation(z):
    # 1 - (1 - exp(-z))/z, evaluated without cancellation at small z.
    out=np.empty_like(z)
    lo=z<2e-3
    y=z[lo]
    out[lo]=y*(.5+y*(-1/6+y*(1/24+y*(-1/120+y/720))))
    out[~lo]=1.+np.expm1(-z[~lo])/z[~lo]
    return out


def transient_fixture():
    # Cross-check the decay exponential and accumulated fill against an ODE
    # on synthetic detailed-balance ion/return coefficients.
    from scipy.integrate import solve_ivp
    for rate in (.02, 3., 170., 2200.):
        for p0, peq in ((.11,.63),(.77,.15),(.44,.44)):
            t=np.array([0.,2e-8,7.2e-5,.002,.02])
            y=solve_ivp(lambda _, y:[rate*(peq-y[0]),y[0]], (0.,t[-1]), [p0,0.],
                        t_eval=t, method='Radau', atol=4e-14, rtol=4e-13)
            assert y.success
            v = p0+(peq-p0)*(-np.expm1(-rate*t))
            g = t*(p0+(peq-p0)*mean_relaxation(rate*t))
            np.testing.assert_allclose(y.y,[v,g],rtol=5e-11,atol=4e-14)


def compare(a, q, order):
    d, endpoints, spin, spectral = residue(a, q, order)
    r = integrate(a, q, order, .90, probe_x=0.)
    H = np.sqrt(np.pi**2*106.75/90)*(86./a['x_end'])/2.435e18
    fac = H*(16*np.pi*KAPPA*.05)**2
    ours = fac*d
    ref = np.c_[r['gross_density_per_reference_H'], r['net_density_per_reference_H']]
    err = abs(ours-ref)/np.maximum(abs(ref), 1e-50)
    assert np.max(err)<8e-10
    return dict(q=q,order=order,incoming_fermion_energy_interval=endpoints,squared_rescaled_spin_residue=spin,
                factorized_gross_and_net_density_over_H=ours.tolist(),
                reference_gross_and_net_density_over_H=ref.tolist(),relative_difference=err.tolist(),
                quasistatic_same_channels=spectral)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid', choices=('fine', 'third'))
    p.add_argument('--order', type=int, default=14)
    p.add_argument('--q', type=float, nargs='+', default=[0., 4.1e-6])
    p.add_argument('--exposure', type=float, nargs='+', default=None,
                   help='Check analytic single-channel parent time evolution in a frozen bath.')
    args = p.parse_args()
    with gzip.open(GRID/f'early-pole-spectrum-{args.grid}-m06-2026-11-26.json.gz', 'rt') as fh:
        data = json.load(fh)
    if args.exposure is not None:
        transient_fixture()
        t=np.asarray(args.exposure)
        if np.any(~np.isfinite(t)) or np.any(t<0.):
            p.error("Use finite nonnegative exposures.")
        cases=[residue(data,q,args.order,t) for q in args.q]
    else:
        cases=[compare(data,q,args.order) for q in args.q]
    print(json.dumps(dict(grid=args.grid,cases=cases),indent=2,allow_nan=False))
