"""Radial principal-value diagnostics with archived bare occupations."""
import gzip
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from early_pole_decay_mode_clock_reference import intervals,Q
from finite_mass_qed_charged_overlap_asymptotic_reference import solve_axes
from finite_mass_qed_thermal_charged_pv_geometry_reference import chart

D = np.longdouble


def angular_kernel(M, m, ms, external_parent_energy, species):
    """Exact azimuth average followed by the polar Cauchy principal value."""
    M, m, ms, E, q = map(D, (M, m, ms, external_parent_energy, Q))
    Es = np.hypot(ms, q)
    z = (E*Es - (M*M + ms*ms - m*m)/2)/q
    ma, w, h, sign = (M, E, z, 1) if species == 'parent' else (m, E-Es, z-q, -1)
    p = np.sqrt((w-ma)*(w+ma))
    c = h/p
    if not -1 < c < 1:
        raise ValueError('This reference requires an interior external Born event.')
    d = (M*M + m*m - ms*ms)/2
    C, K = (M+m)**2/(d+M*m), 4-2*(M+m)*ma/(d+M*m)

    def kernel(delta):
        delta = D(delta)
        t = w+delta
        if delta == 0 or t <= ma:
            raise ValueError('Evaluate only away from the coincidence and at positive charged momentum.')
        r = np.sqrt((t-ma)*(t+ma))
        R = p*r
        # b=(delta_p^2-delta_E^2)/2; factor the difference of squares.
        v = (t+w)/(r+p)
        vminus1 = ma**2*(1/(t+r)+1/(w+p))/(r+p)
        b = delta**2*vminus1*(v+1)/2
        mid = c*(1+b/R)
        width = np.sqrt((1-c)*(1+c)*b*(2*R+b))/R
        # Keep the difference u-mid directly, rather than subtracting angles.
        um = delta*(Es/q-c*v)/r-c*b/R
        u = mid+um
        if abs(u) == 1:
            raise ValueError('Undamped polar endpoint; split the radial integral here.')
        a, l, h = np.arcsinh(np.array([um,-1-mid,1-mid], dtype=D)/width)
        def primitive(y):
            return a+np.log(abs(np.expm1(y-a)))-np.logaddexp(y,-a)
        # Valid with a principal value for |u|<1, ordinary otherwise.
        H = (primitive(h)-primitive(l))/np.hypot(um,width)
        photon = -(h-l)/(2*R)
        phase = np.log(abs((1-u)/(1+u)))
        return (C*phase-2*d*H/R)/(2*sign*q*r)+K*photon/2

    return kernel


def independent_angular_checks(M,m,ms,E,cut,lo,hi):
    fast = angular_kernel(M,m,ms,E,cut)
    direct = chart(M,m,ms,E,cut)
    extracted = solve_axes(M,m,ms,E,cut)[0]
    resolved = [lo+t*(hi-lo) for t in (-.08,.08,.24,.76,.92,1.08)]
    away = max(float(abs(fast(D(z)-D(E))/direct(z)['principal_value_coefficient']-1))
               for z in resolved)
    near = max(float(abs(fast(t)/extracted(t)[0]-1))
               for t in (-1e-5,1e-5,-1e-8,1e-8))
    if max(away,near)>3e-8:
        raise ArithmeticError('Exact angular primitive did not match independent angular integrals.')
    return dict(resolved_and_outside_aperture_relative_defect=away,overlap_relative_defect=near)



def analyze(d,E,cut,is_edge):
    M=d['x_end']
    m=d['mass_ratio']*M
    ms=d['checkpoints'][-1]['scalar_pole_over_T']
    Es=np.hypot(ms,Q)
    ma=M if cut=='parent' else m
    e0=E if cut=='parent' else E-Es
    name='parent' if cut=='parent' else 'fermion'
    edges=np.hypot(d[f'endpoint_{name}_momentum_edges'],ma)
    occ=np.asarray(d[f'endpoint_{name}_occupations'])
    L=3e-4
    def spectrum(e):
        return occ[np.searchsorted(edges,e)-1]
    _,C,ell=solve_axes(M,m,ms,E,cut)
    angular=angular_kernel(M,m,ms,E,cut)
    P=np.sqrt((e0-ma)*(e0+ma))
    # Use a resolved local probe of the two sides to avoid assigning a
    # value exactly at an archived discontinuity.
    jump=spectrum(e0+2e-9)-spectrum(e0-2e-9)
    if bool(jump)!=is_edge:
        raise ValueError('Expected frozen step behavior did not match the chosen event.')
    log_coeff=float(C*P*jump/2)
    def vals(t):
        q=abs(t)
        a,b=angular(q),angular(-q)
        def W(x):
            l=e0+x
            p=np.sqrt((l-ma)*(l+ma))
            return p*spectrum(l)/2.
        p,n=W(q),W(-q)
        raw=(p-n)*(a-b)/2+(p+n)*(a+b)/2
        # Independent reconstruction from the two individually weighted terms.
        if abs(raw-(p*a+n*b))>2e-7*max(1.,abs(raw)):
            raise ValueError('Symmetric radial pairing did not reproduce the directly weighted values.')
        return np.array([raw*q,raw*q-log_coeff])
    lows=(1e-5,1e-6,1e-7,1e-8,1e-9,1e-10)
    def all_values(n):
        t,w=roots_legendre(n)
        cuts=sorted({*lows,L,1e-4,*[x for x in abs(edges-e0) if min(lows)<x<L]})
        v=np.zeros(2)
        out={}
        for a,b in reversed(list(zip(cuts[:-1],cuts[1:],strict=True))):
            l,r=np.log(a),np.log(b)
            q=np.exp((l+r)/2+(r-l)*t/2)
            z=np.array([vals(x) for x in q])
            v+=(r-l)/2*(w@z)
            if a in lows:
                if abs(v[0]-v[1]-log_coeff*np.log(L/a))>6e-12:
                    raise ValueError('Sampled radial logarithm does not match the extracted jump.')
                out[a]=v.copy()
        return out
    rows=[]
    X,Y=all_values(12),all_values(20)
    for low in lows:
        x,y=X[low],Y[low]
        defect=max(abs(x-y))
        if defect>4e-6:
            raise ValueError('The radial angular-PV remainder did not numerically resolve.')
        rows.append(dict(radial_lower_cut=low,subtracted_radial_value=float(y[1]),raw_radial_value=float(y[0]),
                         resolution_error=float(defect)))
    slope=(rows[-1]['raw_radial_value']-rows[-2]['raw_radial_value'])/np.log(10.)
    relative=abs(slope/log_coeff-1) if log_coeff else None
    if is_edge and relative>2e-4:
        raise ArithmeticError(f'Archive jump prediction failed: observed {slope}, analytic {log_coeff}.')
    remainders=np.asarray([v['subtracted_radial_value'] for v in rows])
    observed=np.abs(np.diff(remainders))
    # Odd, unweighted angular singularity is known analytically; after
    # removing the occupation jump the residual at small cutoff may have a
    # log/linear tail, but its last-decade change must shrink.
    if observed[-1]>max(1e-7,.2*observed[-2]):
        raise ValueError(f'The subtracted local remainder did not settle: {observed}.')
    return dict(cut_species=cut,tuned_to_step=is_edge,external_parent_energy=E,
                occupation_jump=float(jump),leading_bare_odd_pv=C,even_angular_log_coefficient=ell,
                predicted_radial_log_coefficient=log_coeff,measured_radial_cutoff_log_slope=float(slope),
                relative_step_slope_defect=relative,subtracted_radial_changes=observed.tolist(),rows=rows)


def full_radial(d,E,cut,lo,hi,is_edge,local_aperture=False):
    """Whole archived support, with separate step and continuous reconstructions."""
    M=d['x_end']
    m=d['mass_ratio']*M
    ms=d['checkpoints'][-1]['scalar_pole_over_T']
    Es=np.hypot(ms,Q)
    ma=M if cut=='parent' else m
    e0=E if cut=='parent' else E-Es
    name='parent' if cut=='parent' else 'fermion'
    p_edges=np.asarray(d[f'endpoint_{name}_momentum_edges'])
    edges=np.hypot(p_edges,ma)-e0
    if is_edge:
        j=np.argmin(abs(edges))
        if not (0<j<len(edges)-1 and abs(edges[j])<1e-14):
            raise ValueError('Failed to pin the declared archive discontinuity.')
        edges[j]=0.
    occ=np.asarray(d[f'endpoint_{name}_occupations'])
    centers=np.r_[edges[0],(edges[1:]+edges[:-1])/2,edges[-1]]
    smooth=np.r_[occ[0],occ,0.]
    kernel=angular_kernel(M,m,ms,E,cut)
    _,C,_=solve_axes(M,m,ms,E,cut)
    r0=np.sqrt((e0-ma)*(e0+ma))
    def occupation(x,kind):
        if kind=='step':
            idx=np.searchsorted(edges,x,side='right')-1
            return occ[idx] if 0<=idx<len(occ) else 0.
        return np.interp(x,centers,smooth,left=0.,right=0.)
    def integrand(x,kind):
        e=e0+x
        return float(kernel(x))*np.sqrt(max(0.,(e-ma)*(e+ma)))*occupation(x,kind)/2
    rho=3e-4
    coefficient=C*r0/2
    def exact_odd(kind):
        if kind=='step':
            cells=zip(edges[:-1],edges[1:],occ,occ,strict=True)
        else:
            cells=zip(centers[:-1],centers[1:],smooth[:-1],smooth[1:],strict=True)
        # Generalized symmetric finite part only for a true step exactly at zero.
        terms=[]
        for a,b,fa,fb in cells:
            slope=(fb-fa)/(b-a)
            constant=fa-slope*a
            terms.append((fb-fa)+constant*np.log((abs(b) or rho)/(abs(a) or rho)))
        return coefficient*fsum(terms)

    def integral(n,L,eps,kind,subtract):
        nodes,weights=roots_legendre(n)
        chi=list(edges if kind=='step' else centers)
        lim=[x for x in (lo-E,hi-E) if chi[0]<x<chi[-1]]
        geometric=[z for r in np.geomspace(L,max(abs(chi[0]),chi[-1]),45)
                   for z in (-r,r) if chi[0]<z<chi[-1]]
        collar=[t+sig*u for t in lim for sig in (-1,1)
                for u in (1e-7,1e-6,1e-5,1e-4,.001,.01)
                if chi[0]<t+sig*u<chi[-1]]
        fixed=sorted({*chi,*lim,*geometric,*collar,-L,L})
        intervals=[(a,b) for a,b in zip(fixed[:-1],fixed[1:],strict=True)
                   if a>=L or b<=-L]
        def f(z):
            return integrand(z,kind)-(coefficient*occupation(z,kind)/z if subtract else 0.)
        tail=[]
        for a,b in intervals:
            # Vanishing endpoint Jacobian resolves an aperture-edge logarithm
            # and the square root at zero internal momentum.
            x=np.pi*(nodes+1)/2
            u=a+(b-a)*np.sin(x/2)**2
            v=[f(z) for z in u]
            tail.append(np.pi*(b-a)/4*np.dot(weights,np.sin(x)*v))
        j=0.
        if is_edge and kind=='step':
            j=coefficient*(occ[np.searchsorted(edges,0.,side='right')-1]
                           -occ[np.searchsorted(edges,0.,side='left')-1])
        local=chi+lim+collar if local_aperture else chi
        bins=sorted({eps,L,*[abs(x) for x in local if eps<abs(x)<L],
                     *[10.**u for u in range(-12,-3) if eps<10.**u<L]})
        paired=[]
        for a,b in zip(bins[:-1],bins[1:],strict=True):
            x,y=np.log(a),np.log(b)
            if local_aperture:
                alpha=np.pi*(nodes+1)/2
                t=np.exp(x+(y-x)*np.sin(alpha/2)**2)
                W=weights*np.sin(alpha)*(y-x)*np.pi/4
            else:
                t=np.exp((x+y+(y-x)*nodes)/2)
                W=None
            v=[z*(f(z)+f(-z))-(0. if subtract else j) for z in t]
            paired.append(np.dot(W,v) if local_aperture else (y-x)/2*np.dot(weights,v))
        return fsum(tail+paired)+ (exact_odd(kind) if subtract else j*np.log(L/rho))

    result={}
    for kind in ('step','linear_energy'):
        a=integral(28,2e-4,2e-11,kind,False)
        b=integral(36,3e-4,1e-12,kind,True)
        c=integral(54,3e-4,1e-12,kind,True)
        defect=abs(c-b)
        if defect>5e-6 or abs(a-c)>2e-5:
            raise ArithmeticError(f'Whole-archive radial replay failed for {cut} {kind}: {a,b,c}.')
        result[kind]=dict(value=c,analytic_odd_full_support=exact_odd(kind),
                          independent_guarded_replay_defect=abs(a-c),quadrature_resolution_defect=defect)
    return result


def run():
    out=[]
    for grid in ('fine','third'):
        p=Path(__file__).parent/f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz'
        with gzip.open(p) as f:d=json.load(f)
        M=d['x_end'];m=d['mass_ratio']*M;ms=d['checkpoints'][-1]['scalar_pole_over_T']
        Es,lo,hi,*_=intervals(M,ms,m,Q)
        E=(lo+hi)/2
        a=[dict(**analyze(d,E,k,False),full_archived_support=full_radial(d,E,k,lo,hi,False))
           for k in ('parent','daughter')]
        cut='daughter' if grid=='fine' else 'parent'
        name='fermion' if cut=='daughter' else 'parent'
        edges=np.hypot(d[f'endpoint_{name}_momentum_edges'][1:-1],m if cut=='daughter' else M)
        ep=edges+(Es if cut=='daughter' else 0.)
        loc=np.flatnonzero((ep>lo+4e-4)&(ep<hi-4e-4))
        idx=loc[np.argmax(abs(np.diff(d[f'endpoint_{name}_occupations'])[loc]))]
        edge=float(ep[idx])
        a.append(dict(**analyze(d,edge,cut,True),
                      full_archived_support=full_radial(d,edge,cut,lo,hi,True)))
        checks={k:independent_angular_checks(M,m,ms,E,k,lo,hi) for k in ('parent','daughter')}
        out.append(dict(grid=grid,controls=a,independent_angular_checks=checks))
    return dict(scope='Paired radial PV diagnostic with frozen bare occupations; no physical thermal cut coefficients or inclusive decay.',
                rows=out)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
