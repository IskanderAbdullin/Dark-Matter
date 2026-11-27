"""A finite closed-path coherence is not automatically a spectral line."""
import json
import numpy as np
from scipy.integrate import quad
from early_pole_differential_phase_separation_reference import missing_corr
from early_pole_finite_wilson_htl_reference import fold
from early_pole_finite_wilson_poles_reference import pole_fold
from early_pole_matched_wilson_geometry_reference import inputs


def run():
    d=inputs(.005)
    # Probe arbitrarily short times in a fixed finite soft momentum domain.
    # The connector closes a loop with area of order t² at fixed velocities.
    v=d['vP']
    dv=d['vB_parallel']-v
    rows=[]
    for t in (.01,.02,.04,.08,.16):
        c=fold('triangle',t,t*dv,v,4,28,16,8.)
        p=pole_fold(t,t*dv,v,'triangle',7,25,8.)
        # Restore t factored out of the normalized weights.
        S=t*(sum(c)+sum(p[:2]))
        rows.append(dict(t=t,covariance_before_constant_prefactor=S,
                         covariance_divided_by_t4=S/(t**4)))
    vals=np.array([r['covariance_before_constant_prefactor'] for r in rows])
    slope=np.log2(vals[1:]/vals[:-1])
    assert np.max(abs(slope[:2]-4))<.035
    # Independent Gaussian-in-space Markov path witness. Finite-time real
    # modulus is valid as a path expectation, but cosine-transforming it
    # as if stationary need not give a positive spectral measure.
    def response(w,upper):
        return quad(lambda u:np.exp(-2*u*missing_corr(.60103446*u,'gaussian')),
                    0.,upper,weight='cos',wvar=w,epsabs=2e-12,limit=500)[0]*2/np.pi
    freqs=np.linspace(0.,8.,161)
    a=[response(w,30.) for w in freqs]
    b=[response(w,60.) for w in freqs]
    assert max(abs(np.array(a)-b))<2e-10 and min(b)<-.01
    return dict(scope='Finite covariance shape and independent counterexample; no physical line inferred.',
                soft_short_time_rows=rows,short_time_log2_slopes=slope.tolist(),
                toy_frequency=freqs.tolist(),naive_toy_cosine_density=b,
                naive_toy_minimum=min(b),frequency_at_min=freqs[np.argmin(b)])


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
