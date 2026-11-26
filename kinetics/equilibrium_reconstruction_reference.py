"""Point-nodal affinity reconstruction control, not a FV occupation closure."""
import json
import numpy as np
from scipy.special import expit
from massive_shell_moments import shell_weights
from pair_projection_reference import bracket_energy


def reconstruct(energy, shell_mean, node_affinity, statistics):
    """Interpolate thermodynamic affinity, then invert to occupation."""
    h = np.asarray(node_affinity, dtype=float)
    if h.shape != np.shape(shell_mean) or not np.all(np.isfinite(h)):
        raise ValueError("finite nodal affinities are required")
    if statistics == "BE" and np.any(h >= 0):
        raise ValueError("BE node affinities must be negative")
    i, j, w = bracket_energy(energy, shell_mean)
    h = (1-w)*h[i]+w*h[j]
    if statistics == "FD":
        return expit(h)
    if statistics == "BE":
        return 1/np.expm1(-h)
    raise ValueError("unknown statistics")


def run():
    rng = np.random.default_rng(241126)
    n = 45000
    M = .005
    m,b = .08*M,.18*M
    E = np.hypot(np.exp(rng.uniform(np.log(.005),np.log(25),n)),M)
    P = np.sqrt((E-M)*(E+M))
    c = rng.uniform(-1,1,n)
    u = (M*M+m*m-b*b)/(2*M)
    k = np.sqrt((u-m)*(u+m))
    eb = (u*E+k*P*c)/M
    ef = E-eb
    # Chemical equilibrium: mu_parent = mu_fermion + mu_boson.
    mb, mf = .995*m, -.02
    f1 = expit(mb+mf-E)
    rows = []
    for cells in (60,120,240):
        q = np.r_[0,np.geomspace(1e-6,40,cells)]
        values = []
        for part, mass, mu, typ in ((eb,m,mb,"BE"),(ef,b,mf,"FD")):
            nn,ee,_,_ = shell_weights(q,mass)
            em = ee/nn
            h = mu-em
            i,j,w = bracket_energy(part,em)
            fn = 1/np.expm1(-h) if typ == "BE" else expit(h)
            values.append((reconstruct(part,em,h,typ), (1-w)*fn[i]+w*fn[j]))
        res=[]
        for col in (0,1):
            fb,ff=values[0][col],values[1][col]
            g=f1*(1-ff)*(1+fb)
            l=(1-f1)*ff*fb
            r=abs(g-l)
            res.append({"weighted_abs_residual":float(r.sum()/(g+l).sum()),
                        "max_per_event_residual":float(np.max(r/(g+l)))})
        assert res[0]["max_per_event_residual"] < 1e-13
        assert res[1]["weighted_abs_residual"] > 1e-7
        # Purely representability / bound check, not an ODE positivity proof.
        assert np.all((values[1][0]>=0)&(values[1][0]<=1))
        rows.append({"cells":cells,"affinity":res[0],"occupation_linear":res[1]})
    assert reconstruct([1.5],[1.,2.],[100.,110.],"FD")[0] <= 1
    try:
        reconstruct([1.5],[1.,2.],[-.1,.3],"BE")
    except ValueError:
        pass
    else:
        raise AssertionError("unphysical BE interpolation cannot be accepted")
    return {"events":n,"seed":241126,"mu_bose_over_mass":.995,"rows":rows,
            "invalid_BE_affinity_rejected":True}


if __name__ == "__main__":
    print(json.dumps(run(),indent=2,allow_nan=False))
