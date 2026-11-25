"""Conservative radial quadrature for discrete classical lattice occupations."""
import numpy as np

def occupation_radial(k2, edges, f, p0, L, cut, fmin=100.):
    e=np.asarray(edges); f=np.where(np.asarray(f)>=fmin,f,0.)
    active=(k2>0)&(k2<=cut**2)
    shell,inv,count=np.unique(np.round(k2[active],12),return_inverse=True,return_counts=True)
    q=p0*np.sqrt(shell)
    assert len(q)>1
    b=np.r_[q[0]/2,(q[:-1]+q[1:])/2,p0*cut]
    assert np.all(np.diff(b)>0)
    cumulative=np.r_[0,np.cumsum(f*np.diff(e**3)/(6*np.pi**2))]
    def integral(q):
        q=np.clip(q,e[0],e[-1])
        j=np.clip(np.searchsorted(e,q,side='right')-1,0,len(f)-1)
        return cumulative[j]+f[j]*(q**3-e[j]**3)/(6*np.pi**2)
    target=np.maximum(np.diff(integral(b)),0.)
    result=np.zeros_like(k2)
    result[active]=(target*L**3/(count*p0**3))[inv]
    return result,dict(radial_lower_q=float(b[0]),radial_lost_number=float(integral(b[0])),
                       exact_cut_number=float(integral(b[-1])),radial_represented_number=float(target.sum()))
