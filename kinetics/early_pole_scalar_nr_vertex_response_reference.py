"""Isotropic nonrelativistic large-N vertex-response diagnostic on IR snapshots."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad


def primitive(k,a):
    # Integral k log|k-a| dk. Its finite limiting term vanishes at k=a.
    y=k-a
    g=np.log(np.where(y!=0,abs(y),1.))
    return y*y*(.5*g-.25)+a*y*(g-1.)


def response(edges,f,P,z,m,g,check=False):
    """Eq. (26), arXiv:1710.11146, with z=2*m*omega/P^2."""
    # Scale all k by P; the logs' constant log(P) cancels.
    l,h=edges[:-1]/P,edges[1:]/P
    a=np.array([(-1-z)/2,(-1+z)/2,(1-z)/2,(1+z)/2])
    c=np.array([1.,1.,-1.,-1.])
    term=(primitive(h[:,None],a)-primitive(l[:,None],a))@c
    re=(term@f)
    left,right=sorted((abs(1-z)/2.,abs(1+z)/2.))
    l1=np.maximum(l,left); h1=np.minimum(h,right)
    imag=np.pi*(np.maximum(h1-l1,0.)*(h1+l1)*.5)@f*np.sign(z)
    norm=m*g*P/(4*np.pi**2)
    qe=0.
    if check:
        total=0.
        for lo,hi,ff in zip(l,h,f,strict=True):
            if ff<=0.:
                continue
            roots=np.unique(a[(a>lo)&(a<hi)]).tolist()
            def func(t):
                return t*(c@np.log(np.maximum(abs(t-a),1e-300)))
            v=quad(func,lo,hi,points=roots,epsabs=2e-11,epsrel=2e-10,limit=180)[0]
            total+=ff*v
        qe=abs(re-total)/max(abs(total),abs(imag),1.)
        assert qe<3e-7,(re,total,qe,P,z)
    return complex(norm*re,norm*imag),qe


def one(s):
    m=s['checkpoints'][-1]['scalar_pole_over_T']
    edges=np.asarray(s['scalar_momentum_edges'])
    p=np.asarray(s['scalar_momentum_nodes'])
    f=np.maximum(np.asarray(s['endpoint_scalar_occupations']),0.)
    nt=np.asarray(s['endpoint_scalar_signed_shell_number_over_T3'])
    # Restrict the input to the ultracold component; the exact hard response
    # needs a relativistic kernel and is not extrapolated from this calculation.
    h=np.searchsorted(edges,1e-5)
    assert edges[h]==1e-5 and edges[h]/m<.02
    K=p[np.argmax(p*p*f)]
    eps=7.65175e-8
    g=eps/(8*m*m)
    result=[]
    for u in (.25,.5,1.,2.,4.):
        for z in (0.,.5,1.,2.):
            val,err=response(edges[:h+1],f[:h],u*K,z,m,g,True)
            neg,_=response(edges[:h+1],f[:h],u*K,-z,m,g)
            assert abs(neg-val.conjugate())<1e-10*max(abs(val),1.)
            result.append(dict(transfer_over_peak=u,frequency_over_P2_div_2m=z,
                               polarization_real=val.real,polarization_imag=val.imag,
                               vertex_square_ratio=1./abs(1+val)**2,
                               quadrature_crosscheck_relative=err))
    return dict(x=s['x_end'],mass_ratio=s['mass_ratio'],peak_q=float(K),m_over_T=m,
                IR_number_fraction=float(nt[:h].sum()/nt.sum()),same_contact_as_lambda_over_8m2=g,
                samples=result)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('spectra',nargs='+')
    a=p.parse_args()
    out={}
    for f in a.spectra:
        ff=Path(f)
        raw=gzip.decompress(ff.read_bytes()) if ff.suffix=='.gz' else ff.read_bytes()
        out[ff.name]=one(json.loads(raw))
    print(json.dumps(out,indent=2,allow_nan=False))
