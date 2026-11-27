"""Finite-part diagnostic for both sequential poles in one actual spacelike slice."""
import gzip,json
from pathlib import Path
import numpy as np
from numpy.polynomial import Polynomial
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import roots_legendre
from early_pole_hard_yukawa_photoproduction_reference import I,G,dot,slash,spinors
from finite_mass_qed_nonnull_photon_ward_poles_reference import setup


def moments(a,b,d,p):
    if not a<0.<b:
        raise ArithmeticError('Diagonal aperture does not contain the root')
    if d==0.:
        m0=1/a-1/b
        x=[m0,np.log(b/(-a))]
        for n in range(2,p+1):x.append((b**(n-1)-a**(n-1))/(n-1))
        return m0,np.array(x),0.
    m0=-(np.arctan(d/b)+np.arctan(d/(-a)))/d
    x=[np.pi/d+m0,.5*np.log((b*b+d*d)/(a*a+d*d))]
    for n in range(2,p+1):x.append((b**(n-1)-a**(n-1))/(n-1)-d*d*x[n-2])
    return m0,np.array(x),d


def one(d,ctrl,rootrow,k,w,c,A=-.35,B=1.45):
    ph=rootrow['azimuth']
    roots=rootrow['found_near_born_intermediate_on_shell_roots']
    if [t['arm'] for t in roots]!=[0,1]:
        raise ArithmeticError('Expected ordered physical poles')
    M,m,q,K,Qp,geom=setup(d,k,w,c)
    e0,e1=(x['parent_E'] for x in roots)
    h=2*dot(q,K)
    par=np.linalg.norm(Qp[1:])
    az=Qp[1:]/par
    azp=np.array([az[2],0.,-az[0]])
    ax=Qp[0]/par
    a0=(M*M+dot(Qp,Qp)-m*m)/2
    low=e0+(e1-e0)*(A-.8)
    high=e0+(e1-e0)*(B+.8)
    def t_of_E(E):
        p,b=geom(E,ph)
        return (dot(K,K)+2*dot(p,K))/h
    def numerator(t):
        e=brentq(lambda z:t_of_E(z)-t,low,high,xtol=2e-15)
        p,b=geom(e,ph)
        pp=(e*Qp[0]-a0)/par
        pt=np.sqrt((e-M)*(e+M)-pp*pp)
        dpz=ax*az[2]+(e-pp*ax)/pt*np.cos(ph)*azp[2]
        J=h/(2*(w-k*dpz))
        if J<=0.:raise ArithmeticError('Nonmonotone pole coordinate')
        u=spinors(p,M)
        v=spinors(b,m).conj().T@G[0]
        S,T,Ks=map(slash,[p,b,K])
        a=[];b0=[]
        for pol in np.eye(4)[1:3]:
            L=slash(pol)
            a.append(v@(S+Ks+M*I)@L@u)
            b0.append(v@L@(T-Ks+m*I)@u)
        a=np.ravel(a);b0=np.ravel(b0)
        c0=a@b0.conj()
        return J*np.array([a@a.conj(),b0@b0.conj(),c0]).astype(complex),e,p,b
    def crop(t):
        y,E,p,b=numerator(t)
        j=[]
        for vec,q0 in ((p,'parent'),(b,'fermion')):
            ed=d[f'endpoint_{q0}_momentum_edges']
            idx=np.searchsorted(ed,np.linalg.norm(vec[1:]))-1
            if idx<0 or idx>=len(ed)-1:raise ArithmeticError('Occupied slice reached outside archive')
            j.append(int(idx))
        return j
    test=np.linspace(A,B,31)
    ref=[crop(t) for t in test]
    if any(t!=ref[0] for t in ref):
        raise ArithmeticError('Selected finite part lands on an occupation step')
    scale=(numerator(0.)[0][0].real+numerator(1.)[0][1].real)/2
    def raw(t):
        x=numerator(t)[0]/scale
        return np.array([x[0].real,x[1].real,x[2].real,x[2].imag])
    checks=np.array([raw(t) for t in test])
    if max(abs(checks[:,3]))>2e-7:
        raise ArithmeticError('Unpolarized mixed numerator was not real')
    def polynomials(degree):
        t=np.cos((np.arange(degree+4)+.5)*np.pi/(degree+4))
        x=(A+B+(B-A)*t)/2
        W=np.array([raw(t) for t in x])
        polys=[Polynomial.fit(x,W[:,i],degree).convert() for i in range(4)]
        for i,tr in enumerate((0.,1.)):
            polys[i]=polys[i]+(raw(tr)[i]-polys[i](tr))
        err=np.max(abs(np.array([p(test) for p in polys]).T-checks))
        if err>2e-7:raise ArithmeticError(('Smooth numerator compression',degree,err))
        return polys,float(err)
    def evaluate(pol,dlt,ratio=1.,limits=None):
        f,g,c0,m0=pol
        a,b=(A,B) if limits is None else limits
        x2=Polynomial([1.,-1.])
        dR=dlt*ratio
        # Partial fractions with different constant diagnostic widths.
        z=complex(1.,dlt-dR)
        c=(c0+1j*m0)/z
        cp=Polynomial(c.coef.real)
        mp=Polynomial(c.coef.imag)
        gy,cy,my=g(x2),cp(x2),mp(x2)
        n=max(len(p.coef) for p in (f,gy,cp,cy,mp,my))+2
        Fx,X,delta=moments(a,b,dlt,n)
        Fy,Y,d2=moments(1-b,1-a,dR,n)
        def s(z,v):return z.coef@v[:len(z.coef)]
        aa=Fx*f(0.)+s(Polynomial(f.coef[1:]),X[1:])
        bb=Fy*g(1.)+s(Polynomial(gy.coef[1:]),Y[1:])
        im=(delta*s(mp,X)-d2*s(my,Y) if dlt else np.pi*(mp(0.)-my(0.)))
        cc=2*(s(cp,X[1:])+s(cy,Y[1:])+im)
        return np.array([aa,bb,cc])
    p,q=polynomials(9),polynomials(13)
    conv=[]
    for delta in (.08,.025,.008,.002,.0005,.000125,.000025,0.):
        result=evaluate(q[0],delta)
        err=max(abs(result-evaluate(p[0],delta)))
        conv.append(dict(relative_width=delta,diagonal_finite_parts_and_interference=result.tolist(),
                         total=float(sum(result)),compression_difference=float(err)))
        if err>1e-6:raise ArithmeticError(('Finite part basis convergence',delta,err))
    def real_numeric(delta,ratio=1.):
        def body(t):
            a,b,x,y=raw(t)
            X=t+1j*delta
            Y=1-t-1j*delta*ratio
            return a/abs(X)**2+b/abs(Y)**2+2*(complex(x,y)/(X*Y)).real
        val=quad(body,A,B,points=[0.,1.],epsabs=1e-8,epsrel=1e-9,limit=140)[0]
        return val-np.pi*(raw(0.)[0]+raw(1.)[1]/ratio)/delta
    unequal=[]
    F=sum(evaluate(q[0],0.))
    for v in (.1,2.,10.):
        a=sum(evaluate(q[0],.034,v))
        err=abs(real_numeric(.034,v)-a)
        seq=[float(sum(evaluate(q[0],d,v))) for d in (.002,.0001,2.5e-6)]
        if err>3e-6 or abs(seq[-1]-F)>4e-5:
            raise ArithmeticError(('Unequal-width exact integration',v,err,seq,F))
        unequal.append(dict(second_to_first_width_ratio=v,finite_width_direct_integral_error=float(err),
                            corrected_integrals_at_first_widths_002_0001_2p5e_minus6=seq))
    narrow=(A+.09,B-.12)
    def regular(t):
        f,g,cr,ci=raw(t)
        return f/(t*t)+g/((1-t)**2)+2*cr/(t*(1-t))
    add=(quad(regular,A,narrow[0],epsabs=1e-10)[0]+
         quad(regular,narrow[1],B,epsabs=1e-10)[0])
    small=sum(evaluate(q[0],0.,limits=narrow))
    if abs(F-(small+add))>4e-8:raise ArithmeticError('Window additivity failed')
    checks_2=[abs(real_numeric(v)-sum(evaluate(q[0],v))) for v in (.055,.022)]
    if max(checks_2)>2e-6:raise ArithmeticError(('Independent normal width integral',checks_2))
    return dict(grid=ctrl,photon_k=k,photon_omega=w,angle_to_scalar=c,azimuth=ph,
                energy_roots=[e0,e1],positive_pole_invariant_gap=h,
                interval_in_two_pole_coordinate=[A,B],fixed_occupation_histogram_indices=ref[0],
                raw_pole_residues_normalized=[raw(0.)[0],raw(1.)[1]],max_numerator_fit_defect=q[1],
                maximum_direct_width_integral_defect=max(checks_2),widths=conv,unequal_width_checks=unequal,
                narrow_wider_window_total_limit_delta=float(small-F),
                exact_window_additivity_abs_defect=float(abs(F-small-add)))


def run():
    root=Path(__file__).parent/'results'
    r=json.loads((root/'finite-mass-qed-nonnull-photon-ward-poles-2026-11-27.json').read_text())
    rows=[]
    for x in r['data']:
        g=x['grid']
        with gzip.open(root/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as f:d=json.load(f)
        for index,c in ((0,0.),(1,-.5)):
            z=x['checks'][index]
            rows.append(one(d,g,z['azimuths'][0],z['k'],z['omega'],c))
    return dict(scope='Real spacelike bare transverse local slice, subtracting leading two sequential diagonal overlaps in a declared common-width convention; not a physical rate or Ward-preserving width.',data=rows)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
