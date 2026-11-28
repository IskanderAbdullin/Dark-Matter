"""All parent energies: exact fixed-q signed real recoil, no virtual finite part."""
import gzip,json
from pathlib import Path
import functools
import numpy as np
from scipy.special import roots_legendre
from scipy.optimize import brentq
from early_pole_decay_mode_clock_reference import intervals,Q

base=Path(__file__).parent/'results'
@functools.lru_cache(maxsize=8)
def gl(n):return roots_legendre(n)


class Source:
    def __init__(self,d):
        self.M=np.longdouble(str(d['x_end']));self.m=self.M*np.longdouble(str(d['mass_ratio']))
        self.s=np.longdouble(str(d['checkpoints'][-1]['scalar_pole_over_T'])); self.q=np.longdouble(str(Q))
        self.Es,self.lo,self.hi,*_=intervals(self.M,self.s,self.m,Q)
        self.pe=np.hypot(np.asarray(d['endpoint_parent_momentum_edges'],dtype=np.longdouble),self.M)
        self.be=np.hypot(np.asarray(d['endpoint_fermion_momentum_edges'],dtype=np.longdouble),self.m)
        self.pf=np.asarray(d['endpoint_parent_occupations'])
        self.bf=np.asarray(d['endpoint_fermion_occupations'])
        dy=np.unique(np.r_[self.pe,self.be+self.Es])
        xx=np.r_[self.lo,dy[(dy>self.lo)&(dy<self.hi)],self.hi]
        mid=(xx[1:]+xx[:-1])/2
        self.stim=np.diff(xx)@(self.pf[np.searchsorted(self.pe,mid)-1]-self.bf[np.searchsorted(self.be,mid-self.Es)-1])
        t=np.arccosh((self.M*self.M+self.m*self.m-self.s*self.s)/(2*self.M*self.m))
        self.B=2*(t/np.tanh(t)-1)

    def signed(self,k,n=30,o=8,az=20):
        k=np.longdouble(k)
        M,m,s,q,Es=self.M,self.m,self.s,self.q,self.Es
        cut=np.unique(np.r_[self.pe,self.be+Es+k,self.lo,self.hi])
        cut=cut[(cut>self.lo-.12)&(cut<self.hi+.12)]
        breaks=[0.,2.]
        if k:
            for ep in cut:
                p=np.sqrt((ep-M)*(ep+M))
                eb=ep-k-Es
                kb=np.sqrt((eb-m)*(eb+m))
                A=(ep-k)*Es-(M*M+s*s-m*m)/2
                low=M*M/(ep+p)
                for v in (q*q-A-q*kb,q*q-A+q*kb):
                    t=(v/k-low)/p
                    if 0.<t<2.:
                        breaks.append(t)
        shift=(m/((self.lo+self.hi)/2))**2/2
        X,W=np.asarray(gl(n),dtype=np.longdouble);x,w=np.asarray(gl(o),dtype=np.longdouble)
        phi=np.cos((np.arange(az,dtype=np.longdouble)+.5)*2*np.pi/az)
        out=np.zeros((2,3),dtype=np.longdouble)
        for thl,thr in zip(sorted(set(breaks))[:-1],sorted(set(breaks))[1:],strict=True):
            phase=np.log([thl+shift,thr+shift])
            zz=(phase[0]+phase[1])/2+np.diff(phase)[0]*X/2
            for t,sup in zip(np.exp(zz)-shift,W*np.diff(phase)[0]*np.exp(zz)/4,strict=True):
                def bounds(E,sign):
                    p=np.sqrt((E-M)*(E+M))
                    el=M*M/(E+p)+p*t
                    zR=np.sqrt((p-k)**2+2*p*k*t)
                    formula=((E-k)*Es-(M*M+s*s-m*m)/2+k*el)-sign*q*zR
                    dP=E/p
                    dL=-(M*M/(E+p)**2)*(1+dP)+t*dP
                    dR=dP*(p-k+k*t)/zR
                    return formula,Es+k*dL-sign*q*dR
                def endpoint(E,sign):
                    for z in range(6):
                        a,b=bounds(E,sign)
                        E-=a/b
                    if abs(bounds(E,sign)[0])>2e-20:
                        raise ValueError('The moving exact support boundary did not close.')
                    return E
                left=endpoint(self.lo,-1.)
                right=endpoint(self.hi,1.)
                xx=np.r_[left,cut[(cut>left)&(cut<right)],right]
                hs=np.diff(xx)/2
                ep=((xx[1:]+xx[:-1])/2)[:,None]+hs[:,None]*x
                p=np.sqrt((ep-M)*(ep+M))
                el=M*M/(ep+p)+p*t
                Qr=np.sqrt((p-k)**2+2*p*k*t)
                cs=((ep-k)*Es-(M*M+s*s-m*m)/2+k*el)/(q*Qr)
                if max(abs(cs).ravel())>1+1e-9: raise ValueError('Unexpected parent support.')
                z=1-t
                zR=(p*z-k)/Qr
                si=p*np.sqrt(t*(2-t))/Qr
                Bn=el[...,None]-Es+q*((cs*zR)[...,None]+(si*np.sqrt(np.maximum(0,1-cs*cs)))[...,None]*phi)
                a,b=-el[...,None],-Bn
                d=(M*M+m*m-s*s)/2+k*(Bn-el[...,None])
                V=(1/a-1/b)/2
                H=M*M/(a*a)+m*m/(b*b)-2*d/(a*b)
                U=b*M*M/a+a*m*m/b-2*d
                A=-4*(d+M*m)*H-8*k*V*U+16*k*k*V*V*a*b
                Wtot=(hs[:,None]*w)*A.mean(axis=-1)*(p/Qr)/(2*((M+m)**2-s*s))*sup
                fp=self.pf[np.searchsorted(self.pe,ep)-1]
                fb=self.bf[np.searchsorted(self.be,ep-k-Es)-1]
                factors=np.asarray([fp*(1-fb),(1-fp)*fb,fp-fb])
                out[0]+=np.sum(factors*Wtot,axis=(1,2))
                outside=(ep<self.lo)|(ep>self.hi)
                out[1]+=np.sum(factors*(Wtot*outside),axis=(1,2))
        return out



def cap_knots(S,limit):
    M,m,s,q,E=S.M,S.m,S.s,S.q,S.Es
    out=[]
    rest=np.r_[S.pe,S.be+E]
    rest=rest[(rest>S.lo-.10)&(rest<S.hi+.10)]
    for edge in rest:
        for moves in (False,True):
            for face in (-1.,1.):
                for z in (-1.,1.):
                    def f(k):
                        ep=edge+k if moves else edge
                        p=np.sqrt((ep-M)*(ep+M))
                        return ((ep-k)*E-(M*M+s*s-m*m)/2+k*(ep-p*z)
                                +face*q*(p-k*z))
                    if f(-limit)*f(limit)<0:
                        v=abs(brentq(f,-limit,limit,xtol=2e-17,rtol=1e-14))
                        if v<limit:
                            out.append(v)
    # Extra knots produced by selecting both edge interpretations are
    # harmless; true support changes must all remain in the list.
    return sorted(set(out))


def integrate_cut(S,n,o,az,outer,tiny):
    max_cap=1e-5
    caps=(2.5e-6,5e-6,max_cap)
    base=S.signed(0.,n,o,az)[0,2]
    T=.9
    def func(k):
        a=S.signed(k,n,o,az)[0]
        b=S.signed(-k,n,o,az)[0]
        v=(a[0]-b[1]-base)/k
        z=(a[2]+b[2]-2*base)/(k*np.expm1(k/T))
        return np.r_[v,z]
    low_coeff=func(tiny)
    knot=cap_knots(S,max_cap)
    wins=np.zeros(2,dtype=np.longdouble)+tiny*low_coeff
    x,w=gl(outer)
    result=[]
    for L,U in zip((tiny,*caps[:-1]),caps,strict=True):
        ends=np.asarray([L,*[q for q in knot if L<q<U],U])
        for a,b in zip(ends[:-1],ends[1:],strict=True):
            r=(a+b)/2+(b-a)*x/2
            wins+=(b-a)/2*(w@np.asarray([func(z) for z in r]))
        result.append(dict(plasma_soft_cap=U,parts=wins.astype(float).tolist(),
                           normalized_before_gauge=float(wins.sum()/S.stim)))
    return result


def compute(src,n,o,az):
    base=src.signed(0.,n,o,az)
    if max(abs(base[0,2]-src.B*src.stim),abs(base[1,2]))>2e-13:
        raise ValueError('Integrated Born support missed the exact leading angular coefficient.')
    rows=[]
    for k in (1e-8,1e-7,3e-7,1e-6,3e-6,1e-5):
        a=src.signed(k,n,o,az)
        b=src.signed(-k,n,o,az)
        vac=a[0,0]-b[0,1]-base[0,2]
        rad=a[0,2]+b[0,2]-2*base[0,2]
        outside=a[1,2]+b[1,2]
        if abs(rad-(a[0,2]-a[1,2]+b[0,2]-b[1,2]-2*base[0,2]+outside))>1e-18:
            raise ValueError('Parent-support decomposition failed.')
        rows.append(dict(photon_energy=k,no_photon_born_context=float(src.stim),
                         Born_times_soft_radiator=float(base[0,2]),
                         vacuum_signed_remainder_per_energy=float(vac/k),
                         thermal_even_remainder_per_energy_squared=float(rad/k**2),
                         signed_real_parent_exterior=float(outside),
                         signed_real_parent_interior_minus_soft=float(rad-outside)))
    return rows


def run():
    rows=[]
    for grid in ('fine','third'):
        with gzip.open(base/f'early-pole-spectrum-{grid}-m06-2026-11-26.json.gz') as f:d=json.load(f)
        src=Source(d)
        a=compute(src,64,12,24)
        b=compute(src,96,18,36)
        for AA,BB in zip(a,b,strict=True):
            k=BB['photon_energy']
            r=lambda x:x['thermal_even_remainder_per_energy_squared']
            err=abs(r(AA)-r(BB))
            if err>max(.035,0.015*abs(r(BB))) or abs(AA['vacuum_signed_remainder_per_energy']-BB['vacuum_signed_remainder_per_energy'])>1e-5:
                raise ValueError(f'Parent-integrated refinement failed: {grid}, {k}, {err}.')
            BB.update(thermal_even_scaled_mesh_defect=err)
        left=integrate_cut(src,64,12,24,12,1e-7)
        right=integrate_cut(src,80,14,28,20,3e-7)
        for l,r in zip(left,right,strict=True):
            er=float(max(abs(np.asarray(l['parts'])-r['parts'])))
            if er>2e-8:
                raise ValueError(f'Signed photon integral mesh or infrared cap failed: {er}.')
            r['photon_mesh_and_low_cap_abs_defect']=er
        rows.append(dict(grid=grid,parent_Born_window=[float(src.lo),float(src.hi)],results=b,integrals=right))
    return dict(scope='Signed soft recoil at one external scalar q integrated over all allowed frozen archive parents; no finite virtual remainder.',rows=rows)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
