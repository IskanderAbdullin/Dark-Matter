"""Regularity of a frozen Born aperture under a signed scalar-energy translation."""
import bisect
import gzip
import json
from pathlib import Path

import mpmath as mp


class Born:
    def __init__(self,data):
        self.M=mp.mpf(str(data['x_end']))
        self.m=self.M*mp.mpf(str(data['mass_ratio']))
        self.ms=mp.mpf(str(data['checkpoints'][-1]['scalar_pole_over_T']))
        self.q=mp.mpf('4.1e-6')
        self.Es=mp.sqrt(self.ms*self.ms+self.q*self.q)
        self.ec=[]
        self.occ=[]
        self.pc=[]
        for name,m in (('parent',self.M),('fermion',self.m)):
            p=list(map(lambda x:mp.mpf(str(x)),data[f'endpoint_{name}_momentum_edges']))
            self.pc.append(p)
            self.ec.append([mp.sqrt(z*z+m*m) for z in p])
            self.occ.append(list(map(lambda x:mp.mpf(str(x)),data[f'endpoint_{name}_occupations'])))
        self.L,self.U=self.aperture(self.Es)

    def aperture(self,E):
        S=(E-self.q)*(E+self.q)
        M,m=self.M,self.m
        A=(M*M+S-m*m)/2
        d2=((M-m)**2-S)*((M+m)**2-S)/4
        T=A*E+self.q*mp.sqrt(d2)
        return (M*M*E*E+d2)/T,T/S

    def evaluate(self,E):
        L,U=self.aperture(E)
        p,b=self.ec
        all_edges=sorted({L,U,*[x for x in p if L<x<U],*[x+E for x in b if L<x+E<U]})
        seed=mp.mpf(0)
        stim=mp.mpf(0)
        for left,right in zip(all_edges[:-1],all_edges[1:],strict=True):
            y=(left+right)/2
            fp=self.occ[0][bisect.bisect_right(p,y)-1]
            fb=self.occ[1][bisect.bisect_right(b,y-E)-1]
            seed+=(right-left)*fp*(1-fb)
            stim+=(right-left)*(fp-fb)
        return seed,stim

    def tune(self):
        # Choose an upward parent occupation edge near the nominal support.
        # Align it exactly with the lower (parallel daughter) vacuum face.
        rows=[]
        for pp,cc,f0,f1 in zip(self.pc[0][1:-1],self.ec[0][1:-1],self.occ[0][:-1],self.occ[0][1:],strict=True):
            if mp.mpf('1.35')<pp<mp.mpf('1.50') and f1>f0:
                rows.append((abs(cc-self.L),pp,f1-f0))
        _,p,jump=min(rows)
        left=mp.sqrt(p*p+self.M*self.M)
        k=p+self.q
        es=left-mp.sqrt(k*k+self.m*self.m)
        if abs(self.aperture(es)[0]-left)>mp.mpf('1e-61'):
            raise ValueError('Tuned energy did not reach the requested shell face.')
        slope=mp.diff(lambda z:self.aperture(z)[0],es)
        return es,p,jump,slope


def run():
    cases=[]
    with mp.workdps(77):
        for grid in ('fine','third'):
            path=Path(__file__).parent/f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz'
            with gzip.open(path) as h: d=json.load(h)
            B=Born(d)
            nominal=B.evaluate(B.Es)
            Ee,p,jump,slope=B.tune()
            tuned=B.evaluate(Ee)
            predicted=jump*slope
            rows=[]
            for h in [mp.mpf(x) for x in ('1e-7','1e-9','1e-11','1e-13','1e-15')]:
                p0,m0=B.evaluate(B.Es+h),B.evaluate(B.Es-h)
                p1,m1=B.evaluate(Ee+h),B.evaluate(Ee-h)
                ns=(p0[1]+m0[1]-2*nominal[1])/(h*h)
                tl=(p1[1]+m1[1]-2*tuned[1])/h
                rows.append(dict(shift=float(h),nominal_quadratic_response=float(ns),
                                 tuned_linear_response=float(tl),
                                 tuned_relative_to_predicted_cusp=float(tl/predicted-1)))
            if abs(rows[-1]['tuned_relative_to_predicted_cusp'])>1e-8:
                raise ValueError('Tuned cusp did not reproduce the local shell jump.')
            if abs(rows[-1]['nominal_quadratic_response']/rows[-2]['nominal_quadratic_response']-1)>1e-8:
                raise ValueError('Nominal aperture did not reach its quadratic limit.')
            resolution=[]
            for w in map(mp.mpf,('1e-7','1e-9','1e-11')):
                dm=mp.diff(lambda e:B.evaluate(e)[1],Ee-w)
                dp=mp.diff(lambda e:B.evaluate(e)[1],Ee+w)
                curv=(dp-dm)/(2*w)
                resolution.append(dict(top_hat_half_width=float(w),finite_curvature=float(curv),
                                       aperture_smearing_scaled_cusp_error=float(abs(2*w*curv/predicted-1))))
            if resolution[-1]['aperture_smearing_scaled_cusp_error']>1e-6:
                raise ValueError('Finite-resolution curvature did not recover the known cusp coefficient.')
            cases.append(dict(grid=grid,q=float(B.q),native_scalar_energy=float(B.Es),
                              native_bare_stimulus=float(nominal[1]),tuned_p=float(p),
                              tuned_scalar_energy=float(Ee),tuned_shell_jump=float(jump),
                              tuned_lower_endpoint_derivative=float(slope),
                              predicted_cusp_coefficient=float(predicted),ladder=rows, finite_resolution=resolution))
    return dict(scope='One-parameter exact Born aperture regularity test; not real-photon recoil or inclusive thermal rates.',rows=cases)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
