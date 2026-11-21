"""Endpoint-adapted log-daughter-energy quadrature for 1<->2 decay.

At fixed mother momentum p and energy E, isotropic decay angles satisfy
q_phi=(E+p*mu)/2, dmu/2=dq_phi/p.
Transform x=log(q_phi/q_lo)/log(q_hi/q_lo) on [0,1].
Then dmu/2=[q_phi/p*log(q_hi/q_lo)] dx.
This is a change of integration variable, not an alteration to the
continuum two-body decay kinematics. Parent thermal bath is still fixed.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.special import roots_legendre, expit
from one_to_two_fixed_bath import FlatThermalDecay


class LogQDaughterKernel(FlatThermalDecay):
    def __init__(self,nq=100,np1=140,nmu=72,M1_T=.1):
        super().__init__(nq=nq,np1=np1,nmu=nmu,M1_T=M1_T)
        xn,xw=roots_legendre(nmu)
        pn,pw=roots_legendre(np1)
        p=13*(pn+1)
        wp=13*pw
        P=p[:,None]
        E=np.sqrt(P*P+M1_T*M1_T)
        q_lo=(E-P)/2
        q_hi=(E+P)/2
        logrange=np.log(q_hi/q_lo)
        x=(xn+1)/2
        wx=xw/2
        q_phi=q_lo*np.exp(logrange*x[None,:])
        q_chi=E-q_phi
        # Legendre weights dx multiplied by dmu/2dx.
        wang=(q_phi*logrange/P)*wx[None,:]
        ang_norm=wang.sum(axis=1)
        phi_mean=(wang*q_phi).sum(axis=1)
        self.angular_normalization_max_error=float(np.max(np.abs(ang_norm-1)))
        self.angular_phi_mean_energy_rel_error=float(np.max(np.abs(phi_mean/(E[:,0]/2)-1)))
        self.min_event_daughter_q=float(min(np.min(q_phi),np.min(q_chi)))
        self.qphi=q_phi.ravel()
        self.qchi=q_chi.ravel()
        self.parentE=np.broadcast_to(E,q_phi.shape).ravel()
        self.parent_fermi=np.broadcast_to(expit(-E),q_phi.shape).ravel()
        self.rateprefactor=(4/(2*math.pi**2)*
            np.broadcast_to(P*P*wp[:,None]*(M1_T/E)*wang,q_phi.shape)).ravel()
        self.ip,self.wp=self._interp_weights(self.qphi)
        self.ic,self.wc=self._interp_weights(self.qchi)
        if self.angular_normalization_max_error>1e-9:
            raise ArithmeticError("Angular change-of-variable normalization failed")
        if self.angular_phi_mean_energy_rel_error>1e-9:
            raise ArithmeticError("Angular daughter-energy mean failed")
