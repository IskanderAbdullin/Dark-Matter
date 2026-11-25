"""Fixed-parameter geometry and smooth-extra-component growth diagnostics.

Includes a covariance-aware, free-ruler distance fit to public DESI DR2 BAO
compression. It does not recompute the drag and visibility redshifts or fit
cosmological parameters other than the common BAO ruler.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.linalg import cho_factor, cho_solve

from direct_activation_consistency import MPL_REDUCED
from direct_subthreshold_friedmann import solve_friedmann
from late_activation_overlap import H, H0, Om
from resonant_mediator_inverse import C, HBAR

MPC_CM = 3.085677581e24
R_BARYON = 3*.02237/(4*2.4728e-5)


class Expansion:
    def __init__(self, time, hubble):
        self.log_h = PchipInterpolator(time, np.log(hubble))
        self.radial_primitive = PchipInterpolator(
            time, C*np.exp(-time)/(MPC_CM*hubble)
        ).antiderivative()

    def h(self, ln_a):
        return np.exp(self.log_h(ln_a))

    def radial(self, z):
        return C/MPC_CM/self.h(-np.log1p(z))

    def transverse(self, z):
        return self.radial_primitive(0.)-self.radial_primitive(-np.log1p(z))

    def sound_horizon(self, z):
        def integrand(log_a):
            a = np.exp(log_a)
            return C/(MPC_CM*np.sqrt(3*(1+R_BARYON*a))*a*self.h(log_a))

        return quad(integrand, -np.log(1e14), -np.log1p(z), epsabs=2e-7, epsrel=2e-9)[0]

    def growth(self, initial_z):
        t = -np.log1p(initial_z)
        slope = self.log_h.derivative()

        def rhs(ln_a, s):
            d, v = s
            matter = Om*(H0/self.h(ln_a))**2*np.exp(-3*ln_a)
            return v, 1.5*matter*d-(2+slope(ln_a))*v

        sol = solve_ivp(rhs, (t, 0.), np.exp(t)*np.ones(2), rtol=2e-9, atol=1e-11,
                        method='DOP853', max_step=.06, dense_output=True)
        if not sol.success:
            raise RuntimeError(sol.message)
        return {str(z): dict(D=float(d), f=float(v/d))
                for z in (0., .5, 1., 2., 5.)
                for d, v in [sol.sol(-np.log1p(z))]}


def observations(directory):
    prefix = directory/'desi_gaussian_bao_ALL_GCcomb'
    mean_file = Path(str(prefix)+'_mean.txt')
    cov_file = Path(str(prefix)+'_cov.txt')
    mean = np.loadtxt(mean_file, dtype=str, skiprows=1)
    cov = np.loadtxt(cov_file)
    signature = {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in (mean_file, cov_file)}
    return mean[:, 0].astype(float), mean[:, 1].astype(float), mean[:, 2], cov, signature


def profile_ruler(data, covariance, distances):
    c = cho_factor(covariance)
    w = cho_solve(c, distances)
    normalization = distances @ w
    scale = data @ w/normalization
    residual = distances*scale-data
    return dict(ruler_Mpc=float(1/scale), ruler_sigma_local_Mpc=float(1/(scale**2*np.sqrt(normalization))),
                chi2=float(residual @ cho_solve(c, residual)),
                predicted_at_best_ruler=(distances*scale).tolist())


def evaluate(expansion, obs):
    zs, data, quantities, covariance, _ = obs
    ds = {float(z): (expansion.radial(z), expansion.transverse(z)) for z in np.unique(zs)}
    values = []
    for z, name in zip(zs, quantities):
        dh, dm = ds[z]
        values.append({'DH_over_rs': dh, 'DM_over_rs': dm, 'DV_over_rs': np.cbrt(z*dm*dm*dh)}[name])
    values = np.asarray(values)
    ruler = expansion.sound_horizon(1060.)
    stellar_ruler = expansion.sound_horizon(1090.)
    dm_cmb = expansion.transverse(1090.)
    fixed = values/ruler-data
    growth = {str(z): expansion.growth(z) for z in (30., 100.)}
    return dict(distance_Mpc=values.tolist(), free_ruler=profile_ruler(data, covariance, values),
                approx_r_drag_fixed_redshift_Mpc=ruler,
                approximate_fixed_ruler_chi2=float(fixed @ cho_solve(cho_factor(covariance), fixed)),
                r_s_fixed_z1090_Mpc=stellar_ruler, D_M_z1090_Mpc=dm_cmb,
                theta_s_fixed_z1090=stellar_ruler/dm_cmb, smooth_extra_component_growth=growth)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bao-dir', type=Path, required=True)
    parser.add_argument('--bao-commit', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    obs = observations(args.bao_dir)
    root = Path(__file__).with_name('results')
    source = json.loads((root/'converged-independent-kernels-2026-11-25.json').read_text(encoding='utf-8-sig'))
    candidates = json.loads((root/'direct-subthreshold-friedmann.json').read_text())
    report = dict(bao_source_commit=args.bao_commit, bao_source_sha256=obs[-1],
                  redshift=obs[0].tolist(), quantities=obs[2].tolist())
    for candidate in candidates:
        mesh = candidate['mesh']
        p = candidate['params']
        row = source[f'parent_angle{mesh}_sparse_k00685'][0]
        settings = candidate['full_friedmann']
        _, history = solve_friedmann(row, p, n=settings['n'], rtol=settings['rtol'],
                                     order=settings['order'], return_history=True)
        t, r = history['log_a'], history['r']
        if mesh == 240:
            base = H(r-1)
            report['LCDM_reference'] = evaluate(Expansion(t, base), obs)
            td = p['Td_on_eV']/(p['z_on']+1)
            rho_eta_rel = p['Neta']*np.pi**2/30*td**4
            eta_h = np.sqrt(base**2 + rho_eta_rel*r**4/(3*(MPL_REDUCED*HBAR)**2))
            report['LCDM_plus_always_relativistic_eta'] = evaluate(Expansion(t, eta_h), obs)
        report[f'direct_pole_{mesh}'] = evaluate(Expansion(t, history['H']), obs)
    report['scope'] = ('Flat background fixed to the archived H0, matter, and radiation densities. '
                       'BAO fit profiles a single ruler, not the cosmological parameters. '
                       'Acoustic redshifts are held fixed. Growth models the extra component as smooth and omits '
                       'early transfer, perturbations, scale dependence, and particle interactions.')
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    for k, v in report.items():
        if isinstance(v, dict) and 'free_ruler' in v:
            print(k, 'free ruler chi2', v['free_ruler']['chi2'],
                  'fixed ruler chi2', v['approximate_fixed_ruler_chi2'],
                  'theta', v['theta_s_fixed_z1090'])


if __name__ == '__main__':
    main()
