"""Audit conservative eta stress and trigger source for a future metric coupling."""
import argparse
import json
from pathlib import Path

import numpy as np

from direct_pole_eta_stability import EtaResponse


def audit(candidate):
    qrule = EtaResponse(32, 1e-12)
    q = qrule.q
    measure = qrule.weight*q*q*qrule.f*candidate['params']['Neta']/(2*np.pi**2)
    p = candidate['params']
    m0, v = p['m_eta_at_vev_eV'], p['vev_eV']
    values = []
    for x, z in ((.03, 8.), (.4, 4.5), (1., 4.), (1.3, 3.7), (.95, 0.)):
        temperature = p['Td_on_eV']/(p['z_on']+1)*(1+z)

        def moments(field, temp):
            yy = m0*field/temp
            e = np.hypot(q, yy)
            r = temp**4*(measure @ e)
            pressure = temp**4/3*(measure @ (q*q/e))
            trace = temp**4*yy**2*(measure @ (1/e))
            source = m0*m0*field*temp*temp/v*(measure @ (1/e))
            return np.array([r, pressure, trace, source])

        rho, pressure, trace, source = moments(x, temperature)
        slope = .12
        eps = 1e-5
        fd = (moments(x+eps, temperature)-moments(x-eps, temperature))/(2*eps)
        fd_time = (moments(x+slope*eps, temperature*np.exp(-eps))[0]
                   -moments(x-slope*eps, temperature*np.exp(eps))[0])/(2*eps)
        direct, _, _ = qrule.coefficients(m0*x, temperature, m0/v, p['Neta'])
        values.append(dict(x=x, z=z,
            stress_trace_error_over_rho=float((rho-3*pressure-trace)/rho),
            mass_exchange_identity_error_over_rho=float(
                (fd_time+3*(rho+pressure)-source*v*slope)/rho),
            potential_source_gradient_relative_error=float(fd[0]/(source*v)-1),
            potential_source_curvature_relative_error=float(fd[3]/(direct*v)-1)))
    return dict(mesh=candidate['mesh'], values=values)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    ref = Path(__file__).with_name('results')/'direct-subthreshold-friedmann.json'
    results = [audit(c) for c in json.loads(ref.read_text())]
    for row in results:
        for point in row['values']:
            for key, amplitude in point.items():
                if 'error' in key:
                    assert abs(amplitude) < 2e-7, (key, amplitude)
    args.output.write_text(json.dumps(results, indent=2, allow_nan=False)+'\n')
    print('Stress and exchange identities pass for 10 snapshots.')


if __name__ == '__main__':
    main()
