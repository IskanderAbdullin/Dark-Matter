"""Positive heavy-forming exposure on passive redshifted, archived supports."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from archived_heavy_recycling_gate_reference import vacuum_upfold
from archived_scalar_pair_reaction_reference import cropped, gross, snapshot_functions
from late_radiation_background_reference import TabulatedRadiation, MPL
from three_species_decay_onsager_reference import m, M, mp


CONTROLS = (
    (4, (3, 8, 8, 3), (4, 12), True),
    (6, (4, 10, 10, 3), (6, 18), True),
    (8, (5, 12, 12, 4), (7, 21), True),
    (10, (7, 16, 24, 16), (9, 27), False),
)


def redshift(occupation, mass, b):
    # Divide lab momentum edges by b and evaluate the original distribution
    # at b times the current momentum; do not refit a two-moment shape.
    def evaluate(energy):
        p_now = np.sqrt(np.maximum((energy - mass) * (energy + mass), 0.))
        return occupation(np.hypot(mass, b * p_now))
    return evaluate


def run_case(table, name, s, normalization):
    edges0 = np.asarray(s["lab_momentum_edges_GeV"])
    T = s["T_visible_GeV"]
    radiation = TabulatedRadiation(table, T)
    n_local = np.sum(s["represented_species_count_a3"]) / s["scale_factor"]**3
    f, _, g, _ = snapshot_functions(s)
    b_phi, b_chi = [edges0[-1] / np.sqrt(M*M - species**2) for species in (mp, m)]

    def sample(b, phi_order, chi_order, split_final):
        H = np.sqrt(radiation.state(b)[1] / 3) / MPL
        edges = edges0 / b
        phi = chi = 0.
        if b < b_phi:
            phi = gross(edges, mp, M, redshift(f, mp, b), np.zeros_like,
                        True, *phi_order, split_final=split_final,
                        exchange_mass=m)[1][1]
        if b < b_chi:
            chi = vacuum_upfold(edges, redshift(g, m, b), *chi_order)
        # Two heavy parents per event; one extra particle per parent if
        # all of them decay promptly into two dark daughters.
        return 2*b**3 / (normalization*H) * np.array([phi, chi])

    panels = np.unique(np.clip(np.log([1., 1.12, 1.35, 1.65, 2., 2.5, 2.8,
                                     max(1., b_phi), max(1., b_chi)]),
                              0., np.log(max(1., b_chi))))
    assert np.array_equal(sample(max(1., b_chi), *CONTROLS[0][1:]), [0., 0.])
    rows = []
    for n, phi_order, chi_order, split in CONTROLS:
        x, w = roots_legendre(n)
        d = np.diff(panels)[:, None] / 2
        centers = (panels[1:, None] + panels[:-1, None]) / 2
        points = np.exp((centers + d*x).ravel())
        rates = np.array([sample(b, phi_order, chi_order, split) for b in points])
        segments = ((d*w)[:, :, None] * rates.reshape(-1, n, 2)).sum(axis=1)
        integrals = segments.sum(axis=0)
        local_rate = sample(1., phi_order, chi_order, split)
        assert np.all(np.isfinite(integrals)) and np.all(integrals >= 0.)
        rows.append(dict(log_b_order=n, phi_order=phi_order, chi_order=chi_order,
                         split_outgoing_panels=split,
                         integrated_phi_and_opposite_stable=integrals.tolist(),
                         anchor_phi_and_opposite_stable_per_log_b=local_rate.tolist(),
                         optimistic_including_same_charge=float(integrals @ [1., 3.]),
                         contribution_by_panel=(segments @ [1., 3.]).tolist()))
    if rows[-1]["optimistic_including_same_charge"] > 0.:
        assert abs(rows[-2]["optimistic_including_same_charge"] /
                   rows[-1]["optimistic_including_same_charge"] - 1.) < .005
    return dict(snapshot=name, anchor_T_GeV=T, endpoint_b=[b_phi, b_chi],
                retained_number_at_anchor_GeV3=float(n_local),
                normalization_at_anchor_GeV3=float(normalization),
                b_panels=np.exp(panels).tolist(), controls=rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos")
    p.add_argument("paths", nargs="+", help="Fixed archived input snapshot JSON")
    p.add_argument("--crop-last", type=float)
    a = p.parse_args()
    rows = []
    for path in a.paths:
        s = json.loads(Path(path).read_text())["snapshots"][-1]
        n = np.sum(s["represented_species_count_a3"]) / s["scale_factor"]**3
        rows.append(run_case(a.eos, Path(path).name, s, n))
    if a.crop_last:
        kept = rows[-1]["controls"][-1]["optimistic_including_same_charge"]
        label = f"{Path(a.paths[-1]).name}: same-state crop to q/T_R={a.crop_last:g}"
        rows.append(run_case(a.eos, label, cropped(s, a.crop_last), n))
        assert rows[-1]["controls"][-1]["optimistic_including_same_charge"] <= kept
    print(json.dumps(dict(scope="Pauli-free positive formation evaluated on passive "
                          "redshifting represented shapes and an external radiation clock. "
                          "Not an evolving signed rate or a bound with redistribution.",
                          cases=rows), indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
