"""Cross-grid and finite-time strong convergence of physical four-reaction
five-population Dirac+antiχ1/2 and real Bose φ on-shell Yukawa kinetics.

This is conditional fixed vacuum Γ and finite parent cutoffs, no SM/FRW,
86-eV late PSD or condensate. Use same finite-dimensional quantum KL.
"""
import json
import math
import numpy as np
from kinetics.phi_bose_five_dirac_full_yukawa_implicit import FullDiracPhysicalMoments
from kinetics.phi_bose_five_dirac_kl_certificate import quantum_KL_certificate


def construct(Nphi=8,P=2.,nrad=4,nrelative=8,ncm=32,naz=4):
    return FullDiracPhysicalMoments(
        P=P,Ppair=P,N1=3,N2=3,Nphi=Nphi,nq=18,
        nrad=nrad,nangle=8,pair_nrad=nrad,
        pair_nrelative=nrelative,pair_ncm=ncm,pair_naz=naz)


def integrate_horizon(model,Tfinal,nominal_steps):
    eta=model.initial_etas()
    start=model.integrate_spectrum(eta)
    start_q=float(sum(start["Q"]))
    start_n=float(start["N"][4])
    initial=eta
    elapsed=0.
    naccepted=0; rejected=0
    step_checks=[]
    while elapsed<Tfinal*(1-1e-12):
        requested=min(Tfinal/nominal_steps,Tfinal-elapsed)
        info=model.adaptive_step(eta,requested)
        elapsed+=info["accepted_dt"]
        eta=info["etas"]
        naccepted+=1
        rejected+=len(info["rejections"])
        step_checks.append((info["entropy_increase"],
                            info["entropy_certificate"],
                            info["energy_defect"],info["U1_defect"]))
        if naccepted>4*nominal_steps+12:
            raise RuntimeError("Unbounded adaptive step rejection")
    final=model.integrate_spectrum(eta)
    cert=quantum_KL_certificate(model,initial,eta,Tfinal)
    return dict(Tfinal_GeV_inverse=Tfinal,steps=naccepted,rejected=rejected,
                delta_Nphi_GeV3=float(final["N"][4]-start_n),
                delta_entropy_GeV3=float(final["S"]-start["S"]),
                absolute_delta_total_energy=float(sum(final["Q"])-start_q),
                absolute_U1_defect=float(final["U1"]-start["U1"]),
                min_Bose_eta=float(min(eta[4])),
                final_source_H=model.source(eta)["H"],
                # This full trajectory is NOT a single backward-Euler step:
                # do not misuse the single-step KL identity for the whole run.
                max_step_energy_defect=float(max(abs(x[2]) for x in step_checks)),
                max_step_U1_defect=float(max(abs(x[3]) for x in step_checks)),
                min_step_entropy_gain=float(min(x[0] for x in step_checks)),
                min_step_entropy_certificate=float(min(x[1] for x in step_checks)),
                note="Global final KL certificate value not used; physical root KL applies to EACH BE step")


def audit():
    rows=[]
    for P in (2.,8.):
        for nphi in (8,12,16):
            m=construct(nphi,P)
            initial=m.source(m.initial_etas())
            records=[integrate_horizon(m,20000.,n) for n in (2,4,8)]
            dx=[v["delta_Nphi_GeV3"] for v in records]
            e1=dx[0]-dx[1]
            e2=dx[1]-dx[2]
            ratio=e1/e2 if e2 else None
            rows.append(dict(P_GeV=P,Nphi=nphi,nvars=m.nvars,events=m.total_events,
                             initial_H=initial["H"],
                             weak_event_charge_defect=initial["exact_event_U1"],
                             weak_event_energy_defect=initial["exact_event_energy"],
                             three_timestep_runs=records,
                             first_order_Richardson_difference_ratio=ratio))
    return dict(schema_version=1,
                scope="nonthermal five-population on-shell 4-reaction physical finite-step strong convergence ONLY; no SM/FRW/Ω/BEC/CLASS-ready",
                rows=rows,class_ready=False)


if __name__=="__main__": print(json.dumps(audit(),indent=2,allow_nan=False))
