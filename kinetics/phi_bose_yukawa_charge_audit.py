"""Symmetry/stoichiometry scope for real-phi off-diagonal Dirac Yukawa.

Interaction: L_int = -y*phi*bar(chi2)*chi1 + h.c., both chi1/chi2 Dirac,
phi REAL. A possible stability Z2 assigns odd chi1,phi and even chi2.

Isolated 1<->2 collision has two positive-count invariants when
particle/antiparticle channels are lumped. But the same Yukawa vertex also
allows chi1+anti(chi1)<->phi+phi via chi2 t/u exchange, and
chi2+anti(chi2)<->phi+phi via chi1 t/u exchange. These order-y^4
processes can break BOTH old positive count invariants; only the
underlying vector-like net fermion U(1) survives as linear continuous
charge. No cross-sections or cosmological relevance are computed here.
"""
import json
import numpy as np

SPECIES=("chi1","anti_chi1","chi2","anti_chi2","phi")
CHANNELS=(
    "chi1 -> chi2 + phi",
    "anti_chi1 -> anti_chi2 + phi",
    "chi1 + anti_chi1 -> phi + phi [virtual chi2]",
    "chi2 + anti_chi2 -> phi + phi [virtual chi1]",
)

def stoichiometry(include_pair=True):
    """Rows=net final minus initial populations, columns=SPECIES."""
    S=np.array([
        [-1, 0,+1, 0,+1],
        [ 0,-1, 0,+1,+1],
        [-1,-1, 0, 0,+2],
        [ 0, 0,-1,-1,+2],
    ],dtype=np.int64)
    return S if include_pair else S[:2]


def analyze():
    S=stoichiometry(True)
    T=stoichiometry(False)
    U1=np.array([1,-1,1,-1,0],dtype=np.int64)
    positive_fermion_count=np.array([1,1,1,1,0],dtype=np.int64)
    old_number_N1_plus_phi=np.array([1,1,0,0,1],dtype=np.int64)
    z2_odd_count=np.array([1,1,0,0,1],dtype=np.int64)
    full_minor=int(round(np.linalg.det(S[:,[0,1,2,4]])))
    return dict(schema_version=1,
                species=SPECIES,channels=CHANNELS,
                decays_only_matrix=T.tolist(),
                full_matrix=S.tolist(),
                decays_only_rank=int(np.linalg.matrix_rank(T)),
                decays_only_nullity=T.shape[1]-int(np.linalg.matrix_rank(T)),
                full_rank=int(np.linalg.matrix_rank(S)),
                full_nullity=S.shape[1]-int(np.linalg.matrix_rank(S)),
                exact_full_rank_four_minor=full_minor,
                fundamental_vector_U1=U1.tolist(),
                full_U1_reaction_changes=(S@U1).tolist(),
                old_total_fermion_reaction_changes=(S@positive_fermion_count).tolist(),
                old_N1_plus_phi_reaction_changes=(S@old_number_N1_plus_phi).tolist(),
                z2_parity_changes_mod2=(S@z2_odd_count%2).tolist(),
                caveat="Allowed additional Yukawa pair topologies only; nonzero amplitudes and actual thermal rates must be calculated separately, no full SM/FRW nor extra channels proved negligible",
                class_ready=False)


if __name__=="__main__":
    print(json.dumps(analyze(),indent=2))
