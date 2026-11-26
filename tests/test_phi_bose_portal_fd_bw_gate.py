from kinetics.phi_bose_portal_fd_bw_gate import ratio,gain

def test_swapped_mb_and_fd():
    a=ratio(.6882487043,.012,False,16)
    b=ratio(.6882487043,.012,True,12)
    c=ratio(.6882487043,.012,True,16)
    assert abs(a["ratio"]-1.0002046815984291)<1e-7
    assert abs(b["ratio"]-c["ratio"])<1e-8
    assert abs(c["ratio"]-1.000193199550118)<1e-7

def test_reduced_fd_gain():
    r=gain(5.,.0001)
    assert 0<r["upper"]<3e-5
