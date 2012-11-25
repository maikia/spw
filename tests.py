import numpy as np
from induc_SPW import count_coincident_ipsps

def test_count_coincident_spikes():
    ipsps_start = [1, 2, 3, 5, 6, 10]
    electrodes = [1, 1, 2, 3, 4, 3]
    
    ipsps = np.rec.fromarrays([ipsps_start, electrodes], names='ipsp_start,electrode')
    n_electr = count_coincident_ipsps(ipsps, 1.)
    expected_n_electr = np.array([2,2,2,2,2, 1])
    print n_electr
    assert (n_electr==expected_n_electr).all()
    
    