import numpy as np
from induc_SPW import *

def test_count_coincident_spikes():
    ipsps_start = [1, 2, 3, 5, 6, 10]
    electrodes = [1, 1, 2, 3, 4, 3]
    
    ipsps = np.rec.fromarrays([ipsps_start, electrodes], names='ipsp_start,electrode')
    n_electr = count_coincident_ipsps(ipsps, 1.)
    expected_n_electr = np.array([2,2,2,2,2, 1])
    assert (n_electr==expected_n_electr).all()
    
def test_count_coincident_spikes_not_sorted():
    ipsps_start = [1, 2, 3, 5, 6, 10]
    electrodes = [1, 1, 2, 3, 4, 3]

    ipsps = np.rec.fromarrays([ipsps_start, electrodes], names='ipsp_start,electrode')
    
    #shuffle randomly and store the order
    idx = np.argsort(np.random.randn(len(ipsps)))
    ipsps = ipsps[idx]
    ipsps_copy = ipsps.copy()
    n_electr = count_coincident_ipsps(ipsps, 1.)
    expected_n_electr = np.array([2,2,2,2,2, 1])
    
    assert (n_electr==expected_n_electr[idx]).all()
    #check whether inputs array was not changed
    assert (ipsps==ipsps_copy).all()

def test_calculate_ipsp_rise_single_electrode():
    electrode = np.ones(4, dtype=np.int32)
    spw_no = np.zeros(4)
    ipsp_start = np.arange(4)
    data = np.vstack([np.zeros(5), np.arange(5)])
    Fs = 1000
    
    ipsps = np.rec.fromarrays([ipsp_start, spw_no, electrode],
                              names='ipsp_start,spw_no,electrode')
    ampls = calculate_ipsp_rise(ipsps, data, Fs)
    
    expected_ampls = np.ones(4)
    expected_ampls[-1] = 0

    assert (ampls==expected_ampls).all(), "%s != %s" % (ampls, expected_ampls)

def test_calculate_ipsp_rise_multi_electrodes():
    electrode = [2, 1, 2, 1, 1]
    spw_no = np.zeros(len(electrode))
    ipsp_start = np.concatenate([np.arange(2), np.arange(3)[::-1]+2])
    data = np.vstack([np.zeros(5), np.arange(5), np.arange(5)])
    Fs = 1000
    
    ipsps = np.rec.fromarrays([ipsp_start, spw_no, electrode],
                              names='ipsp_start,spw_no,electrode')
    ampls = calculate_ipsp_rise(ipsps, data, Fs)
    
    expected_ampls = np.array([1,1,0, 4, 0])

    assert (ampls==expected_ampls).all(), "%s != %s" % (ampls, expected_ampls)
    
def test_calculate_amplitude_of_IPSP():
    #????
    electrode = [1, 2, 1, 2, 1]
    spw_no = np.zeros(len(electrode))
    spw_ipsps_trace = np.concatenate([np.arange(2), np.arange(3)[::-1]+2])
    data = np.vstack([np.arange(5, 10), np.arange(5), np.arange(5)])
    fs = 1000
    
    ipsps = np.rec.fromarrays([ipsp_start, spw_no, electrode],
                              names='ipsp_start,spw_no,electrode')
    
    ampls = calculate_amplitude_of_IPSP(ipsps, data, Fs)
    expected_maxs = np.array([1, 1,0, 4, 0])
    print ampls
    assert  (ampls==expected_ampls).all(), "%s != %s" % (ampls, expected_ampls)