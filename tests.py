import numpy as np
from induc_SPW import *

def test_group_ipsps():
    ipsps_start = [1,7, 2,7,8,2 ,5]
    electrode = np.arange(len(ipsps_start))+1
    ipsps = np.rec.fromarrays([ipsps_start,electrode], names='ipsp_start,electrode')
    group = group_ipsps(ipsps, 1)
    expected_group = np.array([0,1, 0, 1, 1, 0, 2])
    assert (group == expected_group).all(), "%s != %s" % (group, expected_group)

def test_group_ipsps_in_same_electrode():
    ipsps_start = [1, 2, 3]
    electrode =   [1, 1, 1]
    ipsps = np.rec.fromarrays([ipsps_start, electrode], 
                              names='ipsp_start,electrode')
    group = group_ipsps(ipsps, 2)
    expected_group = np.array([0, 1, 2])
    assert (group == expected_group).all(), "%s != %s" % (group, expected_group)

def test_group_ipsps_in_multi_electrode():
    ipsps_start = [2, 1, 3, 4, 5]
    electrode =   [1, 2, 2, 2, 3]
    ipsps = np.rec.fromarrays([ipsps_start, electrode], 
                              names='ipsp_start,electrode')
    group = group_ipsps(ipsps, 2)
    expected_group = np.array([0, 0, 1, 2, 2])
    assert (group == expected_group).all(), "%s != %s" % (group, expected_group)
    
def test_group_ipsps_in_multi_electrode_reverse_order():
    ipsps_start = [2, 3, 1, 4, 5]
    electrode =   [1, 2, 2, 2, 3]
    ipsps = np.rec.fromarrays([ipsps_start, electrode], 
                              names='ipsp_start,electrode')
    group = group_ipsps(ipsps, 2)
    expected_group = np.array([0, 1, 0, 2, 2])
    assert (group == expected_group).all(), "%s != %s" % (group, expected_group)

def test_group_ipsps_in_multi_electrode_unsorted_electrodes():
    ipsps_start = [2, 3, 1, 4, 5]
    electrode =   [2, 3, 3, 3, 1]
    ipsps = np.rec.fromarrays([ipsps_start, electrode], 
                              names='ipsp_start,electrode')
    group = group_ipsps(ipsps, 2)
    expected_group = np.array([1, 2, 1, 0, 0])
    assert (group == expected_group).all(), "%s != %s" % (group, expected_group)

def test_group_ipsps_assign_closest_element_to_existing_group():
    ipsps_start = [2, 2.5, 1]
    electrode =   [1, 2, 2]
    ipsps = np.rec.fromarrays([ipsps_start, electrode], 
                              names='ipsp_start,electrode')
    group = group_ipsps(ipsps, 2)
    expected_group = np.array([0,0,1])
    assert (group == expected_group).all(), "%s != %s" % (group, expected_group)
    
    ipsps_start = [2, 1, 2.5]
    electrode =   [1, 2, 2]
    ipsps = np.rec.fromarrays([ipsps_start, electrode], 
                              names='ipsp_start,electrode')
    group = group_ipsps(ipsps, 2)
    expected_group = np.array([0,1,0])
    assert (group == expected_group).all(), "%s != %s" % (group, expected_group)
    
    ipsps_start = [1, 2, 1, 2.5]
    electrode =   [0, 1, 2, 2]
    ipsps = np.rec.fromarrays([ipsps_start, electrode], 
                              names='ipsp_start,electrode')
    group = group_ipsps(ipsps, 2)
    expected_group = np.array([0, 0,1,0])
    assert (group == expected_group).all(), "%s != %s" % (group, expected_group)
    
    ipsps_start = [1, 2, 2.5, 1]
    electrode =   [0, 1, 2, 2]
    ipsps = np.rec.fromarrays([ipsps_start, electrode], 
                              names='ipsp_start,electrode')
    group = group_ipsps(ipsps, 2)
    expected_group = np.array([0, 0,0,1])
    assert (group == expected_group).all(), "%s != %s" % (group, expected_group)
    
def test_group_ipsps_maximum_one_ipsp_per_electrode():
    ipsps_start = np.cumsum(np.random.rand(1000))
    electrode =   np.random.randint(6, size=1000)
    ipsps = np.rec.fromarrays([ipsps_start, electrode], 
                              names='ipsp_start,electrode')
    group = group_ipsps(ipsps, 1)
    assert np.max(np.bincount(group))==6

def test_group_ipsps_is_stable():
    ipsps_start = np.random.randint(100,size=1000)*1.
    electrode =   np.random.randint(6, size=1000)
    ipsps = np.rec.fromarrays([ipsps_start, electrode], 
                              names='ipsp_start,electrode')
    group1 = group_ipsps(ipsps, 1)
    group2 = group_ipsps(ipsps, 1)
    assert (group1==group2).all()

def test_group_ipsps_check_if_skips_over_some_ipsps():
    ipsps_start = [1, 4, 3, 6]
    electrode =   [0, 0, 1, 0]
    ipsps = np.rec.fromarrays([ipsps_start, electrode], 
                              names='ipsp_start,electrode')
    group = group_ipsps(ipsps, 3)
    expected_group = np.array([0, 1, 1, 2])
    assert (group == expected_group).all(), "%s != %s" % (group, expected_group)

def test_add_rec_field():
    x = np.rec.fromarrays([[1,2], [3,4]], names='a,b')
    c = np.array([5,6])
    y = add_rec_field(x, [c], ['c'])
    assert (y['c']==c).all()

def test_count_coincident_spikes():
    ipsps_start = [1, 2, 3, 5, 6, 10]
    electrodes = [1, 1, 2, 3, 4, 3]
    
    ipsps = np.rec.fromarrays([ipsps_start, electrodes], names='ipsp_start,electrode')
    n_electr = count_coincident_ipsps(ipsps, 1.)
    expected_n_electr = np.array([1,2,2,2,2, 1])
    assert (n_electr==expected_n_electr).all()
    
def test_count_coincident_spikes_not_sorted():
    ipsps_start = [1, 6, 3, 5, 2, 10]
    electrodes = [1, 4, 2, 3, 1, 3]

    ipsps = np.rec.fromarrays([ipsps_start, electrodes], names='ipsp_start,electrode')
    
    #shuffle randomly and store the order
    idx = np.array([0, 3, 5,2,4,1])
    ipsps = ipsps[idx]
    ipsps_copy = ipsps.copy()
    n_electr = count_coincident_ipsps(ipsps, 1.)
    expected_n_electr = np.array([1,2,2,2,2, 1])
    
    assert (n_electr==expected_n_electr[idx]).all(), "%s != %s" % (n_electr, expected_n_electr[idx])
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
    
    expected_ampls = np.array([4,1,0, 0, 1])

    assert (ampls==expected_ampls).all(), "%s != %s" % (ampls, expected_ampls)
    
def test_calculate_amplitude_of_IPSP():
    electrode = [1, 2, 1, 2, 1]
    data = np.vstack([np.zeros(10), np.arange(10), np.arange(10)])
    ipsp_start = [1, 2, 5, 4, 7]
    fs = 1000
    
    ipsps = np.rec.fromarrays([ipsp_start, electrode],
                              names='ipsp_start,electrode')
    
    ampls = calculate_amplitude_of_IPSP(ipsps, data, fs)
    expected_maxs = np.array([3, 1,1, 0, 0])
    assert  (ampls==expected_maxs).all(), "%s != %s" % (ampls, expected_maxs)
    
    
def test_shift_ipsp_to_closest_spike():
    spike_time = [1, 5, 10]
    spike_electrode = [1,1, 1]
    ipsp_start = [2.,3.]
    ipsp_electrode = [1,2]
    ipsp_amplitude = [1,1]
    ipsp_spw_no = [1,1]
    spike_spw_no = [1,1,1]
    group_ipsps = [1,1]
    
    spikes_trace = np.rec.fromarrays([spike_electrode, spike_time, spike_spw_no], 
                                    names='electrode,time,spw_no')
    ipsps_trace = np.rec.fromarrays([ipsp_start, ipsp_electrode, ipsp_amplitude,
                                     ipsp_spw_no, group_ipsps],
                                    names='ipsp_start,electrode,amplitude,spw_no,group')

    shift_spike = 2
    shifted_ipsp = shift_ipsp_start(ipsps_trace, spikes_trace, shift_spike)
    expected_time = np.array([1,1])
    time = shifted_ipsp['ipsp_start']
    assert (time==expected_time).all(), "%s != %s" % (time, expected_time)
    
    
def test_shift_ipsp_to_closest_spike_after():
    spike_time = [4.5, 8, 10]
    spike_electrode = [2,2, 1]
    ipsp_start = [2.,3.]
    ipsp_electrode = [1,2]
    ipsp_amplitude = [1,1]
    ipsp_spw_no = [1,1]
    spike_spw_no = [1,1,1]
    group_ipsps = [1,1]
    
    spikes_trace = np.rec.fromarrays([spike_electrode, spike_time, spike_spw_no], 
                                    names='electrode,time,spw_no')
    ipsps_trace = np.rec.fromarrays([ipsp_start, ipsp_electrode, ipsp_amplitude,
                                     ipsp_spw_no, group_ipsps],
                                names='ipsp_start,electrode,amplitude,spw_no,group')

    shift_spike = 2
    shifted_ipsp = shift_ipsp_start(ipsps_trace, spikes_trace, shift_spike)
    expected_time = np.array([4.5,4.5])
    time = shifted_ipsp['ipsp_start']
    assert (time==expected_time).all(), "%s != %s" % (time, expected_time)