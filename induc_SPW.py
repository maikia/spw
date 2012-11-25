import scipy.signal as signal

from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
import os
import spike_detection as fes
import numpy as np
from matplotlib import pyplot
import pylab as plt
import filtit as filt
#import plotting as plt
import scipy as sc
import detect_waves as dw
from scipy import signal
import data_mang as dat
import folder_manager as fold_mang

def find_envelope(data):
    """ finds envelope of the data"""
    x_hilbert = signal.hilbert(data)
    x_envelope = np.abs(x_hilbert) # creates envelope of the data
    return x_envelope

def ms2pts(ms, fs):
    """ converts number of ms to number of data pts"""
    pts = (fs / 1000.0) * ms
    return pts

def pts2ms(pts, fs):
    """ converts number of pts to number of ms"""
    ms = (pts/fs)*1000.0
    return ms

def detect_spikes(data, thres):
    """ detects all the waves above thres in the data, 
    returns indices of all the found waves """
    above = dw.find_above(data, thres)
    starts, ends = dw.find_startend(above)
    maxs, indcs = dw.max_waves_bartex(data, starts, ends)
    return indcs

def detect_1spike(data, thres, fs, pulse_len = 500):
    """ it detects all the intra spikes and then takes only those which are 
    not closer to each other than pulse_len
    pulse_len is given in ms"""
    spikes = detect_spikes(data, thres)
    pulse_len = ms2pts(pulse_len, fs)
    firsts = []
    if len(spikes) > 0:
        dist = [int(pulse_len+1)]
        dist.extend(np.diff(spikes).tolist())
        firsts = [spikes[i] for i in range(len(dist)) if dist[i] > pulse_len]
    return firsts

def dist_fromSpike(spikes, events, fs):
    """ calculates distance (in points) from
    each event given to last preceding spike"""
    distances = []
    
    for i in range(len(events)):
        ev = events[i]
        sp = np.where(spikes < ev)[0]
        if np.size(sp) > 0:
            sp = spikes[sp[-1]]
            dist = ev - sp
            distances.append(pts2ms(dist, fs))
        #elif np.size(spikes) > 0:
        else:    
            distances.append(float('Inf'))
        
        #print sp
#        if np.size(spikes) == 1:
#            dist = spikes
#            dist = ev-dist
#            distances.append(pts2ms(dist, fs))  
#        elif len(sp) > 0:
#            sp = sp[-1]
#            dist = spikes[sp]
#            dist = ev - dist
#            distances.append(pts2ms(dist, fs))    
#        else:
#            distances.append(float('Inf'))
    return distances

def find_locals(data):
    """ finds all the local maxs and mins"""
    
    change = np.diff(data)
    change[change > 0] = 1
    change[change < 0] = 0

    shift = [0]
    shift.extend(change[:-1])
    subs = change - shift
    
    local_mins = np.where(subs == 1)
    local_maxs = np.where(subs == -1)
    return local_maxs, local_mins

def assign_startend(spw_start, max_idxs):
    """ given list of all the minimums, and list of all the maximums
    find which of the minimums form starts and ends of the waves 
    with peaks at maximum"""
    
    starts = []
    ends = []
    for x in range(len(max_idxs)):
        s = spw_start[0][spw_start[0] < max_idxs[x]]
        e = spw_start[0][spw_start[0] > max_idxs[x]]
        if len(s) > 0 and len(e) > 0:
            starts.append(s[-1])
            ends.append(e[0])
    return starts, ends

def cntRight_rip_in_SPW(starts, ends, ripples, fs, min_dist = 2, max_dist = 80):
    """ calculates proper inter ripple intervals and how many of them 
    are in each wave defined by starts and ends, with ripples being
    all ripple indexes; min_dist and max_dist are given in ms"""

    min_dist = ms2pts(min_dist, fs)
    max_dist = ms2pts(max_dist, fs)
    
    iris = []
    iris_no = []
    for x in range(len(starts)):
        # take tall the ripples for the given wave
        s = starts[x]
        e = ends[x]
        rips = [ripples[o] for o in range(len(ripples)) if (ripples[o] < e and ripples[o]> s)]

        # calculate inter ripple intervals
        if len(rips)> 1:
            iri = np.diff(rips)
            # remove all the ripples which are too far or too close to each other
            iri = [iri[o] for o in range(len(iri)) if iri[o] >= min_dist and iri[o] <= max_dist]
            iris_no.append(len(iri))
            iris.append(iri)
            
        else:
            iris.append(0)
    return iris, iris_no
    
def  find_SPWs(data, fast_data, fs, min_length):
    """ finds all the SPWs in the given data, no smaller than min_length (given in ms)
    it returns all the indexes of the the spws (at their highest amplitude) 
    and their amplitudes"""
    min_length = ms2pts(min_length, fs)
    
    # find SPWs using envelope - hilbert, over fast data and threshold
    fast_databove = fast_data[:]
    fast_databove[fast_databove < 0] = 0
    x_envelope = find_envelope(fast_databove)
    x_envelope = find_envelope(x_envelope)
    
    # find SPWs based on the selected variables and found waves (above)
    thres = np.std(x_envelope) * 2
    waves = dw.find_above(x_envelope, thres) 
    starts, ends = dw.find_startend(waves) # finds all the peaks above threshold
    starts, ends = dw.correct_wave(data, starts, ends, fs) # all too short waves are either concatenated or removed
    maxs, max_idxs = dw.max_waves(data, starts, ends) # finds maxs of each wave

    # take only those waves which are above thres in data_bas (baselined data)
    temp_idxs = np.where(data[np.array(max_idxs)] > 0)
    temp_idxs, = temp_idxs
    max_idxs = np.asarray(max_idxs)
    max_idxs = max_idxs[temp_idxs]
    
    # take only those waves which are long enough
    max_idxs = sorted(max_idxs)
    max_idxs = dw.concat_waves(data, max_idxs, fs)
    lengths = dw.check_width(data, max_idxs, data[max_idxs] * 0.20, fs) # check length of the waves right above 0 baseline
    max_idxs = max_idxs[lengths >= min_length]
    max_ampl = data[max_idxs] 
    return max_idxs, max_ampl



def round_spike_idxs(spike_idxs):
    s_p = []
    for i in range(len(spike_idxs)):
        s_p2 = []
        for j in range(len(spike_idxs[i])):
            s_p1 = []
            for k in range(len(spike_idxs[i][j])):
                s_p1.append(int(round(spike_idxs[i][j][k])))
            s_p2.append(s_p1)
        s_p.append(s_p2)
    return s_p



def update_spw_ripple(starts, ends, spw_idx, save_folder, data_file = 'spw_ripple'):
    """updates the analysis on spw- ripple relation; it saves the indexes of ripples for each SPW, 
    it gives the inter ripple interval if it's not too small and not too big"""
    # find number for ripples in the given spw
    min_d = 2 #ms
    max_d = 80 #ms
    iris_all = []
    iris_no_all = []
    print
    print 'looking for relation between spws and ripples, working on electrode:',         
    for electr in range(np.size(starts, 0)):
        print electr+1,
        i_a_all = []
        i_n_all = []
        for trace in range(len(starts[electr])):
            i_a = []
            i_n_a = []
            for trace in range(len(starts[electr])):
                iris, iris_no = cntRight_rip_in_SPW(starts[electr][trace], ends[electr][trace], spw_idx[electr][trace], min_d, max_d)
                i_a.append(iris)
                i_n_a.append(iris_no)
            i_a_all.append(i_a)
            i_n_all.append(i_n_a)     
        iris_all.append(i_a_all)
        iris_no_all.append(i_n_all)
     
    np.savez(save_folder + data_file, iris_all, iris_no_all) 
    return iris_all, iris_no_all

    
def update_ripples(ripple_data, fs, save_folder, data_file = 'ripples'):
    """ not tested"""
    # find ripples
    rip_idxs = []
    print
    print 'looking for ripples, working on electrode:',  
    for electr in range(len(ripple_data)):
        lmax_idxs_all = []
        print electr+1,
        for trace in range(len(ripple_data[electr])):
            lmax_idxs, lmin_idxs = find_locals(ripple_data[electr][trace])
            #set above std
            rip_thres = np.std(ripple_data[electr][trace])
            lmax_idxs = [lmax_idxs[0][x] for x in range(len(lmax_idxs[0])) if ripple_data[electr][trace][lmax_idxs[0][x]]> rip_thres]            
            lmax_idxs_all.append(lmax_idxs)
        rip_idxs.append(lmax_idxs_all)
        
    np.savez(save_folder + data_file, rip_idxs, fs) 
    return rip_idxs, fs



def update_intraSpikes(save_folder, save_file, load_file, pulse_len = 500):
#(data, fs, save_folder, save_file = "intra_spikes", pulse_len = 500, ):
    """ pulse_len in ms - length of the stimulation pulse in intracellular electrode"""
    npzfile = np.load(save_folder + load_file)
    data = npzfile['data']
    fs = npzfile['fs'] 
    npzfile.close()
    
    
    sp_all = []
    sp_first = []
    
    print 'Detecting intracellular spikes'
    for trace in range(np.size(data,1)):
        spiking_thres = -10 #mV

        # detect only the first spike in the row 
        sp_idx_first = detect_1spike(data[0, trace, :], spiking_thres, fs, pulse_len)
        sp_idx_a = detect_spikes(data[0, trace, :], spiking_thres) # for detecting all the spikes 
        
        
        typ = 'f8'
        spike_idxs_first = pts2ms(sp_idx_first, fs).astype(typ)
        electrodes_first = np.zeros(len(spike_idxs_first), dtype=np.int32)
        traces_first = np.ones(len(spike_idxs_first), dtype=np.int32)*trace
        
        sp_idxs_all = pts2ms(sp_idx_a, fs).astype(typ)
        electrodes_all = np.zeros(len(sp_idxs_all), dtype=np.int32)
        traces_all = np.ones(len(sp_idxs_all), dtype=np.int32)*trace
            
        sp_all.append(np.rec.fromarrays([electrodes_all, traces_all, np.array(sp_idxs_all, dtype=typ)], names='electrode,trace,time'))
        sp_first.append(np.rec.fromarrays([electrodes_first, traces_first, np.array(spike_idxs_first, dtype=typ)], names='electrode,trace,time'))
#        import pdb; pdb.set_trace()
#        part = [0, 200000]
#        data_temp = data[0,0,part[0]:part[1]]
#        
#        t = dat.get_timeline(data_temp, fs, 'ms')
#        plt.plot(t, data_temp)
#        
#        spik = np.concatenate(sp_first)
#        spiks = ms2pts(spik['time'], fs).astype('i4')
#        spiks_temp = spiks[(spiks < part[1]) & (spiks > part[0])]
#        plt.plot(t[spiks_temp - part[0]], data_temp[spiks_temp - part[0]], 'go')
        
        
        
    sp_all = np.concatenate(sp_all)
    sp_first = np.concatenate(sp_first)    
       
    np.savez(save_folder + save_file, spikes_first = sp_first, spikes_all = sp_all, fs = fs) 
    

def load_create(folder_save, filename_save, freq, fs, data, N = 1000):
    # checks if the given folder/file exists and if not calculates the data
    if not os.path.exists(folder_save):
        os.makedirs(folder_save)
    
    try:
        # check if high pass data file already exists for this electrode
        npzfile = np.load(folder_save + filename_save + '.npz')
        data_filt = npzfile['data']
        fs = npzfile['fs']
        npzfile.close()
        
    except IOError as e:
        # not? create it!
        if np.size(freq) == 1:
            data_filt = filt.highPass(freq, fs, data, N)
        else: 
            data_filt = filt.bandPass(freq, fs, data, N) 
        # save it
        
        fold_mang.create_folder(folder_save)
        np.savez(folder_save + filename_save, data = data_filt, fs = fs)
        
    return data_filt, fs


def update_extraspikes(data_load, filter_folder, save_folder, save_file = "ex_spikes", save_filter = 'fast_data_'):
    """ finds and updates the detection of extracellular spikes"""
    
    npzfile = np.load(save_folder + data_load)
    data = npzfile['data']
    fs = npzfile['fs'] 
    npzfile.close()
       
    freq_fast = 750.
    # find extracellular spikes
    folder_name = save_folder + filter_folder
    N = 1000
    print
    print "finding extracellular spikes, working on electrode:",
    idx_all = []
    ampl_all = []
    
    for electr in range(np.size(data,0)):
        print electr,  

        if len(data[electr]) > 1:
            N = 100
        for trace in range(np.size(data,1)):
            data_used = data[electr, trace, :]
            filename_fast = save_filter +str(freq_fast) + '_'+ str(electr) + "_" + str(trace)
            data_fast, fs = load_create(folder_name, filename_fast, freq_fast, fs, data_used, N)
            spike_ampl, spike_idxs = fes.find_extra_spikes(data_used, data_fast, fs) 
            del data_fast
            typ = 'f8'
            spike_idxs = pts2ms(spike_idxs, fs).astype(typ)
            electrodes = np.ones(len(spike_idxs), dtype=np.int32)*electr
            traces = np.ones(len(spike_idxs), dtype=np.int32)*trace    
            idx_all.append(np.rec.fromarrays([electrodes, traces, spike_idxs,  np.array(spike_ampl, dtype=typ)], names='electrode,trace,time, amplitude'))
            del data_used
    
        
    #ampl_all = np.concatenate(ampl_all)         
    idx_all = np.concatenate(idx_all)            
    
    #np.savez(save_folder + save_file, spike_idx = idx_all, spike_ampl = ampl_all, fs = fs)
    np.savez(save_folder + save_file, spike_idx = idx_all, fs = fs)
    del data, idx_all, ampl_all, fs


def update_expikes_params(load_datafile, load_spikefile, save_folder, save_file = "ex_sparamas"):
    """calculate following variables for each spike:
    a =alley to peak
    b =half valley width
    c =half peak width
    
    """
    npzfile = np.load(save_folder + load_datafile)
    data = npzfile['data']
    fs = npzfile['fs'] 
    npzfile.close()
    
    npzfile = np.load(save_folder + load_spikefile)
    spike_idxs = npzfile['spike_idx']
    npzfile.close()

    print
    print "finding parameters of extracellular spikes, working on electrode:", 
    
    f_d = 'filtered/'
    folder_name = save_folder + f_d
    freq_fast = 500.
    freq_smooth = [500, 4000]
    N = 1000
    
    valley_to_peak_norm = []
    half_valley_width_norm = []
    left_most = []
    right_most = []
    norm_factors = []
    ampls = []    
    
    for electr in range(np.size(data,0)):
        print electr,

        if len(data[electr]) > 1:
            N = 100

        for trace in range(np.size(data,1)): 
            
            data_used = data[electr, trace, :]
            spikes_used = spike_idxs[(spike_idxs['electrode'] == electr) & (spike_idxs['trace'] == trace)]['time']
            spikes_used = ms2pts(spikes_used, fs)
            
            filename_fast = 'fast_data' + str(electr) + "_" + str(trace)
            filename_smooth = 'smooth_data' + str(electr) + "_" + str(trace)            
            data_fast, fs = load_create(folder_name, filename_fast, freq_fast, fs, data_used, N)
            data_smooth, fs = load_create(folder_name, filename_smooth, freq_smooth, fs, data_used, N)
            
            val_to_peak, half_val_width, ampl, As, Bs, norm_factor =fes.find_halfampl(data_used, data_fast, data_smooth, fs, spikes_used) #rang = [1, 2], fast_freq = 500.) #, fast_freq = 750.)

            v1, v2 = [], []
            for change in range(len(val_to_peak)):

                v1.append(pts2ms(val_to_peak[change], fs))
                v2.append(pts2ms(half_val_width[change], fs))
            
            
            electrodes = np.ones(len(v1), dtype=np.int32)*electr
            traces = np.ones(len(v1), dtype=np.int32)*trace                 
            
            typ = 'f8'
            valley_to_peak_norm.append(np.rec.fromarrays([electrodes, traces, np.array(v1, dtype=typ)], names='electrode,trace,time'))   
            half_valley_width_norm.append(np.rec.fromarrays([electrodes, traces, np.array(v2, dtype=typ)], names='electrode,trace,time'))         
            ampls.append(np.rec.fromarrays([electrodes, traces, np.array(ampl, dtype=typ)], names='electrode,trace,time'))         
            left_most.append(np.rec.fromarrays([electrodes, traces, np.array(As, dtype=typ)], names='electrode,trace,time'))         
            right_most.append(np.rec.fromarrays([electrodes, traces, np.array(Bs, dtype=typ)], names='electrode,trace,time'))         
            norm_factors.append(np.rec.fromarrays([electrodes, traces, np.array(norm_factor, dtype=typ)], names='electrode,trace,time'))   

                                 
    valley_to_peak_norm = np.concatenate(valley_to_peak_norm)   
    half_valley_width_norm = np.concatenate(half_valley_width_norm)   
    ampls = np.concatenate(ampls)     
    left_most = np.concatenate(left_most)   
    right_most = np.concatenate(right_most)   
    norm_factors = np.concatenate(norm_factors)   
        
    np.savez(save_folder + save_file, spike_idxs = spike_idxs
             , valley2peak = valley_to_peak_norm, halfValleyWidht = half_valley_width_norm
             , left_most = left_most, right_most = right_most, ampls =  ampls, fs = fs)
   
    del data, spike_idxs

def update_filtered(data, fs, save_folder, freq, data_file):
    """ updates filtered data files"""
    
    freq_data = []
    
    print
    print 'filtering the data ',
    print freq,
    print 'working on electrode:',
    
    data_filtered = []
    
    for electr in np.unique(data['electrode']): #range(1, len(data)+1):    
        print electr,    
        for trace in np.unique(data[(data['electrode'] == electr)]['trace']):
            
            #import pdb; pdb.set_trace()
            data_used = data[(data['electrode'] == electr) & (data['trace'] == trace)]['time']
            
            lowcut = 500.0
            highcut = 1250.0
            N = 1000 # number of times it loops around
            # filter the data
            if freq[0] == -1:
                # low pass filter
                freq_dat = filt.lowPass(freq[1], fs, data_used, N)
            elif freq[1] == -1:
                # high pass filter
                freq_dat = filt.highPass(freq[0], fs,data_used, N)
            else:
                # band pass filter
                freq_dat = filt.bandPass(freq, fs, data_used, N)
                
            #import pdb; pdb.set_trace()

            #for trace in np.unique(data[data['electrode'] == electr]['trace']):
            dat = freq_dat
            electrodes = np.ones(len(dat), dtype=np.int32)*electr
            traces = np.ones(len(dat), dtype=np.int32)*trace
        
            data_filtered.append(np.rec.fromarrays([electrodes, traces, dat], names='electrode,trace,time'))
            
            data_filtered = np.concatenate(data_filtered)
        
            plt.figure()
            plt.plot(data[(data['electrode'] == electr) & (data['trace'] == trace)]['time'], 'b')    
            plt.plot(data_filtered[(data_filtered['electrode'] == electr) & (data_filtered['trace'] == trace)]['time'], 'r')
            plt.show()
    
         
    np.savez(save_folder + data_file, data = data_filtered, freq = freq, fs = fs)  
    return freq_data, freq, fs
   
def update_datafile(filename, ex_electr, save_folder, data_file = 'data'):
    """ reads given file and saves the data read into rec array (numpy file format)"""

    data_all = []
    # do the same for every electrode given - read data
    all_data, no_segments = dat.get_data(filename) 
    
    print
    print "reading the data",
    data, fs = dat.get_electrdata(all_data, no_segments, ex_electr)  
    
    np.savez(save_folder + data_file, data = data, fs = fs)
    return data_all, fs
    del data_all


def update_upsample(data, fs, save_folder, uspl = 10, data_file = 'data_uspl_intra'):
    """ upsample the data (mostly used for the spike analysis)"""
    num = np.size(data,1) * uspl
    data_all = []
    for electr in range(np.size(data,0)):
        data_uspl = signal.resample(data[electr], num)
        data_all.append(data_uspl)
    fs2 = fs * uspl
    
    
    np.savez(save_folder + data_file, data_all, fs2)           
    return data_all, fs2    

def update_downsample(data, fs, save_folder, dspl = 2, data_file = 'data_dspl'):
    """downsamples given data
    dspl = downsampling"""
    data_dspl = []
    for electr in range(len(data)):
        data_tmp = [dat[::dspl] for dat in data[electr][:]]
        data_dspl.append(data_tmp)
    fs2 = fs / dspl
    
    np.savez(save_folder + data_file, data_dspl, fs2)           
    return data_dspl, fs2

def find_nearest(array,value):
    idx=(np.abs(array-value)).argmin()
    return idx, array[idx]

def update_highWaves(load_datafile, save_folder, data_file, atten_len = 25):
    """ performs moving average and based on this calculates probable beginning and end of the wave"""
    npzfile = np.load(save_folder + load_datafile)
    data = npzfile['data']
    fs = npzfile['fs']
    npzfile.close()

    
    thres_level = 10
    print
    print "removing averaged baseline and finding possible SPWs:",
    spws_starts = []
    spws_ends = []
    for electr in range(np.size(data,0)): 
        print electr,
        
        for trace in range(np.size(data,1)):

            window = ms2pts(atten_len, fs) # define the size of the window for removing the baseline
            data_used = data[electr, trace, :]
            
            new_dat, moved_avg = filt.remove_baseloc(data_used, window)     
            
            # find beginning and the end of the wave
            di = np.diff(moved_avg)
            di = np.hstack((di, 0))
            di[di>0] = 1
            di[di < 0] = 0

            # find moved_avg above thres_level
            possible_spws = dw.find_above(moved_avg, thres_level)   
            starts, ends = dw.find_startend(possible_spws)
            starts_di, temp = dw.find_startend(di)
            di = np.abs(di-1)
            temp, ends_di = dw.find_startend(di)
            
            
            spw_starts = np.zeros(len(starts))
            spw_ends = np.zeros(len(ends))
            for st_idx, st in enumerate(starts):
                st_used = starts_di[np.where(starts_di <= st)]
                idx, st_value = find_nearest(st_used,st)
                spw_starts[st_idx] = st_value
                
                en = ends[st_idx]
                en_used = ends_di[np.where(ends_di >= en)]
                if len(en_used) == 0:
                    en_used = ends_di
                idx, en_value = find_nearest(en_used,en)
                spw_ends[st_idx] = en_value
                
            # give indexes in ms (not data points!)
            spw_starts = pts2ms(spw_starts.astype(int), fs)
            spw_ends = pts2ms(spw_ends.astype(int), fs)
            
            # trasform everything to be rec arrays
            electrodes = np.ones(len(spw_starts), dtype=np.int32)*electr
            traces = np.ones(len(spw_starts), dtype=np.int32)*trace 
            
            spws_starts.append(np.rec.fromarrays([electrodes, traces, spw_starts], names='electrode,trace,time'))
            #import pdb; pdb.set_trace()
            spws_ends.append(np.rec.fromarrays([electrodes, traces, spw_ends], names='electrode,trace,time'))
                     
    spws_starts = np.concatenate(spws_starts)    
    spws_ends = np.concatenate(spws_ends)    
    
    np.savez(save_folder + data_file, starts = spws_starts, ends = spws_ends, fs = fs)   
    del data

def get_var_used(var, electr, trace):
    var_used = var[(var['electrode'] == electr) & (var['trace'] == trace)]['time']
    return var_used

def update_spikes_ampls(save_folder, save_file, load_spike_file):
    """ find in which electrode spikes have the highest amplitude"""
    
    allow_shift = 0.15
    
    #import pdb; pdb.set_trace() 
    npzfile         = np.load(save_folder + load_spike_file)
    spike_idxs      = npzfile['spike_idx']
    fs              = npzfile['fs']
    npzfile.close()  
    print 'Looking for origin of each spike'
    

    el,tr,sp,am = [],[],[],[]
    # make sure that you are not comparing spikes which appear on different traces
    for trace in np.unique(spike_idxs['trace']):
        # loop trough all the detected spikes sorted timewise
        
        spike_idxs_used = spike_idxs[spike_idxs['trace'] == trace]
        
        sort_idx        = np.argsort(spike_idxs_used['time'])
        posortowane     = spike_idxs_used[sort_idx]
        
        spike_old = -1
        considered_spikes, considered_ampls = [], []
        
        
        # check which spikes are enough close to each other to be considered one
        for (idx, spike) in enumerate(posortowane):
            # find amplitude of this particular spike
            sp_time = spike['time']
            ams = spike['amplitude']
    
            if (np.abs(sp_time - spike_old) < allow_shift) or (len(considered_spikes) == 0):
                # remember this one! It's close enought to the previous one
                considered_spikes.append(spike)
                considered_ampls.append(ams)
            else:
                # it is necessary to sum up previous considered_spikes and start new ones
                highest_ampl = np.argmax(considered_ampls)
                if np.size(highest_ampl) > 1:
                    highest_ampl = highest_ampl[0]
                wining_spike = considered_spikes[highest_ampl] 
                
                el.append(wining_spike['electrode'])
                tr.append(wining_spike['trace'])
                sp.append(wining_spike['time'])
                am.append(considered_ampls[highest_ampl])
                
                considered_ampls = []
                considered_spikes = []
                considered_spikes.append(spike)
                considered_ampls.append(ams)
            
            if idx == len(posortowane):
                # analyse the last spike(s) - repeat of above (might be improved in the future)
                highest_ampl = np.argmax(considered_ampls)
                if np.size(highest_ampl) > 1:
                    highest_ampl = highest_ampl[0]
                wining_spike = considered_spikes[highest_ampl] 

                el.append(wining_spike['electrode'])
                tr.append(wining_spike['trace'])
                sp.append(wining_spike['time'])
                am.append(considered_ampls[highest_ampl])
                  
            spike_old = spike['time']
        
    electrode = np.array(el, dtype='i4')
    trace = np.array(tr, dtype='i4')
    spike = np.array(sp, dtype='f8')    
    amplitude = np.array(am, dtype='f8')   
    large_spikes = np.rec.fromarrays([electrode, trace, spike, amplitude], 
                                               names='electrode, trace, time, amplitude')     
    
    np.savez(save_folder + save_file, fs = fs, spike_idx = large_spikes)
    del spike_idxs, large_spikes    
    
    
    
  
def update_spikes_in_spws(save_folder, save_file, load_spike_file, load_spw_file, win, spw_length = 80):
    """ finds the spikes for each spw"""
    
    npzfile         = np.load(save_folder + load_spike_file)
    #import pdb; pdb.set_trace()  
    spike_idxs      = npzfile['spike_idx']
    npzfile.close()  
    
    npzfile         = np.load(save_folder + load_spw_file)
    ipsps           = npzfile['spw_ipsps']
    npzfile.close()   
    #el, tr, sp, am, spw_no = [], [], [], [], []
    
    print 'Updating spikes in each SPW'
    spw_spike = []
    for spw in np.unique(ipsps['spw_no']):
        ipsp_used = ipsps[ipsps['spw_no'] == spw][0]
        #import pdb; pdb.set_trace()  
        spw_start = ipsp_used['spw_start']
        spw_end = spw_start + spw_length
        trace = ipsp_used['trace']
        #import pdb; pdb.set_trace()  
        spikes = spike_idxs[(spike_idxs['trace'] == trace) & (spike_idxs['time'] >= spw_start + win[0]) & (spike_idxs['time'] <= spw_end + win[1])]
        
        el = spikes['electrode'].astype('i4')
        tr = spikes['trace'].astype('i4')
        sp = spikes['time'].astype('f8')
        am = spikes['amplitude'].astype('f8')
        spw_no = (np.ones(len(spikes)) * spw).astype('i4')
        spw_st = (np.ones(len(spikes)) * spw_start).astype('f8')
        spw_en = (np.ones(len(spikes)) * spw_end).astype('f8')
        
        spw_spike.append(np.rec.fromarrays([el, tr, sp, am, spw_no, spw_st, spw_en],
                                                       names='electrode,trace, time, spike_ampl, spw_no, spw_start, spw_end'))    
        
     
    spw_spike = np.concatenate(spw_spike) 
    #import pdb; pdb.set_trace()  
    np.savez(save_folder + save_file, chosen_spikes = spw_spike)        
    del spike_idxs, spw_spike

     

def update_SPW_spikes_ampl(load_spikefile, load_spwsspike, save_folder, save_name):
    """ finds in which of the electrodes spike has the largest amplitude and possibly plot the traces"""
    
    npzfile         = np.load(save_folder + load_spikefile)
    spike_idxs      = npzfile['spike_idxs']
    ampls           = npzfile['ampls']
    fs              = npzfile['fs']
    npzfile.close()

    # load starts of spws
    npzfile         = np.load(save_folder + load_spwsspike)
    spw_spikes      = npzfile['spw_details']  
    npzfile.close()
    
    # load the parameters of extracellular spikes
    allow_shift = 0.15 # ms
    win         = (-5, 5)
    
    # initiate variables    
    chosen_spikes = []
    spw_len = len(np.unique(spw_spikes['spw_no']))
    el, tr, sp, st, en, sp_no, am = [],[],[],[],[],[],[]
    
    
    for spw in np.unique(spw_spikes['spw_no']):
        print str(spw) + '/' + str(spw_len)
        
        spike_old = -1
        spw_spikes_used = spw_spikes[spw_spikes['spw_no'] == spw]
        
        # sort the data for this trace
        #only_this_trace = spike_idxs[(spike_idxs['trace'] == trace)]
        sort_idx        = np.argsort(spw_spikes_used['spikes'])
        posortowane     = spw_spikes_used[sort_idx]
        
        # use old to not use the same spike twice
        old     = -1 # initiate old to be less then idx
        considered_spikes = []
        considered_ampls = []
        
        # loop trough all the detected spikes sorted timewise
        for (idx, spike) in enumerate(posortowane):
            #import pdb; pdb.set_trace() 
            # find amplitude of this particular spike
            electr = spike['electrode']
            trace = spike['trace']
            spikes = spike_idxs[(spike_idxs['electrode'] == electr) & (spike_idxs['trace'] == trace)]['time']
            id = np.where(spikes == spike['spikes'])[0][0]
            ams = ampls[(ampls['electrode'] == electr) & (ampls['trace'] == trace)]['time']
            
            
            # check if distance to preceeding spike is less than allow_shift
            if (np.abs(spike['spikes'] - spike_old) < allow_shift) or (len(considered_spikes) == 0):
                # remember this one!
                considered_spikes.append(spike)
                considered_ampls.append(ams[id])
            else:
                # it is necessary to sum up previous considered_spikes and start new ones
                highest_ampl = np.argmax(considered_ampls)
                wining_spike = considered_spikes[highest_ampl]
                
                el.append(wining_spike['electrode'])
                tr.append(wining_spike['trace'])
                sp.append(wining_spike['spikes'])
                st.append(wining_spike['spw_start'])
                en.append(wining_spike['spw_end'])
                sp_no.append(wining_spike['spw_no'])
                am.append(considered_ampls[highest_ampl])

                #import pdb; pdb.set_trace() 
                considered_ampls = []
                considered_spikes = []
                considered_spikes.append(spike)
                considered_ampls.append(ams[id])
            spike_old = spike['spikes']
                
    chosen_spikes.append(np.rec.fromarrays([el, tr, sp, am, st, en, sp_no],
                                                       names='electrode,trace, spikes, spike_ampl, spw_start, spw_end, spw_no'))    
    chosen_spikes = np.concatenate(chosen_spikes)
    np.savez(save_folder + save_name, chosen_spikes = chosen_spikes)        
    del spike_idxs, ampls, spw_spikes

def update_induc_spont_spw(save_folder, save_file, load_distances, load_spwfile, max_dist, ext):
    """ checks which spws are initiated and which are sponteneaus"""
    
    npzfile    = np.load(save_folder + load_distances)
    distances      = npzfile['dist_spwspike'] 
    npzfile.close()    
    
    npzfile    = np.load(save_folder + load_spwfile)
    ipsps      = npzfile['ipsps']  
    npzfile.close()   
    #before_pts = ispw.ms2pts(win[0], fs)
    #after_pts = ispw.ms2pts(win[1], fs)
    d = distances['distance']
    initiated_no = distances[(d <= max_dist) & (d>=0)]['spw_no']
    spont_no = distances[(d > max_dist) | (d<0)]['spw_no']
    
    #import pdb; pdb.set_trace() 
    init_set = np.in1d(ipsps['spw_no'], initiated_no, assume_unique= False)
    initiated = ipsps[init_set]
    spont_set = np.in1d(ipsps['spw_no'], spont_no, assume_unique=False)
    spontaneous = ipsps[spont_set]
    
    np.savez(save_folder + save_file, initiated = initiated, spontaneous = spontaneous)
       

def update_dist_SPWfromSpike(save_folder, save_file, load_intrafile, load_spwfile,  max_dist = 3, allowms = 0, spikes = 'first'):
    """ updates all the distances from extracellular spikes to the preceeding it intracellular spike,
    it also returns the indexes of those spikes which are closer than min_dist (in ms) to intracellular spike"""
    
    npzfile         = np.load(save_folder + load_intrafile)
    if spikes == 'first':
        spike_idxs      = npzfile['spikes_first'] # spikes_all
    elif spikes == 'all':
        spike_idxs      = npzfile['spikes_all'] # spikes_all
    fs              = npzfile['fs']
    npzfile.close()
    

    #allow_before = ms2pts(allow_before,fs)
    # load starts of spws
    npzfile         = np.load(save_folder + load_spwfile)
    ipsps      = npzfile['ipsps']  
    npzfile.close()    
    
    spw_no = []
    dist_all = []    
    electrodes = []
    traces = []
    print 'checking distance in trace: '
    #for electr in np.unique(ipsps['electrode']):
    #    print electr,
    #    #min_dist_all = []
    #    ipsps_electr = ipsps[ipsps['electrode'] == electr]
        
    for trace in np.unique(ipsps['trace']):
            #print trace,
            spw_electr_trace = np.unique(ipsps[ipsps['trace'] == trace])
            spw_beginnings = np.unique(spw_electr_trace['spw_start'])
            spike_electr_trace = spike_idxs[spike_idxs['trace'] == trace]['time']

            if len(spike_electr_trace) > 0:
                for spw in spw_beginnings:
                    #make sure that all spw_no are the same
                    spw_numb = spw_electr_trace[spw_electr_trace['spw_start'] == spw]['spw_no']
                    assert len(np.unique(spw_numb))==1
                    spw_numb = spw_numb[0]
                    
                    #import pdb; pdb.set_trace()
                    largest_spikes = np.sort(spike_electr_trace[(spike_electr_trace <= (spw))])
                    if len(largest_spikes) == 0:
                        spike_used = spike_electr_trace[0]
                    else:
                        spike_used = largest_spikes[-1]
                    dist = spw -spike_used
                    
                    dist_all.append(dist)
                    spw_no.append(spw_numb)
                    traces.append(trace)
                    
    
    type = 'f8'
    dists = np.array(dist_all, dtype=type)
    numbers = np.array(spw_no, dtype='i4')
    electrode = np.array(electrodes, dtype='i4')
    traces = np.array(traces, dtype='i4')    
    temp_spw = np.rec.fromarrays([traces, numbers, dists], 
                                               names='trace, spw_no, distance')                
    
    np.savez(save_folder + save_file, dist_spwspike = temp_spw) 
    

def count_coincident_ipsps(spw_ipsps_trace, shift_ipsp):
    """counts number of IPSPs in different electrodes seperated by not more
    than shift_ipsp"""
    i = np.argsort(spw_ipsps_trace['ipsp_start'])
    spw_ipsps_sorted = spw_ipsps_trace[i]
    spw_ipsps_intervals = np.diff(spw_ipsps_sorted['ipsp_start'])
    group_sep = spw_ipsps_intervals>shift_ipsp
    group_sep = np.concatenate([[0], group_sep])
    
    #calculate group IDs by counting seperators
    group_idx = np.cumsum(group_sep)
    electrode = spw_ipsps_sorted['electrode']

    #count number of disitinct electrodes fin each IPSP group
    n_uniq_electrodes_per_group = [len(np.unique(electrode[group_idx==group])) 
                                   for group in range(group_idx.max()+1)]
    
    n_uniq_electrodes_per_group = np.array(n_uniq_electrodes_per_group)
    
    #assign to each IPSP number of coincident IPSPs in other electrodes
    n_electrodes_per_ipsp = n_uniq_electrodes_per_group[group_idx]
    
    #inverse sorting
    inverse_i = np.argsort(i)
    return n_electrodes_per_ipsp[inverse_i]

def calculate_ipsp_rise(spw_ipsps_trace, data_trace, fs):  
    """calculate amplitude rise between IPSPs. for last IPSP in SPW/electrode
    return 0 """
    
    i = np.argsort(spw_ipsps_trace, order=['spw_no', 'electrode', 'ipsp_start'])
    spw_ipsps_sorted = spw_ipsps_trace[i]
    
    time_pts = ms2pts(spw_ipsps_sorted['ipsp_start'], fs).astype(np.uint32)
    ipsp_ampl = data_trace[spw_ipsps_sorted['electrode'], time_pts]
    ampl_rise = np.concatenate([np.diff(ipsp_ampl), [0]])
    
    #set amplitude to zero for last IPSP in each electrode
    last_ipsp = np.append(np.diff(spw_ipsps_sorted['electrode'])!=0, [True])
    ampl_rise[last_ipsp] = 0
    return  ampl_rise
    
    
    
def update_spws_beg(load_datafile, load_spwsipsp, load_spwsspike, save_folder, save_fig, save_file,ext):
    """ checks all the ipsps and corrects them for each spw"""


    
    npzfile         = np.load(save_folder + load_datafile)
    data            = npzfile['data']
    fs              = npzfile['fs']
    npzfile.close()   
    
    npzfile         = np.load(save_folder + load_spwsipsp)
    spw_ipsps       = npzfile['spw_ipsps']
    npzfile.close()   
    
    npzfile         = np.load(save_folder + load_spwsspike)
    spw_spike      = npzfile['chosen_spikes']  
    npzfile.close()     
    import pdb; pdb.set_trace()
    plot_it = True
    
    #import pdb; pdb.set_trace()
    shift_ipsp = 1 # ms
    min_electr = 3 # on how many electrodes IPSP should be detected for the first ipsp (beginning of SPW)
    expected_min_ipsp_ampl = 20 # microV
     
    spw_ipsps_list = []
    all_traces= np.unique(spw_ipsps['trace'])
    
    for trace in all_traces:
        spw_ipsps_trace = spw_ipsps[spw_ipsps['trace']==trace]
        
        n_electrodes_per_ipsp = count_coincident_ipsps(spw_ipsps_trace, shift_ipsp)
        
        ipsp_rise = calculate_ipsp_rise(spw_ipsps_trace, data[:, trace, :], fs)
        
        spw_ipsps_trace = spw_ipsps_trace[(ipsp_rise > expected_min_ipsp_ampl) &
                                          (n_electrodes_per_ipsp >= min_electr)]
        
        spw_ipsps_list.append(spw_ipsps_trace)

        
        
            
        
        
    
    # define variables
    add_it = 100
    
    shift_spike= 0.5 #ms
    proper_ipsps = []

    
    # go through all the spws
    for spw_no in np.unique(spw_ipsps['spw_no']):
        
        # save detected ipsps and spikes for each spw
        print 'updating spw:' + str(spw_no)
        sp_ip_used = spw_ipsps[spw_ipsps['spw_no'] == spw_no] # ipsps
        sp_sp_used = spw_spike[spw_spike['spw_no'] == spw_no] # spikes
        trace = sp_ip_used['trace'][0]                            
        data_temp = data[:,trace,:]           # traces
        sort_ip_idx     = np.argsort(sp_ip_used['ipsp_start'])    
        ip_sorted     = sp_ip_used[sort_ip_idx]
        
        ipsps_temp, spikes_temp, ipsp_ampl_temp, electr_temp, increase_temp, ipsp_ends = [], [], [], [], [], []
        el_ipsp_spw, tr_ipsp_spw, spwno_ipsp_spw, spwstart_isps_spw, spwend_ipsp_spw, ipspno_ipsp_spw = [],[], [], [], [], []
        ipspstrt_ipsp_spw, ampl_ipsp_spw = [], [], []
        
        ipsp_count = 0
        ipsp_old = -1

        #import pdb; pdb.set_trace() 
        spw_start = 0
        for ipsp in ip_sorted:
            # check every detected IPSP
            ipsp_time = ipsp['ipsp_start']
            ipsp_end = ipsp['ipsp_end']
            electr = ipsp['electrode']
        
                            
            if ipsp_time > ipsp_old + shift_ipsp and len(ipsps_temp) != 0:
                # analysis which IPSPs are the same (just detected in different electrodes)
                # check if no of found IPSPS isn't less than min_electr
                
                if len(ipsps_temp) >= min_electr:
                    # next IPSP is different    
                    # find if at least one ipsps causes increase more than expected_min_ipsp_ampl microV
                    if spw_start == 0:
                        # increase must be large only for beginning of SPW - not necessary for beginning of each IPSP
                        idx_temp = [i for i in range(len(increase_temp)) if increase_temp[i] > expected_min_ipsp_ampl]

                    if len(idx_temp) > 0 or spw_start != 0:   # this is possible IPSP  
                        used, spike = 0, np.inf
                        # find which of the electrodes has spike (going from maximum)
                        spike= float('Inf')
                        if sum(spikes_temp) != 0:

                            while spike == float('Inf') and used < len(electr_temp):
                                
                                spike_idx = np.argmin(spikes_temp)
                                if ipsp_ampl_temp[spike_idx] > expected_min_ipsp_ampl:
                                    spike = spikes_temp[spike_idx]
                                used = used+1

                        max_idx = np.argmax(ipsp_ampl_temp)    
                        if spike == float('Inf'): # otherwise use the beginning where the electrode has the highest IPSP
                            spike =  ipsps_temp[max_idx]
                            
                        
                        ipstarts_temp = []
                            
                        # add all this information to temporary lists
                        for el in range(int(min(electr_temp)), int(max(electr_temp))+ 1):
                            # if this electrode had IPSP detected:
                            sp_temp = np.array(spikes_temp)

                            notInf = sp_temp[sp_temp < np.inf]
                            if el in (electr_temp):
                                el_idx = electr_temp.index(el)
                                #if spikes_temp[el_idx] != np.inf:
                                sp = spikes_temp[el_idx]
                                # if amplitude was already found for this electrode
                                ipsp_idx = electr_temp.index(el)
                                ipsp_ampl = ipsp_ampl_temp[ipsp_idx]
                            # if this electrode has no IPSP detected or no spike around
                            if el not in electr_temp or sp == np.inf:
                                # if any spike exists
                                if len(notInf) > 0:
                                    #import pdb; pdb.set_trace() 
                                    data_idxs_temp = ms2pts(notInf, fs).astype(int)
                                    better_spike = np.argmin(data_temp[el,data_idxs_temp])
                                    sp = notInf[better_spike]
                                else:
                                    sp = spike
                                start_pts = ms2pts(spike, fs).astype(int)
                                
                                end_pts =ipsp_ends[max_idx]
                                
                                # for now if IPSPs size is 0, put amplitude to 0, but CHANGE IT!!!!
                                if start_pts >= end_pts:
                                    ipsp_ampl = 0
                                else:
                                    
                                    data_temp_ipsp = data_temp[el, start_pts: end_pts]
                                    #print start_pts
                                    #print end_pts
                                    maxpt = np.argmax(data_temp_ipsp)
                                    ipsp_ampl = data_temp_ipsp[maxpt] - data_temp_ipsp[0]
                                    
                                    #import pdb; pdb.set_trace() 
                                    

                            el_ipsp_spw.append(el)
                            tr_ipsp_spw.append(trace)
                            spwno_ipsp_spw.append(spw_no)
                            
                            ipspno_ipsp_spw.append(ipsp_count)
                            ipspstrt_ipsp_spw.append(sp)
                            ipstarts_temp.append(sp)
                            
                            sp_end = sp_sp_used[sp_sp_used['electrode'] == el]['spw_end']
                            #if len(sp_end) == 0:
                                #if end was not calculated, the longest SPW end will be reused
                            #else:
                            #    spwend_ipsp_spw.append(sp_end[0])
                            
                            ampl_ipsp_spw.append(ipsp_ampl)
                        
                        if spw_start == 0: # if the it's beginning of the SPW
                            #import pdb; pdb.set_trace()  
                            min_idx = np.argmin(ipstarts_temp)
                            spw_start = ipstarts_temp[min_idx]
                            #min_electr = min_electr_other # change to different constraint (for IPSP and not beginning of SPW)   
                            spw_electr_st = el_ipsp_spw[min_idx] 
      
                        spwstart_isps_spw.append(np.ones(len(ipstarts_temp)) * spw_start)
                        spw_duration = 80
                        spwend_ipsp_spw.append(np.ones(len(ipstarts_temp))*(spw_start+spw_duration))
                        ipsp_count = ipsp_count + 1
                    
                ipsps_temp, spikes_temp, ipsp_ampl_temp, electr_temp, increase_temp, ipsp_ends = [], [], [], [], [], []

            # collect all IPSPs
            ipsps_temp.append(ipsp_time)
            # find closest spike which is not further that shift_spike from the detected beginning
            electr_memb = [i in [electr] for i in sp_sp_used['electrode']]
            sp_temp = sp_sp_used[np.where(electr_memb)]['spikes']
            sp_temp = sp_temp[np.where(sp_temp >= ipsp_time - shift_spike)]
            if len(sp_temp) > 0 and sp_temp[0] - shift_spike <= ipsp_time:
                spikes_temp.append(sp_temp[0])
            else:
                spikes_temp.append(float('Inf'))
                
            # calculate amplitude of the ipsp (from the beginning to the highest peak)
            start_pts = ms2pts(ipsp_time, fs).astype(int)
            end_pts = ms2pts(ipsp_end,fs).astype(int)
            if end_pts <= start_pts:
                # in the future don't add this IPSP!!!!
                ipsp_ampl = 0
            else:
                data_temp_ipsp = data_temp[electr,start_pts: end_pts]

                maxpt = np.argmax(data_temp_ipsp)
                ipsp_ampl = data_temp_ipsp[maxpt] - data_temp_ipsp[0]
            ipsp_ampl_temp.append(ipsp_ampl)
            ipsp_ends.append(end_pts)
            
            # check incerase only if still looking for the beginning of the SPW
            if spw_start == 0:
                increase = data_temp_ipsp[-1] - data_temp_ipsp[0]
                increase_temp.append(increase)
                
            electr_temp.append(electr) 
            
            ipsp_old = ipsp_time
        #import pdb; pdb.set_trace()    
        if spw_start != 0:
            # change type of all data, so that it can be added to rec array
            el_ipsp_spw = np.array(el_ipsp_spw, dtype='i4')
            tr_ipsp_spw = np.array(tr_ipsp_spw, dtype='i4')
            spwno_ipsp_spw = np.array(spwno_ipsp_spw, dtype='i4')
            ipspno_ipsp_spw = np.array(ipspno_ipsp_spw, dtype='i4') 
            ipspstrt_ipsp_spw = np.array(ipspstrt_ipsp_spw, dtype='f8')
            ampl_ipsp_spw = np.array(ampl_ipsp_spw, dtype='f8')
            spwend_ipsp_spw = np.concatenate(spwend_ipsp_spw)
            spw_electr_start = np.ones(len(el_ipsp_spw), dtype=np.int32)*spw_electr_st
            spwstart_isps_spw = np.concatenate(spwstart_isps_spw)
            spwstart_isps_spw = np.array(spwstart_isps_spw, dtype='f8')
            temp_spw = np.rec.fromarrays([el_ipsp_spw, tr_ipsp_spw, spwno_ipsp_spw, spw_electr_start, spwstart_isps_spw,
                                                   spwend_ipsp_spw,ipspno_ipsp_spw, ipspstrt_ipsp_spw, ampl_ipsp_spw], 
                                               names='electrode, trace, spw_no, spw_electr_start, spw_start, spw_end, ipsp_no, ipsp_start, ipsp_ampl')
            proper_ipsps.append(temp_spw)
         
        #import pdb; pdb.set_trace()    
        if spw_start != 0 and plot_it and len(spw_electr_start) > 0: 
            fig = plt.figure()   
            spw_min_start = 9000000000
            spw_max_end = -1
            
            for el in np.unique(el_ipsp_spw):
                if el == spw_electr_st:
                    colr = 'r'
                else:
                    colr = 'k'
                idxs, = np.where(el_ipsp_spw == el)
                
                spw_start = spwstart_isps_spw[el_ipsp_spw == el][0]
                spw_end = max(spwend_ipsp_spw[el_ipsp_spw == el])
                spw_st_pts, spw_en_pts = ms2pts(spw_start, fs).astype(int), ms2pts(spw_end, fs).astype(int)
                data_used = data_temp[el,spw_st_pts:spw_en_pts]
                
                spikes = ms2pts(sp_sp_used[sp_sp_used['electrode'] == el]['spikes'], fs).astype(int) - spw_st_pts
                ipsp_start = ms2pts(ipspstrt_ipsp_spw[idxs], fs).astype(int) - spw_st_pts
                
                spw_min_start = min(spw_st_pts, spw_min_start)
                spw_max_end = max(spw_en_pts, spw_max_end)
                
                t = dat.get_timeline(data_used, fs, 'ms') + spw_start
                plt.plot(t,data_used + add_it * el, colr)
                sp_to_plot = spikes[(spikes>0) & (spikes<len(t))]
                plt.plot(t[sp_to_plot],data_used[sp_to_plot] + add_it * el, 'r<')
                ipsp_to_plot = ipsp_start[ipsp_start<len(data_used)]
                plt.plot(t[ipsp_to_plot],data_used[ipsp_to_plot] + add_it * el, 
                         'o', mfc='none', mec='g')

            tit = 'spw: ' + str(spw_no)
            plt.title(tit)
            fig.savefig(save_folder + save_fig + str(spw_no) + ext,dpi=600)
            fig.savefig(save_folder + save_fig + str(spw_no) + '.eps',dpi=600)        
            #plt.show()
            plt.close()
    if len(proper_ipsps) > 0:
        proper_ipsps = np.concatenate(proper_ipsps)

    np.savez(save_folder + save_file, ipsps = proper_ipsps) 
    
    
    
    del spw_ipsps, data



def update_SPW_ipsp_correct(load_datafile, load_spwsipsp, load_spwsspike, save_folder, save_fig, save_file,ext):
    """ checks all the ipsps and corrects them for each spw"""

    npzfile         = np.load(save_folder + load_datafile)
    data            = npzfile['data']
    fs              = npzfile['fs']
    npzfile.close()   
    
    npzfile         = np.load(save_folder + load_spwsipsp)
    spw_ipsps       = npzfile['spw_ipsps']
    npzfile.close()   
    
    npzfile         = np.load(save_folder + load_spwsspike)
    spw_spike      = npzfile['chosen_spikes']  
    npzfile.close()     
    
    plot_it = True

    add_it = 100
    shift_ipsp = 1 # ms
    min_electr_first = 3 # ipsp has to be found in at least that many electrodes
    min_electr_other = 2 # for the first ipsp (also beginning of SPW, check should be stronger
    shift_spike= 0.5 #0.15 # ms
    proper_ipsps = []
    expected_min_ipsp_ampl = 20 # microV
    
    
    
    #import pdb; pdb.set_trace()
    # go through all the spws
    for spw in np.unique(spw_ipsps['spw_no']):
        # save detected ipsps and spikes for each spw
        print 'updating spw:' + str(spw)
        
        sp_ip_used = spw_ipsps[spw_ipsps['spw_no'] == spw] # ipsps
        sp_sp_used = spw_spike[spw_spike['spw_no'] == spw] # spikes
        trace = sp_ip_used['trace'][0]                            
        data_temp = data[:,trace,:]           # traces
        sort_ip_idx     = np.argsort(sp_ip_used['ipsp_start'])    
        ip_sorted     = sp_ip_used[sort_ip_idx]
        spw_no = sp_ip_used['spw_no'][0]
        
        
        #spike_times = sp_sp_used['spikes']
        ipsps_temp, spikes_temp, ipsp_ampl_temp, electr_temp, increase_temp, ipsp_ends = [], [], [], [], [], []
        el_ipsp_spw, tr_ipsp_spw, spwno_ipsp_spw, spwstart_isps_spw, spwend_ipsp_spw, ipspno_ipsp_spw = [],[], [], [], [], []
        ipspstrt_ipsp_spw, ipspend_ipsp_spw, ampl_ipsp_spw = [], [], []
        
        ipsp_count = 0
        ipsp_old = -1
        #beginnings = []
        min_electr = min_electr_first
        #import pdb; pdb.set_trace() 
        spw_start = 0
        for ipsp in ip_sorted:
            #print 'ipsp: ' + str(ipsp_count)
            # check every detected IPSP
            ipsp_time = ipsp['ipsp_start']
            ipsp_end = ipsp['ipsp_end']
            electr = ipsp['electrode']
        
                            
            if ipsp_time > ipsp_old + shift_ipsp and len(ipsps_temp) != 0:
                # analysis which IPSPs are the same (just detected in different electrodes)
                # check if no of found IPSPS isn't less than min_electr
                
                if len(ipsps_temp) >= min_electr:
                    # next IPSP is different    
                    # find if at least one ipsps causes increase more than expected_min_ipsp_ampl microV
                    if spw_start == 0:
                        # increase must be large only for beginning of SPW - not necessary for beginning of each IPSP
                        idx_temp = [i for i in range(len(increase_temp)) if increase_temp[i] > expected_min_ipsp_ampl]

                    if len(idx_temp) > 0 or spw_start != 0:   # this is possible IPSP  
                        used, spike = 0, np.inf
                        # find which of the electrodes has spike (going from maximum)
                        spike= float('Inf')
                        if sum(spikes_temp) != 0:

                            while spike == float('Inf') and used < len(electr_temp):
                                
                                spike_idx = np.argmin(spikes_temp)
                                if ipsp_ampl_temp[spike_idx] > expected_min_ipsp_ampl:
                                    spike = spikes_temp[spike_idx]
                                used = used+1

                        max_idx = np.argmax(ipsp_ampl_temp)    
                        if spike == float('Inf'): # otherwise use the beginning where the electrode has the highest IPSP
                            spike =  ipsps_temp[max_idx]
                            
                        
                        ipstarts_temp = []
                            
                        # add all this information to temporary lists
                        for el in range(int(min(electr_temp)), int(max(electr_temp))+ 1):
                            # if this electrode had IPSP detected:
                            sp_temp = np.array(spikes_temp)

                            notInf = sp_temp[sp_temp < np.inf]
                            if el in (electr_temp):
                                el_idx = electr_temp.index(el)
                                #if spikes_temp[el_idx] != np.inf:
                                sp = spikes_temp[el_idx]
                                # if amplitude was already found for this electrode
                                ipsp_idx = electr_temp.index(el)
                                ipsp_ampl = ipsp_ampl_temp[ipsp_idx]
                            # if this electrode has no IPSP detected or no spike around
                            if el not in electr_temp or sp == np.inf:
                                # if any spike exists
                                if len(notInf) > 0:
                                    #import pdb; pdb.set_trace() 
                                    data_idxs_temp = ms2pts(notInf, fs).astype(int)
                                    better_spike = np.argmin(data_temp[el,data_idxs_temp])
                                    sp = notInf[better_spike]
                                else:
                                    sp = spike
                                start_pts = ms2pts(spike, fs).astype(int)
                                
                                end_pts =ipsp_ends[max_idx]
                                
                                # for now if IPSPs size is 0, put amplitude to 0, but CHANGE IT!!!!
                                if start_pts >= end_pts:
                                    ipsp_ampl = 0
                                else:
                                    
                                    data_temp_ipsp = data_temp[el, start_pts: end_pts]
                                    #print start_pts
                                    #print end_pts
                                    maxpt = np.argmax(data_temp_ipsp)
                                    ipsp_ampl = data_temp_ipsp[maxpt] - data_temp_ipsp[0]
                                    
                                    #import pdb; pdb.set_trace() 
                                    

                            el_ipsp_spw.append(el)
                            tr_ipsp_spw.append(trace)
                            spwno_ipsp_spw.append(spw_no)
                            
                            ipspno_ipsp_spw.append(ipsp_count)
                            ipspstrt_ipsp_spw.append(sp)
                            ipstarts_temp.append(sp)
                            
                            sp_end = sp_sp_used[sp_sp_used['electrode'] == el]['spw_end']
                            #if len(sp_end) == 0:
                                #if end was not calculated, the longest SPW end will be reused
                            #else:
                            #    spwend_ipsp_spw.append(sp_end[0])
                            
                            ampl_ipsp_spw.append(ipsp_ampl)
                        
                        if spw_start == 0: # if the it's beginning of the SPW
                            #import pdb; pdb.set_trace()  
                            min_idx = np.argmin(ipstarts_temp)
                            spw_start = ipstarts_temp[min_idx]
                            #min_electr = min_electr_other # change to different constraint (for IPSP and not beginning of SPW)   
                            spw_electr_st = el_ipsp_spw[min_idx] 
      
                        spwstart_isps_spw.append(np.ones(len(ipstarts_temp)) * spw_start)
                        spw_duration = 80
                        spwend_ipsp_spw.append(np.ones(len(ipstarts_temp))*(spw_start+spw_duration))
                        ipsp_count = ipsp_count + 1
                    
                ipsps_temp, spikes_temp, ipsp_ampl_temp, electr_temp, increase_temp, ipsp_ends = [], [], [], [], [], []

            # collect all IPSPs
            ipsps_temp.append(ipsp_time)
            # find closest spike which is not further that shift_spike from the detected beginning
            electr_memb = [i in [electr] for i in sp_sp_used['electrode']]
            sp_temp = sp_sp_used[np.where(electr_memb)]['spikes']
            sp_temp = sp_temp[np.where(sp_temp >= ipsp_time - shift_spike)]
            if len(sp_temp) > 0 and sp_temp[0] - shift_spike <= ipsp_time:
                spikes_temp.append(sp_temp[0])
            else:
                spikes_temp.append(float('Inf'))
                
            # calculate amplitude of the ipsp (from the beginning to the highest peak)
            start_pts = ms2pts(ipsp_time, fs).astype(int)
            end_pts = ms2pts(ipsp_end,fs).astype(int)
            if end_pts <= start_pts:
                # in the future don't add this IPSP!!!!
                ipsp_ampl = 0
            else:
                data_temp_ipsp = data_temp[electr,start_pts: end_pts]

                maxpt = np.argmax(data_temp_ipsp)
                ipsp_ampl = data_temp_ipsp[maxpt] - data_temp_ipsp[0]
            ipsp_ampl_temp.append(ipsp_ampl)
            ipsp_ends.append(end_pts)
            
            # check incerase only if still looking for the beginning of the SPW
            if spw_start == 0:
                increase = data_temp_ipsp[-1] - data_temp_ipsp[0]
                increase_temp.append(increase)
                
            electr_temp.append(electr) 
            
            ipsp_old = ipsp_time
        #import pdb; pdb.set_trace()    
        if spw_start != 0:
            # change type of all data, so that it can be added to rec array
            el_ipsp_spw = np.array(el_ipsp_spw, dtype='i4')
            tr_ipsp_spw = np.array(tr_ipsp_spw, dtype='i4')
            spwno_ipsp_spw = np.array(spwno_ipsp_spw, dtype='i4')
            ipspno_ipsp_spw = np.array(ipspno_ipsp_spw, dtype='i4') 
            ipspstrt_ipsp_spw = np.array(ipspstrt_ipsp_spw, dtype='f8')
            ampl_ipsp_spw = np.array(ampl_ipsp_spw, dtype='f8')
            spwend_ipsp_spw = np.concatenate(spwend_ipsp_spw)
            spw_electr_start = np.ones(len(el_ipsp_spw), dtype=np.int32)*spw_electr_st
            spwstart_isps_spw = np.concatenate(spwstart_isps_spw)
            spwstart_isps_spw = np.array(spwstart_isps_spw, dtype='f8')
            temp_spw = np.rec.fromarrays([el_ipsp_spw, tr_ipsp_spw, spwno_ipsp_spw, spw_electr_start, spwstart_isps_spw,
                                                   spwend_ipsp_spw,ipspno_ipsp_spw, ipspstrt_ipsp_spw, ampl_ipsp_spw], 
                                               names='electrode, trace, spw_no, spw_electr_start, spw_start, spw_end, ipsp_no, ipsp_start, ipsp_ampl')
            proper_ipsps.append(temp_spw)
         
        #import pdb; pdb.set_trace()    
        if spw_start != 0 and plot_it and len(spw_electr_start) > 0: 
            fig = plt.figure()   
            spw_min_start = 9000000000
            spw_max_end = -1
            
            for el in np.unique(el_ipsp_spw):
                if el == spw_electr_st:
                    colr = 'r'
                else:
                    colr = 'k'
                idxs, = np.where(el_ipsp_spw == el)
                
                spw_start = spwstart_isps_spw[el_ipsp_spw == el][0]
                spw_end = max(spwend_ipsp_spw[el_ipsp_spw == el])
                spw_st_pts, spw_en_pts = ms2pts(spw_start, fs).astype(int), ms2pts(spw_end, fs).astype(int)
                data_used = data_temp[el,spw_st_pts:spw_en_pts]
                
                spikes = ms2pts(sp_sp_used[sp_sp_used['electrode'] == el]['spikes'], fs).astype(int) - spw_st_pts
                ipsp_start = ms2pts(ipspstrt_ipsp_spw[idxs], fs).astype(int) - spw_st_pts
                
                spw_min_start = min(spw_st_pts, spw_min_start)
                spw_max_end = max(spw_en_pts, spw_max_end)
                
                t = dat.get_timeline(data_used, fs, 'ms') + spw_start
                plt.plot(t,data_used + add_it * el, colr)
                sp_to_plot = spikes[(spikes>0) & (spikes<len(t))]
                plt.plot(t[sp_to_plot],data_used[sp_to_plot] + add_it * el, 'r<')
                ipsp_to_plot = ipsp_start[ipsp_start<len(data_used)]
                plt.plot(t[ipsp_to_plot],data_used[ipsp_to_plot] + add_it * el, 
                         'o', mfc='none', mec='g')

            tit = 'spw: ' + str(spw)
            plt.title(tit)
            fig.savefig(save_folder + save_fig + str(spw) + ext,dpi=600)
            fig.savefig(save_folder + save_fig + str(spw) + '.eps',dpi=600)        
            #plt.show()
            plt.close()
    if len(proper_ipsps) > 0:
        proper_ipsps = np.concatenate(proper_ipsps)

    np.savez(save_folder + save_file, ipsps = proper_ipsps) 
    
    
    
    del spw_ipsps, data

def update_SPW_ipsp_ampl(save_folder, save_file, data, fs):
    pass

def update_ipsp_exSpikes(save_folder, save_file):
    pass

def update_SPW_ipsp(load_datafile, load_waves, load_spikes, save_folder, save_file, spw_length = 80):
    # it looks for the ipsps within detected spws - separate for each electrode
    # it take very long to analyse so be patient!
    npzfile         = np.load(save_folder + load_datafile)
    data            = npzfile['data']
    fs              = npzfile['fs']
    npzfile.close()
    #import pdb; pdb.set_trace()
    
    # load starts of spws
    npzfile         = np.load(save_folder + load_spikes)
    spike_idxs      = npzfile['spike_idx']
    npzfile.close()

    npzfile         = np.load(save_folder + load_waves)
    spw_details      = npzfile['spw_details']
    npzfile.close()
    
    spw_len_pts = ms2pts(spw_length, fs)
    print
    print "analyzing SPWs in electrode:",
    plot_it = True
    add_it = 100
    window = 0.5 # ms for calculating moving average
    window = ms2pts(window, fs)
    spw_ipsps = []
    
    
    # take each SPW separately and find ipsps
    spw_len = len(np.unique(spw_details['spw_no']))
    for spw in np.unique(spw_details['spw_no']):
        print str(spw) + '/' + str(spw_len)
        
        spw_used = spw_details[spw_details['spw_no'] == spw]
        trace = int(spw_used['trace'][0])
        spw_no = int(spw_used['spw_no'][0])
        
        # always use minimum detected start for this SPW
        min_start = min(spw_used['spw_start'])
        min_start_pts = ms2pts(min_start, fs).astype(int)
        for electr in range(np.size(data,0)):
            # no matter if wave was detected in this electrode or not, still check for IPSPS
            # as start use the earliest start found in any electrdoe for this spw
            data_used = data[electr, trace, min_start_pts:min_start_pts + spw_len_pts]
                      
            # calculate moving average  
            temp, moved_avg = filt.remove_baseloc(data_used, window)  
            rest_win = len(data_used) - np.floor(len(data_used)%window)
            
            prepare4mean = np.reshape(moved_avg[0:rest_win], [rest_win/window, window])
            meanEach = np.mean(prepare4mean.T, 0)
             
            err_allowed = 1 #mV
            mean_temp   = meanEach[0]
            switch      = 1 # looking for start of rise
            maxs, mins  = [], []
            
            # find ipsps - where the bins start changing their amplitude
            for idx, m in enumerate(meanEach):
                if m * switch > (mean_temp + err_allowed)*switch:
                    # found it (rise or fall)
                    detected = np.argmin(data_used[(idx-1)*window : (idx+1)*window] * switch)
                    if switch >0:
                        maxs.append(detected + (idx-1)*window)
                    else:
                        mins.append(detected + (idx-1)*window)
                    switch = switch * -1 # looking for start of fall
                
                mean_temp = m
                
            if plot_it:
                t = dat.get_timeline(data_used, fs, 'ms') + min_start #spw_used['spw_start'][0]
                plt.plot(t, data_used + add_it*electr)   
                plt.plot(t[maxs], data_used[maxs] + add_it*electr, 'r<')
                plt.plot(t[mins], data_used[mins] + add_it*electr, 'b<')
                
            #typ = 'f8'
            mini = pts2ms(maxs,fs) + min_start #spw_electr_used['spw_start'][0]
            #mini = np.append(mini,spw_electr_used['spw_end'][0])

            ipsp_start = mini.astype('f8')
            #ipsp_end = mini[1:].astype('f8')
            electrodes = np.ones(len(ipsp_start), dtype='i4')*electr
            traces = np.ones(len(ipsp_start), dtype='i4')*trace
            spw_num = np.ones(len(ipsp_start), dtype='i4')*spw_no
            ipsp_no = np.arange(0,len(mini)).astype('i4')
            spw_start = np.ones(len(ipsp_start), dtype = 'f8') * min_start

            #import pdb; pdb.set_trace()  
            spw_ipsps.append(np.rec.fromarrays([electrodes, traces, spw_num, spw_start, 
                                                    ipsp_no, ipsp_start], 
                                                   names='electrode,trace, spw_no, spw_start,ipsp_no,  ipsp_start'))
            
            
            if plot_it: 
                plt.plot(t, moved_avg + add_it*electr)

        if plot_it:
            plt.show()
    spw_ipsps = np.concatenate(spw_ipsps)
    np.savez(save_folder + save_file, spw_ipsps = spw_ipsps) 
    del data, spw_details 
      


def update_highWaves_numb(load_spwsfile, save_folder, data_file):
    # reject all the SPWs where there is less than min_no_wave spws detected, 
    # check which spws belog together and save them
    # it also numbers SPWs
    
    
    # load starts of spws
    npzfile         = np.load(save_folder + load_spwsfile)
    spw_starts      = npzfile['starts']
    #spw_ends        = npzfile['ends']
    npzfile.close()
    
    spw_details = []
    min_no_electr = 0
    max_move = 15 # ms

    spw_details_temp = []
    spw_no = 0
    no_traces = len(np.unique(spw_starts['trace']))
    for trace in np.unique(spw_starts['trace']):
        print 'trace: ' + str(trace + 1) + '/' + str(no_traces)
        
        # get spws for each electrode and check if it's the same
        spw_st_trace = spw_starts[(spw_starts['trace'] == trace)]
        #spw_en_trace = spw_ends[(spw_ends['trace'] == trace)]
        
        sort_idx = np.argsort(spw_st_trace['time'])
        spw_st_sorted = spw_st_trace[:,:,sort_idx]
        #spw_en_sorted = spw_en_trace[:,:,sort_idx]
        
        same = 0
        start_init, end_init = 0, 0 
        
        # analize every single SPW
        for idx, next_spw_st in enumerate(spw_st_sorted):
            #next_spw_en = spw_en_sorted[idx]
            
            # check if start is enclosed between start and end of previous one and beginning is not 
            # further than max_move from new beginning
            st_same = (start_init < next_spw_st['time']) 
            #en_same = (end_init > next_spw_st['time'])
            dist_same = abs(next_spw_st['time'] - start_init) < max_move
            
            if st_same and dist_same:
                same = same + 1
                start_init = min(start_init, next_spw_st['time'])
                #end_init = max(end_init, next_spw_en['time'])              

            elif same >= min_no_electr:
                # checks if this wave was detected in enough electrodes
                same = 0
                # save previous spws   
                spw_details = spw_details + spw_details_temp
                start_init = next_spw_st['time']
                #end_init = next_spw_en['time']        
                spw_details_temp = []
                spw_no = spw_no + 1    
            else:
                # something wrong with previous SPW = remove it
                spw_details_temp = []
                start_init = next_spw_st['time']
                #end_init = next_spw_en['time']
                same = 0

                
            # check the beginning of this SPWs
            electr = next_spw_st['electrode']
            trace = next_spw_st['trace']

            # save temp details of this spw
            typ = 'f8'
            start = next_spw_st['time']
            #end = next_spw_en['time']

            spw_details_temp.append(np.rec.fromarrays([[electr], [trace] ,[start], [spw_no]], names='electrode,trace, spw_start, spw_no'))
    
    #import pdb; pdb.set_trace()          
    # do the same check as before but for the last SPW      
    if same >= min_no_electr - 1:
        # checks if this wave was detected in enough electrodes and save
        spw_details = spw_details + spw_details_temp
    spw_details = np.concatenate(spw_details)   
    np.savez(save_folder + data_file, spw_details = spw_details) 
    del spw_starts

def update_databas(data_load, save_folder, data_file = 'data_bas'):
    """ to stay constant in all the data, it will not only update the data file, but all the
    files which are defined by this data as well
    """
    print
    print "setting the data to baseline, working on electrode:",
    
    npzfile = np.load(save_folder + data_load)
    data = npzfile['data']
    fs = npzfile['fs']
    npzfile.close()

    for electr in range(np.size(data,0)): 
        # remove all mean from each electrode 
        print electr,
        
        electro_data = data[electr, :, :] #data[data['electrode'] == electr]
        mean_datatime = np.mean(electro_data)
        data = data - mean_datatime

    np.savez(save_folder + data_file, data = data, fs = fs)   
    del data
    
def update_spws(data, fast_data = [], spw_data = [], fs = 0, save_folder = '', save_file = 'SPWs', thresh_mult= 2):
    """finds all the spws in the data """
    
    # SPW characteristics 
    min_length = 30# ms
    min_dist = 40
    before = 20 # ms
    after = 80 # ms
    before = ms2pts(before, fs)
    after = ms2pts(after, fs)
    min_len = ms2pts(min_length,fs)
    min_dist = ms2pts(min_dist,fs)  
    
    use_method = 2 # method 1 - hilbert; method 2 - threshold and checking width etc

    ends_spw = []
    
    t = dat.get_timeline(data[0][0], fs, 'ms')
    spw_idxs = []
    spw_maxs = [] 
    starts_spw = []
    lengths_spw= []
    print
    print 'looking for SPWs, working on electrode:',     
    for electr in np.unique(data['electrode']): 
        print electr,
        spw_idxs_trace = []
        ends_spw_trace = []
        spw_maxs_trace = []
        starts_spw_trace = []
        lengths_spw_trace = []
        
        for trace in np.unique(data[data['electrode'] == electr]['trace']):

        # find SPWs
            if use_method == 1:
                # use hilbert approach
                max_idxs, max_ampl = find_SPWs(data[electr][trace], fast_data[electr][trace], fs, min_length)
                # find start of the wave based on the SPW data
                temp, spw_start = find_locals(spw_data[electr][trace])
                starts, ends = assign_startend(spw_start, max_idxs)
                starts_spw.append(starts)
                ends_spw.append(ends)
            elif use_method == 2:
                # use threshold and checking width etc
                noiselevel = np.std(data[electr][trace])
                thres_level = noiselevel * thresh_mult
                mean_data = np.mean(data[electr][trace])
                possible_spws = dw.find_above(data[electr][trace], thres_level)          
                # filter data to remove spikes - below 500
                freq_dat = filt.lowPass(750, fs, data[electr][trace][:])
                # check if the width of the spw is large enough
                starts, ends = dw.find_startend(possible_spws)
                maxs, idxs = dw.max_waves(freq_dat, starts, ends)

                m_dat = mean_data*np.ones((1, len(maxs)))[0]
                lengths, starts, ends = dw.check_width(data[electr][trace], idxs, m_dat, fs)
                maxs, idxs = dw.max_waves(data[electr][trace], starts, ends)

                percent10 = (maxs - mean_data)*0.2
                temp, starts, ends = dw.check_width(data[electr][trace], idxs, percent10, fs)
                lengths, temp, temp = dw.check_width(data[electr][trace], idxs, m_dat, fs)
                len_rightidx = np.where(lengths >= min_len)[0]
                wrong = []

                for l in range(len(len_rightidx)-1):

                    if np.size(len_rightidx) > 1:
                        if ends[len_rightidx[l]] > starts[len_rightidx[l+1]]:
                            # if the wave includes another wave
                            wrong.append(len_rightidx[l+1])
                    else:
                        wrong.append(len_rightidx[l])
                        #else:
                        #    wrong.append(len_rightidx[l])

                
                len_rightidx = dw.remove_from(len_rightidx, wrong)

                # remove those waves which are wrong
                #for l in range(len(wrong)):
                #    len_rightidx = dw.remove_from(len_rightidx, l)
                #print 'still here4'

                idxs =  [idxs[i] for i in len_rightidx]
                lengths = [lengths[i] for i in len_rightidx]
                maxs = [maxs[i] for i in len_rightidx]
                starts = [starts[i] for i in len_rightidx]
                ends = [ends[i] for i in len_rightidx]            
                
                spw_idxs_trace.append(idxs)
                ends_spw_trace.append(ends)
                spw_maxs_trace.append(maxs)
                starts_spw_trace.append(starts)
                lengths_spw_trace.append(lengths)
                percent20 = [percent10[i] for i in len_rightidx]
                
                
    #            
    #            # remove waves which are the same
    #            data_waves = np.zeros((1, len(data[electr])))
    #            for o in range(len(starts)):
    #                data_waves[starts[o]:ends[o]] = data_waves[starts[o]:ends[o]] + 1
    #            print data_waves
    #            print np.unique(data_waves)
                    
                #print lengths[len_right]
                #t = dat.get_timeline(data[electr], fs, 'ms')
        spw_idxs.append(spw_idxs_trace)
        ends_spw.append(ends_spw_trace)
        spw_maxs.append(spw_maxs_trace)
        starts_spw.append(starts_spw_trace)
        lengths_spw.append(lengths_spw_trace)
        
#        plt.figure()
#        plt.plot(t, data[electr][trace])
#        plt.plot(t, freq_dat)
#        #plt.plot(t[])
#        plt.plot(t[idxs], data[electr][trace][idxs], 'ro')
#        plt.plot(t[starts], data[electr][trace][starts], 'go')
#        plt.plot(t[ends], data[electr][trace][ends], 'mo')
##        #plt.plot(t[idxs], percent20, 'rx')
#        plt.hlines(noiselevel, 0, t[-1])
#        plt.hlines(thres_level, 0, t[-1])
#        plt.hlines(mean_data, 0, t[-1])    
##    #                    
##    #            
#        plt.show()
#        

        #spw = dw.cut_waves(spw_data[electr], np.transpose(spw_start).tolist(), fs, before, after)
        #whole_spws.append(spw)
        #all_idxs.append(max_idxs)
        #events_ampl.append(data[electr][max_idxs])
        
    #np.savez(save_folder + save_file, all_idxs, events_ampl, whole_spws, starts, ends, fs) 
    np.savez(save_folder + save_file, spw_idxs, spw_maxs, starts_spw, ends_spw, lengths_spw, fs)
    return spw_idxs, spw_maxs, starts_spw, ends_spw, lengths_spw, fs
 
    

        




#def update_spike_events_dist():    
#    spike_event_dist = "distances"
#    
#    # find the distances from each event to previous intra-spike
#    distance = dist_fromSpike(sp_idx, max_idxs)
#    distance = pts2ms(distance, fs2)
#    distances.append(distance)
#    
#    np.savez(save_folder + spike_event_dist, distances)    
         

## -0000000000000000000000000000000000000000000000000000000 
#use_exectrodes =  [1, 2, 3, 4, 5, 6, 7, 8, 9] 
#save_folder = '/home/maja/PhDProject/SPWs/SPWs/saved_data/cell3/'  
#f_dir = '/home/maja/PhDProject/SPWs/data/induced/Cell 3/'
#f_name = '08112011_0000_gap free.abf'
#filename = f_dir + f_name
#update_datafile(filename, use_exectrodes, save_folder)
##npzfile = np.load(save_folder + "data.npz")
##data = npzfile['arr_0']
##fs = npzfile['arr_1'] 
##npzfile = np.load(save_folder + "data_bas.npz")
##data = npzfile['arr_0']
##fs = npzfile['arr_1'] 
##spw_data = np.load(save_folder + "spw_data.npz")['arr_0']
##fast_data = np.load(save_folder + "fast_data.npz")['arr_0']
#
##
##update_databas(data, fs, save_folder)
##update_filtered(data_bas, fs, save_folder, intra = 1)
##update_events(data, fast_data, spw_data, fs, save_folder, intra = 1, pulse_len = 500)
#npzfile = np.load(save_folder + "spw_data.npz")
#print npzfile.files
#files =  npzfile.files
#for i in range(len(files)):
#    print np.shape(npzfile[files[i]])
#    print npzfile[files[i]]
## -0000000000000000000000000000000000000000000000000000000









#    
#def prepare_data(filename, save_folder): 
#    """assiging the variables ------!!!!"""
#    # set variable
#    scale = 'ms'
#    intra = 1 # no of intracellular electrode, if intra = -1 -> there is no intra electrode
#    use_exectrodes =  [2, 3] # electrodes which are to be used (not including intra electrode)
#            
#    dspl = 1 # downsampling variable
#
#    # how much to cut the data  and after the specified event (SPW or spike)
#    before = 20 # in ms
#    after = 100
#    
#    

#    """ the end of assiging the variables ------!!!!"""
#
#    if intra != -1:
#        use_exectrodes.insert(0, intra)
#    use_electrodes = np.array(use_exectrodes)
#    
#    # read the data (again)
#    data = update_datafile(filename, use_electrodes) 
#               
#    # init variables
#
#    data_bas_all = []
#    SPW_data_all = []
#    ripple_data_all = []
#    fast_data_all = []
#    
#    all_maxs = []
#    max_idxs = []
#    distances = []
#    sp_idx = []
#    events_idx = []
#    events_ampl = []
#    events_idx2 = []
#    events_ampl2 = []
#    spikes = []
#    spikes_a = []    
#    
#
#    
##    for electr in use_electrodes:
##        print electr
##        # reading the data
##        data, fs = dat.get_electrdata(all_data, electr, data_part)
##        data_all.append(data)
#        # downsample the data and find timeline       
#
#        

##        else:
#            # only extracellular data
#            # remove baseline 
#            data_bas = filt.remove_base(data_bas)
#            window = fs2 * 0.8 # define the size of the window for removing the baseline
#            data_bas, temp = filt.remove_baseloc(data_bas, window)
#            

#            

#            

#            
#            # save the indexes and amplitudes; save the whole data
#            #events_idx.append(max_idxs)
#            #events_ampl.append(max_ampl)
#            

#            
#            #idx = range(0,len(data_bas))
#            #events_idx2.append(idx)
#            #events_ampl2.append(data_bas)
#            
#            #ir_more0 =  [o for o in range(len(iris)) if iris[o] > 0]
#
#        data_bas_all.append(data_bas)
