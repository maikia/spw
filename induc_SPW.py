import scipy.signal as signal
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
import os
import spike_detection as fes
from matplotlib import pyplot
import pylab as plt
import filtit as filt
#import plotting as plt
import scipy as sc
import detect_waves as dw
from scipy import signal
import data_mang as dat
import folder_manager as fold_mang
from random import randint

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

def define_first_spikes(spikes, current_pulse):
    """ finds which spikes are the closest to the beginnings of IPSP but no 
    further than shift_spike"""
    
    first_spikes = np.zeros(len(current_pulse))
    
    i_current = np.argsort(current_pulse)
    i_spikes = np.argsort(spikes)
    
    current_sorted = current_pulse[i_current]
    spikes_sorted = spikes[i_spikes]
    
    for idx, t in enumerate(current_sorted):
        #import pdb; pdb.set_trace()
        try:
            idx_spike, nearest_spike = find_nearest(spikes[spikes > t],t)
        except:
            nearest_spike = np.nan
        first_spikes[idx] = nearest_spike
    
    inverse_i = np.argsort(i_current)   
    first_spikes = first_spikes[inverse_i]    
    return first_spikes[~ np.isnan(first_spikes)]

def update_intraSpikes(save_folder, save_file, load_file, pulse_len = 300):
#(data, fs, save_folder, save_file = "intra_spikes", pulse_len = 500, ):
    """ pulse_len in ms - length of the stimulation pulse in intracellular electrode"""
    npzfile = np.load(save_folder + load_file)
    data = npzfile['data']
    fs = npzfile['fs'] 
    npzfile.close()
    way = 1
    
    sp_all = []
    sp_first = []
    thres_mult = 6 #mV
    print 'Detecting intracellular spikes'
    for trace in range(np.size(data,1)):
        data_used =data[0, trace, :]
        spiking_thres = np.std(data_used)* thres_mult

        # detect only the first spike in the row 
        #sp_idx_first = detect_1spike(data_used, spiking_thres, fs, pulse_len)
        sp_idx_a = detect_spikes(data_used, spiking_thres) # for detecting all the spikes 
        if way == 1:
            sp_idx_first = detect_1spike(data[0, trace, :], spiking_thres, fs, pulse_len)
            sp_idx_a = detect_spikes(data[0, trace, :], spiking_thres) # for detecting all the spikes 
            
        else: 
            # detect only the first spike in the row 
            #sp_idx_first = detect_1spike(data_used, spiking_thres, fs, pulse_len)
            sp_idx_a = detect_spikes(data_used, spiking_thres) # for detecting all the spikes 
            
            input_current = detect_spikes(data_used*(-1), 10)
            #import pdb; pdb.set_trace()
            sp_idx_first = define_first_spikes(sp_idx_a, input_current)
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
#        plt.figure()
#        part = [0, 200000]
#        data_temp = data[0,trace,:]#part[0]:part[1]]
#        
#        t = dat.get_timeline(data_temp, fs, 'ms')
#        plt.plot(t, data_temp, 'k')
#        
#        spik = spike_idxs_first
#        spiks = ms2pts(spik, fs).astype('i4')
#        #spiks_temp = spiks[(spiks < part[1]) & (spiks > part[0])]
#        #plt.plot(t[spiks_temp - part[0]], data_temp[spiks_temp - part[0]], 'go')
#        plt.plot(t[spiks], data_temp[spiks], 'ro')
#        #plt.show()
#        
#        plt.show()
#    import pdb; pdb.set_trace()  
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
        return data_filt, fs
        
    except IOError as e:
        # not? create it!
        pass

    if np.size(freq) == 1:
        data_filt = filt.highPass(freq, fs, data, N)
    elif freq[0] == -1:
        data_filt = filt.lowPass(freq[1], fs, data, N)
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
       
    freq_fast = 600.
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
            #plt.show()
    
         
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
    
    thres = .5
    min_length = 15
    min_length_pts = ms2pts(min_length, fs)
    print
    print "removing averaged baseline and finding possible SPWs, electrode:",
    spws_starts = []
    spws_ends = []
    plot_it = False
    
    for electr in range(np.size(data,0)): 
        print electr,
        
        
        thres_data = data[electr, :, :]
        thres_window = np.std(thres_data)
        #print "thres: ", thres_level
        for trace in range(np.size(data,1)):
            #import pdb; pdb.set_trace() 
            
            window = np.floor(ms2pts(thres_window, fs)) # define the size of the window for removing the baseline
            data_used = data[electr, trace, :]
            
            #thres_level = np.std(thres_data) * thres
            #import pdb; pdb.set_trace()
            new_dat, moved_avg = filt.remove_baseloc(data_used, window)     
            thres_level = np.std(moved_avg) * thres
            
            # find beginning and the end of the wave
            di = np.diff(moved_avg)
            di = np.hstack((di, 0))
            di[di>0] = 1
            di[di < 0] = 0

            # find moved_avg above thres_level
            possible_spws = dw.find_above(moved_avg, thres_level)   
            starts, ends = dw.find_startend(possible_spws)
            lengths = ends - starts
            leng_enough = lengths >= min_length_pts
            #import pdb; pdb.set_trace() 
            starts = starts[leng_enough]
            ends = ends[leng_enough]
            starts_di, temp = dw.find_startend(di)
            di = np.abs(di-1)
            temp, ends_di = dw.find_startend(di)
            
            #spw_starts = pts2ms(starts, fs)
            #spw_ends = pts2ms(ends, fs)
            
            #import pdb; pdb.set_trace() 
            spw_starts = np.zeros(len(starts))
            spw_ends = np.zeros(len(ends))
            for st_idx, st in enumerate(starts):
                #import pdb; pdb.set_trace() 
                
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
            spws_ends.append(np.rec.fromarrays([electrodes, traces, spw_ends], names='electrode,trace,time'))
            #import pdb; pdb.set_trace() 
            
            if plot_it:
                
                add_it = 100
                t = dat.get_timeline(data_used, fs, 'ms')
                plt.plot(t,possible_spws* 50 + add_it * electr)
                plt.plot(t, data_used + add_it * electr)
                plt.plot(t, moved_avg + add_it * electr)
                plt.axhline(y=thres_level + add_it * electr, xmin=t[0], xmax=t[-1])
                
                spw_pts = ms2pts(spw_starts, fs).astype('i4')
                #import pdb; pdb.set_trace() 
                plt.plot(t[spw_pts], data_used[spw_pts] + add_it * electr, '*r', ms = 10)

            #plt.show()
             
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
    
    #import pdb; pdb.set_trace() 
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
            

            spike_old = spike['time']
        
    electrode = np.array(el, dtype='i4')
    trace = np.array(tr, dtype='i4')
    spike = np.array(sp, dtype='f8')    
    amplitude = np.array(am, dtype='f8')   
    large_spikes = np.rec.fromarrays([electrode, trace, spike, amplitude], 
                                               names='electrode, trace, time, amplitude')     
    
    np.savez(save_folder + save_file, fs = fs, spike_idx = large_spikes)
    del spike_idxs, large_spikes    
    

     

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

def take_random_elements(arr, no_elements):
    siz = len(arr)
    elements = np.zeros(no_elements)
    #elements = []
    data = arr.copy()
    for numb in range(no_elements):
        siz = siz - 1
        index = randint(0, siz)
        elem = data[index]
        data[index] = data[siz]
        #import pdb; pdb.set_trace()
        elements[numb]= elem
        #elements.append(elem)
    #elements = np.concatenate(elements)
    return elements

def update_add_missing_ipsps(save_folder, save_file, spw_file, data_file):
    """ It goes through every group of IPSPs and checks in other electrodes if this IPSP
    is not missing; if it's missing - it ads it to this group. It also finds the maximum
    of each IPSP in each electrode"""
    npzfile        = np.load(save_folder + spw_file)
    spws = npzfile['spw_ipsps']
    npzfile.close()   
    
    npzfile        = np.load(save_folder + data_file)
    data = npzfile['data']
    fs = npzfile['fs']
    npzfile.close()   
           
#    import pdb; pdb.set_trace()
#    
#    for spw_no in np.unique(spws['spw_no']):
#        spw_used = spws[spws['spw_no'] == spw_no]
#        trace = spw_used['trace'][0]
#        data[:,trace,:]
#        
#        for group in np.unique(spw_used['group']):
#            spw_group = spw_used[spw_used['group'] == group]
#            min()
            
            #
            
    
    # add maximums of IPSPs
    
    
def check_if_ipsp(data_prev_ipsp, min_ampl):
    """ checks what is the maximum of this given trace and if it is larger than both 
    egdes of the trace and if it is higher than the smallest trace than the given min_ampl"""
    #import pdb; pdb.set_trace()
    max_trace = max(data_prev_ipsp)
    left_trace = data_prev_ipsp[0]
    right_trace = data_prev_ipsp[-1]
    if max_trace - max(left_trace, right_trace) > min_ampl:
        return 1
    else:
        return -1


def merge_ipsps(spws_used, group2, group1, data, fs):
    """ merges given ipsp groups into one"""
    spws1 = spws_used[spws_used['group'] == group1]
    spws2 = spws_used[spws_used['group'] == group2]
    electr_used = np.unique(np.concatenate([spws1['electrode'],spws2['electrode']]))
    all_spws = []
    try:
        trace = spws1['trace'][0]
    except:
        import pdb; pdb.set_trace()
    if len(spws1) + len(spws2) == len(electr_used):
        # if the groups are on different electrodes- just merge them to one
        spw_temp = spws_used[spws_used['group'] == group2]
        spw_temp['group'] = np.ones(len(spws_used[spws_used['group'] == group2])) * group1
        spws_used[spws_used['group'] == group2] = spw_temp
    else:    
        # if the groups are on the same electrode choose the IPSP of the smaller amplitude
        new_group = group1
        new_spws = []
        for electr in electr_used:
            spw1_temp = spws1[spws1['electrode'] == electr]
            spw2_temp = spws2[spws2['electrode'] == electr]
            assign = 0
            
            # check if this electrode is in both SPWs
            if len(spw1_temp) > 0 and len(spw2_temp)> 0:
                # choose which one to add
                
                ipsp1_pts = ms2pts(spw1_temp['ipsp_start'], fs).astype('i4')
                ipsp2_pts = ms2pts(spw2_temp['ipsp_start'], fs).astype('i4')
                if data[electr, trace, ipsp1_pts][0] == data[electr, trace, ipsp2_pts][0]:
                    assign = 1
                else:
                    #import pdb; pdb.set_trace()
                    smaller = np.argmin([data[electr, trace, ipsp1_pts], data[electr, trace, ipsp2_pts]])
                    assign = smaller + 1
                
            if len(spw1_temp) == 0 or assign == 2:
                # add only second spw (changing the group)
                spw2_temp['group'] = new_group
                new_spws.append(spw2_temp)
                
            elif len(spw2_temp) == 0 or assign == 1:
                # add only first spw
                new_spws.append(spw1_temp)
            
            #all_spws.append(new_spws) 
        #new_spws = np.concatenate(new_spws)   
        #import pdb; pdb.set_trace()
        spws_used = spws_used[spws_used['group'] != group1]
        spws_used = spws_used[spws_used['group'] != group2]
        new_spws.append(spws_used)
        spws_used = np.concatenate(new_spws)
        #import pdb; pdb.set_trace()
    return spws_used
        
def update_merge_close_groups(save_folder, save_file, spw_file, data_file):  
    """ it merges too close groups of IPSPs, it chooses the one which is lower if there are
    two on the same electrode"""
    # load necessary data
    npzfile        = np.load(save_folder + spw_file)
    spws = npzfile['spw_ipsps']
    npzfile.close()   
    
    npzfile        = np.load(save_folder + data_file)
    data = npzfile['data']
    fs = npzfile['fs']
    npzfile.close()  
    
    min_distance_allowed = 1 # ms
    
    all_new_spws = []
    for spw_no in np.unique(spws['spw_no']):
        spws_used = spws[spws['spw_no'] == spw_no]
        sorted_spws = np.argsort(spws_used[['ipsp_start']])
        spws_used = spws_used[sorted_spws]
        
        prev_group_time = [-5, -5]
        all_groups, gr_idx = np.unique(spws_used['group'], return_index = True)
        order_groups = np.argsort(gr_idx)
        all_groups = all_groups[order_groups]
        for idx, group in enumerate(all_groups):
            
            curr_spw = spws_used[spws_used['group'] == all_groups[idx]]
            curr_group_time = [min(curr_spw['ipsp_start']), max(curr_spw['ipsp_start'])]
            
            if curr_group_time[0] - prev_group_time[1] <= min_distance_allowed:
                #import pdb; pdb.set_trace()
                # the two groups are too close
                spws_used = merge_ipsps(spws_used, all_groups[idx-1], group, data, fs)
                 
            prev_group_time = curr_group_time
                
        all_new_spws.append(spws_used)
    all_new_spws = np.concatenate(all_new_spws)
    #import pdb; pdb.set_trace()        
    # save all_new_spws
    
    np.savez(save_folder + save_file, spw_ipsps = all_new_spws)   
               
#            import pdb; pdb.set_trace()
#        
#    ipsp_len = 10
#    ipsp_size_pts = ms2pts(ipsp_len, fs).astype('i4')
#    win = [-20, 100]
#    ipsp_allowed_error = 0.5
#    min_len_ipsp = 1
#    min_len_ipsp_pts = ms2pts(min_len_ipsp, fs).astype('i4')
#    ipsp_ae_pts = ms2pts(ipsp_allowed_error, fs).astype('i4') # allowed error in points
#    plot_it = False
#    all_new_ipsps = []
#    time_if_no_ipsp = 10 #ms
#    time_if_no_ipsp_pts = ms2pts(time_if_no_ipsp, fs).astype('i4')
#    min_ipsp_height = 15
    
            
    
    
def update_add_missing_electrodes_SPW(save_folder, save_file, spw_file, data_file):
    """ it checks if in neighbouring groups there are no ipsps 
    from the group already detected; it also adds location of maximum
    of each IPSP"""
    npzfile        = np.load(save_folder + spw_file)
    spws = npzfile['spw_ipsps']
    npzfile.close()   
    #min_ampl = 20
    npzfile        = np.load(save_folder + data_file)
    data = npzfile['data']
    fs = npzfile['fs']
    npzfile.close()  
    
    ipsp_len = 10
    ipsp_size_pts = ms2pts(ipsp_len, fs).astype('i4')
    win = [-20, 100]
    ipsp_allowed_error = 0.5
    min_len_ipsp = 1
    min_len_ipsp_pts = ms2pts(min_len_ipsp, fs).astype('i4')
    ipsp_ae_pts = ms2pts(ipsp_allowed_error, fs).astype('i4') # allowed error in points
    plot_it = False
    all_new_ipsps = []
    time_if_no_ipsp = 10 #ms
    time_if_no_ipsp_pts = ms2pts(time_if_no_ipsp, fs).astype('i4')
    min_ipsp_height = 15
    
    for spw_no in np.unique(spws['spw_no']):
        spw_used =  spws[spws['spw_no'] == spw_no]   
        groups = np.unique(spw_used['group'])
        
        # separate data for this SPW
        trace = spw_used['trace'][0]
        spw_start = spw_used['spw_start'][0]
        spw_start_pts = ms2pts(spw_start + win[0], fs).astype('i4')
        spw_end_pts = ms2pts(spw_start + win[1], fs).astype('i4')
        data_trace = data[:,trace,:]
        
        group_times = []
        
        group_times = np.zeros(len(groups))
        # for every group select the average time
        for idx, gr in enumerate(groups):
            all_times = spw_used[spw_used['group'] == gr]['ipsp_start']
            averag = np.mean(all_times)
            group_times[idx] = averag
        
        #import pdb; pdb.set_trace()
        plot_ipsp = []
        # order the IPSPs
        group_order = np.argsort(group_times)
        group_times = group_times[group_order]
        groups = groups[group_order]
        
        all_electrodes = range(len(data_trace))
        
        # analyze every group
        for idx, group in enumerate(groups):
            # check which electrode doesn't detect this group
            group_electr = spw_used[spw_used['group'] == group]['electrode']
            no_group_electr  = np.setdiff1d(all_electrodes, group_electr)
            ipsp_time = spw_used[spw_used['group'] == group]['ipsp_start']
            #group_times[idx]
            
            group_pts_all = ms2pts(ipsp_time, fs).astype('i4')
            
            if idx != 0: # check if not first ipsp, 
                # find the time of the previous
                ipsp_prev = group_times[idx - 1] 
                ipsp_start_previous = ms2pts(ipsp_prev, fs).astype('i4')
            else: 
                ipsp_prev = -1
                ipsp_start_previous = -1

            if idx != len(groups)-1: # check if not last ipsp
                # check the time of the next 
                ipsp_next = group_times[idx + 1]
                #    import pdb; pdb.set_trace() 
                ipsp_end_next = ms2pts(ipsp_next, fs).astype('i4')
            else: 
                ipsp_next = -1 
                ipsp_end_next = -1
                
                              
            for electr in no_group_electr:
                # check what happens between this and previous ipsps
                # uneless it's the first one
                
                # find index of smallest argument between the two IPSPs which also is changing it's sign on derivative
                # check where to look for the min
                search_start = max(ipsp_start_previous, min(group_pts_all)-ipsp_ae_pts)
                search_end = min(ipsp_end_next, max(group_pts_all)+ipsp_ae_pts)
                # check where first derivative is changing its sign
                deriv = np.diff(data_trace[electr, search_start:search_end])
                deriv[deriv < 0] = 0
                deriv[deriv > 0] = 1
                deriv = deriv.tolist()
                try:
                    deriv_1 = np.array(deriv + [deriv[-1]])
                except:
                    import pdb; pdb.set_trace() 
                deriv_2 = np.array([deriv[0]] + deriv)
                change = deriv_1 - deriv_2
                data_temp_1 = data_trace[electr, search_start: search_end].copy()
                data_temp_1[change < 1] = 1000000
                pot_ipsp = np.argmin(data_temp_1) + 1
                
                #if spw_no == 7 and electr == 5 and group == 8:
                #    import pdb; pdb.set_trace() 
                
                group_pts = search_start - ipsp_ae_pts + pot_ipsp
                
                if ipsp_prev == -1:
                    # assign different start - var_ms before this one
                    ipsp_start_pts = max(0, group_pts - time_if_no_ipsp_pts)
                else:
                    ipsp_start_pts = ipsp_start_previous
                if ipsp_next == -1:
                    # assign different end - var_ms before this one
                    ipsp_end_pts = min(group_pts + time_if_no_ipsp_pts, np.size(data_trace, 1))
                else:
                    ipsp_end_pts = ipsp_end_next
                
                #if spw_no == 7 and electr == 5 and group == 8:
                #    import pdb; pdb.set_trace()     
                # find where there is maximum between new ipsp and previous Ipsp
                
                if ipsp_start_pts+1 < group_pts and group_pts < ipsp_end_pts -1:
                    try:
                        x_a = np.argmax(data_trace[electr,ipsp_start_pts:group_pts]) + ipsp_start_pts
                    except:
                        import pdb; pdb.set_trace()    
                    x_b = group_pts
                    x_c = np.argmax(data_trace[electr,group_pts:ipsp_end_pts]) + group_pts
                    
                    y_a = max(data_trace[electr,ipsp_start_pts:group_pts])
                    y_b = data_trace[electr, group_pts]
                    y_c = max(data_trace[electr,group_pts:ipsp_end_pts])
                    
                    # move so that max of first IPSP is start of coordinating system
                    y_I_c = y_c - y_a
                    y_I_b = y_b - y_a
                    y_I_a = y_a - y_a
                    
                    x_I_c = (x_c - x_a) * 1.0 # length of 'triangle'
                    x_I_b = (x_b - x_a) * 1.0
                    x_I_a = (x_a - x_a) * 1.0
                    
                    # calculate the proportion and where is c at x_I_b on the line
                    # crossing between a_I and c_I
                    y_prime_b = (x_I_b * y_I_c) / x_I_c
                    fraction = x_I_b/x_I_c 
                    # check if new x is at smaller than the 0.5 of the expected one
                    allow = False
                    if y_prime_b * 0.7 > y_I_b:
                        # check if one of the maximums does not lay on the x
                        if fraction > 0 and fraction < 1:
                            # check if at least one IPSP (from right or left) is higher than
                            # minimum IPSP height
                            if y_I_c > min_ipsp_height:
                                allow = True
                else:
                    allow = False
                
                #make decision if to add this electrode to this group or not
                #import pdb; pdb.set_trace() 
                # check if the y_I_b is at least twice smaller than y_prime_b
                if allow:
                    # this IPSP can be included
                #if (after_ok + before_ok) > 0:
                    #import pdb; pdb.set_trace()
                    # this IPSP is chosen to add to the given group; find best time
                    data_next_ipsps = data_trace[electr,group_pts:ipsp_end_pts]
                    
                    potential_ipsp = data_trace[electr, group_pts-ipsp_ae_pts: group_pts+ipsp_ae_pts]
                    ipsp_time_pts = np.argmin(potential_ipsp) + group_pts - ipsp_ae_pts
                    new_ipsp_time = pts2ms(ipsp_time_pts, fs).astype('f8')
                    new_ipsp = spw_used[spw_used['group'] == group][0]
                    new_ipsp['electrode'] = electr.astype('i4')
                    new_ipsp['ipsp_start'] = new_ipsp_time
                    new_amplitude = max(max(data_next_ipsps) - np.min(potential_ipsp), 0)
                    new_ipsp['amplitude'] = new_amplitude
                    
                    all_new_ipsps.append(new_ipsp)
                    if plot_it:
                        plot_ipsp.append(new_ipsp)
                # if add find the correct time
                
            # find the max of each ipsp in each electrode
        
        if plot_it and spw_no == 7:
            data_used = data_trace[:,spw_start_pts:spw_end_pts]
            plt.figure()
            add_it = 150          
            
            t = dat.get_timeline(data_used[0, :], fs, 'ms')
            for electr in range(len(data_used)):
                plt.plot(t, data_used[electr] + electr * add_it)

                new_ones = np.concatenate([spw_used, plot_ipsp])
                ipsps_used = new_ones[new_ones['electrode'] == electr]['ipsp_start']
                ipsps_to_plot = ms2pts(ipsps_used - spw_start - win[0], fs).astype('i4')
                plt.plot(t[ipsps_to_plot], data_used[electr, ipsps_to_plot] + electr * add_it, 'ro', ms = 3)

                    
                ipsps_used = spw_used[spw_used['electrode'] == electr]['ipsp_start']
                ipsps_to_plot = ms2pts(ipsps_used - spw_start - win[0], fs).astype('i4')
                plt.plot(t[ipsps_to_plot], data_used[electr, ipsps_to_plot] + electr * add_it, 'go', ms = 3)
                #import pdb; pdb.set_trace()

            plt.show()
    new_ones = np.concatenate([spws, all_new_ipsps])           
    
    np.savez(save_folder + save_file, spw_ipsps = new_ones)          
           

def update_remove_with_to_few_ipsps(save_folder, save_file, spw_file, to_remove):
    npzfile        = np.load(save_folder + spw_file)
    spws = npzfile['spw_ipsps']
    npzfile.close()    
    #import pdb; pdb.set_trace()
    assert len(np.unique(spws[['trace','spw_start']])) == len(np.unique(spws['spw_no']))
    chosen_spws = []
    for spw_no in np.unique(spws['spw_no']):
        spw_used = spws[spws['spw_no'] == spw_no]
        no_groups = len(np.unique(spw_used['group']))
        if no_groups >= to_remove:
            chosen_spws.append(spw_used)
    
    chosen_spws = np.concatenate(chosen_spws)
    np.savez(save_folder + save_file, spw_ipsps = chosen_spws)



def fill_gap_in_all_groups(spws):
    """ it takes the data given and fills all the places where 
    the groups are not in neighbouring electrodes"""
    
    #import pdb; pdb.set_trace()
    new_ipsps = []

    for spw_no in np.unique(spws['spw_no']):
        spw_used = spws[spws['spw_no'] == spw_no]
        for group in np.unique(spw_used['group']):
            group_used = spw_used[spw_used['group'] == group]
            electrodes = group_used['electrode']
            
            # check if all the electrodes in the row are defined for this group
            luki = np.diff(electrodes)
            #import pdb; pdb.set_trace()
            if np.any(luki > 1):
                #import pdb; pdb.set_trace()
                luki = np.diff(electrodes)
                new_electrodes = np.arange(electrodes[0], electrodes[-1] + 1)
                missing_electrodes = np.setdiff1d(new_electrodes, electrodes).astype('i4')
                mis_no = len(missing_electrodes)
                
                # detect all the variables for this new IPSPs
                groups = (np.ones(mis_no) * group).astype('i4')
                trace = group_used['trace'][0]
                
                traces = (np.ones(mis_no) * trace).astype('i4')
                
                spw_no = group_used['spw_no'][0]
                spw_nos = (np.ones(mis_no) * spw_no).astype('i4')
                
                spw_start = group_used['spw_start'][0]
                spw_starts = (np.ones(mis_no) * spw_start).astype('f8')
                
                #spw_end = group_used['spw_end'][0]
                #spw_ends = (np.ones(mis_no) * spw_end).astype('f8')
                
                ipsps_no = max(group_used['ipsp_no'])
                ipsps_nos = (np.ones(mis_no) * ipsps_no).astype('i4')
                
                ipsp_start = min(group_used['ipsp_start'])
                ipsp_starts = (np.ones(mis_no) * ipsp_start).astype('f8')
                
                temp_ampls = (np.zeros(mis_no) * ipsps_no).astype('f8')
                new_ipsps.append(np.rec.fromarrays([missing_electrodes, traces, 
                                                        spw_nos, spw_starts,
                                                        ipsps_nos, ipsp_starts, temp_ampls, groups],
                                                           names='electrode,trace, spw_no, spw_start, ipsp_no, ipsp_start, amplitude, group'))   
                
    # make sure that all the variables used are of the same type
    new_dtype = np.dtype([('electrode', '<i4'), 
                          ('trace', '<i4'),
                          ('spw_no', '<i4'), 
                          ('spw_start', '<f8'),
                          ('ipsp_no', '<i4'), 
                          ('ipsp_start', '<f8'),
                          ('amplitude', '<f8'),
                          ('group', '<i4')])
   
    new_ipsps.append(spws.astype(new_dtype))
    all_ipsps = np.concatenate(new_ipsps)
    return all_ipsps

    
def update_fill_gap_between_ipsp_groups(save_folder, save_file, spw_file, data_file):
    """ takes all the groups of the spws and IPSPs and fills the electrode which
    did not have ipsp detected but electrodes from both sides have this group"""
    npzfile        = np.load(save_folder + spw_file)
    spws = npzfile['spw_ipsps']
    npzfile.close()   
    

    #import pdb; pdb.set_trace()
    assert len(np.unique(spws[['trace','spw_start']])) == len(np.unique(spws['spw_no']))    
    
    npzfile        = np.load(save_folder + data_file)
    data = npzfile['data']
    fs = npzfile['fs']
    npzfile.close()       
    all_ipsps = fill_gap_in_all_groups(spws)

    #import pdb; pdb.set_trace()
    for trace in range(np.size(data,1)):
        data_trace = data[:,trace,:]

        new_ipsps = all_ipsps[all_ipsps['trace'] == trace]
        amplitude = calculate_amplitude_of_IPSP(new_ipsps, data_trace, fs)
        all_ipsps['amplitude'] [all_ipsps['trace'] == trace]= amplitude
        #import pdb; pdb.set_trace()
        
        
    assert len(np.unique(spws[['trace','spw_start']])) == len(np.unique(spws['spw_no'])) 
    np.savez(save_folder + save_file, spw_ipsps = all_ipsps)
    
    
def update_equalize_number_spws(save_folder, save_file, induc_spont, load_distances):
    """ it takes the same number of both induc and spontaneous spws"""
    npzfile    = np.load(save_folder + load_distances)
    distances      = npzfile['dist_spwspike'] 
    npzfile.close()     
    
    npzfile    = np.load(save_folder + induc_spont)
    spont = npzfile['spontaneous']
    init = npzfile['initiated']
    npzfile.close()     
    
    way = 'random' # 'random', 'closest', 
    print 'setting the same number or spontaneous and induced spws'
    if way == 'random':
        # take randomly from both groups the number of elements equal smaller group
        
        no_elements = min(len(np.unique(spont['spw_no'])), len(np.unique(init['spw_no'])))
        #import pdb; pdb.set_trace()
        # choose init elements
        
        if len(np.unique(spont['spw_no'])) == no_elements:
            # there is less spontaneous events
            to_correct = init            
        else:
            # there is less induced events
            to_correct = spont
            
        chosen_ones = take_random_elements(np.unique(to_correct['spw_no']), no_elements)
        chosen_ones = np.sort(chosen_ones)
        selected_ones = []
        #import pdb; pdb.set_trace()
        for spw_no in chosen_ones:
            selected_ones.append(to_correct[to_correct['spw_no'] == spw_no])
            
        if len(np.unique(spont['spw_no'])) == no_elements:
            try:
                selected_init = np.concatenate(selected_ones)
            except:
                import pdb; pdb.set_trace()
            selected_spont = spont
        else:
            # there is less induced events         
            selected_spont = np.concatenate(selected_ones)
            selected_init = init
#    elif way == 'closest':
        
    
    np.savez(save_folder + save_file, initiated = selected_init, spontaneous = selected_spont)
    
def update_induc_spont_spw(save_folder, save_file, load_distances, load_spwfile, max_dist, ext):
    """ checks which spws are initiated and which are sponteneaus"""
    npzfile    = np.load(save_folder + load_distances)
    
    error_allowed = max_dist[0]
    distances      = npzfile['dist_spwspike'] 
    npzfile.close()    
    npzfile    = np.load(save_folder + load_spwfile)
    #ipsps      = npzfile['ipsps']  
    ipsps = npzfile['spw_ipsps']
    #import pdb; pdb.set_trace()
    npzfile.close()   
    #before_pts = ispw.ms2pts(win[0], fs)
    #after_pts = ispw.ms2pts(win[1], fs)
    dist = distances['distance']
    print "checking which spws are initiated and which are sponteneaus"
    initiated_no = distances[(dist <= max_dist[1]) & (dist>=error_allowed)]['spw_no']
    #import pdb; pdb.set_trace()
    spont_no = distances[(dist > max_dist[1]) | (dist<error_allowed)]['spw_no']
    #unsure = unsure_spws
     
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
    
    error_allowed = -0.5
    #allow_before = ms2pts(allow_before,fs)
    # load starts of spws
    npzfile         = np.load(save_folder + load_spwfile)
    #import pdb; pdb.set_trace()
    ipsps      = npzfile['spw_ipsps']  
    npzfile.close()    
    
    spw_no = []
    dist_all = []    
    electrodes = []
    traces = []
    print 'checking distance from intracellular spike '
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
                    try:
                        assert len(np.unique(spw_numb))==1
                    except:
                        import pdb; pdb.set_trace()
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

def group_ipsps(spw_ipsps_trace, shift_ipsp):
    """Assign coincident spikes (closer than shift_ipsp) to one group
    and return the group's index"""
    
    def _find_repeated_group_ids(ids):
        sorted_ids = np.sort(ids)
        repeated = (np.diff(sorted_ids)==0)
        repeated_ids = np.unique(sorted_ids[repeated])
        return repeated_ids[repeated_ids>=0]
    
    
    order = np.argsort(spw_ipsps_trace['ipsp_start'])
    spw_ipsps_sorted = spw_ipsps_trace[order]
    
    group_ids = np.empty(len(spw_ipsps_sorted), dtype=np.int)
    group_ids.fill(-1)
    electrodes = spw_ipsps_sorted['electrode']
    ipsp_start = spw_ipsps_sorted['ipsp_start']
    
    
    first_el = np.min(electrodes)

    group_ids[electrodes==first_el] = np.arange(np.sum(electrodes==first_el))
    #print np.min(electrodes)
    #print electrodes
    #assert (electrodes>0).all()
    max_group = group_ids.max()
    for el in np.unique(electrodes)[1:]:
        ipsps_to_assign, = np.where((electrodes==el))
        ipsps_assigned, = np.where(group_ids>-1)
        while len(ipsps_to_assign) > 0:
            max_group = group_ids.max()

            if len(ipsps_assigned)==0:
                group_ids[ipsps_to_assign] = np.arange(len(ipsps_to_assign))+max_group+1
                break
             
            time_assigned = ipsp_start[ipsps_assigned]
            time_to_assign = ipsp_start[ipsps_to_assign]
            
            assert (np.diff(time_assigned.argsort(kind='mergesort'))==1).all()
            i = np.searchsorted(time_assigned, time_to_assign)
            #compare times with ipsp to the left and right
            left = time_assigned[np.maximum(0,i-1)]
            right = time_assigned[np.minimum(i, len(time_assigned)-1)]
            closer_to_left = (time_to_assign-left)< (right-time_to_assign)
            i[closer_to_left] = i[closer_to_left]-1
            i = np.minimum(np.maximum(i,0), len(time_assigned)-1)
            dist = np.abs((time_assigned[i]-time_to_assign))
            
            new_group_ids = -1*np.ones(len(ipsps_to_assign), dtype=int) 
            
            new_group_ids[dist<=shift_ipsp] = group_ids[ipsps_assigned[i[dist<=shift_ipsp]]]
            new_group_ids[dist>shift_ipsp] = np.arange(np.sum(dist>shift_ipsp))+max_group+1
            
            #assure that only one ipsp is assigned to any given group
            rep_gids = _find_repeated_group_ids(new_group_ids)
            for gid in rep_gids:
                idx_in_group, =  np.where(new_group_ids==gid)
                idx_in_electrode, = np.where(group_ids==gid)
                closest_electrode = idx_in_electrode[np.argmax(electrodes[idx_in_electrode])]
            
                dist_to_electrode = np.abs(time_to_assign[idx_in_group]-ipsp_start[closest_electrode])
                i_closest_ipsp = np.argmin(dist_to_electrode)
                new_group_ids[idx_in_group] = -1
                new_group_ids[idx_in_group[i_closest_ipsp]] = gid
                
            group_ids[ipsps_to_assign] = new_group_ids
           
            #remove ipsps from assign list that were already used for this electrode
            was_used  = np.in1d(group_ids[ipsps_assigned], new_group_ids)
            ipsps_assigned = ipsps_assigned[~was_used] 
            ipsps_to_assign = ipsps_to_assign[group_ids[ipsps_to_assign]<0]
            #if closest ipsps is closer thas shift_ipsp assign the new ipsps to the same group
           
    #print group_ids
    assert (group_ids>-1).all()
    
    #sort_order = np.argsort(spw_ipsps_trace['ipsp_start'])
    #spw_ipsps_sorted = spw_ipsps_trace[sort_order]
    #spw_ipsps_intervals = np.diff(spw_ipsps_sorted['ipsp_start'])
    inverse_i = np.argsort(order)
    return group_ids[inverse_i]

def count_coincident_ipsps(spw_ipsps_trace, shift_ipsp):
    """counts number of IPSPs in different electrodes seperated by not more
    than shift_ipsp"""
    #if spw_ipsps_trace['spw_start'][0] == 149.3:
    #    import pdb; pdb.set_trace()
    i = np.argsort(spw_ipsps_trace['ipsp_start'])
    spw_ipsps_sorted = spw_ipsps_trace[i]
        
    group_idx = group_ipsps(spw_ipsps_sorted, shift_ipsp)

        
    electrode = spw_ipsps_sorted['electrode']

    #count number of disitinct electrodes in each IPSP group
    n_uniq_electrodes_per_group = [len(np.unique(electrode[group_idx==group])) 
                                   for group in range(group_idx.max()+1)]
    
    n_uniq_electrodes_per_group = np.array(n_uniq_electrodes_per_group)
    
    #assign to each IPSP number of coincident IPSPs in other electrodes
    n_electrodes_per_ipsp = n_uniq_electrodes_per_group[group_idx]

    #inverse sorting
    inverse_i = np.argsort(i)
    return n_electrodes_per_ipsp[inverse_i]


def calculate_ipsp_amplitude(ipsp_trace, data_trace, fs):
    """ checks how 'high' is each of the beginning of ipsps"""
    #import pdb; pdb.set_trace()
    amplitudes = np.zeros(len(ipsp_trace))
    for idx, time in enumerate(ipsp_trace['ipsp_start']):
        time_pts = ms2pts(time, fs).astype('i4')
        electr = ipsp_trace[idx]['electrode']
        amplitudes[idx] = data_trace[electr, time_pts]
    return amplitudes

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
    
    #inverse sorting
    inverse_i = np.argsort(i)
    return  ampl_rise[inverse_i]
       
    

def calculate_amplitude_of_IPSP(spw_ipsps_trace, data_trace, fs):
    """calculate the maximum amplitude of each IPSP, for last IPSP return 0"""
   
    sort_order = np.argsort(spw_ipsps_trace, order=['electrode', 'ipsp_start'])
    spw_ipsps_sorted = spw_ipsps_trace[sort_order]
    
    time_pts = ms2pts(spw_ipsps_sorted['ipsp_start'], fs).astype(np.uint32)
    
    ipsps_starts = time_pts[:-1]
    ipsps_ends = time_pts[1:]
    
    last_ipsp = np.append(np.diff(spw_ipsps_sorted['electrode'])!=0, [True])
    
    ipsp_ampl = np.zeros(len(spw_ipsps_sorted))
    idxs, = np.where(~last_ipsp)
    for i in idxs:
        first, last = ipsps_starts[i], ipsps_ends[i]
        electrode = spw_ipsps_sorted['electrode'][i]
        single_ipsp = data_trace[electrode, first:last]
        try:
            #import pdb; pdb.set_trace()
            ipsp_ampl[i] = np.max(single_ipsp)-data_trace[electrode, first]
        except:
            ipsp_ampl[i] = np.nan
    
    # return to the original order
    inverse_i = np.argsort(sort_order)       
    return  ipsp_ampl[inverse_i]
    
def define_spikes_closest2IPSP_starts(spikes, ipsps):
    """ finds which spikes are the closest to the beginnings of IPSP but no 
    further than shift_spike"""
    
    i_ipsp = np.argsort(ipsps, order=['spw_no', 'electrode', 'ipsp_start'])
    ipsps_sorted = ipsps[i_ipsp]
    #i_spike = np.argsort(spikes, order=['spw_no', 'electrode', 'time'])
    #spikes_sorted = spikes[i_spike]    
    
    closest_events = []
    for t in ipsps_sorted:
        #if t['ipsp_start'] == 424.25:
        #    import pdb; pdb.set_trace()
            
        electr_index = np.arange(0,len(spikes))
        min_idx = spikes['electrode']==t['electrode']
        possible_electr = electr_index[min_idx]
        if len(possible_electr) > 0:
            i = np.argsort(np.abs(t['ipsp_start']-spikes[possible_electr]['time']))
            i = possible_electr[i]
        
            #&             (spikes['spw_no']==t['spw_no'])]
            closest = i[0] if len(i)>0 else np.nan
        else:
            closest = np.nan
        closest_events.append(closest)
        
    closest_events = np.array(closest_events)
    #s = take_element_or_nan(ipsps_sorted, np.array(closest_events))
    #spikes_sorted[closest_events]
    inverse_i = np.argsort(i_ipsp)       
    return closest_events[inverse_i]

def take_element_or_nan(x, i):
    y = np.empty(len(i), dtype='f8')
    y.fill(np.nan)
    isnan = np.isnan(i)
    not_nan_i = i[~isnan].astype('i4')
    #import pdb; pdb.set_trace()
    y[~isnan]=x[not_nan_i]
    return y

def add_rec_field(recarray, arrs, field_names):
    #import pdb; pdb.set_trace()
    names = list(recarray.dtype.names)
    field_data = [recarray[d] for d in names]
    if np.ndim(arrs) > 1:
        field_data += arrs
    else:
        field_data.append(arrs)
    new_names = ','.join(names+field_names)
    #import pdb; pdb.set_trace()
    return np.rec.fromarrays(field_data, names=new_names)

def shift_ipsp_start(ipsps_trace, spikes_trace, shift_spike):
    """shift ipsp either to closest spike or two largerst IPSP in the group.
    returns new rec array with IPSPs shifted in time"""
    
    group_ids = ipsps_trace['group']
    
    # check which spikes are the closest to found beginnings of IPSPs 
    closest_spike_idx = define_spikes_closest2IPSP_starts(spikes_trace, 
                                                      ipsps_trace)
    
    closest_spike_idx = closest_spike_idx.astype('f4')
    closest_spike_time = take_element_or_nan(spikes_trace['time'], closest_spike_idx)

    dist_to_spike = np.abs(ipsps_trace['ipsp_start'] - closest_spike_time) 
    closest_spike_idx[dist_to_spike>shift_spike] = np.nan
    
    amplitudes = ipsps_trace['amplitude']
    ipsp_start = ipsps_trace['ipsp_start']
    
    new_ipsp_start = np.empty(len(ipsps_trace))
    new_ipsp_start.fill(np.nan)
    
    for gid in np.unique(group_ids):
        ipsps_idx, = np.where(group_ids==gid)
        d = closest_spike_idx[ipsps_idx]
        
        if np.isnan(d).all():
            #did not find close spike, find IPSP with largest amplitude
            i = np.argmax(amplitudes[ipsps_idx])
            t = ipsp_start[ipsps_idx[i]] 
        else:
            #found at least one spike, shift to spike closest to largest IPSP
            ipsps_close_to_spikes = ipsps_idx[~np.isnan(d)]
            i = np.argmax(amplitudes[ipsps_close_to_spikes])
            t = closest_spike_time[ipsps_close_to_spikes[i]]
        new_ipsp_start[ipsps_idx] = t 
        
    assert ~np.isnan(new_ipsp_start).any()
    new_ipsps_trace = ipsps_trace.copy()
    new_ipsps_trace['ipsp_start'] = new_ipsp_start
    
    return new_ipsps_trace

def shift_spw_to_first_ipsp(ipsps_trace, min_distance_between):
    ipsps_trace = ipsps_trace.copy()
    #import pdb; pdb.set_trace()
    spw_no_used = 0
    spw_start_old = -min_distance_between -1
    for spw_no in np.unique(ipsps_trace['spw_no']):
        spw_start = ipsps_trace[ipsps_trace['spw_no'] == spw_no]['ipsp_start'].min()
        spw_starts = np.ones(len(ipsps_trace[ipsps_trace['spw_no'] == spw_no]['ipsp_start']))
        spw_starts = spw_starts*spw_start
        if (spw_starts[0] - spw_start_old) < min_distance_between:
            spw_starts = spw_starts * 0 + spw_start_old
            ipsps_trace['spw_start'][ipsps_trace['spw_no'] == spw_no] = spw_starts.astype('f8')
            ipsps_trace['spw_no'][ipsps_trace['spw_no'] == spw_no] = spw_no_used.astype('i4')
            #import pdb; pdb.set_trace()
        else:
            spw_no_used = spw_no.astype('i4')
            ipsps_trace['spw_start'][ipsps_trace['spw_no'] == spw_no] = spw_starts.astype('f8')
        spw_start_old = spw_starts[0]
        #ipsps_trace[ipsps_trace['spw_no'] == spw_no]['spw_start'] = spw_starts
        
        #import pdb; pdb.set_trace()
    return ipsps_trace

def calculate_max_in_given_patch(data, points, distanse_from_point, fs):
    """ calculates the amplitude of the IPSPs given"""
    index = np.argsort(points, order=['electrode', 'ipsp_start'])
    points_ordered = points[index]
    electr = points_ordered['electrode'][0]
    points_pts = ms2pts(points_ordered['ipsp_start'], fs).astype('i4')
    distance_pts = ms2pts(distanse_from_point, fs)
    ipsp_maxs = np.zeros(len(points_ordered))
    
    for idx, pt_idx in enumerate(points_pts):
        # check what is closer - next IPSP or given number of points
        remove_ipsp = False
        electr = points_ordered['electrode'][idx]
        if idx < len(points_pts) - 1 and electr == points_ordered['electrode'][idx+1]:
            # not the last IPSP and it's the same electrode as the next one
            next_ipsp = points_pts[idx+1] - pt_idx
            if next_ipsp <= 1:
                remove_ipsp = True
                #import pdb; pdb.set_trace()
            #assert next_ipsp > 1 # it's the same ipsp as previously, this IPSP should be removed?
        else:
            next_ipsp = distance_pts + 1
        #import pdb; pdb.set_trace()
        if not remove_ipsp:
            ipsp_end =  min(distance_pts, next_ipsp)
            
            ipsp_maxs[idx] = np.max(data[electr, pt_idx:pt_idx + ipsp_end]) - data[electr, pt_idx]
            
            
        if ipsp_maxs[idx] < 0 or remove_ipsp: #the maximum amplitude of this IPSP is smaller then 0
            # which would give very occured IPSP - it must be removed! (give -1)
            ipsp_maxs[idx] = -1
        
    index_reversed = np.argsort(index)
    return ipsp_maxs[index_reversed]

def calc_distance_between(spw_points, min_dist):
    """ calculates distance between given IPSPs and removes boolean structure
    saying which are long enough and which are not"""
    
    index = np.argsort(spw_points, order=['electrode', 'ipsp_start'])
    points_ordered = spw_points[index]    
    
    dist_electr = []
    last_idx = 0
    for electr in np.unique(points_ordered['electrode']):
        #import pdb; pdb.set_trace()
        dist = np.diff(points_ordered[points_ordered['electrode'] == electr]['ipsp_start']).tolist()

        dist_electr = dist_electr + dist + [100]

    dist_electr = np.array(dist_electr)
    distances_alright = dist_electr > min_dist

    assert (dist_electr>=0).all()
    index_reversed = np.argsort(index)
    return distances_alright[index_reversed]

def remove_too_short(distance_too_short, ipsp_amplitudes):
    """ if the distance is too short, checks which amplitude is higher (from
    right or from the left and returns only this for further use"""
    ipsp_used = distance_too_short
    
    left_ampl = ipsp_amplitudes[:-1][distance_too_short[:-1]]
    right_ampl = ipsp_amplitudes[1:][distance_too_short[:-1]]
    
    comparision = left_ampl - right_ampl
    # it's negative - right is larger, if it's positive, left is larger
    
    #ipsp_used[]
    
    comparision[comparision < 0] = left_ampl[comparision < 0]
    comparision[comparision >= 0] = right_ampl
    
    import pdb; pdb.set_trace()

def update_spws_first_max(save_folder, spws, datafile, save_file, window = [-3, 3]):
    npzfile         = np.load(save_folder + datafile)
    data            = npzfile['data']
    fs              = npzfile['fs']
    npzfile.close()   
    
    npzfile         = np.load(save_folder + spws)
    spws      = npzfile['spw_ipsps']
    npzfile.close()       
    #win0 = ms2pts(window[0], fs)
    
    
#    traces = []
    for spw_no in np.unique(spws['spw_no']):
        print 
        print spws['spw_start'][spws['spw_no'] == spw_no]
        #import pdb; pdb.set_trace() 
        spw_used = spws[spws['spw_no'] == spw_no]
        spw_used = np.sort(spw_used, order = 'electrode')
        trace = spws['trace'][0]

        spw_start = spw_used['spw_start'][0] + window[0]
        spw_end = spw_used['spw_start'][0] + window[1]
        spw_start_pts = ms2pts(spw_start, fs).astype('i4')
        spw_end_pts = ms2pts(spw_end, fs).astype('i4')
        
        data_used = data[:,trace,spw_start_pts:spw_end_pts].copy()

        electr_max = np.argmax(np.max(data_used, axis = 1))
        #import pdb; pdb.set_trace()
        peak = spw_start_pts + np.argmax(data_used[electr_max, :])
        new_start = pts2ms(peak, fs).astype('f8')

        spws['spw_start'][spws['spw_no'] == spw_no] = np.ones(len( spws['spw_start'][spws['spw_no'] == spw_no])) * new_start
        print spws['spw_start'][spws['spw_no'] == spw_no]
    np.savez(save_folder + save_file, spw_ipsps = spws) 

def update_remove_too_small_spws(load_datafile, load_spwsipsp, min_ampl, save_folder, save_file, ext):
    """ checks which SPWs in it's highest point are smaller than min_ampl and removes those"""
    # load all the necessary data
    npzfile         = np.load(save_folder + load_datafile)
    data            = npzfile['data'] # data
    fs              = npzfile['fs'] # sampling frequency
    npzfile.close()      

    npzfile         = np.load(save_folder + load_spwsipsp)
    spws       = npzfile['spw_ipsps'] # detected so far spws and IPSPs
    npzfile.close()  

    # order spws according first to trace, and then to their start
    i = np.argsort(spws, order=['trace', 'spw_start'])
    spws = spws[i]

    ipsp_length = 10 # ms
    spw_no_all = np.unique(spws['spw_no'])
    
    plot_it = False

    print "Removing SPWs of amplitude smaller than: " + str(min_ampl) + "; working on SPW number: "
    
    new_spw_no = 0
    new_spws = []
    
    
    for idx, spw_no in enumerate(spw_no_all):
        print spw_no
        
        # get parameters for this spw
        spw_used = spws[spws['spw_no'] == spw_no]
        trace = spw_used['trace'][0]
        spw_start = spw_used['spw_start'][0]
        
        # check when is next SPW
        if idx < len(spw_no_all) - 1:
            # not the end of all the SPWs yet
            
            next_spw_used = spws[spws['spw_no'] == spw_no_all[idx+1]]
            if trace == next_spw_used['trace'][0]:
                # not the end of this trace yet
                # check distance between last IPSP of this and before the next IPSP
                distance = next_spw_used['spw_start'][0] - max(spw_used['ipsp_start'] )
                assert distance > 0
                
                # if distance is smaller from 0 it means that new SPW start before the old SPW ends - WRONG!
                if distance < ipsp_length:
                    add_to_end = distance 
                else:
                    # if next SPW is too far away, just use ipsp_length    
                    add_to_end = ipsp_length
            else:
                # new trace, use ipsp_length
                add_to_end = ipsp_length
        else:
            # use ipsp_length
            add_to_end = ipsp_length
                
        spw_end = max(spw_used['ipsp_start']) + add_to_end        
        
        spw_start_pts = ms2pts(spw_start, fs).astype('i4')
        spw_end_pts = ms2pts(spw_end, fs).astype('i4')
        
        data_trace = data[:, trace, spw_start_pts: spw_end_pts]
        
        max_spw = np.max(data_trace,1)
        max_arg = np.argmax(data_trace,1)
        max_idx = np.argmax(max_spw)
        max_spw = max_spw[max_idx]
        max_arg = [max_idx, max_arg[max_idx]]
        #import pdb; pdb.set_trace()
        
        if max_spw >= min_ampl:
            # save this spw and assign new spw_no to it
            spw_used['spw_no'] = np.ones(len(spw_used['spw_no'])) * new_spw_no
            new_spws.append(spw_used)
            new_spw_no = new_spw_no + 1
        
        # plot it if necessary
        if plot_it and max_spw >= min_ampl:
            plt.figure()
            add_it = 150
            t = dat.get_timeline(data_trace[0, :], fs, 'ms')
            for electr in range(len(data_trace)):
                plt.plot(t, data_trace[electr, :] + electr * add_it)
                
                if electr == max_arg[0]:
                    #import pdb; pdb.set_trace()
                    plt.plot(t[max_arg[1]], data_trace[electr, max_arg[1]] + electr * add_it, 'ro', ms = 8)
                    plt.text(t[max_arg[1]], data_trace[electr, max_arg[1]] + electr * add_it, str(max_spw))
                
            plt.show()
        
    new_spws = np.concatenate(new_spws)
    np.savez(save_folder + save_file, spw_ipsps = new_spws)    
    

def corect_ipsps(load_datafile, load_spwsipsp, load_spwsspike, save_folder, save_fig, save_file,ext):
    """ checks all the chosen preliminary ipsps and chooses the ones to use """
    
    # load all the necessary data
    npzfile         = np.load(save_folder + load_datafile)
    data            = npzfile['data'] # data
    fs              = npzfile['fs'] # sampling frequency
    npzfile.close()   
    
    npzfile         = np.load(save_folder + load_spwsipsp)
    spw_ipsps       = npzfile['spw_ipsps'] # detected so far spws and IPSPs
    npzfile.close()   

    npzfile         = np.load(save_folder + load_spwsspike)
    spw_spike      = npzfile['spike_idx']  # detected spikes
    npzfile.close()     
    
    plot_it = False
    
    shift_spike= 1      # ms - if there is spike close by, beginning of IPSP will be shifted to it
    min_length_ipsp = 2 # ipsp cannot be shorter than this
    
    print "correcting IPSPs, trace:"
    all_traces= np.unique(spw_ipsps['trace'])
    spw_ipsps_all = []
    
    # treat all the IPSPS
    for trace in all_traces:
        print str(trace) + ', '
        spw_ipsps_trace = spw_ipsps[spw_ipsps['trace']==trace]
        spikes_trace = spw_spike[spw_spike['trace']==trace] # spikes are found properly
        
        # check which spikes are the closest to each ipsp
        closest_spike_idx = define_spikes_closest2IPSP_starts(spikes_trace, 
                                                      spw_ipsps_trace) 
        
        #import pdb; pdb.set_trace()
        closest_spike_times = take_element_or_nan(spikes_trace['time'], closest_spike_idx)
        spikes_close_enough = np.abs(closest_spike_times - spw_ipsps_trace['ipsp_start']) < shift_spike
        spw_ipsps_trace['ipsp_start'][spikes_close_enough] = closest_spike_times[spikes_close_enough]
        
        #import pdb; pdb.set_trace()
        distance_between = calc_distance_between(spw_ipsps_trace[['electrode','ipsp_start']], min_length_ipsp)
        spw_ipsps_trace = spw_ipsps_trace[distance_between]

        spw_ipsps_all.append(spw_ipsps_trace)

    spw_ipsps_all = np.concatenate(spw_ipsps_all)
    
    # for testing purposes it is possible to plot the chosen trace of the data
    if plot_it:
        plt.figure()
        trace = 2
        data_used = data[:,trace,:]
        add_it = 150
        t = dat.get_timeline(data_used[0,:], fs, 'ms')
        for electr in range(len(data_used)):
            #import pdb; pdb.set_trace()
            ipsp_new_temp = spw_ipsps_all[(spw_ipsps_all['electrode'] == electr) & (spw_ipsps_all['trace'] == trace)]
            ipsp_new = ipsp_new_temp['ipsp_start']
            ipsp_new_pts = ms2pts(ipsp_new, fs).astype('i4')
            
            ipsp_old = spw_ipsps[(spw_ipsps['electrode'] == electr) & (spw_ipsps['trace'] == trace)]['ipsp_start']
            ipsp_old_pts = ms2pts(ipsp_old, fs).astype('i4')
            
            plt.plot(t, data_used[electr, :] + electr * add_it)
            #import pdb; pdb.set_trace()
            
            plt.plot(t[ipsp_new_pts], data_used[electr, ipsp_new_pts]+ electr * add_it, 'go', alpha = 0.6, markersize=8)
            plt.plot(t[ipsp_old_pts], data_used[electr, ipsp_old_pts]+ electr * add_it, 'r<', alpha = 0.6)
            
            spikes_trace = spw_spike[spw_spike['trace']==trace]
            spike_times = spikes_trace[(spikes_trace['electrode'] == electr)]['time']

            spike_pts = ms2pts(spike_times, fs).astype('i4')
            plt.plot(t[spike_pts], data_used[electr, spike_pts]+ electr * add_it, 'y*') 
    plt.show()     
    #import pdb; pdb.set_trace()
    np.savez(save_folder + save_file, spw_ipsps = spw_ipsps_all)    


def update_spws_beg(load_datafile, load_spwsipsp, load_spwsspike, save_folder, save_fig, save_file,ext, win = [-20, 80]):
    """ checks all the ipsps and corrects them for each spw"""
    
    # load the necessary data
    npzfile         = np.load(save_folder + load_datafile)
    data            = npzfile['data'] # data
    fs              = npzfile['fs'] # sampling frequency
    npzfile.close()   
    
    npzfile         = np.load(save_folder + load_spwsipsp)
    spw_ipsps       = npzfile['spw_ipsps'] # preliminary IPSPs and SPWs
    npzfile.close()   
    
    npzfile         = np.load(save_folder + load_spwsspike)
    spw_spike      = npzfile['spike_idx'] # detected spikes
    npzfile.close()     
    
    # define parameters
    distanse_from_point = 5 # ms - length of IPSP
    shift_ipsp = 2.5 # ms
    min_electr_first = 2 # on how many electrodes IPSP should be detected for the first ipsp (beginning of SPW)
    expected_min_ipsp_ampl = 25 # microV
    shift_spike= 1 #ms
    plot_it = False
    all_ipsps = []
    
    print "correcting beginning of IPSPs"
    
    all_traces= np.unique(spw_ipsps['trace'])

    # treat all the IPSPS
    for trace in all_traces:
        spw_ipsps_trace = spw_ipsps[spw_ipsps['trace']==trace]
        data_trace = data[:,trace,:] # use only data for this trace
        #import pdb; pdb.set_trace()
        old_last_ipsp = -1
        # analize every spw separately
        old_start = -1
        for spw_no in np.unique(spw_ipsps_trace['spw_no']):
            
            # calculate the amplitude for each IPSP (from start to distance_from_point after)
            spw_ipsps_spw_no = spw_ipsps_trace[spw_ipsps_trace['spw_no'] == spw_no]
            max_ampls = calculate_max_in_given_patch(data_trace, spw_ipsps_spw_no[['electrode','ipsp_start']], distanse_from_point, fs)
            
            # check which IPSPs are of too low amplitude
            spw_ipsps_first = spw_ipsps_spw_no[max_ampls >= expected_min_ipsp_ampl]
                 
            # check in how many electrodes each IPSP appears
            if len(spw_ipsps_first) > 0:
                n_electrodes_per_ipsp = count_coincident_ipsps(spw_ipsps_first, shift_ipsp)
                spw_ipsps_first = spw_ipsps_first[(n_electrodes_per_ipsp >= min_electr_first)]
            
            # is there any candidate for beginning of SPW?
            if len(spw_ipsps_first) > 0: #& min_ipsp_groups:
                #import pdb; pdb.set_trace()
                
                # move the first IPSPs to the closest spike or to the spike which is first
                group_ids = group_ipsps(spw_ipsps_first, shift_ipsp).astype('i4')
                
                #import pdb; pdb.set_trace()
                ipsp_amplitudes = calculate_max_in_given_patch(data_trace, spw_ipsps_first[['electrode','ipsp_start']], distanse_from_point, fs)
                spw_ipsps_first = add_rec_field(spw_ipsps_first, [ipsp_amplitudes, group_ids], ['amplitude', 'group'])                
                spikes_trace = spw_spike[spw_spike['trace']==trace]
                
                # save group which is first
                group_selected = spw_ipsps_first[spw_ipsps_first['ipsp_start'] == min(spw_ipsps_first['ipsp_start'])]['group'][0]
                values = spw_ipsps_first[spw_ipsps_first['group'] == group_selected]
                
                # move the ipsps to the beginning 
                spw_ipsps_first = shift_ipsp_start(spw_ipsps_first, spikes_trace, shift_spike)
                
                # set the old chosen first group to the new values
                for val in values:
                    val_ipsp_no = val['ipsp_no']
                    val_electr = val['electrode']
                    
                    spw_ipsps_spw_no[(spw_ipsps_spw_no['electrode'] == val['electrode']) & 
                                         (spw_ipsps_spw_no['ipsp_no'] == 
                                          val_ipsp_no)] = spw_ipsps_first[(spw_ipsps_first['electrode'] == val_electr) & 
                                         (spw_ipsps_first['ipsp_no'] == val_ipsp_no)]
                                          
                # set spw start to be the smalles ipsp_start from the chosen groups
                spw_ipsps_spw_no['spw_start'] = np.ones(len(spw_ipsps_spw_no)) * min(spw_ipsps_first['ipsp_start'])
                
                assert len(np.unique(spw_ipsps_first['spw_start'])) == 1
                #import pdb; pdb.set_trace()
                
                # remove all the IPSPs for this SPW which are before detected beginning
                ipsps_selected = spw_ipsps_spw_no[spw_ipsps_spw_no['ipsp_start'] >= spw_ipsps_spw_no['spw_start'][0]]
                proper_start = ipsps_selected['spw_start'][0]
                
                assert np.all(ipsps_selected['ipsp_start'] >= ipsps_selected['spw_start'][0])

                #import pdb; pdb.set_trace()
                all_ipsps.append(ipsps_selected)
                assert proper_start > old_start # make sure that this spw is after previous one
                try:
                    assert proper_start > old_last_ipsp # make sure that this spw is after last IPSP of previous one
                except:
                    import pdb; pdb.set_trace()
                    
                old_last_ipsp = max(ipsps_selected['ipsp_start'])
                old_start = proper_start

                if plot_it:
                    
                    
                    plt.figure()
                    add_it = 150
                    starts_pts = ms2pts(proper_start, fs)
                    win = [-20, 100]
                    win_pts = [ms2pts(win[0], fs), ms2pts(win[1], fs)]
                    data_used = data_trace[:, max(0,starts_pts + win_pts[0]): min(starts_pts + win_pts[1], np.size(data_trace, 1))]
                    
                    t = dat.get_timeline(data_used[0, :], fs, 'ms')
                    for electr in range(len(data_trace)):
                        plt.plot(t, data_used[electr] + electr * add_it)
                        try:
                            pass
                            #data_trace_filt2 = filt.bandPass([3, 300], fs, data_trace[electr,:])
                            #data_used_filt2 = data_trace_filt2[max(0,starts_pts + win_pts[0]): min(starts_pts + win_pts[1], np.size(data_trace, 1))]
                            #plt.plot(t, data_used_filt2+ electr * add_it)
#                            
                        except:
                            import pdb; pdb.set_trace()
                        spw_used = ipsps_selected[ipsps_selected['electrode'] == electr]
                        ipsps_all = ipsps_selected[ipsps_selected['electrode'] == electr]
                        ipsp_old = spw_ipsps_trace[spw_ipsps_trace['spw_no'] == spw_no]
                        ipsp_old = ipsp_old[ipsp_old['electrode'] == electr]
                        #spw_old = spw_ipsps[]
                        if len(spw_used) > 0:
                            start_used_all = ms2pts(spw_used['spw_start'], fs).astype('i4')
                            start_used = ms2pts(ipsps_all['ipsp_start'], fs).astype('i4')
                            start_old = ms2pts(ipsp_old['ipsp_start'], fs).astype('i4')
                            #import pdb; pdb.set_trace()
                            try:
                                # plot all the old IPSPs
                                
                                all_old = (start_old - win_pts[0] - starts_pts).astype('i4')
                                plt.plot(t[all_old], data_used[electr, all_old] + electr * add_it, 'mx', markersize = 8)
                                
                                # plot all the new IPSPs
                                all_used = (start_used_all - win_pts[0] - starts_pts).astype('i4')
                                plt.plot(t[all_used], data_used[electr, all_used] + electr * add_it, 'go', markersize = 8)
                                                                
                                # plot all the possible starts of SPWs selected
                                used = (start_used - win_pts[0] - starts_pts).astype('i4')
                                plt.plot(t[used], data_used[electr, used] + electr * add_it, 'ro')
                                
  
                            except:
                                import pdb; pdb.set_trace()
                    plt.show()           
         
    all_ipsps = np.concatenate(all_ipsps)
    
    # update the numbers of SPWs 
    old_no, proper_no = np.unique(all_ipsps['spw_no'], return_inverse=True)
    all_ipsps['spw_no'] = proper_no
   
    np.savez(save_folder + save_file, spw_ipsps = all_ipsps) 

def move_closest_group_to_spw_start(ipsps):
    """ moves first ipsps to spw start """
    
    # find group which is the closest to start
    for spw_no in np.unique(ipsps['spw_no']):
        ipsps_used = ipsps[ipsps['spw_no'] == spw_no]
        min_ipsp = min(ipsps_used['ipsp_start'])
        min_group = ipsps_used[ipsps_used['ipsp_start'] == min_ipsp]['group'][0]
        ipsps_used['ipsp_start'][ipsps_used['group'] == min_group] = ipsps_used[ipsps_used['group'] == min_group]['spw_start']
        ipsps[ipsps['spw_no'] == spw_no] = ipsps_used
    
    return ipsps

def remove_alone_ipsps(ipsps_trace, min_electr, shift_ipsp):
    """ calculates how many ipsps are there in any group and
    returns true if it should be removed, false if it shouldn't
    (saves always the first ipsp in spw"""
    n_electrodes_per_ipsp = count_coincident_ipsps(ipsps_trace, shift_ipsp)
    # find group which is the closest to start
    to_save = []
    #if ipsps_trace['spw_start'][0] == 149.3:
    #    import pdb; pdb.set_trace()
    
    for idx, spw  in enumerate(ipsps_trace):
        #import pdb; pdb.set_trace()
        if (spw['spw_start'] == spw['ipsp_start']) or n_electrodes_per_ipsp[idx] >= min_electr:
            to_save.append(True)
        else:
            to_save.append(False) 
    return np.array(to_save)    

def separate_groups(ipsps, max_length_ipsp_pts, init_spw_no = 0):
    """ checks how far are ipsps from each other, if further than max_length_ipsp_pts
    it separates them to two SPWs"""
    spw_no = init_spw_no
    new_ipsps = []
    for spw_no in np.unique(ipsps['spw_no']):
        spw_used = ipsps[ipsps['spw_no'] == spw_no]
        all_min_groups = np.zeros(len(np.unique(spw_used['group'])))
        
        #group_no = []
        # take the minimum IPSP in each group
        for idx, group in enumerate(np.unique(spw_used['group'])):
            all_min_groups[idx] = min(spw_used[spw_used['group'] == group]['ipsp_start'])
            #group_no.append(group) # remember which group it belongs to

        idx_sort = np.argsort(all_min_groups)
        differences = np.hstack([max_length_ipsp_pts + 1, np.diff(all_min_groups[idx_sort])])
        #import pdb; pdb.set_trace()
        differences = differences > max_length_ipsp_pts
        #if len(differences) > 1:
        
        separate_spws = np.cumsum(differences) - 1
        idx_reverse = np.argsort(idx_sort)
        separate_spws = separate_spws[idx_reverse]
        
        group_nos = separate_spws + init_spw_no
        init_spw_no = max(group_nos) + 1
        #import pdb; pdb.set_trace()
        
        for idx, group in enumerate(np.unique(spw_used['group'])):
            use_spw_no = group_nos[idx]
            # change spw no
            spw_used['spw_no'][spw_used['group'] == group] = (np.ones(len(spw_used[spw_used['group'] == group])) * use_spw_no).astype('i4')

            
        # change spw_start
        for spw_no in np.unique(spw_used['spw_no']):
            new_spw_start = min(spw_used[spw_used['spw_no'] == spw_no]['ipsp_start'])
            spw_used['spw_start'][spw_used['spw_no'] == spw_no] = np.ones(len(spw_used[spw_used['spw_no'] == spw_no])) * new_spw_start
        
        #print group_nos
        #if len(np.unique(group_nos)) > 1:
        #    import pdb; pdb.set_trace()
        new_ipsps.append(spw_used)
    new_ipsps = np.concatenate(new_ipsps)
    #import pdb; pdb.set_trace()
    return new_ipsps, init_spw_no

def separate_if_increase_too_low(ipsps_trace, ipsp_rise, expected_min_ipsp_rise, ipsp_amplitude, min_ipsp_ampl, init_spw_no):
    """ it checks each IPSP in the given spws, if it's not giving any increase to the next
    IPSP, the SPW is separated; it also checks if this new SPW does not become too short
    if it does it totally removes the previous IPSP"""
    all_new_ipsps = []
    
    
    init_spw_no = init_spw_no - 1
    init_spw_start = -1
    for spw_no in np.unique(ipsps_trace['spw_no']):
        
        init_spw_no = init_spw_no + 1 # new spw no for next spw
        spw_used = ipsps_trace[ipsps_trace['spw_no'] == spw_no]
        rise_used = ipsp_rise[ipsps_trace['spw_no'] == spw_no]
        ampl_used = ipsp_amplitude[ipsps_trace['spw_no'] == spw_no]
        #print 'new'
        all_groups = np.unique(spw_used['group'])
        if init_spw_start == spw_used['spw_start'][0] or spw_used['spw_start'][0] < init_spw_start:
            init_spw_no = init_spw_no -1
        else:
            init_spw_start = spw_used['spw_start'][0]
        
        #print init_spw_start
        all_group_times = np.zeros(len(all_groups))
        
        #import pdb; pdb.set_trace()
        # check the times for each group
        for idx, group in enumerate(all_groups):
            ipsp_time = min(spw_used[spw_used['group'] == group]['ipsp_start'])
            all_group_times[idx] = ipsp_time
        #import pdb; pdb.set_trace()    
        
        # order each group according to their starting times
        time_order = np.argsort(all_group_times)
        all_group_times = all_group_times[time_order]
        all_groups = all_groups[time_order]
        
        # go through every possible group
        for idx, group in enumerate(all_groups):
            #import pdb; pdb.set_trace() 

            # check if there is any IPSP after it
            if idx < len(all_groups)-1:
                # check which electrode are used by next IPSP
                next_electrodes = spw_used[spw_used['group'] == all_groups[idx + 1]]['electrode']
                
                range_next = np.arange(min(next_electrodes), max(next_electrodes)+1)
                current_electrodes = spw_used[spw_used['group'] == group]['electrode']
                range_current = np.arange(min(current_electrodes), max(current_electrodes)+1)
                #import pdb; pdb.set_trace() 
                #used_electr = np.intersect1d(range_next,range_current)
                
                # use only those electrodes
                group_used = (spw_used['group'] == group) #& np.in1d(spw_used['electrode'], used_electr)
                rise_used_group = rise_used[group_used]
                ampl_used_group = ampl_used[group_used]
                    
                if len(rise_used_group) == 0 or (max(rise_used_group) < expected_min_ipsp_rise) and max(ampl_used_group) < min_ipsp_ampl: 
                    #import pdb; pdb.set_trace()
                    # the rise is not large enough - new spw must be formed!
                    if init_spw_start == min(spw_used[spw_used['group'] == all_groups[idx+1]]['ipsp_start']):
                        update_new_start = init_spw_start
                        update_spw_no = init_spw_no                        
                    else:
                        update_spw_no = init_spw_no + 1
                        update_new_start = min(spw_used[spw_used['group'] == all_groups[idx+1]]['ipsp_start']) 
                
                else:
                    update_new_start = init_spw_start
                    update_spw_no = init_spw_no
                    print 2

            else:
                update_spw_no = init_spw_no
                update_new_start = init_spw_start #spw_used[spw_used['group'] == group]['spw_start']
                print 3

            print init_spw_start
            print init_spw_no
            print
            #else:
            #    # if it is last one just update the spw no, spw_start etc
            #    pass
            #import pdb; pdb.set_trace() 
            new_ipsp = spw_used[spw_used['group'] == group]
            
            
            new_ipsp['spw_no'] = (np.ones(len(new_ipsp)) * init_spw_no).astype('i4') 
            new_ipsp['spw_start'] = (np.ones(len(new_ipsp)) * init_spw_start).astype('f8') 
            all_new_ipsps.append(new_ipsp)

            init_spw_no = update_spw_no
            init_spw_start = update_new_start
            

    all_new_ipsps = np.concatenate(all_new_ipsps)  
    if len(np.unique(all_new_ipsps['spw_no'])) != len(np.unique(all_new_ipsps['spw_start'])):
        import pdb; pdb.set_trace()

    return all_new_ipsps, init_spw_no + 1


def separate_if_increase_too_low_only_first(ipsps_trace, ipsp_rise, expected_min_ipsp_rise, init_spw_no):
    """ it checks each first IPSP in the given spws, if it's not giving any increase to the next
    IPSP, the SPW is separated"""
    new_ipsps = []
    for spw_no in np.unique(ipsps_trace['spw_no']):
        spw_used = ipsps_trace[ipsps_trace['spw_no'] == spw_no]
        
        min_group = spw_used[spw_used['ipsp_start'] == min(spw_used['ipsp_start'])]['group'][0]
        rise_used = ipsp_rise[ipsps_trace['spw_no'] == spw_no]
        rise_used_group = rise_used[spw_used['group'] == min_group]
        # check if there is more than one IPSP group in this SPW
        if (len(np.unique(spw_used['group'])) > 1) and max(rise_used_group) < expected_min_ipsp_rise:
            # check if at least one rise in the first group is > expected_min_ipsp_rise
            # rise is too low, need of separation
            # set new spw_no
            spw_used['spw_no'][spw_used['group'] == min_group] = (np.ones(len(spw_used[spw_used['group'] == min_group])) * init_spw_no).astype('i4')            
            
            spw_used['spw_no'][spw_used['group'] != min_group]= (np.ones(len(spw_used[spw_used['group'] != min_group])) * init_spw_no + 1).astype('i4')
            
            new_spw_start = min(spw_used[spw_used['spw_no'] == init_spw_no + 1]['ipsp_start'])
            spw_used['spw_start'][spw_used['spw_no'] ==  init_spw_no + 1] = np.ones(len(spw_used[spw_used['spw_no'] ==  init_spw_no + 1])) * new_spw_start

        else:
            # everything alright; rise is fine; only update the spw_no
            spw_used['spw_no'] = (np.ones(len(spw_used)) * init_spw_no).astype('i4')
            
        add_init = len(np.unique(spw_used['spw_no']))
        #if add_init > 1:
        #    import pdb; pdb.set_trace()
        init_spw_no = init_spw_no + add_init
        new_ipsps.append(spw_used)
    new_ipsps = np.concatenate(new_ipsps)    
    return new_ipsps, init_spw_no

def update_ipsps_groups(save_folder, ipsps_groups, load_spwsipsp, load_datafile, save_file):
    
    npzfile         = np.load(save_folder + load_datafile)
    data            = npzfile['data']
    fs              = npzfile['fs']
    npzfile.close()   
    
    npzfile         = np.load(save_folder + load_spwsipsp)
    spw_ipsps       = npzfile['spw_ipsps']
    npzfile.close()           

    shift_ipsp = 1.2 # ms
    new_ipsps = []
    all_traces= np.unique(spw_ipsps['trace'])
    distance_from_point = 10 #ms
    assert len(np.unique(spw_ipsps[['trace','spw_start']])) == len(np.unique(spw_ipsps['spw_no']))    

    for trace in all_traces:    
        ipsps_trace = spw_ipsps[spw_ipsps['trace']==trace]
        
        # divide the IPSPs to groups
        #ipsp_amplitudes = calculate_amplitude_of_IPSP(ipsps_trace, 
        #                                          data[:, trace, :], fs)
        
        ipsp_amplitudes = calculate_max_in_given_patch(data[:, trace, :], ipsps_trace[['electrode','ipsp_start']], distance_from_point, fs)
      
        group_ids = group_ipsps(ipsps_trace, shift_ipsp)
        ipsps_trace = add_rec_field(ipsps_trace, [ipsp_amplitudes, group_ids],
                                     ['amplitude', 'group'])    
        new_ipsps.append(ipsps_trace)
    new_ipsps = np.concatenate(new_ipsps)
    #import pdb; pdb.set_trace()
    assert len(np.unique(new_ipsps[['trace','spw_start']])) == len(np.unique(new_ipsps['spw_no']))  

    np.savez(save_folder + save_file, spw_ipsps = new_ipsps)     

def remove_not_chosen_groups(spws, shorten_spws):
    """ removes only the groups of ipsps which were not used in the second 
    given rec array, it leaves the ipsps. So IPSPs must follow given earlier
    rule at least in one electrode, not necessary in all"""
    # each group have to be separate in each spw
    #assert len(np.unique(shorten_spws['group'])) == len(np.unique(shorten_spws[['group', 'spw_no']]))
    #assert len(np.unique(spws['group'])) == len(np.unique(spws[['group', 'spw_no']]))
               
#    if len(np.unique(shorten_spws['group'])) != len(np.unique(shorten_spws[['group', 'spw_no']])):
#        import pdb; pdb.set_trace()
        
    all_groups_left = np.unique(shorten_spws[['group', 'spw_no']])
    spws_left = []
    for spw in all_groups_left:
        #import pdb; pdb.set_trace()
        group = spw[0]
        spw_no = spw[1]
        spws_left.append(spws[(spws['group'] == group) & (spws['spw_no'] == spw_no)])
    if len(spws_left) > 0:
        return np.concatenate(spws_left)
    #import pdb; pdb.set_trace()
    else:
        return spws_left 


def update_spws_ipsp_beg(load_datafile, filter_folder, load_spwsipsp, load_spwsspike, save_folder, save_fig, save_file,ext, win = [-20, 80], save_filter = 'ipsp_filt_'):
    """ checks all the ipsps and corrects them for each spw"""

    npzfile         = np.load(save_folder + load_datafile)
    data            = npzfile['data']
    fs              = npzfile['fs']
    npzfile.close()   
    
    npzfile         = np.load(save_folder + load_spwsipsp)
    spw_ipsps       = npzfile['spw_ipsps']
    npzfile.close()   
    
    npzfile         = np.load(save_folder + load_spwsspike)
    spw_spike      = npzfile['spike_idx']  
    npzfile.close()     
    
    plot_it = False
    distanse_from_point = 10 # ms - used for calculating the amplitude
    #import pdb; pdb.set_trace()
    shift_ipsp = 1.2 # maximal possible shift of IPSP
    min_electr = 2 # on how many electrodes IPSP should be detected for the first ipsp (beginning of SPW)
    expected_min_ipsp_ampl = 30 # microV, min allowed IPSP amplitude
    expected_total_min = 3
    min_length_ipsp = 2 # minimum length of iPSP
    max_length_ipsp = 13 # maximum length of IPSP
    expected_min_ipsp_rise = 15 # minimim rise of ipsp (compared to the height of next IPSP start
    all_traces= np.unique(spw_ipsps['trace'])
    
    n= 5. # factor by which to downsample the data
    new_fs = fs*1./n
    data_downsampled = data[:,:,::n]
 
    # treat all the IPSPS
    print 'improving beginning of IPSP in trace: '
    folder_name = save_folder + filter_folder
    freq_slow = 300 # frequency to use for filtering the data
    base_window = 5 # window for moving average
    base_window = ms2pts(base_window, fs)
    min_ipsp_ampl = 5 # micro V hight of each first IPSP - used to separate SPWs
    all_ipsps = []
    init_spw_no = 0
    if np.size(data,2) < 100000:
        N_samples = 256
    else:
        N_samples = 1024

    #import pdb; pdb.set_trace()
    assert len(np.unique(spw_ipsps[['trace','spw_start']])) == len(np.unique(spw_ipsps['spw_no']))    

    for trace in all_traces:
        
        spw_ipsps_trace = spw_ipsps[spw_ipsps['trace']==trace]       
        print str(trace) + ' / ' + str(max(all_traces))
        # use data only for this trace
        data_trace = data[:,trace,:]
        data_downsampled_trace =  data_downsampled[:,trace,:]
        
        # use downsampled data only for this trace
        data_trace_filt = np.zeros([np.size(data_trace, 0),np.size(data_trace, 1)] )
        #data_trace_down_filt = np.zeros([np.size(data_downsampled_trace, 0),np.size(data_downsampled_trace, 1)] )
        
        # filter the data in all electrodes
        for electr in range(len(data_trace)):
            filename_slow = save_filter +'_' + str(freq_slow) + '_'+ str(electr) + "_" + str(trace)
            data_trace_filt[electr,:], temp = load_create(folder_name, filename_slow, [-1, freq_slow], fs, data_trace[electr,:], N = N_samples)
            data_trace_filt[electr,:], temp = filt.remove_baseloc(data_trace_filt[electr,:], base_window)  
            filename_down_slow = save_filter +'_' + str(freq_slow) + '_'+ str(electr) + "_" + str(trace) + "_down"        
            #data_trace_down_filt[electr,:], temp = load_create(folder_name, filename_down_slow, [-1, freq_slow], new_fs, data_downsampled_trace[electr,:], N = N_samples)      
            #data_trace_down_filt[electr,:], temp = filt.remove_baseloc(data_trace_down_filt[electr,:], base_window/n)  
        
          
        # check which IPSPs are of too low amplitude (on filtered data so that amplitude is not checked on spikes
        max_ampls = calculate_max_in_given_patch(data_trace, spw_ipsps_trace[['electrode','ipsp_start']], distanse_from_point, fs)  
        #max_ampls[max_ampls >= expected_min_ipsp_ampl]
        
        #if len(np.unique(spw_ipsps_trace['group'])) != len(np.unique(spw_ipsps_trace[['group', 'spw_no']])):
        #    import pdb; pdb.set_trace()
        ipsps_trace_temp = spw_ipsps_trace[max_ampls >= expected_min_ipsp_ampl]
        

        ipsps_trace = remove_not_chosen_groups(spw_ipsps_trace, ipsps_trace_temp)
        #import pdb; pdb.set_trace()
        
        if len(ipsps_trace) > 0:
            max_ampls = calculate_max_in_given_patch(data_trace_filt, ipsps_trace[['electrode','ipsp_start']], distanse_from_point, fs)
            ipsps_trace = ipsps_trace[max_ampls >= expected_total_min]
        
        
        if len(ipsps_trace) > 0:
            # remove IPSPs which are too close from each other
            #import pdb; pdb.set_trace()
            distance_between = calc_distance_between(ipsps_trace[['electrode','ipsp_start']], min_length_ipsp) #, max_ampls)      
            ipsps_trace = ipsps_trace[distance_between]
            
#            

            #import pdb; pdb.set_trace()
#            
#            # fill all the IPSP groups 
#            # import pdb; pdb.set_trace()
#            
#            ipsps_trace = fill_gap_in_all_groups(ipsps_trace)
#
#            
#            # remove the group which has not enough IPSPs
#            # (save the group which is the first from the correction even if it has too few IPSPs)
#
#            to_save = remove_alone_ipsps(ipsps_trace, min_electr, shift_ipsp)
#            ipsps_trace = ipsps_trace[to_save]
#
#            
        if len(ipsps_trace) > 0:
            # separate all the groups of ipsps which are further apart than max_length_ipsp_pts 
            ipsps_trace, init_spw_no = separate_groups(ipsps_trace, max_length_ipsp, init_spw_no)

        # if there are no groups don't do anything
        if len(ipsps_trace) > 0:
            # assigns time of closest group as spw beginning
            
            ipsps_trace = move_closest_group_to_spw_start(ipsps_trace)
            
#            # calculate the rise between the beginning of IPSP and next IPSP
#            
#            #import pdb; pdb.set_trace()
#            ipsp_rise = calculate_ipsp_rise(ipsps_trace, data_trace_down_filt, new_fs)
#            ipsp_amplitude = calculate_ipsp_amplitude(ipsps_trace, data_trace_down_filt, new_fs)
#
#            ipsps_trace, init_spw_no = separate_if_increase_too_low(ipsps_trace, ipsp_rise, expected_min_ipsp_rise, ipsp_amplitude, min_ipsp_ampl, init_spw_no)
#            #if len(np.unique(ipsps_trace[['trace','spw_start']])) == len(np.unique(ipsps_trace['spw_no'])):
#            #    import pdb; pdb.set_trace()
#        
        all_ipsps.append(ipsps_trace)
        temp = np.unique(np.concatenate(all_ipsps))
        if len(np.unique(temp[['trace','spw_start']])) != len(np.unique(temp['spw_no'])):
            import pdb; pdb.set_trace()
    all_ipsps = np.concatenate(all_ipsps)
    
    assert len(np.unique(all_ipsps[['trace','spw_start']])) == len(np.unique(all_ipsps['spw_no']))   
    np.savez(save_folder + save_file, spw_ipsps = all_ipsps)
    
    # -----------------------------------------------
    if plot_it:
        add_it = 350
        
        before = ms2pts(win[0], fs).astype('i4')
        after = ms2pts(win[1], fs).astype('i4')
        before_down = ms2pts(win[0], new_fs).astype('i4')
        after_down = ms2pts(win[1], new_fs).astype('i4')


        #     go through all the spws
        for spw_no in [0]: #range(10, 17): #np.unique(all_ipsps['spw_no']):
            print spw_no
            #import pdb; pdb.set_trace()
            fig = plt.figure()   
            #spw_min_start = 9000000000
            #spw_max_end = -1
            
            spw_used = all_ipsps[all_ipsps['spw_no'] == spw_no]
            #ipsps_old = spw_ipsps[spw_ipsps['trace']==trace]
              
            
            #import pdb; pdb.set_trace()
            trace = spw_used['trace'][0]
            #spikes_used = all_spikes[all_spikes['trace'] == trace]
            
            ipsps_old = spw_ipsps[spw_ipsps['trace'] == trace]
            
            
            spw_start = spw_used['spw_start'][0]
            #import pdb; pdb.set_trace()
            spw_start_pts = ms2pts(spw_start, fs).astype('i4')
            
            plot_start = max(0,spw_start_pts + before)
            plot_end = min(spw_start_pts + after,  np.size(data,2))
            spw_start_pt = spw_start_pts - plot_start
            
            spw_start_pts_down = ms2pts(spw_start, new_fs).astype('i4')
            
            plot_start_down = max(0,spw_start_pts_down + before_down)
            plot_end_down = min(spw_start_pts_down + after_down,  np.size(data,2)/n)
            spw_start_pt_down = spw_start_pts_down - plot_start_down
                     
            for electr in range(np.size(data,0)):
                plot_add = add_it * electr
                
                ipsps_used = spw_used[spw_used['electrode'] == electr]
                ipsps_used_pts = ms2pts(ipsps_used['ipsp_start'], fs).astype('i4')
                #ipsps_used_pts_down = ms2pts(ipsps_used['ipsp_start'], new_fs).astype('i4')
                
#                sp_used = spikes_used[spikes_used['electrode'] == electr]
#                sp_used_pts = ms2pts(sp_used['electrode'], fs).astype('i4')
                
               
                data_used = data[electr,trace, plot_start: plot_end]
                #import pdb; pdb.set_trace()
                filename_down_slow = save_filter +'_' + str(freq_slow) + '_'+ str(electr) + "_" + str(trace) + "_down"  
                data_filt, temp = load_create(folder_name, filename_down_slow, [-1, freq_slow], new_fs, data_downsampled[electr,trace, :], N = N_samples)
                #data_filt, temp = filt.remove_baseloc(data_filt, base_window/n)  

                #filename_fast = save_filter +str(freq_fast) + '_'+ str(electr) + "_" + str(trace)
                #data_filt, fs = load_create(folder_name, filename_fast, freq_fast, fs, data_filt)

                data_filt = data_filt[plot_start_down: plot_end_down]
                data_down = data_downsampled[electr,trace, plot_start_down: plot_end_down]
                #data_filt = data_trace_down_filt[electr,trace, plot_start_down: plot_end_down]
                #import pdb; pdb.set_trace()
                t = dat.get_timeline(data_used, fs, 'ms')
                t_down = dat.get_timeline(data_filt, new_fs, 'ms')
                
                #import pdb; pdb.set_trace()
                plt.plot(t, data_used + plot_add, 'b')
                plt.plot(t_down, np.correlate(data_filt,data_down) + plot_add, 'r')
                plt.plot(t_down, data_down + plot_add, 'g')
                ipsps_to_plot = ipsps_used_pts - plot_start
                
                # plot new ipsps
                try:
                    plt.plot(t[ipsps_to_plot], data_used[ipsps_to_plot] + plot_add, 'ko', ms=6)
                except:
                    import pdb; pdb.set_trace()
                # plot old ipsps
                #ipsps_old = spw_ipsps[spw_ipsps['trace']==trace]
                ipsps_old_used = ipsps_old[(ipsps_old['electrode'] == electr) & (ipsps_old['ipsp_start'] >= pts2ms(plot_start, fs)) & (ipsps_old['ipsp_start'] <= pts2ms(plot_end, fs))]
                ipsps_old_pts = ms2pts(ipsps_old_used['ipsp_start'], fs).astype('i4')
                
                ipsps_to_old = ipsps_old_pts - plot_start
                #import pdb; pdb.set_trace()
                plt.plot(t[ipsps_to_old], data_used[ipsps_to_old] + plot_add, 'ko', mfc='none', ms=8)
                
                # plot spikes
#                spikes_to_plot = sp_used_pts - plot_start
#                spikes_to_plot = spikes_to_plot[spikes_to_plot < len(t)]
                #plt.plot(t[spikes_to_plot], data_used[spikes_to_plot] + plot_add, 'kx', ms=7)
                plt.axvline(t[spw_start_pt])

            plt.show()
            #import pdb; pdb.set_trace()
            tit = 'spw: ' + str(spw_no)
            plt.title(tit)
            fig.savefig(save_folder + save_fig + str(spw_no) + ext,dpi=600)
            fig.savefig(save_folder + save_fig + str(spw_no) + '.eps',dpi=600)        

            plt.close()


    
    del spw_ipsps, data, all_ipsps



def update_spws_beg_backup(load_datafile, load_spwsipsp, load_spwsspike, save_folder, save_fig, save_file,ext, win = [-20, 80]):
    """ checks all the ipsps and corrects them for each spw"""

    npzfile         = np.load(save_folder + load_datafile)
    data            = npzfile['data']
    fs              = npzfile['fs']
    npzfile.close()   
    
    npzfile         = np.load(save_folder + load_spwsipsp)
    spw_ipsps       = npzfile['spw_ipsps']
    npzfile.close()   
    
    npzfile         = np.load(save_folder + load_spwsspike)
    #import pdb; pdb.set_trace()
    spw_spike      = npzfile['spike_idx']  
    npzfile.close()     
    
    plot_it = False
    distanse_from_point = 5 # ms
    #import pdb; pdb.set_trace()
    shift_ipsp = 1 # ms
    min_electr = 2 # on how many electrodes IPSP should be detected for the first ipsp (beginning of SPW)
    expected_min_ipsp_ampl = 30 # microV
    shift_spike= 1 #ms
    min_length_ipsp = 3
    
     
    spw_ipsps_list = []
    all_traces= np.unique(spw_ipsps['trace'])

# treat all the IPSPS



# use only for the beginning of SPW (stricter rules

# use for all other IPSPs
    
    
    
    for trace in all_traces:
        spw_ipsps_trace = spw_ipsps[spw_ipsps['trace']==trace]

        # calculate the amplitude for each IPSP (from start to distanse_from_point after)
        data_trace = data[:,trace,:]
        
        distance_between = calc_distance_between(spw_ipsps_trace[['electrode','ipsp_start']])
        spw_ipsps_trace = spw_ipsps_trace[distance_between >= min_length_ipsp]

        # check which IPSPs are of too low amplitude
        max_ampls = calculate_max_in_given_patch(data_trace, spw_ipsps_trace[['electrode','ipsp_start']], distanse_from_point, fs)
        spw_ipsps_trace = spw_ipsps_trace[max_ampls >= expected_min_ipsp_ampl]

        # calculate the rise between the beginning of IPSP and next IPSP
        #ipsp_rise = calculate_ipsp_rise(spw_ipsps_trace, data[:, trace, :], fs)
        
        #spw_ipsps_trace_rised = spw_ipsps_trace[(np.abs(ipsp_rise) > expected_min_ipsp_ampl)]        
        # check in how many electrodes each IPSP appears
        n_electrodes_per_ipsp = count_coincident_ipsps(spw_ipsps_trace, shift_ipsp)
        #
        spw_ipsps_trace = spw_ipsps_trace[(n_electrodes_per_ipsp >= min_electr)]
        #spw_ipsps_trace = spw_ipsps_trace_rised
        #group_ids = group_ipsps(spw_ipsps_trace, shift_ipsp)
        if len(spw_ipsps_trace) > 0:
            cipsps = count_coincident_ipsps(spw_ipsps_trace, shift_ipsp)
            while np.min(cipsps)<min_electr:
                print np.sum(cipsps<min_electr)
                spw_ipsps_trace = spw_ipsps_trace[cipsps>=min_electr]
                cipsps = count_coincident_ipsps(spw_ipsps_trace, shift_ipsp)
        spw_ipsps_list.append(spw_ipsps_trace)
    spw_selected = np.concatenate(spw_ipsps_list)
    
    
    all_traces_update = np.unique(spw_selected['trace'])
    all_ipsps = []
    all_spikes = []
    for trace in all_traces_update:
        #import pdb; pdb.set_trace()
        ipsps_trace = spw_selected[spw_selected['trace']==trace]
        spikes_trace = spw_spike[spw_spike['trace']==trace]
        all_spikes.append(spikes_trace)

        ipsp_amplitudes = calculate_amplitude_of_IPSP(ipsps_trace, 
                                                  data[:, trace, :], fs)
        group_ids = group_ipsps(ipsps_trace, shift_ipsp)
        ipsps_trace = add_rec_field(ipsps_trace, [ipsp_amplitudes, group_ids],
                                     ['amplitude', 'group'])
        

        if len(ipsps_trace) > 0:
            #import pdb; pdb.set_trace()
            ipsps_trace = shift_ipsp_start(ipsps_trace, spikes_trace, 
                                       shift_spike)
        
            #shift spw to first ipsp
            ipsps_trace = shift_spw_to_first_ipsp(ipsps_trace)
             
        all_ipsps.append(ipsps_trace)
    all_ipsps = np.concatenate(all_ipsps)
    all_spikes = np.concatenate(all_spikes)
        

        # plot 
    if plot_it:
        add_it = 100
        before = ms2pts(win[0], fs).astype('i4')
        after = ms2pts(win[1], fs).astype('i4')
        
        #     go through all the spws
        for spw_no in np.unique(all_ipsps['spw_no']):
        #import pdb; pdb.set_trace()    
            fig = plt.figure()   
            spw_min_start = 9000000000
            spw_max_end = -1
            
            spw_used = all_ipsps[all_ipsps['spw_no'] == spw_no]
            spikes_used = all_spikes[all_spikes['spw_no'] == spw_no]
            trace = spw_used['trace'][0]
            spw_start = spw_used['spw_start'][0]
            spw_start_pts = ms2pts(spw_start, fs).astype('i4')
            
            plot_start = max(0,spw_start_pts + before)
            plot_end = min(spw_start_pts + after,  np.size(data,2))
            
            spw_start_pt = spw_start_pts - plot_start

            
            for electr in range(np.size(data,0)):
                plot_add = add_it * electr
                
                ipsps_used = spw_used[spw_used['electrode'] == electr]
                ipsps_used_pts = ms2pts(ipsps_used['ipsp_start'], fs).astype('i4')
                
                sp_used = spikes_used[spikes_used['electrode'] == electr]
                sp_used_pts = ms2pts(sp_used['time'], fs).astype('i4')
                
                data_used = data[electr,trace, plot_start: plot_end]
                t = dat.get_timeline(data_used, fs, 'ms')
                
                
                plt.plot(t, data_used + plot_add)
                ipsps_to_plot = ipsps_used_pts - plot_start
                #try:
                plt.plot(t[ipsps_to_plot], data_used[ipsps_to_plot] + plot_add, 'ko', mfc='none', ms=7)
                spikes_to_plot = sp_used_pts - plot_start
                spikes_to_plot = spikes_to_plot[spikes_to_plot < len(t)]
                plt.plot(t[spikes_to_plot], data_used[spikes_to_plot] + plot_add, 'kx', ms=7)
                plt.axvline(t[spw_start_pt])

            #plt.show()
            tit = 'spw: ' + str(spw_no)
            plt.title(tit)
            fig.savefig(save_folder + save_fig + str(spw_no) + ext,dpi=600)
            fig.savefig(save_folder + save_fig + str(spw_no) + '.eps',dpi=600)        

            plt.close()

    np.savez(save_folder + save_file, spw_ipsps = all_ipsps, spikes = all_spikes) 
    
    del spw_ipsps, data, all_ipsps


def update_SPW_ipsp(load_datafile, filter_folder, load_waves, load_spikes, save_folder, save_file, spw_length = 80, save_filter = 'ipsp_filt_'):
    # it looks for the ipsps within detected spws - separate for each electrode
    # it take very long to analyse so be patient!
    
    npzfile         = np.load(save_folder + load_datafile)
    data            = npzfile['data']
    fs              = npzfile['fs']
    npzfile.close()
    
    # load starts of spws
    npzfile         = np.load(save_folder + load_spikes)
    spike_idxs      = npzfile['spike_idx']
    npzfile.close()

    npzfile         = np.load(save_folder + load_waves)
    spw_details      = npzfile['spw_details']
    npzfile.close()
    #import pdb; pdb.set_trace()

    # order spws according first to trace, and then to their start
    i = np.argsort(spw_details, order=['trace', 'time'])
    spw_details = spw_details[i]
    
    folder_name = save_folder + filter_folder
    
    print
    print "analyzing SPW no:",
    plot_it = False
    add_it = 100
    window = 1.0 #0.5 # ms for calculating moving average
    
    freq_fast = 300
    freq_slow = 100
    window = ms2pts(window, fs)
    spw_ipsps = []
    
    data_length_ms = pts2ms(np.size(data,2), fs)
    
    # take each SPW separately and find ipsps
    spw_len = len(np.unique(spw_details['spw_no']))
        
    all_spw_numbers = np.unique(spw_details['spw_no'])
    
    for idx_spw, spw in enumerate(all_spw_numbers):
        
        spw_used = spw_details[spw_details['spw_no'] == spw]
        trace = int(spw_used['trace'][0])
        spw_no = int(spw_used['spw_no'][0])
        
        #if trace == 2:
        #    import pdb; pdb.set_trace()

        # always use minimum detected start for this SPW
        min_start = min(spw_used['time'])
        min_start_pts = ms2pts(min_start, fs).astype(int)
        
        if idx_spw < len(all_spw_numbers) - 1:
            # there exists spw after this one
            spw_next = spw_details[spw_details['spw_no'] == all_spw_numbers[idx_spw + 1]]
            trace_next = spw_next['trace'][0]

            if trace_next == trace:
                # the two SPWs are on the same trace
                time_next = min(spw_next['time'])
                distance = time_next - min_start
                assert distance > 0
                if distance < spw_length:
                    add_length = distance
                else:
                    add_length = spw_length
            else:
                add_length = spw_length
        else:
            add_length = spw_length
        #import pdb; pdb.set_trace()
        spw_max_len = min_start + add_length

        #----
        
#        temp_trace = 0
#        all_time = []
#        for sp in np.unique(spw_details['spw_no']):
#            sp_us = spw_details[spw_details['spw_no'] == sp]
#            tim = min(sp_us['time'])
#            trace = sp_us['trace'][0]
#            
#            if trace != temp_trace:
#                sol =  np.diff(all_time)
#                print trace
#                print sol
#                print len(sol < 0)
#                print 
#                all_time = []
#            else:
#                all_time.append(tim)    
#            temp_trace = trace
#            #print trace
#                
                
        
#        # find the end of the IPSP
#        if idx_spw < len(all_spw_numbers)-1 and spw_details[spw_details['spw_no'] == all_spw_numbers[idx_spw + 1]]['trace'][-1] == trace:
#            # check if it's the same trace
#            next_idx = all_spw_numbers[idx_spw + 1]
#            next_spw_time = spw_details[spw_details['spw_no'] == next_idx]['time'][-1]
#            spw_next = min(next_spw_time, spw_max_len)
#        else:
#            spw_next = min(spw_max_len, data_length_ms)

        spw_end_pts = ms2pts(spw_max_len, fs)

        print str(spw) + '/' + str(spw_len)

        for electr in range(np.size(data,0)):
                
            data_used = data[electr, trace, min_start_pts:spw_end_pts]
            data_trace = data[:,trace,:]        
            
            """---- """  
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
                        #import pdb; pdb.set_trace() 
                        mins.append(detected + (idx-1)*window) # <- normally maxs
                    else:
                        pass
                        #mins.append(detected + (idx-1)*window)
                    switch = switch * -1 # looking for start of fall
                
                mean_temp = m
            """ --- """
            
            """ +++ """
            # Richard's idea on checking the filtered data's second derivative 0 point
            #import pdb; pdb.set_trace()

            filename_fast = save_filter +str(freq_fast) + '_'+ str(electr) + "_" + str(trace)
            filename_slow = save_filter +'_' + str(freq_slow) + '_'+ str(electr) + "_" + str(trace)

            data_trace_filt, fs = load_create(folder_name, filename_fast, freq_fast, fs, data_trace[electr,:])
            data_trace_filt, fs = load_create(folder_name, filename_slow, [-1, freq_slow], fs, data_trace[electr,:])
            data_used_filt = data_trace_filt[min_start_pts: spw_end_pts]
            
            # calculate derivative
            first_deriv = np.diff(data_used_filt)
            second_deriv = np.diff(first_deriv)
            second_deriv[second_deriv > 0] = 1
            second_deriv[second_deriv < 0] = 0
            second_deriv2 = np.hstack((second_deriv[1:], second_deriv[-1]))
            peaks = second_deriv - second_deriv2

            peaks2nd = np.hstack((peaks, [0, 0]))
            peaks_init, = np.where(peaks2nd == -1)

            peaks_idx = peaks_init + min_start_pts

            for star in range(len(peaks_idx) -1):
                ipsp_start = np.argmin(data_trace[electr, peaks_idx[star]:peaks_idx[star+1]]) + peaks_init[star]
                ipsp_max = np.argmax(data_trace[electr, peaks_idx[star]:peaks_idx[star+1]]) + peaks_init[star]
                maxs.append(ipsp_start)
                
            """ +++ """                            
            #import pdb; pdb.set_trace()
            
            # combine the two methods together
            maxs = np.unique(maxs+mins).astype('i4') # join both functions together
            if plot_it and trace == 2:
                t = dat.get_timeline(data_used, fs, 'ms') + min_start #spw_used['spw_start'][0]
                plt.plot(t, data_used + add_it*electr)   
                plt.plot(t[maxs], data_used[maxs] + add_it*electr, 'r>', alpha = 0.6)

            mini = pts2ms(maxs,fs) + min_start #spw_electr_used['spw_start'][0]

            ipsp_start = mini.astype('f8')
            electrodes = np.ones(len(ipsp_start), dtype='i4')*electr
            traces = np.ones(len(ipsp_start), dtype='i4')*trace
            spw_num = np.ones(len(ipsp_start), dtype='i4')*spw_no
            ipsp_no = np.arange(0,len(mini)).astype('i4')
            spw_start = np.ones(len(ipsp_start), dtype = 'f8') * min_start
            #spw_end = np.ones(len(ipsp_start), dtype = 'f8') * spw_next
            #import pdb; pdb.set_trace()  
            spw_ipsps.append(np.rec.fromarrays([electrodes, traces, spw_num, spw_start,
                                                    ipsp_no, ipsp_start], 
                                                names='electrode, trace, spw_no, spw_start, ipsp_no, ipsp_start'))
        # plot it if necessary
        if plot_it and trace == 2:
            plt.show()

    spw_ipsps = np.concatenate(spw_ipsps)
    np.savez(save_folder + save_file, spw_ipsps = spw_ipsps) 
    del data, spw_details 
      

def check_SPW_length(spw_starts, spw_ends, min_length):
    # checks which SPWs are of correct lenghts and returns only those correct ones
    #import pdb; pdb.set_trace()
    long_enough = (spw_ends['time'] - spw_starts['time']) >= min_length
    return spw_starts[long_enough], spw_ends[long_enough]
    
    
def group_spws(spw_start, max_move):
    """Assign spws to groups"""
    #import pdb; pdb.set_trace() 
    order = np.argsort(spw_start['time'])
    
    spw_st_sorted = spw_start[order]
    
    group_ids = np.empty(len(spw_st_sorted), dtype=np.int)
    group_ids.fill(-1)
    
    #electrodes = spw_st_sorted['electrode']
    spw_start_times = spw_st_sorted['time']
    
    distances = np.diff(spw_start_times)
    same_wave = distances > max_move
    group_ids =np.concatenate([[True], same_wave])
    group_ids = np.cumsum(group_ids)
    
    inverse_i = np.argsort(order)
    return group_ids[inverse_i]

def unique2d(a):
    order = np.lexsort(a.T)
    a = a[order]
    diff = np.diff(a, axis=0)
    ui = np.ones(len(a), 'bool')
    ui[1:] = (diff != 0).any(axis=1)
    inverse_i = np.argsort(order)
    a = a[inverse_i]
    ui = ui[inverse_i]
    return a[ui], ui

def update_highWaves_numb(load_spwsfile, save_folder, data_file):
    # reject all the SPWs where there is less than min_no_wave spws detected, 
    # check which spws belog together and save them
    # it also numbers SPWs
    
    # load starts of spws
    npzfile         = np.load(save_folder + load_spwsfile)
    spw_starts      = npzfile['starts']
    spw_ends        = npzfile['ends']
    npzfile.close()
    
    #import pdb; pdb.set_trace() 
    min_no_electr = 1
    max_move = 15 # ms
    min_spw_length = 5 #ms
    #min_dist_between_starts = 20
    all_spw_traces = []
    #spw_details_temp = []
    spw_no = 0
    a = 0
    no_traces = len(np.unique(spw_starts['trace']))
    
    # eliminate spws which are too short
    spw_starts, spw_ends = check_SPW_length(spw_starts, spw_ends, min_spw_length)
    last_group = 0
    for trace in np.unique(spw_starts['trace']):
        #import pdb; pdb.set_trace() 
        print 'trace: ' + str(trace + 1) + '/' + str(no_traces)
        # get spws for each electrode and check if it's the same
        spw_st_trace = spw_starts[(spw_starts['trace'] == trace)]               
        
        sort_idx = np.argsort(spw_st_trace['time'])
        spw_st_sorted = spw_st_trace[:,:,sort_idx]
        
        # remove those starts of spws which are too close from each other in one electrode
        for electr in np.unique(spw_starts['electrode']):
            spw_st_electr = spw_st_sorted[(spw_st_sorted['electrode'] == electr)]
            
            distance = np.diff(spw_st_electr['time'])
            too_close = distance < max_move
            if np.sum(too_close) > 0:
                print '2 spws are trying to be too close, go to: update_highWaves_numb'
                import pdb; pdb.set_trace()
        
        group_ids = group_spws(spw_st_sorted, max_move)
        #import pdb; pdb.set_trace()
        bins = np.bincount(group_ids)
        number_in_bin = bins[group_ids]
        found_in_enough_electr = number_in_bin >= min_no_electr
        
        # check if electrodes in each of the groups are unique
        electrodes_used = spw_st_sorted[found_in_enough_electr]
        group_ids = group_ids[found_in_enough_electr]
        no_unique = np.zeros([2, len(group_ids)])
        no_unique[0,:] = group_ids
        no_unique[1, :] = electrodes_used['electrode']
        #import pdb; pdb.set_trace()
        unique_values, unique_index = unique2d(no_unique.T)
        electrodes_used = electrodes_used[unique_index]
        
        if len(unique_values) != len(electrodes_used):
            print 'SPW on the same electrode was given the same group, go to: update_highWaves_numb'
            import pdb; pdb.set_trace()
        
        selected_spw_trace = electrodes_used

        #import pdb; pdb.set_trace() 
        spw_no = group_ids[unique_index]
        move_to_ordered = np.unique(spw_no)
        move_to_ordered = move_to_ordered.searchsorted(spw_no)
        
        groups_ordered = move_to_ordered + last_group
        assert len(groups_ordered) != 0

        last_group = max(groups_ordered)+1
        
        # add spw no to already existing 'electrode', 'trace', and 'time' 
        # of a SPW
        new_spw_trace = add_rec_field(selected_spw_trace, groups_ordered, ['spw_no'])
        
        
        all_spw_traces.append(new_spw_trace)
        #if trace == 2:

            
    spw_details = np.concatenate(all_spw_traces)   
    
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
        data[electr, :, :] = data[electr, :, :] - mean_datatime



    print
    print "setting the data to local baseline, working on electrode:",
    local_length = 5000
    # remove by more local baseline
    for electr in range(np.size(data, 0)):
        print electr,
        for trace in range(np.size(data, 1)):
            new_dat, moved_avg = filt.remove_baseloc(data[electr, trace, :], local_length)     
            #import pdb; pdb.set_trace() 
            #electro_data = data[electr, :, :] #data[data['electrode'] == electr]
            #mean_datatime = np.mean(electro_data)
            data[electr, trace, :] = data[electr, trace, :] - moved_avg

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
