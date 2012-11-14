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
    maxs, indcs = dw.max_waves(data, starts, ends)
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

def update_dist_fromSpike(sp_idx_intr, sp_idx_extr, fs, save_folder, data_file = 'dist_fromSpike', max_dist = 3, intra = 0, allowms = 0):
    """ updates all the distances from extracellular spikes to the preceeding it intracellular spike,
    it also returns the indexes of those spikes which are closer than min_dist (in ms) to intracellular spike"""
    distances = []
    min_distances = []
    # allow is variable which alows to count allow ms before the spike
    allow = ms2pts(allowms, fs)
    
    
    for electr in range(np.size(sp_idx_extr,0)):
        dist_all = []
        min_dist_all = []
        for trace in range(np.size(sp_idx_extr[electr],0)):
          

            #print np.size(sp_idx_intr[0][0])
            if intra == 1:
            #if np.size(sp_idx_intr[0][0]) > 10:
                
                dist = dist_fromSpike(sp_idx_intr[trace], sp_idx_extr[electr][trace], fs)
               
#                plt.plot(sp_idx_intr[trace], 'ro')
#                plt.plot(sp_idx_extr[electr][trace], 'go')
#                plt.show()                
            else:
                #print sp_idx_intr[trace]
                #a = b
#                if allow > 0:
#                    sp_chaged = []
#                    for i in range(len(sp_idx_intr[electr][trace])):
#
#                        sp_chaged.append(sp_idx_intr[electr][trace][i] - allow)
#                        #[sp_idx_intr[trace][i] - allow for )]
#                sp_idx_intr[electr][trace] = sp_chaged
                #print sp_idx_intr[trace] 

                
                dist = dist_fromSpike(sp_idx_intr[electr][trace] - allow, sp_idx_extr[electr][trace], fs)
            if allowms > 0:
                dist = [dist[i]-allowms for i in range(len(dist))]
            #print sp_idx_extr[electr][trace]
            #print dist

            inter = [sp_idx_extr[electr][trace][i] for i in range(len(dist)) if dist[i] <=max_dist]
            dist_all.append(dist)
            min_dist_all.append(inter)
        distances.append(dist_all)
        min_distances.append(min_dist_all)

    np.savez(save_folder + data_file, distances, min_distances, fs, max_dist) 
    return distances, min_distances, fs, max_dist

    
#
#def update_inter_extr_spike():
#    """ not tested,
#    it finds and saves the distances of all of the extracellular spikes to closest, proceeding 
#    intracellular spike"""
#    pass

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




def update_intraSpikes(data, fs, save_folder, save_file = "intra_spikes", pulse_len = 500, ):
    """ pulse_len in ms - length of the stimulation pulse in intracellular electrode"""
#    print 'data'
#    print np.shape(data)
    #print np.shape(data)
    sp_idx_first_all = []
    sp_idx_all = []
    for trace in range(len(data)):
        spiking_thres = -10 #mV
        #pulse_len 
        # detect only the first spike in the row 
        sp_idx_first = detect_1spike(data[trace][:], spiking_thres, fs, pulse_len)
        sp_idx_a = detect_spikes(data[trace][:], spiking_thres) # for detecting all the spikes 
        #print np.shape(sp_idx_first)
        #print np.shape(sp_idx_a)
        sp_idx_first_all.append(sp_idx_first)
        sp_idx_all.append(sp_idx_a)
#    print 'all'
#    print sp_idx_first_all
#    print np.size(sp_idx_first_all,0)
#    print sp_idx_all
#    print np.size(sp_idx_all,0)
    
       
    np.savez(save_folder + save_file, sp_idx_first_all, sp_idx_all, fs) 
    return sp_idx_first_all, sp_idx_all, fs

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
        np.savez(folder_save + filename_save, data = data_filt, fs = fs)
        
    return data_filt, fs


def update_extraspikes(data_load, save_folder, save_file = "ex_spikes"):
    """ finds and updates the detection of extracellular spikes"""
    
    npzfile = np.load(save_folder + data_load)
    data = npzfile['data']
    fs = npzfile['fs'] 
    npzfile.close()
       
    freq_fast = 500.
    # find extracellular spikes
    f_d = 'filtered/'
    folder_name = save_folder + f_d
    N = 1000
    print
    print "finding extracellular spikes, working on electrode:",
    idx_all = []
    ampl_all = []
    for electr in np.unique(data['electrode']):
        print electr,  

        if len(data[electr]) > 1:
            N = 100
        for trace in np.unique(data[(data['electrode'] == electr)]['trace']):
            data_used = data[(data['electrode'] == electr) & (data['trace'] == trace)]['time']
            filename_fast = 'fast_data' + str(electr) + "_" + str(trace)
            data_fast, fs = load_create(folder_name, filename_fast, freq_fast, fs, data_used, N)
            spike_ampl, spike_idxs = fes.find_extra_spikes(data_used, data_fast, fs) #(data[electr][trace], fs)
            del data_fast
            typ = 'f8'
            spike_idxs = pts2ms(spike_idxs, fs).astype(typ)
            electrodes = np.ones(len(spike_idxs), dtype=np.int32)*electr
            traces = np.ones(len(spike_idxs), dtype=np.int32)*trace    
            idx_all.append(np.rec.fromarrays([electrodes, traces, spike_idxs], names='electrode,trace,time'))
            ampl_all.append(np.rec.fromarrays([electrodes, traces, np.array(spike_ampl, dtype=typ)], names='electrode,trace,time'))
            del data_used
    
    #import pdb; pdb.set_trace()
    ampl_all = np.concatenate(ampl_all)         
    idx_all = np.concatenate(idx_all)            
    
    np.savez(save_folder + save_file, spike_idx = idx_all, spike_ampl = ampl_all, fs = fs)
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
    
    for electr in np.unique(data['electrode']):
        print electr,

        if len(data[electr]) > 1:
            N = 100

        for trace in np.unique(data[(data['electrode'] == electr)]['trace']):
            
            data_used = data[(data['electrode'] == electr) & (data['trace'] == trace)]['time']
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
   
def update_datafile(filename, ex_electr, save_folder, data_file = 'data', data_part = 'all'):
    """ reads given file and saves the data read into rec array (numpy file format)"""
    
    data_all = []
    
    # do the same for every electrode given - read data
    all_data, no_segments = dat.get_data(filename) 
    
    print
    print "reading the data, working on electrode:",
    
    for electr in ex_electr:
        print electr,
        data, fs = dat.get_electrdata(all_data, no_segments, electr, data_part)
        data_all.append(data)
    data_all = np.concatenate(data_all)
    np.savez(save_folder + data_file, data = data_all, fs = fs)
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

    
    thres_level = 20
    print
    print "removing averaged baseline:",
    spws_starts = []
    spws_ends = []
    for electr in np.unique(data['electrode']): 
        print electr,
        
        for trace in np.unique(data[data['electrode'] == electr]['trace']):

            window = ms2pts(atten_len, fs) # define the size of the window for removing the baseline
            data_used = data[(data['electrode'] == electr) & (data['trace'] == trace)]['time']
            
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
            spws_ends.append(np.rec.fromarrays([electrodes, traces, spw_ends], names='electrode,trace,time'))
                     
    spws_starts = np.concatenate(spws_starts)    
    spws_ends = np.concatenate(spws_ends)    
    
    np.savez(save_folder + data_file, starts = spws_starts, ends = spws_ends, fs = fs)   
    del data

def get_var_used(var, electr, trace):
    var_used = var[(var['electrode'] == electr) & (var['trace'] == trace)]['time']
    return var_used


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
                #chosen_spikes.append(np.rec.fromarrays([el, tr, sp, am, st, en, sp_no],
                #                                       names='electrode,trace, spikes, spike_ampl, spw_start, spw_end, spw_no'))                
                
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
    
    plot_it = False
    add_it = 100
    shift_ipsp = 1 # ms
    in_min_electrodes = 3 # ipsp has to be found in at least that many electrodes
    shift_spike= 0.15 # ms
    proper_ipsps = []
    expected_min_ipsp_ampl = 15 # microV
    # go through all the spws
    for spw in np.unique(spw_ipsps['spw_no']):
        # save detected ipsps and spikes for each spw
        print 'updating spw:' + str(spw)
        sp_ip_used = spw_ipsps[spw_ipsps['spw_no'] == spw]
        sp_sp_used = spw_spike[spw_spike['spw_no'] == spw]
        trace = sp_ip_used['trace'][0]
        data_temp = data[data['trace'] == trace]
        sort_ip_idx     = np.argsort(sp_ip_used['ipsp_start'])
        ip_sorted     = sp_ip_used[sort_ip_idx]
        spw_no = sp_ip_used['spw_no'][0]
        spike_times = sp_sp_used['spikes']
        ipsps_temp, spikes_temp, ipsp_ampl_temp, electr_temp, increase_temp = [], [], [], [], []
        ipsp_old = -1
        beginnings = []
        fig = plt.figure()
        for ipsp in ip_sorted:
            ipsp_time = ipsp['ipsp_start']
            ipsp_end = ipsp['ipsp_end']
            electr = ipsp['electrode']
                            
            if ipsp_time > ipsp_old + shift_ipsp and len(ipsps_temp) != 0:
                # analysis which IPSP to consider
                # check if no of found IPSPS isn't less than in_min_electrodes
                if len(ipsps_temp) >= in_min_electrodes:
                    
                    # find those ipsps where it causes increase more than expected_min_ipsp_ampl microV
                    idx_temp = [i for i in range(len(increase_temp)) if increase_temp[i] > expected_min_ipsp_ampl]
                    
                    electr_temp = np.array(electr_temp)[idx_temp]
                    ipsps_temp = np.array(ipsps_temp)[idx_temp]
                    spikes_temp = np.array(spikes_temp)[idx_temp]
                    ipsp_ampl_temp = np.array(ipsp_ampl_temp)[idx_temp]

                    
                    used, spike = 0, 0
                    # find which of the electrodes has spike (going from maximum)
                    if sum(spikes_temp) != 0:
                        while spike == 0 or used < len(electr_temp):
                            max_idx = np.argmax(ipsp_ampl_temp)
                            spike = spikes_temp[max_idx]
                            used = used+1
                    else:
                        max_idx = np.argmax(ipsp_ampl_temp)
                        spike =  ipsps_temp[max_idx]
                        
                    electr = electr_temp[max_idx]
                    spw_start = spike
                    
                    #used_beginning = sp_sp_used[sp_sp_used[] == spike] 
                    #spikes_temp.append([sp_sp_used[idx_sp]])
                    #spikes_temp = np.rec.fromrecords(spikes_temp, dtype=sp_sp_used.dtype)
                    #beginnings.append(spikes_temp)
                    #beginnings.append(np.rec.fromarrays([electr, trace, spw_no], 
                    #                       names='electrode,trace, ipsp_start, ipsp_end, amplitude, swp_no'))
                    
                    #import pdb; pdb.set_trace()    

                ipsps_temp, spikes_temp, ipsp_ampl_temp, electr_temp, increase_temp = [], [], [], [], []

            # collect all potentially the same IPSPs
            ipsps_temp.append(ipsp_time)
            # find closest spike which is not further that shift_spike from the detected beginning

            electr_memb = [i in [electr-1, electr, electr+1] for i in sp_sp_used['electrode']]
            sp_temp = sp_sp_used[np.where(electr_memb)]['spikes']
            sp_temp = sp_temp[np.where(sp_temp >= ipsp_time - shift_spike)]
            if len(sp_temp) > 0 and sp_temp[0] - shift_spike <= ipsp_time:
                spikes_temp.append(sp_temp[0])
            else:
                spikes_temp.append(0)
            # calculate amplitude of the ipsp (from the beginning to the highest peak)
            start_pts = ms2pts(ipsp_time, fs).astype(int)
            end_pts = ms2pts(ipsp_end,fs).astype(int)
            data_temp_ipsp = data_temp[data_temp['electrode'] == electr]['time'][start_pts: end_pts]
            maxpt = np.argmax(data_temp_ipsp)
            ipsp_ampl = data_temp_ipsp[maxpt] - data_temp_ipsp[0]
            ipsp_ampl_temp.append(ipsp_ampl)
            
            increase = data_temp_ipsp[-1] - data_temp_ipsp[0]
            increase_temp.append(increase)
            electr_temp.append(electr) 
            
            
            ipsp_old = ipsp_time
            
        import pdb; pdb.set_trace() 
        proper_ipsps.append(np.rec.fromarrays([beginnings['electrode'], beginnings['trace'], beginnings['spw_no'],
                                            beginnings['spikes']], 
                                           names='electrode,trace, spw_no, ipsp_start'))
        
        
        #import pdb; pdb.set_trace()
        if len(beginnings) > 0:
            #if mean()
            beginnings = np.concatenate(beginnings)
            #beginnings = np.unique(beginnings)
            begs_pts = beginnings['spikes']
            begs_pts = ms2pts(begs_pts, fs).astype(int)
            #import pdb; pdb.set_trace()
            
            proper_ipsps.append(np.rec.fromarrays([beginnings['electrode'], beginnings['trace'], beginnings['spw_no'],
                                                        beginnings['spikes']], 
                                                       names='electrode,trace, spw_no, ipsp_start'))
            
        if plot_it:    
            spw_min_start = 9000000000
            spw_max_end = -1
            #import pdb; pdb.set_trace() 
            for electr in np.unique(sp_sp_used['electrode']):
    
                spw_start = sp_sp_used[sp_sp_used['electrode'] == electr]['spw_start'][0]
                spw_end = sp_sp_used[sp_sp_used['electrode'] == electr]['spw_end'][0]
                spw_st_pts, spw_en_pts = ms2pts(spw_start, fs).astype(int), ms2pts(spw_end, fs).astype(int)
                data_used = data_temp[data_temp['electrode'] == electr]['time'][spw_st_pts:spw_en_pts]
                spikes = ms2pts(sp_sp_used[sp_sp_used['electrode'] == electr]['spikes'], fs).astype(int) - spw_st_pts
                ipsp_start = ms2pts(sp_ip_used[sp_ip_used['electrode'] == electr]['ipsp_start'], fs).astype(int) - spw_st_pts
                
                spw_min_start = min(spw_st_pts, spw_min_start)
                spw_max_end = max(spw_en_pts, spw_max_end)
                
                t = dat.get_timeline(data_used, fs, 'ms') + spw_start
                plt.plot(t,data_used + add_it * electr)
                plt.plot(t[spikes],data_used[spikes] + add_it * electr, 'ro')
                plt.plot(t[ipsp_start],data_used[ipsp_start] + add_it * electr, 'gx')
            
            data_used = data_temp[data_temp['electrode'] == electr]['time'][spw_min_start:spw_max_end]
            t = dat.get_timeline(data_used, fs, 'ms') + pts2ms(spw_min_start, fs)
            if len(beginnings) > 0:
                begs_pts = begs_pts - spw_min_start
                #t_beg =
                for a in t[begs_pts]:
                    plt.vlines(a, -200, 1200)
                #plt.plot(t[begs_pts], data_used[begs_pts], 'bo')
                #import pdb; pdb.set_trace() 
            tit = 'spw: ' + str(spw) 
            plt.title(tit)
            fig.savefig(save_folder + save_fig + str(spw) + ext,dpi=600)
            fig.savefig(save_folder + save_fig + str(spw) + '.png',dpi=600)        
            #plt.show()
            plt.clf()
    proper_ipsps = np.concatenate(proper_ipsps)
    np.savez(save_folder + save_file, ipsps = proper_ipsps) 
    del spw_ipsps, spw_spikes, data

def update_SPW_ipsp_ampl(save_folder, save_file, data, fs):
    pass

def update_ipsp_exSpikes(save_folder, save_file):
    pass

def update_SPW_ipsp(load_datafile, load_spwsspike, save_folder, save_file):
    # it looks for the ipsps within detected spws - separate for each electrode
    # it take very long to analyse so be patient!
    npzfile         = np.load(save_folder + load_datafile)
    data            = npzfile['data']
    fs              = npzfile['fs']
    npzfile.close()
    
    # load starts of spws
    npzfile         = np.load(save_folder + load_spwsspike)
    spw_spikes      = npzfile['spw_details']
    npzfile.close()
    
    print
    print "analyzing SPWs in electrode:",
    plot_it = False
    add_it = 100
    window = 0.5 # ms for calculating moving average
    window = ms2pts(window, fs)
    spw_ipsps = []
    # take each SPW separately and find ipsps
    #import pdb; pdb.set_trace() 
    spw_len = len(np.unique(spw_spikes['spw_no']))
    for spw in np.unique(spw_spikes['spw_no']):
        print str(spw) + '/' + str(spw_len)
        
        spw_spike_used = spw_spikes[spw_spikes['spw_no'] == spw]
        for electr in np.unique(data['electrode']):
            # check if in this electrode spw was recorded

            exists = electr in np.unique(spw_spike_used['electrode'])                
            
            # this electrode was recorded
            if exists:
                spw_electr_used = spw_spike_used[spw_spike_used['electrode'] == electr]
                trace = int(spw_electr_used['trace'][0])
                start_ms = spw_electr_used['spw_start'][0]
                end_ms = spw_electr_used['spw_end'][0]
                start  = ms2pts(start_ms, fs).astype(int)
                end = ms2pts(end_ms, fs).astype(int)
                spw_no = spw_electr_used['spw_no'][0]
                
                data_used = data[(data['electrode'] == electr) & (data['trace'] == trace)]['time'][start:end]
                          
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
                    t = dat.get_timeline(data_used, fs, 'ms') + spw_electr_used['spw_start'][0]
                    plt.plot(t, data_used + add_it*electr)   
                    plt.plot(t[maxs], data_used[maxs] + add_it*electr, 'r<')
                    plt.plot(t[mins], data_used[mins] + add_it*electr, 'b<')
                    
                typ = 'f8'
                mini = pts2ms(maxs,fs) + spw_electr_used['spw_start'][0]
                mini = np.append(mini,spw_electr_used['spw_end'][0])

                ipsp_start = mini[:-1]
                ipsp_end = mini[1:]
                electrodes = np.ones(len(ipsp_start), dtype=typ)*electr
                traces = np.ones(len(ipsp_start), dtype=typ)*trace
                spw_num = np.ones(len(ipsp_start), dtype=typ)*spw_no
                ipsp_no = range(1,len(mini))


                spw_ipsps.append(np.rec.fromarrays([electrodes, traces, spw_num,
                                                        ipsp_start, ipsp_end, ipsp_no], 
                                                       names='electrode,trace, spw_no, ipsp_start, ipsp_end, ipsp_no'))
                
                
                if plot_it: 
                    plt.plot(t, moved_avg + add_it*electr)
        #import pdb; pdb.set_trace() 

        if plot_it:
            plt.show()
    spw_ipsps = np.concatenate(spw_ipsps)
    np.savez(save_folder + save_file, spw_ipsps = spw_ipsps) 
    del data, spw_spikes 
      


def update_SPWspikes(load_spikefile, load_spwsfile, save_folder, data_file):
    # reject all the SPWs where there is less than min_no_wave spws detected, 
    # check which spws belog together and save them
    npzfile         = np.load(save_folder + load_spikefile)
    spike_idxs      = npzfile['spike_idxs']
    fs              = npzfile['fs']
    npzfile.close()
    
    # load starts of spws
    npzfile         = np.load(save_folder + load_spwsfile)
    spw_starts      = npzfile['starts']
    spw_ends        = npzfile['ends']
    npzfile.close()
    
    
    print
    print "analyzing SPWs in electrode:",
    spw_details = []
    min_no_wav = 3
    max_move = 30 # ms

    spw_details_temp = []
    spw_no = 0
    no_traces = len(np.unique(spw_starts['trace']))
    for trace in np.unique(spw_starts['trace']):
        print 'trace: ' + str(trace + 1) + '/' + str(no_traces)
        # get spws for each electrode and check if it's the same
        spw_st_trace = spw_starts[(spw_starts['trace'] == trace)]
        spw_en_trace = spw_ends[(spw_ends['trace'] == trace)]
        
        sort_idx = np.argsort(spw_st_trace['time'])
        spw_st_sorted = spw_st_trace[:,:,sort_idx]
        spw_en_sorted = spw_en_trace[:,:,sort_idx]
        
        same = 0
        start_init, end_init = 0, 0 
        
        # analize every single SPW
        for idx, next_spw_st in enumerate(spw_st_sorted):
            next_spw_en = spw_en_sorted[idx]
            
            # check if start is enclosed between start and end of previous one and beginning is not 
            # further than max_move from new beginning
            st_same = (start_init < next_spw_st['time']) 
            en_same = (end_init > next_spw_st['time'])
            dist_same = abs(next_spw_st['time'] - start_init) < max_move
            
            if st_same and en_same and dist_same:
                same = same + 1
                start_init = min(start_init, next_spw_st['time'])
                end_init = max(end_init, next_spw_en['time'])              

            elif same >= min_no_wav:
                # checks if this wave was detected in enough electrodes
                same = 0
                # save previous spws   
                spw_details = spw_details + spw_details_temp
                start_init = next_spw_st['time']
                end_init = next_spw_en['time']        
                spw_details_temp = []
                spw_no = spw_no + 1    
            else:
                # something wrong with previous SPW = remove it
                spw_details_temp = []
                start_init = next_spw_st['time']
                end_init = next_spw_en['time']
                same = 0

                
            # check the beginning of this SPWs
            
            electr = next_spw_st['electrode']
            trace = next_spw_st['trace']
            
            # read all the spikes for this electrode between start and end
            spikes_used  = get_var_used(spike_idx, electr, trace)
            spikes_usedMs = spikes_used[(spikes_used < next_spw_en['time'])&(spikes_used > next_spw_st['time'])]

            # save temp details of this spw
            #spiki = spikes_usedMs
            typ = 'f8'
            electrodes = np.ones(len(spikes_usedMs), dtype=typ)*electr
            start = np.ones(len(spikes_usedMs), dtype=typ)*next_spw_st['time']
            end = np.ones(len(spikes_usedMs), dtype=typ)*next_spw_en['time']
            traces = np.ones(len(spikes_usedMs), dtype=typ)*trace
            spw_num = np.ones(len(spikes_usedMs), dtype=typ)*spw_no

            
            #start = np.array(next_spw_st['time'], dtype=typ)
            #end = np.array(next_spw_en['time'], dtype=typ)
            #electr
            #for sp in spikes_usedMs:
            #
            #valley_to_peak_norm.append(np.rec.fromarrays([electrodes, traces, np.array(v1, dtype=typ)], names='electrode,trace,time'))  
            spw_details_temp.append(np.rec.fromarrays([electrodes, traces, spikes_usedMs ,start  , end, spw_num], names='electrode,trace, spikes, spw_start, spw_end, spw_no'))
    
    #import pdb; pdb.set_trace()          
    # do the same check as before but for the last SPW      
    if same >= min_no_wav - 1:
        # checks if this wave was detected in enough electrodes and save
        spw_details = spw_details + spw_details_temp
    spw_details = np.concatenate(spw_details)   
    np.savez(save_folder + data_file, spw_details = spw_details) 
    del spike_idxs, spw_starts, spw_ends

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
    
    data_bas = []
    for electr in np.unique(data['electrode']): 
        # remove all mean from each electrode 
        print electr,
        
        electro_data = data[data['electrode'] == electr]
        mean_datatime = np.mean(electro_data['time'])
        
        for trace in np.unique(data[data['electrode'] == electr]['trace']):
            data[(data['electrode'] == electr) & (data['trace'] == trace)]['time'] = data[(data['electrode'] == electr) & (data['trace'] == trace)]['time'] - mean_datatime
  
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
