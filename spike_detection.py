import scipy.signal as signal

from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle

import numpy as np
from matplotlib import pyplot
import pylab as py
import filtit as filt
import plotting as plt
import detect_waves as dw
from scipy import signal
import data_mang as dat
#import math as math
import induc_SPW as ispw

def max_leftright(data, pt_idx):
    """ find local max from left and right to pt_idx"""
    change = np.diff(data)
    left = change[:pt_idx-2]
    if len(np.where(left >= 0)[0]) > 0:
        mleft = np.where(left >= 0)[0][-1] + 1
    else:
        mleft = 1
        
    right = change[pt_idx+2:]        
    if len(np.where(right <=0)[0]) > 0:
        mright = np.where(right <=0)[0][0] + pt_idx+1 
    else:
        mright = len(right) + pt_idx+1
    return mleft, mright

def ms2pts(ms, fs):
    """ converts number of ms to number of data pts"""
    pts = (fs / 1000.0) * ms
    return pts

def pts2ms(pts, fs):
    """ converts number of pts to number of ms"""
    ms = (pts*1./fs)*1000.0
    return ms
   
def find_extra_spikes(data, fast_data, fs): #, fast_freq = 750.):
    """ finds spikes in the extracellular data, of fs sampling frequency
    it will first highpass filter the data above vfast_freq to detect the spikes
    """

    allow_error = 0.2 # ms
    alerr_pts = ms2pts(allow_error,fs)
    
    # find spikes 
    spike_data = fast_data * (-1)
    spike_data[spike_data < 0] = 0
    thres = np.std(spike_data)*5 # use thres as 6 * standard deviation
    
    waves = dw.find_above(spike_data, thres)
    starts, ends = dw.find_startend(waves) # finds all the peaks above threshold
    
    mins, min_idxs = dw.max_waves(data*(-1), starts-int(alerr_pts), ends+int(alerr_pts)) # finds mins of each wave    
    scale = 'ms'

#    t = dat.get_timeline(data, fs, scale)    
#    py.plot(t, data)
#    py.plot(t, spike_data)
#    py.plot(t[min_idxs], spike_data[min_idxs], 'o')
#    py.show()
    return mins, min_idxs
    
#
def cut_spikes(data, fs, min_idxs, rang = [1, 3]):
    """ cut out spikes from the given data 
    between min_idxs-rang[0]: min_idxs+rang[1]
    rang is in ms
    """
    waves, min_idxs = dw.cut_waves(data, min_idxs, fs, rang[0], rang[1])
    return waves, min_idxs


def find_halfampl(data, data_stat, data_smooth, fs, min_idxs, rang = [1, 3]): #, fast_freq = 500., smooth_freq = [500, 2000]):
    """ detects values such as:
    a  valley to peak
    b  half valley width
    data filtered to fast_freq     """ 
    
    # cut the dataparts to given range
    waves, min_idxs2 = cut_spikes(data_stat, fs, min_idxs) # fast data
    waves_smooth, min_idxs2 = cut_spikes(data_smooth, fs, min_idxs) # slow data
    norm_waves, min_idxs2 = cut_spikes(data, fs, min_idxs) # slow data
    
    # initialize variables to find
    half_bs, half_as, ampls, relative_ampl = [], [], [], []
    As, Bs, norm_factor = [], [], []
    peak = ms2pts(rang[0],fs) # place is always peak (in pts)
    
    # go through every possible wave (spike)
    for i in range(len(waves)):

        
        if len(waves[i]) > 0:
            # finding the maxs from left and right
            #mleft, mright = max_leftright(waves[i], peak)
            mleft_smooth, mright_smooth = max_leftright(waves_smooth[i], peak)
            half_a = mright_smooth - peak
            
            # finding the end of the spike
            #temp, spadek = max_leftright(waves_smooth[i] * (-1), mright_smooth)

            # get the amplitude of the spike before normalization (calculated on the original data)
            maxRight = norm_waves[i][mright_smooth]
            maxLeft = norm_waves[i][mleft_smooth]
            top = max(maxRight, maxLeft)
            temp_ampl = top - norm_waves[i][peak]
#            if temp_ampl < 0:
#                t = dat.get_timeline(waves_smooth[i], fs, 'ms')
#                py.figure()
#                py.plot(t, norm_waves[i], 'r')
#                py.plot(t, waves[i], 'g')
#                #py.plot(t, norm_waves2, 'b')
#    #            
#                #py.plot(t[mleft], norm_waves[i][mleft], 'ob')
#                #py.plot(t[mright], norm_waves[i][mright], 'ob')
#                py.plot(t[mleft_smooth], norm_waves[i][mleft_smooth], 'og')
#                py.plot(t[mright_smooth], norm_waves[i][mright_smooth], 'og')   
#                #y.plot(t[find_start], norm_waves[i][find_start], '<r')          
#                
#                #py.vlines(t[peak], norm_waves[i][peak], norm_waves[i][peak]+temp_ampl)
#                py.show()
            # and relative amplitude: dist(maxleft to maxright)/dist(max_right to peak)
            dist1 = norm_waves[i][mright_smooth] - norm_waves[i][mleft_smooth]
            dist2 = np.abs(norm_waves[i][mleft_smooth] - norm_waves[i][peak])
            temp_rel_amp = dist1/dist2
            
            # normalize
            # move the waves to the baseline
            based = norm_waves[i][mleft_smooth]
            norm_waves[i] = norm_waves[i] - based
            max_pt = norm_waves[i][mleft_smooth]
            min_pt = norm_waves[i][peak]
            temp_norm_factor = np.abs(norm_waves[i][peak])
            
            if min_pt == 0 or temp_ampl <= 0 or temp_norm_factor < 0: # or min_pt_smooth == 0:
                #print 'appending 0 because amplitude is equal 0'
                half_bs.append(0)
                half_as.append(0)
                #half_cs.append(0)
                ampls.append(0)
                relative_ampl.append(0)
                As.append(0)
                Bs.append(0)
                norm_factor.append(0)
                continue              
            else:
                
                
                norm_waves[i] = norm_waves[i]/temp_norm_factor
            
            # calculate half height as half between base and min (peak of the spike)
            min_pt = -1.
            half_pt = np.mean([max_pt, min_pt])
            
            # find which part belongs to the spike below the half hight and calculate the width
            all_below = np.zeros(np.size(norm_waves[i]))
            all_below[waves[i] > half_pt] = 1
            
            # calculate half_c (detection of D is not perfect so better not to use)
#            C_c = norm_waves[i][mright_smooth]
#            D_c = norm_waves[i][spadek]
#            half_hC = (C_c + D_c)/2.0
#    
#            all_below_c = np.ones(np.size(waves_smooth[i]))
#            all_below_c[waves_smooth[i] > half_hC] = 0    

#            t = dat.get_timeline(waves_smooth[i], fs, 'ms')
#            py.figure()
#            py.plot(t, norm_waves[i], 'r')
#            py.plot(t, waves[i], 'g')
#            py.plot(t, norm_waves2, 'b')
##            
#            #py.plot(t[mleft], norm_waves[i][mleft], 'ob')
#            #py.plot(t[mright], norm_waves[i][mright], 'ob')
#            py.plot(t[mleft_smooth], norm_waves[i][mleft_smooth], 'og')
#            py.plot(t[mright_smooth], norm_waves[i][mright_smooth], 'og')   
            #y.plot(t[find_start], norm_waves[i][find_start], '<r')          
            
            #py.vlines(t[peak], norm_waves[i][peak], norm_waves[i][peak]+temp_ampl)
            #py.show()                
            try: 
#                find_start_c = np.nonzero(all_below_c[:mright_smooth])[0][-1]
                #print len(np.nonzero(all_below[:peak])[0])

                if len(np.nonzero(all_below[:peak])[0]) == 0:
                    find_start = mleft_smooth
                else: 
                    find_start = np.nonzero(all_below[:peak])[0][-1]
                    if find_start < mleft_smooth:
                        find_start = mleft_smooth
                
                find_end = np.nonzero(all_below[peak:])[0][0]+peak 

#                find_end_c = np.nonzero(all_below_c[mright_smooth:])[0][0]+mright_smooth            
#                if find_end > mright:
#                    find_end = mright
#    
            except IndexError:


                
                #print 'appending 0 because of half widths'
                half_bs.append(0)
                half_as.append(0)
                relative_ampl.append(0)
                ampls.append(0)
                As.append(0)
                Bs.append(0)
                norm_factor.append(0)
#                continue 
              
            half = find_end - find_start
            half_bs.append(half)
            half_as.append(half_a)
            ampls.append(temp_ampl)
            
            relative_ampl.append(temp_rel_amp)
            pk = ms2pts(rang[0], fs)

            As.append(mleft_smooth - pk)
            Bs.append(mright_smooth - pk)
            norm_factor.append(temp_norm_factor)
    
            t = dat.get_timeline(waves[i], fs, 'ms')
        else:        
            
            #print 'appending 0 because wave is too short'
            half_bs.append(0)
            half_as.append(0)
            #half_cs.append(0)  
            ampls.append(0)
            relative_ampl.append(0)
            As.append(0)
            Bs.append(0)
            norm_factor.append(0)
            
    return half_as, half_bs, ampls, As, Bs, norm_factor
