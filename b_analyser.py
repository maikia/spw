#!/usr/bin/env python
#coding=utf-8

import numpy as np
import matplotlib.pyplot as plt
import os

import induc_SPW as ispw

def find_max_corr(x, y):
    
    nomean_x = x - np.mean(x)
    nomean_y = y - np.mean(y)

    cxy = np.correlate(nomean_x,nomean_y,mode='full')
    i = cxy.argmax()
    lag = i - len(cxy)/2
    return lag
    

def extract_spws(x, indices, fs, win):

    l = win[1]-win[0]
    ind = []
    spws = np.zeros((l, len(indices)))
    wrong = []
    for i, idx in enumerate(indices):
        electrode, trace, time = idx
        if time+win[0] > 0 and time+win[1] < np.size(x,2):
            spws[:,i] = x[electrode, trace, time+win[0]:time+win[1]]
        else:
            wrong.append(i)
            #indices.remove(i)
            #indices.remove(indices['time'] == i)
            #print 'error!'
    if len(wrong) > 0:
        np.delete(indices,wrong)
        np.delete(spws,wrong)
    return spws, indices


def align_spws(data, idx, fs, win = (-40, 60)):


    spws, idx = extract_spws(data, idx, fs, win)

    spw_mean = spws.mean(1)
    
    lags = [find_max_corr(spw_mean, spws[:,i]) 
            for i in range(spws.shape[1])]
    idx = idx.copy()

    idx['time'] = idx['time']-lags

    return idx

def correlation_matrix(spws, subtract_mean=True):
    spws = spws.copy()
    if subtract_mean:
        spws = spws - spws.mean(1)[:, np.newaxis]
    spws = spws - spws.mean(0)
    
    return np.corrcoef(spws.T)


def load_data(folder):

    npz_data = np.load(os.path.join(folder, 'data_bas.npz'))
    fs = npz_data['arr_1']
    data = np.array(npz_data['arr_0'])
    
    npz_data.close()

    npz_data = np.load(os.path.join(folder, 'SPWs.npz'))
    idx = idx2recarr(npz_data['arr_2'])
    
    npz_data.close()

    return data, fs, idx

def idx2recarr(idx, time_dtype='i4'):
    n_el,n_traces = idx.shape
    arrs = []
    for el in range(n_el):
        for tr in range(n_traces):
            indices = np.array(idx[el][tr], dtype=time_dtype)
            electrodes = np.ones(len(indices), dtype=np.int32)*el
            traces = np.ones(len(indices), dtype=np.int32)*tr
            temp_arr = np.rec.fromarrays([electrodes, traces, indices],
                                       names='electrode,trace,time')
            arrs.append(temp_arr)
    rec_arr = np.concatenate(arrs)
    return rec_arr


if __name__=='__main__':

    #path = '/Users/bartosz/Desktop/maja/saved_data/cell6/gap_free' 
    path = '/home/maja/PhDProject/SPWs/SPWs/saved_data/cell6/gap_free'
    data, fs, idx = load_data(path)
    
    win = (-10, 30)
    win_pts = (ispw.ms2pts(win[0], fs), ispw.ms2pts(win[1], fs))
    
    #calculate and plot mean spw in each electrode
    spws, idx = extract_spws(data, idx, fs, win_pts)
    
    spw_mean = np.array([spws[:, idx['electrode']==i].mean(1) 
                         for i in np.unique(idx['electrode'])]).T
    plt.figure()
    #plt.plot(spw_mean)
    lines = plt.plot(spw_mean)
    labels = ['El. %d' % i for i in np.unique(idx['electrode'])]
    plt.legend(lines, labels)
    
    # correlation matrix in a single electrode
    
    el_id = 0
    
    idx_el  = idx[idx['electrode']==el_id]

    idx_aligned = align_spws(data, idx_el, fs, win_pts)
    spws_el, idx_aligned = extract_spws(data, idx_aligned, fs, win_pts)
    cxy = correlation_matrix(spws_el)
    
    plt.figure()
    plt.imshow(cxy, interpolation='nearest', origin='lower')
    
    plt.figure()
    plt.subplot(211)
    plt.plot(spws_el)
    plt.subplot(212)
    plt.plot(spws[:,idx['electrode']==el_id]) #plot all traces in
                                              # electrode el_id
    
    # idx[idx['electrode']==1]['trace'] # all traces for el. 1
    
    # find sws coincidences
    bin = 5.
    n_el = np.max(idx['electrode'])+1
    events = np.arange(idx['time'].min(), idx['time'].max(), bin*fs/1000.)
    bin_id = np.searchsorted(events, idx['time'])
    ev_id = np.searchsorted(np.unique(bin_id), bin_id)
    
    spw_locs = np.empty((ev_id[-1]+1, n_el))
    spw_locs.fill(np.inf)

    spw_locs[ev_id, idx['electrode']] = idx['time']
    
    plt.figure()
    plt.imshow(~np.isinf(spw_locs), aspect='auto', interpolation='nearest')
    
    # fill in spws in "empty" electrodes
    min_spws_per_bin=3
    has_enough_spws = np.sum(~np.isinf(spw_locs),1)>min_spws_per_bin
    spw_locs_selected = spw_locs[has_enough_spws,:]
    i,j = np.where(np.isinf(spw_locs_selected))
    fill_values = spw_locs_selected.min(1)
    spw_locs_selected[i,j] = fill_values[j]
    
    
    plt.show()
    







