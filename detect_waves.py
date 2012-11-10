from scipy import fftpack, signal
import numpy as np
import math as math
import pylab as py
import datetime

def points_fromMs(points, fs):
    """ calculates number of points from number of ms"""
    return (fs / 1000.) * points

def ms2pts(ms, fs):
    """ converts number of ms to number of data pts"""
    pts = (fs / 1000.0) * ms
    return pts

def pts2ms(pts, fs):
    """ converts number of pts to number of ms"""
    ms = (pts/fs)*1000.0
    return ms


def find_above(data, thres):
    """ gives True for each index above and False
    for each index below thres; returns either:
    above - in bulean
    waves - in 01 format"""
    above = (data > thres)
    waves = np.zeros(len(data))
    waves[above != True] = 0
    waves[above == True] = 1 
    return waves

def find_startend(above):
    """ above must be give in 0 for no event and 1 for event 
    format"""
    shifted = np.array(above)
    temp = shifted[0]
    shifted[:-1] = above[1:]
    shifted[-1] = temp
    st_end = shifted - above # beginnings and ends of waves
    starts = np.where(st_end == 1)
    ends = np.where(st_end == -1)
    starts = starts[0]
    ends = ends[0]
    # remove waves if wave have already started or didn't finish in the course
    # of the recordings

    if len(starts) > 1:
        if starts[0] > ends[0]:
            ends = np.setdiff1d(ends, ends[0])
  
        if ends[-1] < starts[-1]:
            starts = np.setdiff1d(starts, starts[-1])
    return starts, ends+1

def max_waves(data, starts, ends):
    """ returns the maximums and its indexes of each wave 
    defined by indexes: starts and ends in the data"""
    maxs = []
    idxs = []
    allow_err = 0 # mV
    for wave in range(len(starts)):
        loc_min = 0
        st = starts[wave]
        ed = ends[wave]
        while loc_min == 0:
            if st >= ed:
                loc_min = -1   
            else:
                #loc_min = np.diff(data[starts[wave]:ends[wave]])
                # check if the start and end are not less or more than the format of the data
                if st < 0:
                    st = 0
                if ed >= len(data):
                    ed = len(data)-1
                dat_temp = data[st:ed]
                m_ax = max(dat_temp)
                
                idx = st + np.where(dat_temp == m_ax)[0][0]
                # check if truely data min 
                if data[idx - 1] > data[idx] + allow_err:
                    # there are larger points from the left
                    #print 'removing left'
    #                py.figure()
    #                t = range(0,len(data))
    #                py.plot(t, data)
    #                py.plot(t[st:ed], data[st:ed], 'r')
    #                py.plot(t[idx], data[idx], 'ro')
    #                py.xlim([st-50,ed+50])              
    #                py.show()  
                    st = st+1                              
                    
                elif data[idx + 1] > data[idx] + allow_err:
                    # there are larger points from the right
                    #print 'removing right'
    #                py.figure()
    #                t = range(0,len(data))
    #                py.plot(t, data)
    #                py.plot(t[st:ed], data[st:ed], 'r')
    #                py.plot(t[idx], data[idx], 'ro')
    #                py.xlim([st-50,ed+50])    
    #                py.show()
                    ed = ed-1    
      
                else:
                    loc_min = 1
                if ed - st < 2:
                    loc_min = -1                          
        
#        py.plot()
#        t = range(0,len(data))
#        py.plot(t, data)
#        py.plot(t[st:ed], data[st:ed], 'r')
#        py.plot(t[idx], data[idx], 'ro')
#        py.xlim([st-20,ed+20])  
        if loc_min == 1:       
            #t = range(0,len(data))
            maxs.append(m_ax)
            idxs.append(idx)
    return maxs, idxs

def length_waves(starts, ends):
    """ calculates the length of the waves"""
    length = ends - starts
    return length

def dist_waves(starts, ends):
    """ calculates the distance between the waves; output array is one shorter
    than each of the input arrays"""
    s = starts[1:]
    e = ends[:-1]
    dist = length_waves(e, s)
    return dist

def dist_small(dist, dist_min):
    """ returns too small indexes """
    small_idx = np.where(dist < dist_min)
    return small_idx[0][:]

def remove_from(dist, rem_idx):
    """ removes all the dist which are smaller than dist_min """
    dist = np.setdiff1d(dist, rem_idx)
    return dist
    
def find_locals(data):
    """ finds local minimas and maximas
     it will not work on very noisy data """
    maxs = np.r_[True, data[1:] > data[:-1]] & np.r_[data[:-1] > data[1:], True]
    mins = np.r_[True, data[1:] < data[:-1]] & np.r_[data[:-1] < data[1:], True]
    return maxs, mins

def correct_wave(data, starts, ends, fs, min_dist=0, min_length=0):
    """ if the two waves are too close to each other, they are checked 
    for the correct length and height; if either length or hight is too small they are 
    concatenated to become one wave; if any of the waves is too small or too short and 
    there is no other wave close by - it is removed from the list
    min_dist [ms], min_length [ms]"""
    min_len = points_fromMs(fs, min_length)
    min_dist = points_fromMs(fs, min_dist)

    maxs, max_idx = max_waves(data, starts, ends)
    dist = dist_waves(starts, ends)
    small_idx = dist_small(dist, min_dist)
    s_no = []
    e_no = []
    l_waves = length_waves(starts, ends)
    # concatenate too short waves and update starts and ends
    for idx in small_idx:
        s = starts[idx + 1]
        e = ends[idx]
        #lenOne = ends[idx] - starts[idx]
        #lenTwo = ends[idx + 1] - starts[idx + 1]
        #if lenOne < min_len or lenTwo < min_len:
        s_no.append(s)
        e_no.append(e)
          
    starts = np.setdiff1d(starts, s_no)
    ends = np.setdiff1d(ends, e_no)
    # check for minimum length and remove events which are too short
    lengths = ends - starts
    too_short = np.where(lengths < min_len)[0]
    starts = np.setdiff1d(starts, starts[too_short])
    ends = np.setdiff1d(ends, ends[too_short])
    
    s_no = []
    e_no = []
    max_min = 0.7 #0.3 # var for setting how much of max must be smallest part of the wave to still form one wave
    
    # concatenate the wave which are one
    for idx in range(len(starts)):
        if len(starts) != (idx + 1) and ends[idx] < starts[idx + 1]:
            min_between = min(data[ends[idx]:starts[idx + 1]])
            max1 = max(data[starts[idx]:ends[idx]])
            max2 = max(data[starts[idx + 1]:ends[idx + 1]])
            
            if max1 * max_min < min_between or  max2 * max_min < min_between:
                s = starts[idx + 1]
                e = ends[idx]
                s_no.append(s)
                e_no.append(e)

    starts = np.setdiff1d(starts, s_no).tolist()
    ends = np.setdiff1d(ends, e_no).tolist()



    
    # divide the waves which are not one                                  
#    for idx in range(len(starts)):
#        part1 = data[starts[idx]:((ends[idx] - starts[idx]) / 2) + starts[idx]].tolist()
#        part2 = data[((ends[idx] - starts[idx]) / 2) + starts[idx] + 1: ends[idx]].tolist()
#        between = part1[len(part1) / 2:] + part2[:len(part2) / 2]
#        
#        if max(part1) * max_min > min(between) or max(part2) * max_min > min(between):
#            end_2wave = data[starts[idx]:ends[idx]].tolist()
#            start_2wave = data[starts[idx]:ends[idx]].tolist()
#            starts.append(start_2wave.index(min(start_2wave)) + starts[idx])
#            ends.append(ends[idx])
#            ends[idx] = between.index(min(between)) + starts[idx] + (len(part1) / 2)
    starts.sort()
    ends.sort()
        
    return starts, ends

def check_width(data, locs, heights, fs):
    """ at the heights given, at the locations (locs) given check what is the length
    of the event, it also returns begining and end of each wave"""
    lengths = []
    ends = []
    starts = []
    before = 100 # ms
    after = 100 # ms
    before_pts = ms2pts(before, fs)
    after_pts = ms2pts(after, fs)
    #t0 = datetime.datetime.now()
    # now do something that consumes time
    for idx in range(len(locs)):
        #if idx - before_pts > 0 and idx + after_pts < len(data): # <-check this!!!
        # check where the wave starts
        #print 'new'
        #print  datetime.datetime.now() - t0 
        #t0 = datetime.datetime.now()
        temp = np.where(data[locs[idx] - before_pts:locs[idx]] < heights[idx]) + locs[idx] - before_pts
        #print temp
        #temp = [data[:locs[idx]] for ]
        
        #temp = 
        #print 'where 1'
        #print datetime.datetime.now() - t0 
        #t0 = datetime.datetime.now()
        if np.size(temp) > 0:
            str = temp[0][-1]
            starts.append(str)
            # check where the wave ends
            temp = np.where(data[locs[idx] + 1:locs[idx] + after_pts] < heights[idx]) + locs[idx]
            #print temp
#            temp = 
            #print '1'
            #print datetime.datetime.now() - t0 
            #t0 = datetime.datetime.now()
        else:
            starts.append(0)
            #print '2'
            #print datetime.datetime.now() - t0 
            #t0 = datetime.datetime.now()
        if np.size(temp) > 0:
            en = temp[0][0]
            lengths.append(en - str)
            ends.append(en)
            #print '3'
            #print datetime.datetime.now() - t0 
            #t0 = datetime.datetime.now()      
        else:
            lengths.append(0)
            ends.append(0)
            #print '4'
            #print datetime.datetime.now() - t0 
            #t0 = datetime.datetime.now()                     
    return lengths, starts, ends

def concat_waves(data, maxs, fs, min_dist=10):
    # given the different maximums, and min length of the wave, and the data
    # if the two maximums are too close together - it take out the smaller out
    min_dist = (fs / 1000) * min_dist

    m_no = []
    small_idx = [x - maxs[i - 1] for i, x in enumerate(maxs)][1:]

    small_idx = (small_idx < min_dist)
    small_idx = np.where(small_idx == True)

    for idx in small_idx[0]:
        if data[maxs[idx]] < data[maxs[idx + 1]]:
            m_no.append(maxs[idx])
        elif data[maxs[idx]] > data[maxs[idx + 1]]:
            m_no.append(maxs[idx + 1])
    maxs = np.setdiff1d(maxs, m_no)
    return maxs       

def cut_waves(data, point_idx, fs, before, after):
    """before and after are given in ms"""
    wave = []
    before = points_fromMs(fs, before)
    after = points_fromMs(fs, after)
    #print np.shape(data)
    for idx in point_idx:

        if (idx - before < 0):         
            point_idx = point_idx[1:]
            wave.append([])
        elif (idx + after > len(data)):
        # if wave finishes or ends outside the data set, remove the wave
            #print 'got in'
            #point_idx.remove(idx)
            point_idx = point_idx[:-1]
            wave.append([])
            #print 'removed'
        else:
            spw = data[int(idx - before): int(idx + after)]
            wave.append(spw)
    #print np.shape(point_idx)
    return wave, point_idx

if __name__ == '__main__':
    pass
