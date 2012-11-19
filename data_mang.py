'''
Created on 25 juin 2012

@author: Maja
'''

import numpy as np
import neo.io
import math

def get_timeline(data, fs, scale='ms'):
    """ returns the timeline of the given data 
    in prefered scale: 'sec' (defoult), 'min' or 'ms'
    """
    divider = {'min': 1.0 / 60,
               'sec': 1.0,
               'ms' : 1000,
               }
    t = len(data) / fs * divider[scale]
    #return np.arange(0, t, (1. / fs) * divider[scale]) #(1. / fs) * divider[scale])
    return np.linspace(0, t, len(data))


def get_electrdata(data, no_segments, electrode=[1], data_part = 'all'):
    """ for given analogdata get the array of data for only one electrode and
    it's sampling rate"""
    electr = []
    #electrode = electrode - 1
    #import pdb; pdb.set_trace()
    dat = np.zeros((len(electrode), no_segments, np.size(data,2)))
    for electr in electrode:
        for trace in range(no_segments):
            el = electr - 1
            dat[el, trace, :] = data[trace][el].magnitude
        
        #electrodes = np.ones(len(dat), dtype=np.int32)*electrode
        #traces = np.ones(len(dat), dtype=np.int32)*i
        #electr.append(np.rec.fromarrays([electrodes, traces, dat], names='electrode,trace,time'))
    #electr = np.concatenate(electr)
    
    fs = data[0][0].sampling_rate
    fs = fs.magnitude

    #if data_part != 'all':
    #    # if data_part is 'all' then the whole data will be used
    #    electr = electr[:][data_part[0]: data_part[1]]
    return dat, fs
                                                                                                                                   
def generate_dummy(sim_time=60, fs=20):
    nsamps = sim_time * fs
    t = np.linspace(0, sim_time, nsamps)
    freqs = [0.1, 0.5, 1, 4]
    data = 0
    for i in range(len(freqs)):
        data += np.cos(2 * math.pi * freqs[i] * t)
    return data, t

def get_data(fname, ifsaved = 0):
    """ reads and passes the .abf data
    as for now only the first segment is considered
    ifsaved """
    print "Patience please, processing ", fname, " ..."
    reader = neo.io.AxonIO(filename=fname)
    block = reader.read_block()
    data = []
    for i in range(len(block.segments)):
        seg = block.segments[i]
        data.append(seg.analogsignals)
    return data, len(block.segments)

def get_data_oscillo(fname, ifsaved = 0):
    """ reads and passes the .abf data
    as for now only the first segment is considered
    ifsaved """
    print "Patience please, processing ", fname, " ..."
    reader = neo.io.AxonIO(filename=fname)
    block = reader.read_block() #(cascade=True, lazy = False)
    seg = block.segments[0]
    data = seg.analogsignals
    return data

#>>> r = io.AxonIO(filename='File_axon_1.abf')
#>>> bl = r.read_block(lazy=False, cascade=True)
#>>> print bl.segments
#[<neo.core.segment.Segment object at 0x105516fd0>]
#>>> print bl.segments[0].analogsignals
#[<AnalogSignal(array([ 2.18811035,  2.19726562,  2.21252441, ...,  1.33056641,
#        1.3458252 ,  1.3671875 ], dtype=float32) * pA, [0.0 s, 191.2832 s], sampling rate: 10000.0 Hz)>]
#>>> print bl.segments[0].eventarrays
#[]

if __name__ == '__main__':
    pass

##f_dir = '/home/maja/PhDProject/SPWs/data/induced/Cell 1/25102011_0023_stim oscillo_ inducing cell.abf'
#f_dir = '/home/maja/PhDProject/SPWs/data/induced/Cell 5/17112011_0002_stim oscillo.abf'
#f_dir = '/home/maja/PhDProject/SPWs/data/induced/Cell 1/25102011_0019_gap free.abf'
#f_dir = '/home/maja/PhDProject/SPWs/data/dat.abf'
#get_data_oscillo(f_dir)
