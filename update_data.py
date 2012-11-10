import induc_SPW as ispw
import numpy as np
#import read_data as reader
import pylab as py
import data_mang as dat
import folder_manager as fold_mng

def up_datafile(filename, save_folder, save_file = 'data.npz', ext_electrodes = [1], intr_electrode = 1, data_part = 'all', reanalize = False):
    """ updates only the datafile for the given values """
    # check if folder already exists
    fold_mng.create_folder(save_folder)
    
    # check if this file already exists
    exists = fold_mng.file_exists(save_folder, save_file)
    if reanalize or not exists:
        data_all, fs = ispw.update_datafile(filename, ext_electrodes, save_folder, data_file = save_file, data_part = data_part)
    else:
        print 'raw data file already exists'


def up_databas(save_folder, save_file = "data_dspl.npz", load_file = 'data.npz', reanalize = False):
    """ it downsamples the data taken from the given file and saves it in another file"""
    # check if folder already exists
    fold_mng.create_folder(save_folder)
    
    # check if this file already exists
    exists = fold_mng.file_exists(save_folder, save_file)
    if reanalize or not exists:
        # load the data
        npzfile = np.load(save_folder + load_file)
        data = npzfile['data']
        fs = npzfile['fs']
        npzfile.close()
        
        ispw.update_databas(data, fs, save_folder, data_file = save_file)
    else:
        print 'raw data was already moved to the baseline'    

def up_highWaves(save_folder, save_file = "data_movavg.npz", load_datafile = 'data.npz', reanalize = False):
    """ it subtracts moving average from the data"""
    # check if folder already exists
    fold_mng.create_folder(save_folder)
    
    # check if this file already exists
    exists = fold_mng.file_exists(save_folder, save_file)
    if reanalize or not exists:
        # load the data
        npzfile = np.load(save_folder + load_datafile)
        data = npzfile['data']
        fs = npzfile['fs']
        npzfile.close()
        
        ispw.update_highWaves(data, fs, save_folder, data_file = save_file, atten_len = 25)
    else:
        print 'raw data was already moved to the baseline' 

def up_spws_spikes(save_folder, save_file = 'spws_params.npz', load_spwsfile = 'spws_potential', load_spikefile = 'spikes_params.npz', reanalize = False):
    """ it finds the characteristics of each spw"""
    # check if folder already exists
    fold_mng.create_folder(save_folder)
    
    # check if this file already exists
    exists = fold_mng.file_exists(save_folder, save_file)
    if reanalize or not exists:   
        # load spike params
        npzfile         = np.load(save_folder + load_spikefile)
        spike_idxs      = npzfile['spike_idxs']
        fs              = npzfile['fs']
        npzfile.close()
        
        # load starts of spws
        npzfile         = np.load(save_folder + load_spwsfile)
        spw_starts      = npzfile['starts']
        spw_ends        = npzfile['ends']
        npzfile.close()
        
        ispw.update_SPWspikes(fs, save_folder, save_file, spw_starts, spw_ends, spike_idxs)
    else:
        print 'spws were already analysed'    

def up_spws_spikes_ampl(save_folder, save_file = 'data.npz', load_spwsspike = 'SPWs_spikes.npz', load_spikefile = 'spikes_params.npz', reanalize = False):
    """ Finds which of the spikes detected is of the same amplitude as coresponding highest spike - it returns only the highest amplitude spikes"""
    
    fold_mng.create_folder(save_folder)
    
    # check if this file already exists
    exists = fold_mng.file_exists(save_folder, save_file)
    if reanalize or not exists:
        # load spike params
        npzfile         = np.load(save_folder + load_spikefile)
        spike_idxs      = npzfile['spike_idxs']
        ampls           = npzfile['ampls']
        fs              = npzfile['fs']
        npzfile.close()

        # load starts of spws
        npzfile         = np.load(save_folder + load_spwsspike)
        spw_spikes      = npzfile['spw_details']  
        npzfile.close()
        
        ispw.update_SPW_spikes_ampl(save_folder, save_file, spw_spikes, spike_idxs, ampls, fs)  
    else:
        print 'spws were already analysed'        
    
def up_spws_ipsp_beg(save_folder, save_file = 'save_it.npz', load_datafile = 'data.npz', load_spwsipsp = 'spws.npz', load_spwsspike = 'spw_spike.npz', reanalize = False):       
    """analyse the ipsps in each SPWs - finds the beginnings, and removes those which are not correct"""
    # check if folder already exists
    fold_mng.create_folder(save_folder)
    
    # check if this file already exists
    exists = fold_mng.file_exists(save_folder, save_file)
    if reanalize or not exists:
        # load the data
        npzfile         = np.load(save_folder + load_datafile)
        data            = npzfile['data']
        fs              = npzfile['fs']
        npzfile.close()   
        
        npzfile         = np.load(save_folder + load_spwsipsp)
        spw_ipsps       = npzfile['spw_ipsps']
        npzfile.close()   
        
        npzfile         = np.load(save_folder + load_spwsspike)
        spw_spikes      = npzfile['chosen_spikes']  
                
        ispw.update_SPW_ipsp_correct(save_folder, save_file, data, fs, spw_ipsps, spw_spikes)  
        
def up_spws_ipsp_ampl(save_folder, save_file = 'save_it.npz', load_datafile = 'data.npz', load_spwsipsp = 'spws.npz', reanalize = False):       
    """analyse the ipsps and checks in which electrode it has the highest amplitude"""
    # check if folder already exists
    fold_mng.create_folder(save_folder)
    
    # check if this file already exists
    exists = fold_mng.file_exists(save_folder, save_file)
    if reanalize or not exists:
        # load the data
        npzfile         = np.load(save_folder + load_datafile)
        data            = npzfile['data']
        fs              = npzfile['fs']
        npzfile.close()   
        
        npzfile         = np.load(save_folder + load_spwsipsp)
        
        import pdb; pdb.set_trace()
             
        ispw.update_SPW_ipsp_ampl(save_folder, save_file, data, fs) 
         
def update_ipsp_exSpikes(save_folder, save_file):
    pass

        
def up_SPW_ipsp(save_folder, save_file = 'spws_params.npz', load_datafile = "data_movavg.npz", load_spwsspike = 'spws_potential', reanalize = False):
    """ it finds the characteristics of each spw"""
    # check if folder already exists
    fold_mng.create_folder(save_folder)
    
    # check if this file already exists
    exists = fold_mng.file_exists(save_folder, save_file)
    if reanalize or not exists:
        # load the data
        npzfile         = np.load(save_folder + load_datafile)
        data            = npzfile['data']
        fs              = npzfile['fs']
        npzfile.close()
        
        # load starts of spws
        npzfile         = np.load(save_folder + load_spwsspike)
        spw_spikes      = npzfile['spw_details']

        npzfile.close()
        
        ispw.update_SPW_ipsp(data, fs, save_folder, save_file, spw_spikes)
    else:
        print 'spws were already analysed'    



def up_filtered(save_folder, save_file = 'spw_data.npz', load_file = "data_dspl.npz", freq = [1.5, 500.0], reanalize = False):
    """ filteres the data from the given file to given frequencies"""
    # check if folder already exists
    fold_mng.create_folder(save_folder)
    
    # check if this file already exists
    exists = fold_mng.file_exists(save_folder, save_file)
    if reanalize or not exists:
        npzfile = np.load(save_folder + load_file)
        data = npzfile['data']
        fs = npzfile['fs']        

        data_filt, freq, fs_data = ispw.update_filtered(data, fs, save_folder, freq, data_file = save_file)   
    else:
        print 'raw data was already filtered' 

def up_extraspikes(save_folder, save_file = "ex_spikes", load_file = "data_dspl.npz", reanalize = False):
    """ finding extracellular spikes in the data """
    # check if folder already exists
    fold_mng.create_folder(save_folder)
    
    # check if this file already exists
    exists = fold_mng.file_exists(save_folder, save_file)
    if reanalize or not exists:
        npzfile = np.load(save_folder + load_file)
        data = npzfile['data']
        fs = npzfile['fs'] 
        spike_idxs, spikes_ampl, fs_espikes = ispw.update_extraspikes(data, fs, save_folder, save_file = save_file)
    else:
        print 'spikes were already found'
    return spike_idxs, spikes_ampl, fs_espikes

def up_spws(save_folder, save_file = 'spw_data.npz', load_file = 'spw_data.npz', reanalize = False):
    """ updates details of the spws"""
    # check if folder already exists
    fold_mng.create_folder(save_folder)
    
    # check if this file already exists
    exists = fold_mng.file_exists(save_folder, save_file)
    if reanalize or not exists:
        npzfile = np.load(save_folder + load_file)
        data = npzfile['data']
        fs = npzfile['fs']        

        spw_idxs, spw_maxs, starts_spw, ends_spw, lengths_spw, fs_spws = ispw.update_spws(data, fs = fs, save_folder = save_folder, save_file = save_file) 
    else:
        print 'raw data was already filtered' 

def up_expikes_params(save_folder, save_file = 'spw_data.npz', load_datafile = 'spw_data.npz', load_spikefile = 'spikefile.npz', reanalize = False):
    """ finds different parameters of the spike and returns them in ms"""
    # check if folder already exists
    fold_mng.create_folder(save_folder)
    
    # check if this file already exists
    exists = fold_mng.file_exists(save_folder, save_file)
    if reanalize or not exists:
        npzfile = np.load(save_folder + load_datafile)
        data = npzfile['data']
        fs = npzfile['fs'] 
        
        npzfile = np.load(save_folder + load_spikefile)
        spike_idx = npzfile['spike_idx']
        
        ispw.update_expikes_params(data, fs, save_folder, spike_idx, save_file = save_file)  
    else:
        print 'spikes parameters were already calculated'







def up_downsample(save_folder, dspl = 2, data_file = 'data_dspl'):
    """ updates only downsampling data"""
    data_all, fs = reader.read_datafile(save_folder) # read existing data
    data_dspl, fs_dspl = ispw.update_downsample(data_all, fs, save_folder, dspl, data_file = 'data_dspl')
    return data_dspl, fs_dspl









def up_ipsp(save_folder):
    """ updates details of the spws"""
    data_bas, fs = reader.read_databas(save_folder)
    spw_idxs, spw_maxs, starts_spw, ends_spw, lengths_spw, fs_spws = ispw.update_spws(data_bas, fs = fs, save_folder = save_folder, save_file = 'IPSPs', thresh_mult= 1.5)
    
    #spw_idxs, spw_maxs, starts_spw, ends_spw, lengths_spw, fs = reader.read_SPWs(save_folder)
    #print spw_idxs
    return spw_idxs, spw_maxs, starts_spw, ends_spw, lengths_spw, fs_spws




def up_ripples(save_folder): 
    data_ripple, freq, fs = reader.read_filtered_ripple(save_folder, save_file = "ripple_data.npz")
    rip_idxs, fs_ripple = ispw.update_ripples(data_ripple, fs, save_folder, data_file = 'ripples')
    return rip_idxs, fs_ripple

def up_spw_ripple(save_folder):
    spw_idxs, spw_maxs, starts_spw, ends_spw, lengths_spw, fs = reader.read_SPWs(save_folder)
    ris_all, iris_no_all = ispw.update_spw_ripple(starts_spw, ends_spw, spw_idxs, save_folder, data_file = 'spw_ripple')
    return ris_all, iris_no_all

def update_fs(fs_1, fs_2, values = []):
    # all the values will be updated to the fs_2
    if fs_1 != fs_2:
        fs_div = fs_2/ fs_1  # if the sampling rate is different  
        # the sampling rates in the two data sets are different        
        
        for val in range(len(values)):
            #print np.size(spw_idxs)
            for electr in range(len(values[val])):
                for trace in range(len(values[val][electr])):
                    for individual in range(len(values[val][electr][trace])):
                        values[val][electr][trace][individual] = values[val][electr][trace][individual] * fs_div
            values[val] = ispw.round_spike_idxs(values[val])
        fs = fs_1 * fs_div               
    else:
        fs = fs_2
    return values, fs
        

def up_dist_fromSPW(save_folder, intra = 0):
    spw_idxs, spw_maxs, starts_spw, ends_spw, lengths_spw, fs_spw = reader.read_SPWs(save_folder)
    spike_idxs, spikes_ampl, fs_spike  = reader.read_extraspikes(save_folder)

    values, fs = update_fs(fs_spw, fs_spike, values = [spw_idxs, starts_spw, ends_spw])
    spw_idxs = values[0]
    starts_spw = values[1]
    ends_spw = values[2]
    distances, min_distances, fs, max_dist = ispw.update_dist_fromSpike(starts_spw, spike_idxs, fs, save_folder, data_file = 'dist_fromSPW', allowms = 5, intra = 0)
    

def up_intrafile(filename, save_folder, data_part = 'all'):
    """ read intracellular data"""
    data_intra, fs = ispw.update_datafile(filename, [1], save_folder, data_file = 'intra', data_part = data_part)
    return data_intra, fs
    
def up_downsample_intra(save_folder, dspl = 2):
    """ updates only downsampling data"""
    data_intra, fs = reader.read_datafile(save_folder, 'intra.npz') # read existing data
    data_dspl_intra, fs_new_down = ispw.update_downsample(data_intra, fs, save_folder, dspl = 2, data_file = 'data_dspl_intra')
    return data_dspl_intra, fs_new_down
        
def up_intraSpikes(save_folder):
    data_intra, fs = reader.read_datafile(save_folder, save_file = "intra.npz")
    sp_idx_first, sp_idx_all, fs = ispw.update_intraSpikes(data_intra[0], fs, save_folder, save_file = "intra_spikes", pulse_len = 500)
    return sp_idx_first, sp_idx_all, fs
    
def up_dist_fromSpike(save_folder, intra = 0):
    
    sp_intra_first, sp_intra_all, fs_intra = reader.read_intra_spikes(save_folder, save_file = "intra_spikes.npz")
    spike_extra, spikes_ampl, fs_extra = reader.read_extraspikes(save_folder, save_file = "ex_spikes.npz")
    
    values, fs = update_fs(fs_extra, fs_intra, values = [spike_extra])
    spike_extra = values[0]

    distances, min_distances, fs, max_dist = ispw.update_dist_fromSpike(sp_intra_all, spike_extra, fs_intra, save_folder, data_file = 'dist_fromSpike', max_dist = 3, intra = intra)

    return distances, min_distances, fs, max_dist
        
def up_dist_SpwfromSpike(save_folder):
    sp_intra_first, sp_intra_all, fs_intra = reader.read_intra_spikes(save_folder)
    spw_idxs, spw_maxs, starts_spw, ends_spw, lengths_spw, fs_spw = reader.read_SPWs(save_folder)
    print 'spws'
    print len(spw_idxs[0][0])
    # 

    data_intra, fs = reader.read_datafile_intra(save_folder, save_file = "intra.npz")
#    print sp_intra_first
#    print sp_intra_all

    values, fs = update_fs(fs_spw, fs_intra, values = [starts_spw])
    starts_spw = values[0]

    #py.plot(t[sp_intra_first[0]], data_intra[0][0][sp_intra_first[0]], 'go')

#    t = dat.get_timeline(data_intra[0][0], fs, 'ms')
#    py.plot(t, data_intra[0][0])
#
#    py.plot(t[sp_intra_first[0]], data_intra[0][0][sp_intra_first[0]], 'go')
#    py.figure()
#   
#    #
#    py.show()
    # check if fs are the same!!!
    distances, min_distances, fs, max_dist = ispw.update_dist_fromSpike(sp_intra_first, starts_spw, fs, save_folder, data_file = 'dist_SpwfromSpike', max_dist = 10, intra = 1)
    
def update_all(filename, ext_electrodes, save_folder, intr_electrode = 1, data_part = 'all'):
    """ should be done each time the data is to be run from the beginning to the end again"""
    # read all the data
    data_all, fs = ispw.update_datafile(filename, ext_electrodes, save_folder, data_file = 'data', data_part = data_part)
    data_dspl, fs_dspl = ispw.update_downsample(data_all, fs, save_folder, data_file = 'data_dspl')
    data_bas, fs_bas = ispw.update_databas(data_dspl, fs_dspl, save_folder, data_file = 'data_bas') # move to baseline
    
    # update filtered data
    SPW_freq = [1.5, 500.0]
    ripple_freq = [100.0, 300.0]
    fast_freq = [750.0, -1]
    data_spw, freq, fs_data = ispw.update_filtered(data_bas, fs_bas, save_folder, SPW_freq, 'spw_data')
    data_ripple, freq, fs_ripple = ispw.update_filtered(data_bas, fs_bas, save_folder, ripple_freq, "ripple_data")
    data_fast, freq, fs_fast = ispw.update_filtered(data_bas, fs_bas, save_folder, fast_freq, "fast_data")
    #uspl = 1
    #data_uspl_extra, fs_new_up = update_upsample(data_all, fs, save_folder, uspl = uspl, data_file = 'data_uspl_extra')
    spike_idxs, spikes_ampl, fs_espikes = ispw.update_extraspikes(data_all, fs, save_folder, save_file = "ex_spikes")
    spike_idxs, all_valley_to_peak, all_half_valley_width, all_half_peak_width, fs_new = ispw.update_expikes_params(data_all, fs_espikes, save_folder, spike_idxs, save_file = "ex_sparamas")
    spw_idxs, spw_maxs, starts_spw, ends_spw, lengths_spw, fs_spws = ispw.update_spws(data_bas, data_fast, data_spw, fs_bas, save_folder, save_file = 'SPWs')
    rip_idxs, fs_ripple = ispw.update_ripples(data_ripple, fs_ripple, save_folder, data_file = 'ripples')
    ris_all, iris_no_all = ispw.update_spw_ripple(starts_spw, ends_spw, spw_idxs, save_folder, data_file = 'spw_ripple')
    
    if fs_spws != fs_espikes:
        fs_div = fs_espikes/ fs_spws  # if the sampling rate is different  
        print np.shape(spw_idxs[0])
        for i in range(len(spw_idxs)):
            for g in range(len(spw_idxs[i])):
                spw_idxs[i][g] = spw_idxs[i][g] * fs_div
                starts_spw[i][g] = starts_spw[i][g] * fs_div
                ends_spw[i][g] = ends_spw[i][g] * fs_div
   
        spw_idxs = ispw.round_spike_idxs(spw_idxs)
        starts_spw = ispw.round_spike_idxs(starts_spw)
        ends_spw = ispw.round_spike_idxs(ends_spw)
        fs_spws = fs_spws * fs_div
        
    distances, min_distances, fs, max_dist = ispw.update_dist_fromSpike(starts_spw, spike_idxs, fs_spws, save_folder, data_file = 'dist_fromSPW')

    if intr_electrode != -1:
        data_intra, fs = ispw.update_datafile(filename, [intr_electrode], save_folder, data_file = 'intra', data_part = data_part)
        data_dspl_intra, fs_new_down = ispw.update_downsample(data_intra, fs, save_folder, dspl = 2, data_file = 'data_dspl_intra')
        sp_idx_first, sp_idx_all, fs = ispw.update_intraSpikes(data_intra[0], fs, save_folder, save_file = "intra_spikes", pulse_len = 500)
        
        # compare extra and intra data
        
        distances, min_distances, fs, max_dist = ispw.update_dist_fromSpike(sp_idx_all, spike_idxs, fs, save_folder, data_file = 'dist_fromSpike', max_dist = 3)
        distances, min_distances, fs, max_dist = ispw.update_dist_fromSpike(sp_idx_first, starts_spw, fs, save_folder, data_file = 'dist_SpwfromSpike', max_dist = 10)