import induc_SPW as ispw
import numpy as np
#import read_data as reader
import pylab as py
import data_mang as dat
import folder_manager as fold_mng
import gc #garbage collector


def up_intraSpikes(save_folder, save_file = 'intra_spikes.npz', load_file = 'data_intra.npz', reanalize = False):
    """ finds the spikes in the intracellular data"""
    # check if folder already exists
    fold_mng.create_folder(save_folder)
    # check if this file already exists
    exists = fold_mng.file_exists(save_folder, save_file)
    if reanalize or not exists:
        # load the data
        ispw.update_intraSpikes(save_folder, save_file = save_file, load_file = load_file, pulse_len = 500)
    else:
        print 'intracellular spikes were already found previously'    
    gc.collect()    

def up_intrafile(filename, save_folder, save_file = 'data_intra.npz', int_electrodes = [1], reanalize = False):
    """ read intracellular data"""
    # check if folder already exists
    fold_mng.create_folder(save_folder)
    
    # check if this file already exists
    exists = fold_mng.file_exists(save_folder, save_file)
    if reanalize or not exists:
        # load the data
        ispw.update_datafile(filename, int_electrodes, save_folder, data_file = save_file)
    else:
        print 'raw intracellular data was already loaded'    
    gc.collect()
    

def up_datafile(filename, save_folder, save_file = 'data.npz', ext_electrodes = [1], intr_electrode = 1, reanalize = False):
    """ updates only the datafile for the given values """
    # check if folder already exists
    fold_mng.create_folder(save_folder)
    
    # check if this file already exists
    exists = fold_mng.file_exists(save_folder, save_file)
    if reanalize or not exists:
        ispw.update_datafile(filename, ext_electrodes, save_folder, data_file = save_file)
    else:
        print 'raw data file already exists'
    gc.collect()

def up_databas(save_folder, save_file = "data_dspl.npz", load_file = 'data.npz', reanalize = False):
    """ it downsamples the data taken from the given file and saves it in another file"""
    # check if folder already exists
    fold_mng.create_folder(save_folder)
    
    # check if this file already exists
    exists = fold_mng.file_exists(save_folder, save_file)
    if reanalize or not exists:
        # load the data
        ispw.update_databas(data_load = load_file, save_folder = save_folder, data_file = save_file)
    else:
        print 'raw data was already moved to the baseline'    
    gc.collect()

def up_highWaves(save_folder, save_file = "data_movavg.npz", load_datafile = 'data.npz', reanalize = False):
    """ it subtracts moving average from the data"""
    # check if folder already exists
    fold_mng.create_folder(save_folder)
    
    # check if this file already exists
    exists = fold_mng.file_exists(save_folder, save_file)
    if reanalize or not exists:
        # load the data       
        ispw.update_highWaves(load_datafile, save_folder, data_file = save_file, atten_len = 25)
    else:
        print 'raw data was already moved to the baseline' 
    gc.collect()
    
def up_highWaves_numb(save_folder, save_file = 'spws_params.npz', load_spwsfile = 'spws_potential', reanalize = False):
    """ it finds the characteristics of each spw"""
    # check if folder already exists
    fold_mng.create_folder(save_folder)
    
    # check if this file already exists
    exists = fold_mng.file_exists(save_folder, save_file)
    if reanalize or not exists:   
        # load spike params
       
        ispw.update_highWaves_numb(load_spwsfile, save_folder, save_file)
    else:
        print 'spws were already analysed'    
    gc.collect()
        
def up_spws_spikes_ampl(save_folder, save_file = 'data.npz', load_spwsspike = 'SPWs_spikes.npz', load_spikefile = 'spikes_params.npz', reanalize = False):
    """ Finds which of the spikes detected is of the same amplitude as coresponding highest spike - it returns only the highest amplitude spikes"""
    
    fold_mng.create_folder(save_folder)
    
    # check if this file already exists
    exists = fold_mng.file_exists(save_folder, save_file)
    if reanalize or not exists:      
        ispw.update_SPW_spikes_ampl(load_spikefile, load_spwsspike, save_folder, save_file)  
    else:
        print 'spws were already analysed'        
    gc.collect()


def up_spikes_in_spw(save_folder, save_file ='spikes_in_spws', load_spike_file = 'spikes.npz', load_spw_file = 'spws.npz', reanalize = False, win = [-20, 20]):
    """ finds all the spikes in a given SPW"""
    # check if folder already exists
    fold_mng.create_folder(save_folder)
    
    # check if this file already exists
    exists = fold_mng.file_exists(save_folder, save_file)
    if reanalize or not exists:
        # load the data   
        ispw.update_spikes_in_spws(save_folder, save_file, load_spike_file, load_spw_file, win)
    else:
        print 'spikes were already found in the SPW'
    gc.collect()        
    
    

def up_spikes_ampl(save_folder, save_file ='spikes_in_spws', load_spike_file = 'spikes.npz', reanalize = False):
    """ finds all the spikes for the given spw"""
    # check if folder already exists
    fold_mng.create_folder(save_folder)
    
    # check if this file already exists
    exists = fold_mng.file_exists(save_folder, save_file)
    if reanalize or not exists:
        # load the data   
        ispw.update_spikes_ampls(save_folder, save_file, load_spike_file)
    else:
        print 'origins of spikes were already calculated'
    gc.collect()    

def up_correct_ipsps(save_folder, save_fig = 'spw_ipsp', save_file = 'save_it.npz', load_datafile = 'data.npz', load_spwsipsp = 'spws.npz', load_spwsspike = 'spw_spike.npz', reanalize = False, ext = '.pdf'):
    """analyse all the ipsps and correct them"""
    # check if folder already exists
    fold_mng.create_folder(save_folder)
    
    # check if this file already exists
    exists = fold_mng.file_exists(save_folder, save_file)
    if reanalize or not exists:
        fig_fold_name = 'SPW_IPSPs/'
        fold_mng.create_folder(save_folder + fig_fold_name)
          
        ispw.corect_ipsps(load_datafile, load_spwsipsp, load_spwsspike, save_folder, fig_fold_name + save_fig, save_file, ext)  
    gc.collect()    

def up_spws_first_max(save_folder, save_file, spws, datafile, reanalize = False):
    """ alignes spws on the first maximum within the window"""
 # check if folder already exists
    fold_mng.create_folder(save_folder)
    
    # check if this file already exists
    exists = fold_mng.file_exists(save_folder, save_file)
    if reanalize or not exists:
        fig_fold_name = 'SPW_IPSPs/'
        fold_mng.create_folder(save_folder + fig_fold_name)
          
        ispw.update_spws_first_max(save_folder, spws, datafile, save_file, window = [-1, 3])
    gc.collect()

def up_divide_to_groups(load_datafile, load_spwsipsp, save_folder, save_file, reanalize): 
    """analyse the beginning of each SPW - finds the beginnings - time and location"""
    # check if folder already exists
    fold_mng.create_folder(save_folder)
    
    # check if this file already exists
    exists = fold_mng.file_exists(save_folder, save_file)
    if reanalize or not exists:
        fig_fold_name = 'SPW_IPSPs/'
        fold_mng.create_folder(save_folder + fig_fold_name)
          
        ispw.divide_to_groups(load_datafile, load_spwsipsp, save_folder, save_file)  
    gc.collect()
  
def up_spws_beg(save_folder, save_fig = 'spw_ipsp', save_file = 'save_it.npz', load_datafile = 'data.npz', load_spwsipsp = 'spws.npz', load_spwsspike = 'spw_spike.npz', reanalize = False, ext = '.pdf'): 
    """analyse the beginning of each SPW - finds the beginnings - time and location"""
    # check if folder already exists
    fold_mng.create_folder(save_folder)
    
    # check if this file already exists
    exists = fold_mng.file_exists(save_folder, save_file)
    if reanalize or not exists:
        fig_fold_name = 'SPW_IPSPs/'
        fold_mng.create_folder(save_folder + fig_fold_name)
          
        ispw.update_spws_beg(load_datafile, load_spwsipsp, load_spwsspike, save_folder, fig_fold_name + save_fig, save_file, ext)  
    gc.collect()
        
        
def up_spws_ipsp_beg(save_folder, save_fig = 'spw_ipsp', save_file = 'save_it.npz', load_datafile = 'data.npz', load_spwsipsp = 'spws.npz', load_spwsspike = 'spw_spike.npz', reanalize = False, ext = '.pdf'):       
    """analyse the ipsps in each SPWs - finds the beginnings, and removes those which are not correct"""
    # check if folder already exists
    fold_mng.create_folder(save_folder)
    
    # check if this file already exists
    exists = fold_mng.file_exists(save_folder, save_file)
    if reanalize or not exists:
        fig_fold_name = 'SPW_IPSPs/'
        fold_mng.create_folder(save_folder + fig_fold_name)
        # load the data   
        ispw.update_SPW_ipsp_correct(load_datafile, load_spwsipsp, load_spwsspike, save_folder, fig_fold_name + save_fig, save_file, ext)  
    gc.collect()
            
#def up_spws_ipsp_ampl(save_folder, save_file = 'save_it.npz', load_datafile = 'data.npz', load_spwsipsp = 'spws.npz', reanalize = False):       
#    """analyse the ipsps and checks in which electrode it has the highest amplitude"""
#    # check if folder already exists
#    fold_mng.create_folder(save_folder)
#    
#    # check if this file already exists
#    exists = fold_mng.file_exists(save_folder, save_file)
#    if reanalize or not exists:
#        # load the data
#        npzfile         = np.load(save_folder + load_datafile)
#        data            = npzfile['data']
#        fs              = npzfile['fs']
#        npzfile.close()   
#        
#        npzfile         = np.load(save_folder + load_spwsipsp)
#        npzfile.close()
#        import pdb; pdb.set_trace()
#             
#        ispw.update_SPW_ipsp_ampl(save_folder, save_file, data, fs) 
#    gc.collect()
             
def update_ipsp_exSpikes(save_folder, save_file):
    pass

def equalize_number_spws(save_folder, save_file, induc_spont, load_distances, reanalize):
    fold_mng.create_folder(save_folder)
    
    # check if this file already exists
    exists = fold_mng.file_exists(save_folder, save_file)

    if reanalize or not exists:
        # load the data        
        ispw.update_equalize_number_spws(save_folder = save_folder, save_file = save_file, induc_spont = induc_spont, load_distances = load_distances)
    else:
        print 'Initiated SPWs were already saved'    
    gc.collect()      


def up_fill_gap_between_ipsp_groups(save_folder, save_file, spw_file, data_file, reanalize = False):
    fold_mng.create_folder(save_folder)
    
    # check if this file already exists
    exists = fold_mng.file_exists(save_folder, save_file)

    if reanalize or not exists:
        # load the data        
        ispw.update_fill_gap_between_ipsp_groups(save_folder = save_folder, save_file = save_file, spw_file = spw_file, data_file = data_file)
    else:
        print 'Initiated SPWs were already saved'    
    gc.collect()  
    
def up_induc_spont_spw(save_folder, save_file = 'i_s_spws', load_distances = 'distances.npz', load_spwfile = 'spws.npz', max_init_dist = 10, reanalize = False, ext = '.pdf'):
    """ it finds which spws are initiated and which are sponteneaus"""
    fold_mng.create_folder(save_folder)
    
    # check if this file already exists
    exists = fold_mng.file_exists(save_folder, save_file)

    if reanalize or not exists:
        # load the data        
        ispw.update_induc_spont_spw(save_folder = save_folder, save_file = save_file, load_distances = load_distances, load_spwfile = load_spwfile, max_dist =max_init_dist, ext = ext)
    else:
        print 'Initiated SPWs were already saved'    
    gc.collect()       
    
    

def up_dist_SpwfromSpike(save_folder, save_file = 'spw_dist.npz', load_intrafile = 'intra_data.npz', load_spwfile = 'spw_data.npz', spikes = 'all', reanalize = False):
    """ it finds the distance intracellular spike to each spw"""
    # check if folder already exists
    fold_mng.create_folder(save_folder)
    
    # check if this file already exists
    exists = fold_mng.file_exists(save_folder, save_file)

    if reanalize or not exists:
        # load the data        
        ispw.update_dist_SPWfromSpike(save_folder = save_folder, save_file = save_file, load_intrafile = load_intrafile, load_spwfile = load_spwfile, max_dist = 15, spikes = spikes)
    else:
        print 'distances of spws to intracellular spikes were already calculated'    
    gc.collect()    
    
    
def up_SPW_ipsp(save_folder, filter_folder, save_file = 'spws_params.npz', load_datafile = "data_movavg.npz", load_waves = 'spws.npz', load_spikes = 'spws_potential', induced_dist = 7, reanalize = False):
    """ it finds the characteristics of each spw"""
    # check if folder already exists
    fold_mng.create_folder(save_folder)
    
    # check if this file already exists
    exists = fold_mng.file_exists(save_folder, save_file)
    if reanalize or not exists:
        # load the data        
        ispw.update_SPW_ipsp(load_datafile, filter_folder, load_waves, load_spikes, save_folder, save_file)
    else:
        print 'spws were already analysed'    
    gc.collect()


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
        npzfile.close()
        
        data_filt, freq, fs_data = ispw.update_filtered(data, fs, save_folder, freq, data_file = save_file)   
    else:
        print 'raw data was already filtered' 
    gc.collect()
    
def up_extraspikes(save_folder, filter_folder,  save_file = "ex_spikes", load_file = "data_dspl.npz", spikes_filter = 'filter_', reanalize = False):
    """ finding extracellular spikes in the data """
    # check if folder already exists
    fold_mng.create_folder(save_folder)
    
    # check if this file already exists
    exists = fold_mng.file_exists(save_folder, save_file)
    if reanalize or not exists:  
        ispw.update_extraspikes(data_load = load_file, filter_folder = filter_folder, save_folder = save_folder, save_file = save_file, save_filter = spikes_filter)
    else:
        print 'spikes were already found'
    gc.collect()

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
        npzfile.close()
        
        spw_idxs, spw_maxs, starts_spw, ends_spw, lengths_spw, fs_spws = ispw.update_spws(data, fs = fs, save_folder = save_folder, save_file = save_file) 
    else:
        print 'raw data was already filtered' 
    gc.collect()
    
def up_expikes_params(save_folder, save_file = 'spw_data.npz', load_datafile = 'spw_data.npz', load_spikefile = 'spikefile.npz', reanalize = False):
    """ finds different parameters of the spike and returns them in ms"""
    # check if folder already exists
    fold_mng.create_folder(save_folder)
    
    # check if this file already exists
    exists = fold_mng.file_exists(save_folder, save_file)
    if reanalize or not exists:       
        ispw.update_expikes_params(load_datafile, load_spikefile, save_folder, save_file = save_file)  
    else:
        print 'spikes parameters were already calculated'
    gc.collect()






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
    


    
def up_downsample_intra(save_folder, dspl = 2):
    """ updates only downsampling data"""
    data_intra, fs = reader.read_datafile(save_folder, 'intra.npz') # read existing data
    data_dspl_intra, fs_new_down = ispw.update_downsample(data_intra, fs, save_folder, dspl = 2, data_file = 'data_dspl_intra')
    return data_dspl_intra, fs_new_down
        

    
def up_dist_fromSpike(save_folder, intra = 0):
    
    sp_intra_first, sp_intra_all, fs_intra = reader.read_intra_spikes(save_folder, save_file = "intra_spikes.npz")
    spike_extra, spikes_ampl, fs_extra = reader.read_extraspikes(save_folder, save_file = "ex_spikes.npz")
    
    values, fs = update_fs(fs_extra, fs_intra, values = [spike_extra])
    spike_extra = values[0]

    distances, min_distances, fs, max_dist = ispw.update_dist_fromSpike(sp_intra_all, spike_extra, fs_intra, save_folder, data_file = 'dist_fromSpike', max_dist = 3, intra = intra)

    return distances, min_distances, fs, max_dist
        

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