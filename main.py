import scipy.signal as signal
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pylab as py
from scipy import signal
from matplotlib.ticker import NullFormatter
import update_data as updater
import analyser as analyser
import folder_manager as fold_mng
from configure_path import find_folders
from convert_data_bas import rec2array


# if new settings applied - it will rerun everything again
# to select the settings go to induc_SPW
def work_on_all(filename, save_folder, ext_electrodes = [1, 2, 3, 4, 5, 6, 7], intr_electrode = 1, data_part = 'all'):
    """ work on given data file - proceed all the variables and save them """
    ext = '.pdf'
    
    print 'working on: ' +  filename
    reanalize = True # set to true if to analyze the data no matter if it was already analysed or not
    
    delete_old = False #!!!!
    
    run_all_functions = False
    if delete_old:
        fold_mng.remove_folder(save_folder)
    win = [0, 0]
    
    
    raw_data        = 'data.npz'
    if run_all_functions:
        updater.up_datafile(filename, save_folder = save_folder, save_file = raw_data, ext_electrodes = ext_electrodes, intr_electrode = 1, reanalize = reanalize)
    
    plot_folder = 'full_data/'
    save_plots = 'data'
    if run_all_functions:
        analyser.display_data(save_folder, plot_folder, save_plots, raw_data, trace = 0, part = [0, 10000], ext = ext)
    
    raw_baselined   = "data_bas.npz"
    if run_all_functions:
        updater.up_databas(save_folder, save_file = raw_baselined, load_file = raw_data, reanalize = reanalize)
    
    spikes_raw      = 'spikes.npz'
    spikes_filter = 'fast_data_'
    filter_folder = 'filtered/'
    if run_all_functions:
        updater.up_extraspikes(save_folder, filter_folder, save_file = spikes_raw, load_file = raw_baselined, spikes_filter =  spikes_filter, reanalize = reanalize)
    
    spikes_params   = 'spikes_params.npz'
    if run_all_functions:
        updater.up_expikes_params(save_folder, save_file = spikes_params, load_datafile = raw_baselined, load_spikefile = spikes_raw, reanalize = reanalize)
    
    spikes_largest = 'spikes_largest.npz'
    if run_all_functions:
        updater.up_spikes_ampl(save_folder, save_file =spikes_largest, load_spike_file = spikes_raw, reanalize = reanalize)
     
    SPWs_potential  = 'potential_SPWs.npz'
    if run_all_functions:
        updater.up_highWaves(save_folder, save_file = SPWs_potential, load_datafile = raw_baselined,reanalize = reanalize)
    
    SPWs_spikes     = 'spws_spikes.npz'
    if run_all_functions:
        #does not use previously selected largest spikes
        updater.up_spws_spikes(save_folder, save_file = SPWs_spikes, load_spwsfile = SPWs_potential, load_spikefile = spikes_params, reanalize = reanalize)
    
    SPWs_spikes_ampl= 'spw_spikes_ampl.npz'
    if run_all_functions:
        updater.up_spws_spikes_ampl(save_folder, save_file = SPWs_spikes_ampl, load_spwsspike = SPWs_spikes, load_spikefile = spikes_params, reanalize = reanalize)
    
    SPWs_ipsps      = 'spws_params.npz' #'spws_ipsps.npz'
    if run_all_functions:
        updater.up_SPW_ipsp(save_folder, save_file = SPWs_ipsps, load_datafile = raw_baselined, load_spwsspike = SPWs_spikes, reanalize = reanalize)
    
    SPWs_ipsps_beg  = 'spw_ipsps_beg.npz'
    if run_all_functions:
        updater.up_spws_ipsp_beg(save_folder,  save_fig = 'spw_ipsp', save_file = SPWs_ipsps_beg, load_datafile = raw_baselined, load_spwsipsp = SPWs_ipsps, load_spwsspike = SPWs_spikes_ampl, reanalize = reanalize, ext = ext)
    
    spikes_inSPWs = 'spikes_inSpw.npz'
    if run_all_functions:
        # uses previously selected largest spikes
        updater.up_spikes_in_spw(save_folder, save_file =spikes_inSPWs, load_spike_file = spikes_largest, load_spw_file = SPWs_ipsps_beg, reanalize = reanalize, win = win)
    
    #print intr_electrode
    if intr_electrode == 1:
        
        data_intra = 'data_intra.npz'
        if run_all_functions:
            updater.up_intrafile(filename, save_folder, save_file = data_intra, int_electrodes = [0], reanalize = reanalize)
        
        plot_folder = 'full_data/'
        save_plots = 'data_intra'
        if run_all_functions:
            analyser.display_data(save_folder, plot_folder, save_plots, data_intra, trace = 0, part = [0, 100000], ext = ext)
        
        intra_spikes = 'intra_spikes.npz'
        if run_all_functions:
            updater.up_intraSpikes(save_folder, save_file = intra_spikes, load_file = data_intra, reanalize = reanalize)
        
        dist_spw_inspikes = 'spw_dist2first.npz'
        if run_all_functions:
            updater.up_dist_SpwfromSpike(save_folder, save_file = dist_spw_inspikes, load_intrafile = intra_spikes, load_spwfile = SPWs_ipsps_beg, spikes = 'first', reanalize = reanalize)
        
        dist_spw_inspikes2all = 'spw_dist2all.npz'
        if run_all_functions:
            updater.up_dist_SpwfromSpike(save_folder, save_file = dist_spw_inspikes2all, load_intrafile = intra_spikes, load_spwfile = SPWs_ipsps_beg, spikes = 'all', reanalize = reanalize)
        
        induc_spont_spw = 'induc_spont_spw.npz'
        max_dist = 10 # ms
        if run_all_functions:
            updater.up_induc_spont_spw(save_folder, save_file = induc_spont_spw, load_distances = dist_spw_inspikes, load_spwfile = SPWs_ipsps_beg, max_init_dist = max_dist, reanalize = reanalize, ext = ext)
    
    ##### - analyser - #####
        solutions_folder = 'plots/'
        
        
        numIpsp2distance = 'numIPSP_distance'
        if run_all_functions:
            analyser.plot_noIpsps2distance(save_folder, solutions_folder+numIpsp2distance + '/', save_plots = numIpsp2distance, spw_file = SPWs_ipsps_beg, dist_file = dist_spw_inspikes, ext = ext)
        
        dist_spw2psike = 'dist_spw2spike'
        if run_all_functions:
            analyser.plot_dist_spw2spike(save_folder, solutions_folder+dist_spw2psike + '/', save_plots = dist_spw2psike, dist_file = dist_spw_inspikes, ext = ext)
        
        alignedSPWs = 'aligned_SPWs'
        if run_all_functions:

            analyser.plot_alignedSPW(save_folder, solutions_folder+alignedSPWs + '/', save_plots = alignedSPWs, data_file = raw_baselined, intra_data_file = data_intra, spike_file = induc_spont_spw, intra_spikes = intra_spikes, ext = ext)
        
#        alignedSPWs_2all = 'aligned_SPWs2allSpikes'
#        if run_all_functions:
#            analyser.plot_alignedSPW(save_folder, solutions_folder+alignedSPWs + '/', save_plots = alignedSPWs_2all, data_file = raw_baselined, intra_data_file = data_intra, spw_file = induc_spont_spw, dist_file = dist_spw_inspikes2all, ext = ext)
        
        spikes_inSPWs_plot = 'spikes_inSPWs'
        if run_all_functions:
            analyser.plot_spikes4spw(save_folder, solutions_folder+spikes_inSPWs_plot + '/', 
                                 save_plots = spikes_inSPWs_plot, data_file = raw_baselined, 
                                 spike_data = spikes_inSPWs, spw_data = induc_spont_spw, 
                                 spikes_filter = [], ext = ext, win = win)
        
    
        spikes_inSPWs_plot_fig3a = 'spikes_inSPWs_fig3a'
        if run_all_functions:
            analyser.plot_spikes4spw(save_folder, solutions_folder+spikes_inSPWs_plot_fig3a + '/', 
                                 save_plots = spikes_inSPWs_plot_fig3a, data_file = raw_baselined, 
                                 spike_data = spikes_inSPWs, spw_data = induc_spont_spw, 
                                 spikes_filter = filter_folder + spikes_filter, ext = ext, win = win)
        
        
        spikePerElectrode = 'spike_per_electrode'
        if not run_all_functions:
            analyser.plot_spike(save_folder, solutions_folder + spikePerElectrode + '/', save_plots = spikePerElectrode, 
                            spike_data = spikes_inSPWs, spw_data = induc_spont_spw, 
                            ext = ext, win = win)
    
    
    #ipsp_exSpikes = 'ipsp_exSpikes.npz'
    #updater.update_ipsp_exSpikes(save_folder, save_file = ipsp_exSpikes)
#    def plot_extra_spike_distribut(save_folder, solutions_fold+spikes_inSPWs + '/', save_plots = spikes_inSPWs, spike_data = spikes_params, spw_data = SPWs_ipsps_beg):
#    """ it takes the spikes (max amplitude) for each SPW """

    


#    clean_filter = [-1, 2] # filter the data for using for further analysis
#    raw_clean_data = 'clean_data.npz'
#    updater.up_filtered(save_folder, save_file = raw_clean_data, load_file = raw_baselined, freq = clean_filter, reanalize = reanalize)
#    
#    spw_data = "spws.npz"
#    updater.up_spws(save_folder, save_file = spw_data, load_file = raw_clean_data, reanalize = reanalize)
    
    
    
    
#    updater.up_downsample(save_folder)
    
#    SPW_freq = [1.5, 500.0]
#    ripple_freq = [100.0, 300.0]
    fast_freq = [500.0, -1]
    smooth_freq = [500, 2000]
#    updater.up_filtered(save_folder, SPW_freq, save_file = 'spw_data')
#    updater.up_filtered(save_folder, ripple_freq, save_file = "ripple_data")
    ##updater.up_filtered(save_folder, fast_freq, save_file = 'fast_data')
    #updater.up_filtered(save_folder, smooth_freq, save_file = 'smooth_data')
    
    #updater.up_extraspikes(save_folder)
    #updater.up_expikes_params(save_folder)
##    
#    updater.up_spws(save_folder)
#    updater.up_ipsp(save_folder)
#    updater.up_ripples(save_folder)
#    updater.up_spw_ripple(save_folder)

    #updater.up_dist_fromSPW(save_folder, intra = intr_electrode)
#    
##    
#    if intr_electrode:
#        updater.up_intrafile(filename, save_folder)
        #updater.up_intraSpikes(save_folder)
        #updater.up_dist_fromSpike(save_folder, intr_electrode)
        #updater.up_dist_SpwfromSpike(save_folder)
        
#    if intr_electrode == 1:
#        # intracellular electrode was recorded and should be processed
#        updater.up_intrafile(filename, save_folder, data_part = data_part)
#        updater.up_downsample(save_folder, dspl = 2)
    
def analyse_all(filename, save_folder, intra = 1): 
    ext = '.pdf'  
    #ext = '.png'  
    print save_folder 
#    analyser.plot_data_all(save_folder)
#    analyser.plot_data_bas(save_folder)
#    analyser.plot_spw_data(save_folder)
#    analyser.plot_ripple_data(save_folder)
#    analyser.plot_fast_data(save_folder)
    #analyser.correlate_spw(save_folder)
    #analyser.plot_extraspikes(save_folder)
    #analyser.work_on_swp_spikes(save_folder)
    
    
    #analyser.plot_mean_waves_elecrwise_from_saved(save_folder, save_ext = ext)
    #analyser.calculate_corr_matrix_electrowise(save_folder, save_ext = ext, spws = 1, win = (-10, 30), save_name = 'corr_electr_spws') 
    #analyser.cluster_corrs(save_folder, read_corr = 'corr_electr_spws.npz', save_clust = 'clust_electrSPW.npz', n_clusters = 3)
    #if intra == 1:
        #analyser.plot_corr_with_all(save_folder, ext, read_clust = 'clust_electrSPW.npz')
        #analyser.group_mean(save_folder, ext, read_clust = 'clust_electrSPW.npz', read_idx = 'corr_electr_spws.npz')
        #analyser.plot_data_spikes_spws(save_folder, ext)
    #analyser.plot_extra_params(save_folder, ext)
    
    if intra == 1:
        #analyser.plot_extra_params_novel(save_folder, save_name = 'chosen_spikes.npz', ext = ext) # <-- popraw!
        #analyser.detect_cell4spikes(save_folder, save_name = 'chosen_spikes.npz')
        #analyser.save_chosen(save_folder, save_name = 'chosen_spikes.npz', save_ext = '.png')
        #analyser.max_between_spikes(save_folder,  load_chosen  = 'chosen_idx_val.npz', save_name = 'spike_maxs.npz', save_ext = ext)
        
        analyser.spw_spike_plotter(save_folder, save_name = 'main_spw', save_ext = ext)
        #analyser.find_ipsps(save_folder, load_file = 'IPSPs.npz', save_name = 'isps_clean', save_ext = ext)
        #analyser.plot_distSPWfromSpike(save_folder, save_ext = ext)
   
    #analyser.exspikes_clust_separate(save_folder, load_chosen = 'chosen_idx_val.npz', save_clust = 'exclust_sep', save_ext = ext)
    #analyser.plot_clust(save_folder, load_clust = 'exclust_sep.npz', save_ext = ext)
    #if intra == 1:
        #analyser.find_clust_extraspikes(save_folder, n_clusters = 3, save_name = 'extra_clusters', save_clust = 'extra_clusters', save_ext = ext)
        #analyser.plot_extra_clusters(save_folder, save_name = 'extra_spikes_clust', save_clust = 'extra_clusters.npz', save_ext = ext, win = [-1, 2])
        #analyser.plot_hist_dist_ex2intra(save_folder, ext, read_clust = 'extra_clusters.npz')
        #analyser.plot_hist_dist_ex2SPW(save_folder, ext, read_clust = 'extra_clusters.npz')


#
#f_name = ''

  
#save_folder = '/home/maja/PhDProject/SPWs/SPWs/saved_data/cell5/part1_test/'
#save_folder = '/home/maja/PhDProject/SPWs/SPWs/saved_data/cell6/test'


if __name__=='__main__':
    all = ((1, 1, 0), (1, 1, 1), (1, 2, 1), (3, 1, 0), (3, 1, 1), (4, 1, 1),
           (5, 1, 0), (5, 2, 0), (5, 3, 0), (5, 4, 0), (5, 5, 0), (5, 1, 1),
           (6, 1, 0), (6, 1, 1), (7, 1, 0), (8, 1, 1), (9, 1, 1), (10, 1, 1), 
           (11, 1, 1))
    
    update = 1
    analyse = 0
    
    if update == 1:
        #for nex in [15]:
        #for nex in range(len(all)):
        for nex in [15, 17]: #range(1, 15):
            filename, save_folder, intra  = find_folders(all[nex][0], all[nex][1], all[nex][2])
            
            ex_electr = range(intra, 7+intra)
            print 'intra ' + str(intra)
            work_on_all(filename, save_folder, ex_electr, intra)
    
    if analyse == 1:
        #for nex in range(12, len(all)):
        #for nex in [5]: #range(12, len(all)):
        #for nex in [11]:
        for nex in range(15):
            filename, save_folder, intra  = find_folders(all[nex][0], all[nex][1], all[nex][2])
            analyse_all(filename, save_folder, intra)
        

    
