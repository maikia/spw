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
from configure_path import get_save_folder
from convert_data_bas import rec2array
import logging

# if new settings applied - it will rerun everything again
# to select the settings go to induc_SPW
def work_on_all(filename, save_folder, ext_electrodes = [1, 2, 3, 4, 5, 6, 7], intr_electrode = 1, data_part = 'all'):
    """ work on given data file - proceed all the variables and save them """
    ext = '.png'
    
    print 'working on: ' +  filename
    reanalize = True # set to true if to analyze the data no matter if it was already analysed or not
    
    delete_old = False #!!!!
    
    run_all_functions = True
    
    if delete_old: # !!!!!!
        fold_mng.remove_folder(save_folder)
    win = [0, 0]
    
    
    raw_data        = 'data.npz'
    if not run_all_functions:
        # reads the raw data from the file
        updater.up_datafile(filename, save_folder = save_folder, save_file = raw_data, ext_electrodes = ext_electrodes, intr_electrode = 1, reanalize = reanalize)

# - slvaro
#    s_f = '/home/maja/PhDProject/dokumenty/Alvaro/'
#    s_file = 'data.npz'
#    if not run_all_functions:
#        # reads the raw data from the file
#        updater.up_datafile(filename, save_folder = s_f, save_file = s_file, ext_electrodes = [3], intr_electrode = 0, reanalize = reanalize)    
#    if run_all_functions:
#        analyser.display_data_no_electrode(s_f, '', 'trace1', s_file, part = 'all', ext = ext)
    
    plot_folder = 'full_data/'
    save_plots = 'data'
    if not run_all_functions:
        # saves part of the data 
        analyser.display_data(save_folder, plot_folder, save_plots, raw_data, trace = 0, part = [0, 10000], ext = ext)
    
    raw_baselined   = "data_bas.npz"
    if not run_all_functions:
        # removes the baseline
        updater.up_databas(save_folder, save_file = raw_baselined, load_file = raw_data, reanalize = reanalize)
    
    spikes_raw      = 'spikes.npz'
    spikes_filter = 'fast_data_'
    filter_folder = 'filtered/'
    if not run_all_functions:
        # finds all the extracellular spikes
        updater.up_extraspikes(save_folder, filter_folder, save_file = spikes_raw, load_file = raw_baselined, spikes_filter =  spikes_filter, reanalize = reanalize)
    
    # this funtion might be used for future spike clustering, but shall not be used for other reasons so far
    #spikes_params   = 'spikes_params.npz'
    #if run_all_functions:
    #    updater.up_expikes_params(save_folder, save_file = spikes_params, load_datafile = raw_baselined, load_spikefile = spikes_raw, reanalize = reanalize)
    
    spikes_largest = 'spikes_largest.npz'
    if not run_all_functions:
        # checks which spikes have the highest amplitude (if the same spike is detected in multiple electrodes)
        #     only the spike of the highest amplitude will be kept - structrue stays the same
        updater.up_spikes_ampl(save_folder, save_file =spikes_largest, load_spike_file = spikes_raw, reanalize = reanalize)
     
    SPWs_potential  = 'potential_SPWs.npz'
    if not run_all_functions:
        updater.up_highWaves(save_folder, save_file = SPWs_potential, load_datafile = raw_baselined,reanalize = reanalize)
 
    SPWs_potential_numb  = 'potential_SPWs_numbered.npz'
    if not run_all_functions:
        # it numbers which wave is the same SPWs and assigns number to them
        updater.up_highWaves_numb(save_folder, save_file = SPWs_potential_numb, load_spwsfile = SPWs_potential, reanalize = reanalize)
    
    SPWs_ipsps      = 'spws_params.npz' #'spws_ipsps.npz'
    if not run_all_functions:
        # it finds the preliminary IPSPs for each of the detected waves
        updater.up_SPW_ipsp(save_folder, filter_folder, save_file = SPWs_ipsps, load_datafile = raw_baselined, load_waves = SPWs_potential_numb, load_spikes = spikes_largest, reanalize = reanalize)
    
    ipsps_corrected = 'ipsps_corrected.npz'
    if not run_all_functions:
        # correct the IPSPs (no matter if used for SPW start or for later
        updater.up_correct_ipsps(save_folder, save_fig = 'spw_ipsp', save_file = ipsps_corrected, load_datafile = raw_baselined, load_spwsipsp = SPWs_ipsps, load_spwsspike = spikes_largest, reanalize = reanalize, ext = ext)
    
    spws_large_enough = 'spw_large_enough.npz'
    min_amplitude_of_spw = 40 #microV SPW in any point, in any electrode has to be at least this amplitude
    if not run_all_functions:
        updater.up_remove_too_small_spws(save_folder, save_file = spws_large_enough, load_datafile = raw_baselined, load_spwsipsp = ipsps_corrected, min_ampl = min_amplitude_of_spw, reanalize = reanalize, ext = ext)
    
    SPWs_ipsps_beg  = 'SPWs_ipsps_beg.npz'
    if not run_all_functions:
        # finding properly each of the IPSP
        # it combines information on Waves/Ipsps and spikes to find the beginning of the SPW 
        updater.up_spws_beg(save_folder, save_fig = 'spw_ipsp', save_file = SPWs_ipsps_beg, load_datafile = raw_baselined, load_spwsipsp = ipsps_corrected, load_spwsspike = spikes_largest, reanalize = reanalize, ext = ext)

    ipsps_groups = 'ipsps_grouped.npz'
    if not run_all_functions:
        # put IPSPs to groups
        updater.up_group_ipsps(save_folder, ipsps_groups, SPWs_ipsps_beg, raw_baselined, save_file = ipsps_groups, reanalize = reanalize)    

    SPWs_ipsps_corrected = 'SPWs_ipsps_corrected.npz'
    if not run_all_functions:
        updater.up_fill_gap_between_ipsp_groups(save_folder, SPWs_ipsps_corrected, ipsps_groups, data_file = raw_baselined, reanalize = reanalize)    
    
    SPWs_all_IPSPs = 'SPWs_all_ipsps.npz'
    if not run_all_functions:
        updater.up_spws_ipsp_beg(save_folder, filter_folder, save_fig = 'spw_ipsp', save_file = SPWs_all_IPSPs, load_datafile = raw_baselined, load_spwsipsp = SPWs_ipsps_corrected, load_spwsspike = spikes_largest, reanalize = reanalize, ext = ext)
    
#    SPWs_ipsps_first_max = 'SPWs_ipsps_first_max.npz'
#    if run_all_functions:
#        updater.up_spws_first_max(save_folder, save_file = SPWs_ipsps_first_max, spws = SPWs_ipsps_beg, datafile = raw_baselined, reanalize = reanalize)

    #display_group_data(spws, spw_used[subgroups == sub], data, fs, tit = group_name)
    

    
    # not done yet!
#    if not run_all_functions:
#        updater.up_add_missing_ipsps(save_folder, SPWs_ipsps_corrected, SPWs_ipsps_corrected, data_file = raw_baselined, reanalize = reanalize)
#    spikes_inSPWs = 'spikes_inSpw.npz'
#    if run_all_functions:
#        # uses previously selected largest spikes
#        updater.up_spikes_in_spw(save_folder, save_file =spikes_inSPWs, load_spike_file = spikes_largest, load_spw_file = SPWs_ipsps_beg, reanalize = reanalize, win = win)




#    # ----> check if the group does not exist on other electrodes
    SPWs_missing_link = 'SPWs_missing_link.npz'
    if not run_all_functions:
        updater.up_add_missing_electrodes_SPW(save_folder, SPWs_missing_link, SPWs_all_IPSPs, data_file = raw_baselined, reanalize = reanalize)
    
    SPWs_merged = 'SPWs_merged.npz'
    if not run_all_functions:
        updater.up_merge_close_groups(save_folder, SPWs_merged, SPWs_missing_link, data_file = raw_baselined, reanalize = reanalize)
        
    SPWs_ipsps_corrected2 = 'SPWs_ipsps_corrected.npz'
    if run_all_functions:
        updater.up_fill_gap_between_ipsp_groups(save_folder, SPWs_ipsps_corrected2, SPWs_merged, data_file = raw_baselined, reanalize = reanalize)
#        
#    #if run_all_functions:
#    #    updater.separate
#    spws_large_enough2 = 'spw_large_enough2.npz'
#    min_amplitude_of_spw = 40 #microV SPW in any point, in any electrode has to be at least this amplitude
#    if not run_all_functions:
#        updater.up_remove_too_small_spws(save_folder, save_file = spws_large_enough2, load_datafile = raw_baselined, load_spwsipsp = SPWs_missing_link, min_ampl = min_amplitude_of_spw, reanalize = reanalize, ext = ext)

    SPWs_ipsps_final = 'SPWs_ipsps_final.npz'
    min_no_ipsps = 3
    if run_all_functions:
        updater.up_remove_with_to_few_ipsps(save_folder, SPWs_ipsps_final, SPWs_ipsps_corrected2, to_remove = min_no_ipsps, reanalize = reanalize)

    #print intr_electrode
    if intr_electrode == 1:
        data_intra = 'data_intra.npz'
        if not run_all_functions:
            updater.up_intrafile(filename, save_folder, save_file = data_intra, int_electrodes = [0], reanalize = reanalize)
        
        data_intra_base = 'data_baseintra.npz'
        if not run_all_functions:
            # removes the baseline
            updater.up_databas(save_folder, save_file = data_intra_base, load_file = data_intra, reanalize = reanalize)
            
        plot_folder = 'full_data/'
        save_plots = 'data_intra'
        if not run_all_functions:
            # plots part of the data given (in part)
            analyser.display_data(save_folder, plot_folder, save_plots, data_intra_base, trace = 0, part = [0, 100000], ext = ext)
        
        intra_spikes = 'intra_spikes.npz'
        if not run_all_functions:
            # detects intracellular spikes
            updater.up_intraSpikes(save_folder, save_file = intra_spikes, load_file = data_intra_base, reanalize = reanalize)
        
        ##SPWs_ipsps_corrected2 = SPWs_ipsps_final 
        if run_all_functions:
            # it makes the plot to exactly analyse each SPW
            analyser.plot_data_interactive(save_folder, load_datafile = raw_baselined, load_spw_ipsps = SPWs_ipsps_final  , 
                                           load_spikefile = spikes_largest, load_spikesall = spikes_raw, 
                                           load_ipspsOld =  SPWs_merged, spw_base = SPWs_potential_numb,
                                           load_dataintrafile = data_intra_base, load_intraSpikes = intra_spikes)
#    
        dist_spw_inspikes = 'spw_dist2first.npz'
        if run_all_functions:
            # finds the closest distance spw to the proceeding intracellular spike
            updater.up_dist_SpwfromSpike(save_folder, save_file = dist_spw_inspikes, load_intrafile = intra_spikes, load_spwfile = SPWs_ipsps_final, spikes = 'first', reanalize = reanalize)
        
        induc_spont_spw = 'induc_spont_spw.npz'
        max_dist = [-1.0, 7] # ms
        if run_all_functions:
            # checks which SPWs are induced and which are spontaneous (if it's further than max_dist[1] it is spontaneous)
            # if any error is being allowed it should be given in max_idst[0], e.g. -0.5 (half milisecond before intra spike
            updater.up_induc_spont_spw(save_folder, save_file = induc_spont_spw, load_distances = dist_spw_inspikes, load_spwfile = SPWs_ipsps_final, max_init_dist = max_dist, reanalize = reanalize, ext = ext)
#    
        induc_spont_equal = 'induc_spont_equal.npz'
        if run_all_functions:
            # counts spontaneous and initiated SPWs and it randomly choses set of SPWs from the bigger set so that there is equal number in both sets
            updater.equalize_number_spws(save_folder, save_file = induc_spont_equal, induc_spont = induc_spont_spw, load_distances = dist_spw_inspikes, reanalize = reanalize)

##    ##### - analyser - #####
#        solutions_folder = 'plots/'
#        
#        numIpsp2distance = 'numIPSP_distance'
#        if not run_all_functions:
#            analyser.plot_noIpsps2distance(save_folder, solutions_folder+numIpsp2distance + '/', save_plots = numIpsp2distance, spw_file = SPWs_ipsps_corrected, dist_file = dist_spw_inspikes, ext = ext)
#        
#        dist_spw2psike = 'dist_spw2spike'
#        if not run_all_functions:
#            analyser.plot_dist_spw2spike(save_folder, solutions_folder+dist_spw2psike + '/', save_plots = dist_spw2psike, dist_file = dist_spw_inspikes, ext = ext)
#        
#        alignedSPWs = 'aligned_SPWs'
#        if not run_all_functions:
#
#            analyser.plot_alignedSPW(save_folder, solutions_folder+alignedSPWs + '/', save_plots = alignedSPWs, data_file = raw_baselined, intra_data_file = data_intra_base, spike_file = induc_spont_spw, intra_spikes = intra_spikes, ext = ext)
#        
##        alignedSPWs_2all = 'aligned_SPWs2allSpikes'
##        if run_all_functions:
##            analyser.plot_alignedSPW(save_folder, solutions_folder+alignedSPWs + '/', save_plots = alignedSPWs_2all, data_file = raw_baselined, intra_data_file = data_intra, spw_file = induc_spont_spw, dist_file = dist_spw_inspikes2all, ext = ext)
#        
##        spikes_inSPWs_plot = 'spikes_inSPWs'
##        if run_all_functions:
##            analyser.plot_spikes4spw(save_folder, solutions_folder+spikes_inSPWs_plot + '/', 
##                                 save_plots = spikes_inSPWs_plot, data_file = raw_baselined, 
##                                 spike_data = spikes_largest, spw_data = induc_spont_spw, 
##                                 spikes_filter = [], ext = ext, win = win)
#        
#        
#        
#        spikePerElectrode = 'spike_per_electrode'
#        hist_spike_bins = 'all_dists_hist.npz'
#        if not run_all_functions: #->
#            analyser.plot_spike(save_folder, solutions_folder + spikePerElectrode + '/', save_plots = spikePerElectrode, 
#                            save_file = hist_spike_bins, spike_data = spikes_largest, spw_data = induc_spont_equal, 
#                            ext = ext, win = win)
#            
#    
#        alignmend_spws = 'alignmend_spws'
#        if not run_all_functions:
#            analyser.plot_spw_amplitude(save_folder, solutions_folder + spikePerElectrode + '/', save_plots = alignmend_spws, 
#                            data_file = raw_baselined, spw_data = induc_spont_spw, ext = ext)
#        
#        group_per_isps = 'group_per_isps.npz'
#        if not run_all_functions:
#            analyser.plot_spw_ipsps_no_groups(save_folder, save_file = group_per_isps, data_file = raw_baselined, 
#                                              spw_data = SPWs_ipsps_corrected, ext = ext)
#        
#        spw_more_ipsps = 'spw_more_ipsps.npz'
#        if not run_all_functions:
#            analyser.remove_with_less_ipsps(save_folder, save_file = spw_more_ipsps, 
#                                            spw_data = SPWs_ipsps_corrected,min_ipsps_group = [3])
#        
#        
#        spikes_inSPWs_plot_fig3a = 'spikes_inSPWs_fig3a'
#        if not run_all_functions:
#            analyser.plot_spikes4spw(save_folder, solutions_folder+spikes_inSPWs_plot_fig3a + '/', 
#                                 save_plots = spikes_inSPWs_plot_fig3a, data_file = raw_baselined, 
#                                 spike_data = spikes_largest, spw_data = induc_spont_spw, 
#                                 spikes_filter = filter_folder + spikes_filter, ext = ext, win = win)
#            
#            
#        used_spw_data = induc_spont_equal
#        group_per_isps_all = 'group_per_isps_all.npz'
#        if run_all_functions:
#            analyser.plot_spw_ipsps_no_groups_all(save_folder, save_file = group_per_isps_all, data_file = raw_baselined, 
#                                              spw_data = used_spw_data, ext = ext)#
#
#                                              #spw_data = spw_more_ipsps, ext = ext)
#        
#        final_results = 'final_results'
##        dendogram = 'dendogram'
##        if not run_all_functions:
##            analyser.plot_dendograms(save_folder, plot_folder = solutions_folder + final_results + '/', 
##                                      plot_file = dendogram, data_file = raw_baselined, 
##                                      spw_groups = group_per_isps_all, spw_details = used_spw_data,
##                                      spike_data = spikes_raw , ext = ext, win = win)
#        
#        
#        groups_w_firing_rate = 'groups_w_firing_rate'
#        if run_all_functions:
#            analyser.plot_groups_w_fr(save_folder, plot_folder = solutions_folder + final_results + '/', 
#                                      plot_file = groups_w_firing_rate, data_file = raw_baselined, 
#                                      spw_groups = group_per_isps_all, spw_details = used_spw_data,
#                                      spike_data = spikes_raw , ext = ext, win = win)
#                                    # spikes_largest
#         
#        
#        plot_ampl_synch = 'ampl_synchrony'
#        save_file = 'ampl_sync_dat'
#        if not run_all_functions: #->
#            analyser.plot_amplitude_vs_synchrony(save_folder, save_file, plot_folder = solutions_folder + final_results + '/', 
#                                                 plot_file = plot_ampl_synch, data_file = raw_baselined,
#                                                 spw_groups = group_per_isps_all,spw_details = SPWs_ipsps_final, ext = ext) 
        

    
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

        analyser.spw_spike_plotter(save_folder, save_name = 'main_spw', save_ext = ext)



if __name__=='__main__':
    all = ((1, 1, 0), (1, 1, 1), (1, 2, 1), (3, 1, 0), (3, 1, 1), (4, 1, 1),
           (5, 1, 0), (5, 2, 0), (5, 3, 0), (5, 4, 0), (5, 5, 0), (5, 1, 1),
           (6, 1, 0), (6, 1, 1), (7, 1, 1), (8, 1, 1), (9, 1, 1), (10, 1, 1), 
           (11, 1, 1))
    
    electr_placement = {'1': (5, 6), '3': (2, 3), '4': (2, 3),
                        '5': (5, 6), '6': (5, 6), '7': (5, 6),
                        '8': (3, 4), '9': (2, 3), '10': (2, 3),
                        '11': (6, 7)}
    # (cell_no, between_electr, and_electr) 
    
    update = 1
    analyse = 0
    sum_up_all = 0
    
    logging.basicConfig(level=logging.DEBUG)
    all_figures_folder = solutions_folder = 'plots/'
    if update == 1:

        #for nex in [15]: [1, 2, 4, 5, 11, 13, 14, 15, 16, 17, 18]
        for nex in [18]: #[5, 11, 13, 14, 15, 16, 17, 18]: #range(len(all)): #range(18, len(all)): # - 2, len(all)): #[5]: #range(12, len(all)):

        #t = importOdfSpreadsheet(file, sheet)
        #for nex in [15, 17]: #range(1, 15):
            filename, save_folder, intra  = find_folders(all[nex][0], all[nex][1], all[nex][2])
            #import pdb; pdb.set_trace()
            ex_electr = range(intra, 8+intra)
            print 'intra ' + str(intra)
            work_on_all(filename, save_folder, ex_electr, intra)
    
    if sum_up_all == 1:
        spike = False
        ampl_synch = True
        solutions_folder = get_save_folder() + 'solutions/'
        fold_mng.create_folder(solutions_folder)
        
        if spike:
            # plot gathared number of spikes on the beginning (imshow)
            file_name = 'max_electr_origin.npz'
            distance_sponts = []
            distance_inits = []
            distance_spont_diffs = []
            distance_inits_diffs = []
            
            
            for nex in range(len(all)):
                
                filename, save_folder, intra  = find_folders(all[nex][0], all[nex][1], all[nex][2])
                # check if it exists:
                exists = fold_mng.file_exists(save_folder, file_name)
                if exists:
                    cell = all[nex][0]
                    electr_place = np.mean(electr_placement[str(cell)])
                    
                    npzfile = np.load(save_folder + file_name)
                    init = npzfile['init']
                    spont = npzfile['spont']
                    init_diff = npzfile['init_diff']
                    spont_diff = npzfile['spont_diff']
                    npzfile.close()
                    
                    # calculate the distances of the major activity from the intra electrode
                    distance_sponts.append(np.abs(spont - electr_place))
                    distance_inits.append(np.abs(init - electr_place))
                    distance_spont_diffs.append(np.abs(spont_diff - electr_place))
                    distance_inits_diffs.append(np.abs(init_diff - electr_place))
            
            dists = [distance_sponts, distance_inits]
            dists_diff = [distance_spont_diffs, distance_inits_diffs]
            
            
            plt.figure()
            plt.hist(dists, normed = True)
            plt.title('normal, blue = spont, green = init')
            
            plt.figure()
            plt.hist(dists_diff, normed = True)
            plt.title('diff, blue = spont, green = init')
            
            plt.show()
        
        if ampl_synch:
            
            all_ampls = []
            all_syncs = []
            all_cells = []
            for nex in range(len(all)):
                filename, save_folder, intra  = find_folders(all[nex][0], all[nex][1], all[nex][2])
                file_name = 'ampl_sync_dat.npz'
                exists = fold_mng.file_exists(save_folder, file_name)
                
                if exists:
                    npzfile = np.load(save_folder + file_name)
                    groups = npzfile['group']
                    ampls = npzfile['all_ampls']
                    syncs = npzfile['all_syncs']
                    npzfile.close()
                    
                    if len(groups) == 1:
                        # there is only one group
                        all_cells.append(all[nex][0])
                        all_ampls.append(ampls)
                        all_syncs.append(syncs)
                
            plot_ampl_synch = 'ampl_synchrony'

            analyser.plot_amplitude_vs_synchrony_all(plot_folder = solutions_folder, 
                                                 plot_file = plot_ampl_synch, cells = all_cells, 
                                                 amplitudes = all_ampls, synchronise= all_syncs, 
                                                 ext = '.pdf')

    
    if analyse == 1:
        #for nex in range(12, len(all)):
        #for nex in [5]: #range(12, len(all)):
        #for nex in [11]:
        for nex in range(5):
            filename, save_folder, intra  = find_folders(all[nex][0], all[nex][1], all[nex][2])
            analyse_all(filename, save_folder, intra)
        

    
