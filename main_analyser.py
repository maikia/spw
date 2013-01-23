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
def update_all_plots_one_cell(filename, save_folder, ext_electrodes = [1, 2, 3, 4, 5, 6, 7], intr_electrode = 1, data_part = 'all'):
    """ work on given data file"""
    ext = '.png'
    #================files which were previously analysed =======================
    spws_file = 'SPWs_ipsps_final.npz'
    distances = 'spw_dist2first.npz'
    equal_init_spont = 'induc_spont_equal.npz'
    raw_data = "data_bas.npz"
    intra_data = 'data_baseintra.npz'
    spike_file = 'spikes_largest.npz'
    intra_spike_file = 'intra_spikes.npz'
    #==============================================================
    #import pdb; pdb.set_trace()
    plots_folder = 'plots/'
    print 'working on: ' +  filename
    
    #solutions_folder = 'plots/'
    delete_old = False #!!!! it will delete all the previously saved plots so be careful!
    run_all_functions = True
    
    if delete_old: # !!!!!!
        fold_mng.remove_folder(save_folder + plots_folder)
        
    win = [0, 0]
    
    if intr_electrode == 1:     
        
        plot_name_ipsps_2_dist = 'numIPSP_distance'
        save_plot_in = plots_folder+plot_name_ipsps_2_dist + '/'
        
        if not run_all_functions:
            """ plots relation between distance from the spike and number of ipsp groups in a SPW """
            fold_mng.create_folder(save_folder + save_plot_in)
            analyser.plot_noIpsps2distance(save_folder,save_plot_in , save_plots = plot_name_ipsps_2_dist, spw_file = spws_file, dist_file = distances, ext = ext)
        
        dist_spw2psike = 'dist_spw2spike'
        save_plot_in = plots_folder+dist_spw2psike + '/'
        if not run_all_functions:
            """ it plots histogram showing number of SPWs distant from spike"""
            fold_mng.create_folder(save_folder + save_plot_in)
            analyser.plot_dist_spw2spike(save_folder, save_plot_in, save_plots = dist_spw2psike, dist_file = distances, ext = ext)
        
        alignedSPWs = 'aligned_SPWs'
        save_plot_in = plots_folder+ alignedSPWs + '/'
        if run_all_functions:
            """ plots few SPWs - not really cool function"""
            fold_mng.create_folder(save_folder + save_plot_in)
            analyser.plot_alignedSPW(save_folder, save_plot_in, save_plots = alignedSPWs, data_file = raw_data, intra_data_file = intra_data, induc_spont = equal_init_spont, intra_spikes = intra_spike_file, ext = ext)
        import pdb; pdb.set_trace()
        
#        alignedSPWs_2all = 'aligned_SPWs2allSpikes'
#        if run_all_functions:
#            analyser.plot_alignedSPW(save_folder, solutions_folder+alignedSPWs + '/', save_plots = alignedSPWs_2all, data_file = raw_baselined, intra_data_file = data_intra, spw_file = induc_spont_spw, dist_file = dist_spw_inspikes2all, ext = ext)
        
#        spikes_inSPWs_plot = 'spikes_inSPWs'
#        if run_all_functions:
#            analyser.plot_spikes4spw(save_folder, solutions_folder+spikes_inSPWs_plot + '/', 
#                                 save_plots = spikes_inSPWs_plot, data_file = raw_baselined, 
#                                 spike_data = spikes_largest, spw_data = induc_spont_spw, 
#                                 spikes_filter = [], ext = ext, win = win)
        
        
        
        spikePerElectrode = 'spike_per_electrode'
        hist_spike_bins = 'all_dists_hist.npz'
        if not run_all_functions: #->
            analyser.plot_spike(save_folder, solutions_folder + spikePerElectrode + '/', save_plots = spikePerElectrode, 
                            save_file = hist_spike_bins, spike_data = spikes_largest, spw_data = induc_spont_equal, 
                            ext = ext, win = win)
            
    
        alignmend_spws = 'alignmend_spws'
        if not run_all_functions:
            analyser.plot_spw_amplitude(save_folder, solutions_folder + spikePerElectrode + '/', save_plots = alignmend_spws, 
                            data_file = raw_baselined, spw_data = induc_spont_spw, ext = ext)
        
        group_per_isps = 'group_per_isps.npz'
        if not run_all_functions:
            analyser.plot_spw_ipsps_no_groups(save_folder, save_file = group_per_isps, data_file = raw_baselined, 
                                              spw_data = SPWs_ipsps_corrected, ext = ext)
        
        spw_more_ipsps = 'spw_more_ipsps.npz'
        if not run_all_functions:
            analyser.remove_with_less_ipsps(save_folder, save_file = spw_more_ipsps, 
                                            spw_data = SPWs_ipsps_corrected,min_ipsps_group = [3])
        
        
        spikes_inSPWs_plot_fig3a = 'spikes_inSPWs_fig3a'
        if not run_all_functions:
            analyser.plot_spikes4spw(save_folder, solutions_folder+spikes_inSPWs_plot_fig3a + '/', 
                                 save_plots = spikes_inSPWs_plot_fig3a, data_file = raw_baselined, 
                                 spike_data = spikes_largest, spw_data = induc_spont_spw, 
                                 spikes_filter = filter_folder + spikes_filter, ext = ext, win = win)
            
            
        used_spw_data = induc_spont_equal
        group_per_isps_all = 'group_per_isps_all.npz'
        if run_all_functions:
            analyser.plot_spw_ipsps_no_groups_all(save_folder, save_file = group_per_isps_all, data_file = raw_baselined, 
                                              spw_data = used_spw_data, ext = ext)#

                                              #spw_data = spw_more_ipsps, ext = ext)
        
        final_results = 'final_results'
#        dendogram = 'dendogram'
#        if not run_all_functions:
#            analyser.plot_dendograms(save_folder, plot_folder = solutions_folder + final_results + '/', 
#                                      plot_file = dendogram, data_file = raw_baselined, 
#                                      spw_groups = group_per_isps_all, spw_details = used_spw_data,
#                                      spike_data = spikes_raw , ext = ext, win = win)
        
        
        groups_w_firing_rate = 'groups_w_firing_rate'
        if run_all_functions:
            analyser.plot_groups_w_fr(save_folder, plot_folder = solutions_folder + final_results + '/', 
                                      plot_file = groups_w_firing_rate, data_file = raw_baselined, 
                                      spw_groups = group_per_isps_all, spw_details = used_spw_data,
                                      spike_data = spikes_raw , ext = ext, win = win)
                                    # spikes_largest
         
        
        plot_ampl_synch = 'ampl_synchrony'
        save_file = 'ampl_sync_dat'
        if not run_all_functions: #->
            analyser.plot_amplitude_vs_synchrony(save_folder, save_file, plot_folder = solutions_folder + final_results + '/', 
                                                 plot_file = plot_ampl_synch, data_file = raw_baselined,
                                                 spw_groups = group_per_isps_all,spw_details = SPWs_ipsps_final, ext = ext) 


def update_all_plots_all_cells(filename, save_folder, ext_electrodes = [1, 2, 3, 4, 5, 6, 7], intr_electrode = 1, data_part = 'all'):
    pass

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
            update_all_plots_one_cell(filename, save_folder, ex_electr, intra)
    
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