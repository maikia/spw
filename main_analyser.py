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

#import bottleneck as bn

# if new settings applied - it will rerun everything again
# to select the settings go to induc_SPW
def update_all_plots_one_cell(filename, save_folder, ext_electrodes = [1, 2, 3, 4, 5, 6, 7], intr_electrode = 1, data_part = 'all'):
    """ work on given data file"""
    ext = '.eps'
    #================files which were previously analysed =======================
    names = ['max_2_', 'min_3_', 'all_', 'min_2_', 'max_1_'] # depending on max number of IPSPs used it should be added before
     # name of the file: spws_file, distances, equal_init_spont
    name_used = names[2]
    
    spws_file = name_used + 'SPWs_ipsps_final.npz'
    distances = name_used + 'spw_dist2first.npz'
    equal_init_spont = name_used + 'induc_spont_equal.npz'
    raw_data = "data_bas.npz"
    intra_data = 'data_baseintra.npz'
    spike_file = 'spikes_largest.npz'
    intra_spike_file = 'intra_spikes.npz'
    all_spikes = 'spikes.npz'
    #==============================================================
    #import pdb; pdb.set_trace()
    plots_folder = 'plots/' + name_used
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
            """ plots relation between distance from the spike and number of ipsp groups in a SPW"""
            fold_mng.create_folder(save_folder + save_plot_in)
            analyser.plot_dist_spw2spike(save_folder, save_plot_in, save_plots = dist_spw2psike, dist_file = distances, ext = ext)

        firing_rate = name_used + 'firing_rate'
        save_plot_in = plots_folder+ firing_rate + '/'
        if not run_all_functions:
            fold_mng.create_folder(save_folder + save_plot_in)
            analyser.plot_fr_after_spike(save_folder, plot_folder = save_plot_in, 
                                      plot_file = firing_rate, intra_spikes = intra_spike_file,
                                      spike_data = spike_file , ext = ext, win = win)

        firing_rate_and_spws = name_used + 'firing_rate_and_spws'
        save_plot_in = plots_folder+ firing_rate_and_spws + '/'
        if not run_all_functions:
            fold_mng.create_folder(save_folder + save_plot_in)
            analyser.plot_fr_after_spike_and_distances_after_spike(save_folder, plot_folder = save_plot_in,  
                                      plot_file = firing_rate_and_spws, intra_spikes = intra_spike_file, dist_file = distances,
                                      spike_data = spike_file , ext = ext)

        separateSPWs = 'separate_SPWs'
        if not run_all_functions:
            analyser.plot_different_SPWs(save_folder, save_plot_in, save_plots = separateSPWs, data_file = raw_data, intra_data_file = intra_data, induc_spont = equal_init_spont, intra_spikes = intra_spike_file, ext = '.png')
        
        alignedSPWs = 'aligned_SPWs'
        save_plot_in = plots_folder+ alignedSPWs + '/'
        if not run_all_functions:
            fold_mng.create_folder(save_folder + save_plot_in)
            analyser.plot_alignedSPW(save_folder, save_plot_in, save_plots = alignedSPWs, data_file = raw_data, intra_data_file = intra_data, induc_spont = equal_init_spont, intra_spikes = intra_spike_file, ext = '.eps')

        
        spikePerElectrode = name_used + 'spike_per_electrode'
        hist_spike_bins = name_used + 'all_dists_hist.npz'
        save_name_max_electr = name_used + 'max_electr_origin.npz'
        save_plot_in = plots_folder+ spikePerElectrode + '/'
        if not run_all_functions: 
            fold_mng.create_folder(save_folder + save_plot_in)
            analyser.plot_spike(save_folder, save_plot_in, save_plots = spikePerElectrode, 
                            save_file = hist_spike_bins, save_name_max_electr = save_name_max_electr, 
                            spike_data = spike_file, spw_data = equal_init_spont, 
                            ext = ext, win = win)
            
        #import pdb; pdb.set_trace()
        group_per_isps_all = name_used + 'group_per_isps_all.npz'
        if not run_all_functions: # ok
            analyser.plot_spw_ipsps_no_groups_all(save_folder, save_file = group_per_isps_all, data_file = raw_data, 
                                              spw_data = equal_init_spont, ext = ext)#



        
        groups_w_firing_rate = name_used + 'groups_w_firing_rate'
        save_plot_in = plots_folder+ groups_w_firing_rate + '/'
        if not run_all_functions:
            fold_mng.create_folder(save_folder + save_plot_in)
            analyser.plot_groups_w_fr(save_folder, plot_folder = save_plot_in, 
                                      plot_file = groups_w_firing_rate, data_file = raw_data, 
                                      spw_groups = group_per_isps_all, spw_details = equal_init_spont,
                                      spike_data = all_spikes , ext = '.eps', win = win)
                                    # spikes_largest           
        
        cumulative_plot = 'cumulative_plot'
        save_plot_in = plots_folder+ cumulative_plot + '/'
        save_file = name_used + 'cum_change_variance.npz'
        if not run_all_functions:
            fold_mng.create_folder(save_folder + save_plot_in)
            analyser.cum_distribution_funct(save_folder, save_file, plot_folder = save_plot_in, plot_file = cumulative_plot, data_file = raw_data, 
                                      spw_details = equal_init_spont,
                                      ext = ext, win = win)
        
        #import pdb; pdb.set_trace() 
#
        
        plot_ampl_synch = 'ampl_synchrony'
        save_file = name_used + 'ampl_sync_dat'
        save_plot_in = plots_folder+ plot_ampl_synch + '/'
        if run_all_functions: # ok
            fold_mng.create_folder(save_folder + save_plot_in)
            analyser.plot_amplitude_vs_synchrony(save_folder, save_file, plot_folder = save_plot_in, 
                                                 plot_file = plot_ampl_synch, data_file = raw_data,
                                                 spw_groups = group_per_isps_all, spw_details = equal_init_spont, ext = ext) 


def update_all_plots_all_cells(filename, save_folder, ext_electrodes = [1, 2, 3, 4, 5, 6, 7], intr_electrode = 1, data_part = 'all'):
    pass

if __name__=='__main__':
#    all = ((1, 1, 0), (1, 1, 1), (1, 2, 1), (3, 1, 0), (3, 1, 1), (4, 1, 1),
#           (5, 1, 0), (5, 2, 0), (5, 3, 0), (5, 4, 0), (5, 5, 0), (5, 1, 1),
#           (6, 1, 0), (6, 1, 1), (7, 1, 1), (8, 1, 1), (9, 1, 1), (10, 1, 1), 
#           (11, 1, 1))
    
    all = ((1, 2, 1), (3, 1, 0), (3, 1, 1), (4, 1, 1),
           (5, 1, 0), (5, 2, 0), (5, 3, 0), (5, 4, 0), (5, 5, 0), (5, 1, 1),
           (6, 1, 0), (7, 1, 1), (8, 1, 1), (9, 1, 1), (10, 1, 1), 
           (11, 1, 1))    
    electr_placement = {'1': (5, 6), '3': (2, 3), '4': (2, 3),
                        '5': (5, 6), '6': (5, 6), '7': (5, 6),
                        '8': (3, 4), '9': (2, 3), '10': (2, 3),
                        '11': (1, 2)}
    # (cell_no, between_electr, and_electr) 
    
    update = 1
    sum_up_all = 1
    
    logging.basicConfig(level=logging.DEBUG)
    all_figures_folder = solutions_folder = 'plots/'
    if update == 1:

        for nex in range(13, len(all)):
        #for nex in [15]: #range(len(all)): #[18]: #[5, 11, 13, 14, 15, 16, 17, 18]: #range(len(all)): #range(18, len(all)): # - 2, len(all)): #[5]: #range(12, len(all)):

        #t = importOdfSpreadsheet(file, sheet)
        #for nex in [15, 17]: #range(1, 15):
            filename, save_folder, intra  = find_folders(all[nex][0], all[nex][1], all[nex][2])
            
            ex_electr = range(intra, 8+intra)
            print 'intra ' + str(intra)
            update_all_plots_one_cell(filename, save_folder, ex_electr, intra)
    
    if sum_up_all == 1:
        names = ['max_2_', 'min_3_', 'all_', 'min_2_', 'max_1_'] # depending on max number of IPSPs used it should be added before
        # name of the file: spws_file, distances, equal_init_spont
        name_used = names[2]
        

        spike = False
        ampl_synch = True
        cum_change_var = False
        
        solutions_folder = get_save_folder() + 'solutions/'
        fold_mng.create_folder(solutions_folder)
        
        if cum_change_var:
            # gather all the cumulative change of variance and plot them on one plot
            file_name = name_used + 'cum_change_variance.npz'
            
            all_var_spont = []
            all_var_init = []
            
            for nex in range(len(all)):
                #import pdb; pdb.set_trace()
                filename, save_folder, intra  = find_folders(all[nex][0], all[nex][1], all[nex][2])
                # check if it exists:       
                exists = fold_mng.file_exists(save_folder, file_name)

                if exists:  
                    print all[nex][0]               
                    npzfile = np.load(save_folder + file_name)
                    
                    init_temp = npzfile['cum_change_init']
                    spont_temp = npzfile['cum_change_spont']
                    timeline = npzfile['timeline']
                    fs = npzfile['fs']
                    npzfile.close()
                    all_var_spont.append(spont_temp.tolist())
                    all_var_init.append(init_temp.tolist())
            fs = fs.tolist()
            fs = int(fs)      
            all_var_spont = np.array(all_var_spont)
            all_var_init = np.array(all_var_init)
            #import pdb; pdb.set_trace()
            plot_cum = 'cumulative_distribution'
            analyser.plot_all_cum_change_var(plot_folder = solutions_folder, 
                                                 plot_file = plot_cum, all_var_spont = all_var_spont, 
                                                 all_var_init = all_var_init, timeline = timeline, fs = fs,
                                                 ext = '.png')
                    
                    
        if spike:
            # plot gathared number of spikes on the beginning (imshow)
            file_name = name_used + 'max_electr_origin.npz'
            distance_sponts = []
            distance_inits = []
            distance_spont_diffs = []
            distance_inits_diffs = []
            
            all_spont_spikes = []
            all_induc_spikes = []
            all_numbers = []
            
            for nex in range(len(all)):
                #import pdb; pdb.set_trace()
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
                    spont_spikes = npzfile['spont_all']
                    induc_spikes = npzfile['induc_all']
                    spw_numb = npzfile['spw_number']
                    npzfile.close()
                    #print cell
                    spont_spikes =  np.transpose(spont_spikes).tolist()
                    induc_spikes =  induc_spikes.tolist()
                    print spw_numb
                    #import pdb; pdb.set_trace()
                    all_spont_spikes.append(np.array(spont_spikes))
                    all_induc_spikes.append(np.array(induc_spikes))
                    all_numbers.append(int(spw_numb))
                    
                    
                    # calculate the distances of the major activity from the intra electrode
                    distance_sponts.append(np.abs(spont - electr_place))
                    distance_inits.append(np.abs(init - electr_place))
                    distance_spont_diffs.append(np.abs(spont_diff - electr_place))
                    distance_inits_diffs.append(np.abs(init_diff - electr_place))
            #import pdb; pdb.set_trace()    
            np.savetxt(solutions_folder + 'test.txt', np.transpose(all_induc_spikes), delimiter='\t',  fmt='%i') 
                    
            
            dists = [distance_sponts, distance_inits]
            dists_diff = [distance_spont_diffs, distance_inits_diffs]
            
            
            analyser.plot_spikes(save_folder = solutions_folder, save_name = 'init_distance', distances = dists, ext = '.png')
            analyser.plot_spikes(save_folder = solutions_folder, save_name = 'init_distance normed', distances = dists_diff, ext = '.png')

            
        
        if ampl_synch:
            
            all_ampls1 = []
            all_ampls2 = []
            all_syncs1 = []
            all_syncs2 = []
            all_cells = []
            all_groups1 = []
            all_groups2 = []
            for nex in range(len(all)):
                #if all[nex][0] == 1 or all[nex][0] == 6:
                #   import pdb; pdb.set_trace()
                filename, save_folder, intra  = find_folders(all[nex][0], all[nex][1], all[nex][2])
                file_name = name_used + 'ampl_sync_dat.npz'
                exists = fold_mng.file_exists(save_folder, file_name)
                print filename
                if exists:
                    
                    npzfile = np.load(save_folder + file_name)
                    groups = npzfile['group']
                    ampls = npzfile['all_ampls']
                    syncs = npzfile['all_syncs']
                    group_nos = npzfile['groups_ipsp']
                    
                    npzfile.close()
                    all_ampls1.append(ampls[0])
                    all_ampls2.append(ampls[1])
                    all_syncs1.append(syncs[0])
                    all_syncs2.append(syncs[1])
                    all_groups1.append(group_nos[0])
                    all_groups2.append(group_nos[1])
                    all_cells.append(all[nex][0])                                      
                    
                else:
                    print save_folder + file_name
                    print 'does not exist'   
            #import pdb; pdb.set_trace()    
            plot_folder = solutions_folder + name_used + '/'
            fold_mng.create_folder(plot_folder)
            plot_ampl_synch = 'ampl_synchrony'
            
            analyser.plot_amplitude_vs_synchrony_all(plot_folder = plot_folder, 
                                                 plot_file = plot_ampl_synch, cells = all_cells, 
                                                 amplitudes = [all_ampls1, all_ampls2], synchronise= [all_syncs1, all_syncs2], 
                                                 group_nos = [all_groups1, all_groups2], names = groups, ext = '.eps')

