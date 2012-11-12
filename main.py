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


# if new settings applied - it will rerun everything again
# to select the settings go to induc_SPW
def work_on_all(filename, save_folder, ext_electrodes = [1, 2, 3, 4, 5, 6,7], intr_electrode = 1, data_part = 'all'):
    """ work on given data file - proceed all the variables and save them """

    print 'working on: ' +  filename
    reanalize = True # set to true if to analyze the data no matter if it was already analysed or not
    
    raw_data        = 'data.npz'
    #updater.up_datafile(filename, save_folder = save_folder, save_file = raw_data, ext_electrodes = ext_electrodes, intr_electrode = 1, data_part = 'all', reanalize = reanalize)

    raw_baselined   = "data_bas.npz"
    #updater.up_databas(save_folder, save_file = raw_baselined, load_file = raw_data, reanalize = reanalize)
    
    spikes_raw      = 'spikes.npz'
    #updater.up_extraspikes(save_folder, save_file = spikes_raw, load_file = raw_baselined, reanalize = reanalize)
    
    spikes_params   = 'spikes_params.npz'
    #updater.up_expikes_params(save_folder, save_file = spikes_params, load_datafile = raw_baselined, load_spikefile = spikes_raw, reanalize = reanalize)
     
    SPWs_potential  = 'potential_SPWs.npz'
    #updater.up_highWaves(save_folder, save_file = SPWs_potential, load_datafile = raw_baselined,reanalize = reanalize)
    
    SPWs_spikes     = 'spws_spikes.npz'
    #updater.up_spws_spikes(save_folder, save_file = SPWs_spikes, load_spwsfile = SPWs_potential, load_spikefile = spikes_params, reanalize = reanalize)
    
    SPWs_spikes_ampl= 'spw_spikes_ampl.npz'
    #updater.up_spws_spikes_ampl(save_folder, save_file = SPWs_spikes_ampl, load_spwsspike = SPWs_spikes, load_spikefile = spikes_params, reanalize = reanalize)
    
    SPWs_ipsps      = 'spws_params.npz' #'spws_ipsps.npz'
    #updater.up_SPW_ipsp(save_folder, save_file = SPWs_ipsps, load_datafile = raw_baselined, load_spwsspike = SPWs_spikes, reanalize = reanalize)
    
    SPWs_ipsps_beg  = 'spw_ipsps_beg.npz'
    updater.up_spws_ipsp_beg(save_folder, save_file = SPWs_ipsps_beg, load_datafile = raw_baselined, load_spwsipsp = SPWs_ipsps, load_spwsspike = SPWs_spikes_ampl, reanalize = reanalize)
 
    #SPWs_ipsps_ampl = 'spw_ipsps_ampl.npz'
    #updater.up_spws_ipsp_ampl(save_folder, save_file = SPWs_ipsps_ampl, load_datafile = raw_baselined, load_spwsipsp = SPWs_ipsps, reanalize = reanalize)
        
    #ipsp_exSpikes = 'ipsp_exSpikes.npz'
    #updater.update_ipsp_exSpikes(save_folder, save_file = ipsp_exSpikes)
    

    


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
        
def find_folders(cell_no = 1, part = 1, oscil = 0):
    """ assigns appropriate folders and filenames for different data number given
    possible combinations:
    1. (1, 1, 0), (1, 1, 1), (1, 2, 1)
    3. (3, 1, 0), (3, 1, 1)
    4. (4, 1, 1)
    5. (5, 1, 0), (5, 2, 0), (5, 3, 0), (5, 4, 0), (5, 5, 0), (5, 1, 1), 
    6. (6, 1, 0), (6, 1, 1), 
    7. (7, 1, 0) , 
    8. (8, 1, 1), 
    9. (9, 1, 1), 
    10. (10, ),
    11. (11, )
    """
    
    save_path = '/home/maja/PhDProject/SPWs/data/induced/'
    save_data = '/home/maja/PhDProject/SPWs/SPWs/saved_data/'
    save_data = '/home/maja/phdProject/analysis/spw 09 11/'
    #save_data = 'C:/phd2/phd/saved_data/'
    #save_data = 'C:/PhDProject/saved_data/'
    
    if cell_no == 1:
        f_dir = 'Cell 1/'
        if oscil == 0:
            save_folder = 'cell1/gap_free/'
            f_name = '25102011_0019_gap free.abf' # cell1
            intra = 1
        else:
            if part == 1:
                save_folder = 'cell1/oscil_1/'
                f_name = '25102011_0020_stim oscillo_same slice but other cell not inducing SPWs.abf'
                intra = 1
            elif part == 2:
                save_folder = 'cell1/oscil_2/'
                f_name = '25102011_0023_stim oscillo_ inducing cell.abf'
                intra = 1
        
    elif cell_no == 3:
        f_dir = 'Cell 3/'
        if oscil == 0:
            save_folder = 'cell3/gap_free/'
            f_name = '08112011_0000_gap free.abf' # cell3 (no intra)
            intra = 0
        else:
            f_name = '08112011_0017_stim oscillo.abf'
            save_folder = 'cell3/oscil/'
            intra = 1
    elif cell_no == 4:
        f_dir = 'Cell 4/'
        f_name = '08112011_0020_stim oscillo_same slice as cell 2.abf'
        save_folder = 'cell4/'
        intra = 1
    elif cell_no == 5:
        f_dir = 'Cell 5/'
        if oscil == 0:
            intra = 0
            if part == 1:
                f_name = '17112011_0000_gap free_p1.abf' # cell5 p1
                save_folder = 'cell5/part1/'
            elif part == 2:
                f_name = '17112011_0000_gap free_p2.abf'
                save_folder = 'cell5/part2/'
            elif part == 3:
                f_name = '17112011_0000_gap free_p3.abf'
                save_folder = 'cell5/part3/'
            elif part == 4:
                f_name = '17112011_0000_gap free_p4.abf'
                save_folder = 'cell5/part4/'
            elif part == 5:
                f_name = '17112011_0000_gap free_p5.abf'
                save_folder = 'cell5/part5/'
        else:
            f_name = '17112011_0002_stim oscillo.abf'
            save_folder = 'cell5/oscil/'
            intra = 1
    elif cell_no == 6:
        f_dir = 'Cell 6/'
        if oscil == 0:
            f_name = '02122011_0000_gap free stim.abf'# cell6
            save_folder = 'cell6/gap_free/'
            intra = 1
        else:
            f_name = '02122011_0001_oscillo stim.abf'
            save_folder = 'cell6/oscil/'
            intra = 1
    elif cell_no == 7:
        f_dir = 'Cell 7/'
        f_name = '02122011_0006_gap free stim.abf'
        save_folder = 'cell7/'
        intra = 1
    elif cell_no == 8:
        f_dir = 'Cell 8/'
        f_name = '08122011_0012_stim oscillo_as control take gap free recording from cell 5.abf' # cell8
        save_folder = 'cell8/'
        intra = 1
    elif cell_no == 9:
        f_dir = 'Cell 9/'        
        f_name = '08122011_0013_stim gap free_same slice as cell 4.abf'
        save_folder = 'cell9/'
        intra = 1
    elif cell_no == 10:
        f_dir = 'Cell 10/'
        f_name = '09122011_0003.abf'
        save_folder = 'cell10/'
        intra = 1
    elif cell_no == 11:
        f_dir = 'Cell 11/'
        f_name = '09122011_0004.abf'
        save_folder = 'cell11/'
        intra = 1
    read_name = save_path + f_dir + f_name
    save_folder = save_data + save_folder
    return read_name, save_folder, intra  
        
#test_correctness()

#

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
        for nex in [14, 15]: #range(5, 16): #[14, 15]: #range(7, len(all)):
            filename, save_folder, intra  = find_folders(all[nex][0], all[nex][1], all[nex][2])
            
            ex_electr = range(1+intra, 8+intra)
            work_on_all(filename, save_folder, ex_electr, intra)
    
    if analyse == 1:
        #for nex in range(12, len(all)):
        #for nex in [5]: #range(12, len(all)):
        #for nex in [11]:
        for nex in range(14, 15):
            filename, save_folder, intra  = find_folders(all[nex][0], all[nex][1], all[nex][2])
            analyse_all(filename, save_folder, intra)
        
#import pdb; pdb.set_trace()



#ispw.update_all(read_name, ext_electrodes, save_folder, int_electrode, data_part)
#test_corect(save_folder)
#work_on_spws(save_folder)

#group_spws(save_folder)
#work_on_swp_spikes(save_folder)
#corr_spws_electrwise(save_folder)

#corr_spws_electrwise(save_folder, data_file = 'spw_corr')
#plot_distance_meancorrSPW(save_folder)

#all_corr, all_move = correlate_spw(save_folder, 1)
#group_corr(save_folder)




#f_name = ['25102011_0019_gap free.abf', '08112011_0000_gap free.abf', '17112011_0000_gap free_p1.abf', '17112011_0000_gap free_p2.abf', '17112011_0000_gap free_p3.abf', '17112011_0000_gap free_p4.abf', '17112011_0000_gap free_p5.abf', '02122011_0000_gap free stim.abf', '02122011_0006_gap free stim.abf', '08122011_0013_stim gap free_same slice as cell 4.abf'] 
## cell1, cell3, cell5_p1, cell5_p2, cell5_p3, cell5_p4, cell5_p5, cell6, cell7,  cell9
#f_dir = '/home/maja/PhDProject/SPWs/data/induced/'
#f_dir_cell = ['Cell 1/', 'Cell 3/', 'Cell 5/', 'Cell 5/', 'Cell 5/', 'Cell 5/', 'Cell 5/', 'Cell 6/', 'Cell 7/', 'Cell 9/']
#save_folder = '/home/maja/PhDProject/SPWs/SPWs/saved_data/'
#fave_folder_cell = ['cell1/', 'cell3/', 'cell5/part1/', 'cell5/part2/', 'cell5/part3/', 'cell5/part4/', 'cell5/part5/', 'cell6/', 'cell7/', 'cell1/', 'cell1/']
#for i in range(len(f_name)):
#
#read_name = f_dir + f_name  
#
#print len(f_name)
#print len(f_dir_cell)
#
#
#
#ispw.update_all(read_name, ext_electrodes, save_folder, int_electrode, data_part)
#
#
#
#
#int_electrode = 1
#f_dir = '/home/maja/PhDProject/SPWs/data/induced/Cell 7/'
#f_name = '02122011_0006_gap free stim.abf' # cell7 - error check!
#save_folder = '/home/maja/PhDProject/SPWs/SPWs/saved_data/cell7/'
#read_name = f_dir + f_name  
#ispw.update_all(read_name, ext_electrodes, save_folder, int_electrode, data_part)
















#data_uspl_extra, fs_new_up = read_upsample(save_folder)
#distances, min_distances, fs = read_dist_fromSpike(save_folder)
#spike_idxs, spikes_ampl, whole_spikes, fs_new = read_extraspikes(save_folder)
#if fs_new_up != fs:
#    fs_div = fs_new_up / fs # if the sampling rate is different
#    spike_idxs = [spike_idxs[i] / fs_div for i in range(len(spike_idxs))]
#data_uspl_extra, fs_new_up = ispw.update_upsample(data_all, fs, save_folder, uspl = 5, data_file = 'data_uspl_extra')

#
#work_on_extraspikes(save_folder)
#plotting(save_folder)
#work_on_spws(save_folder)

#sp_idx_first, sp_idx_all, fs = read_intra_spikes(save_folder)
#spike_idxs, spikes_ampl, whole_spikes, fs_new = read_extraspikes(save_folder)
#
#if fs_new != fs:
#    fs_div = fs_new / fs # if the sampling rate is different
#    spike_idxs = [spike_idxs[i] / fs_div for i in range(len(spike_idxs))]
#
#
#ispw.update_dist_fromSpike(sp_idx_all, spike_idxs, fs, save_folder)


#work_on_all(save_folder, read_name)
#work_on_spws(save_folder)
#work_on_extraspikes(save_folder)
#work_on_spikes(save_folder)
# run again all the settings

#work_on_extraspikes(save_folder)   

#work_on_all(save_folder, f_cell)

#work_on_spikes(save_folder, f_cell)
#py.show()

#ispw.update_extraspikes(data, fs, save_name)
#npzfile = np.load(save_name + "distances.npz")
#print npzfile.files
#files =  npzfile.files
#for i in range(len(files)):
#    #print npzfile[files[i]]
#    print np.shape(npzfile[files[i]])
#    print npzfile[files[i]]
    #print len(npzfile[files[i]])
    #print npzfile[files[i]][0]
    #print npzfile[files[i]][1]            
    #print np.size(data)
    #pass
    


#Nevertheless, if you calculate for exemple the mean half duration of extracellular spikes in function of latency from intracellular spike, 
#I think that you will find at least for the first 3ms a mean value that is less than what you will find after. That means that at this point mostly interneurons are actives. What you can do is also to take as interneurons only spikes that have a half duration less than ... something. And exclude others. You will have in this way the timing of inteurneurons discharges in respect with intracellular spikes and space (spikes taken independently from each electrode)
#It is clear that all electrodes do not record spikes with the same quality, so do what you can, not more.
#I will discuss it better later. In conclusion, do what you can, but since it is not the main purpose for now, do not spend too much time on it.
#
#To be continued
#
#Michael
    
