import numpy as np
import data_mang as dat
import pylab as plt
#import read_data as reader
import update_data as updater
import induc_SPW as ispw
import scipy.signal as signal
import b_analyser as ba
#import cluster as clust
from matplotlib.ticker import NullFormatter
from mpl_toolkits.mplot3d import Axes3D

def define_colors(no_colors = 8):

    #RGB_tuples = [(x*1.0/no_colors, x*0.5/no_colors,x*0.5/no_colors) for x in range(no_colors)]
    #RGB_tuples = np.array(RGB_tuples)
    #RGB_tuples = [(1, 0, 0), (1, 0.5, 0), (1, 1, 0), (0.5, 1, 0), (0, 1, 0.5), (0, 1, 1), (0, 0.5, 1), (0, 0, 1)]
    RGB_tuples = [(0.5, 1, 0), (0, 1, 0.5), (0, 1, 1), (0, 0.5, 1), (0, 0, 1), (0.5, 0, 1), (1, 0, 1), (1, 0, 0.5)]
    return RGB_tuples

def micro():
    return "$\mu$"

def use_waves(save_folder, spws =1):    
    if spws == 1:
        # reads the data necessary for analysing spws 
        data_bas, fs = reader.read_databas(save_folder)
        data_bas = np.array(data_bas)
        spw_idxs, spw_maxs, starts_spw, ends_spw, lengths_spw, fs = reader.read_SPWs(save_folder)
        starts_wave = ba.idx2recarr(starts_spw)
    else:
        # extracellular spikes
        data_bas, fs = reader.read_datafile(save_folder)
        spike_idxs, spikes_ampl, fs_new  = reader.read_extraspikes(save_folder)
        if fs != fs_new:
            print 'uncompatible sampling rates!!!'
        starts_wave = ba.idx2recarr(spike_idxs)
    return data_bas, fs, starts_wave

def plot_mean_waves_elecrwise_from_saved(save_folder, save_ext):
    # plot mean SPWs on all of the electrodes
    data_bas, fs, starts_spw = use_waves(save_folder, spws = 1)
    plot_mean_wave_elecrwise(data_bas, fs, starts_spw, save_folder, save_ext = save_ext)
    
    # plot mean spikes on all of the electrodes
    data_bas, fs, starts_spw = use_waves(save_folder, spws = 0)
    plot_mean_wave_elecrwise(data_bas, fs, starts_spw, save_folder, save_ext = save_ext, save_name = 'avgSpikes', title = 'mean spikes for each electrode (baseline was not removed)', win = [-2, 5])
    




def calculate_corr_matrix_electrowise(save_folder, save_name = 'corr_electr_spws', save_ext = '.png', title = 'correlation, el: ', win = (-10, 30), spws = 1):
    add = 'wyniki/'
    scale = 'ms'
    data, fs, starts_waves = use_waves(save_folder, spws)

    win_pts = (ispw.ms2pts(win[0], fs), ispw.ms2pts(win[1], fs))
                 
    waves, idx = ba.extract_spws(data, starts_waves, fs, win_pts)             
    cxy_all = []    
    t = dat.get_timeline(waves, fs, scale) 

    # remove all the waves which are outside data
    idx_aligned_all = []
    
    for electr in range(len(data)):
        idx_el  = idx[idx['electrode']==electr]
        idx_aligned = ba.align_spws(data, idx_el, fs, win_pts)
        spws_el, idx_aligned = ba.extract_spws(data, idx_aligned, fs, win_pts)
        cxy = ba.correlation_matrix(spws_el)

        if np.size(cxy) > 1:
            cxy_all.append(cxy)
            idx_aligned_all.append(idx_aligned)
            
            # save plot of correlation matrices
            #fig = plt.figure(figsize=(8.0, 5.0))    

#            plt.imshow(cxy, interpolation='nearest', origin='lower')0
#            tit = title + str(electr)
#            plt.title(tit)
#            plt.xlabel('wave no.')
#            plt.ylabel('wave no.')
#            save_nom = save_folder + add + save_name + str(electr) + save_ext
#            fig.savefig(save_nom,dpi=600)        
#            plt.close()
            
            # save all the plot for this electrode original and moved SPWs
            for with_mean in range(2):
                
                
                
                fig = plt.figure(figsize=(8.0, 5.0))
                plt.subplot(211)
                if with_mean == 0:
                    plt.plot(t, spws_el)
                elif with_mean == 1:
                    plt.plot(t, spws_el, 'b')
                    plt.plot(t, np.mean(spws_el, 1), 'r')
                tit = 'originally found spws, electr:' + str(electr)
                plt.title(tit)
                xlab = 'time [' + scale + ']'
                ylab = 'voltage [' + micro() + 'V]'
                #plt.xlabel(xlab)
                plt.ylabel(ylab)
                plt.xlim([0, t[-1]])
                plt.subplot(212)
                if with_mean == 0:
                    plt.plot(t, waves[:,idx['electrode']==electr]) #plot all traces in
                elif with_mean == 1:
                    plt.plot(t, waves[:,idx['electrode']==electr], 'b')
                    plt.plot(t, np.mean(waves[:,idx['electrode']==electr], 1), 'r')
                
                plt.title('spws corrected by correlation')
                xlab = 'time [' + scale + ']'
                ylab = 'voltage [' + micro() + 'V]'
                plt.xlabel(xlab)
                plt.ylabel(ylab)      
                save_nom = save_folder + add + 'compare_shift' + str(with_mean) + save_name + str(electr) + save_ext
                fig.savefig(save_nom,dpi=600)
                plt.clf()
                plt.close(fig)
#        elif np.size(cxy) == 1:
#              cxy_all.append(np.zeros(1))
#              idx_aligned_all.append(np.zeros(1))
        else:
            cxy_all.append(np.zeros(np.size(cxy)))
            idx_aligned_all.append(np.zeros(np.size(cxy)))
                
    np.savez(save_folder + save_name, cxy_all, idx_aligned_all)


def find_clust_extraspikes(save_folder, n_clusters = 3, save_name = 'extra_clusters', save_clust = 'extra_clusters', save_ext = '.png'):
    #n_pts = 100
    
    spike_idxs, all_valley_to_peak, all_half_valley_width, all_half_peak_width, fs_new = reader.read_ex_sparamas(save_folder)
    a_s = ba.idx2recarr(all_valley_to_peak, time_dtype=np.float32)
    c_s = ba.idx2recarr(all_half_peak_width, time_dtype=np.float32)
    spike_id = ba.idx2recarr(all_half_peak_width)
    spik_clust = []
    
    for electr in range(len(spike_idxs)):
#        spik = []
#        for trace in range(len(spike_idxs[electr])):
#            for val in range(len(spike_idxs[electr][trace])):
        a = a_s[a_s['electrode']==electr]['time']
        c = c_s[c_s['electrode']==electr]['time']

        X = np.vstack((a,c)).T
        X = X/np.std(X,0);
        X = X-X.mean(0)
        X = X/np.std(X,0)

        #X[:n_pts/2,:]+=5
    
        cluster = clust.kmeans(n_clusters, X)
    
        symbols = ['>', '<', 'o']
        cols = ['r', 'b', 'k']
        fig = plt.figure()
        for i in range(n_clusters):
            plt.plot(X[cluster==i,0], X[cluster==i,1], symbols[i], color=cols[i], alpha = 0.03)

        tit = 'clustering extracellular spikes, electr: ' + str(electr)
        plt.title(tit)
        plt.xlabel('all valley to peak (a)')
        plt.ylabel('half pkeas (c)')
        nom = save_folder + 'wyniki/' + save_name + str(electr) + save_ext
        fig.savefig(nom,dpi=600)
        plt.close()
            
        spik_clust.append(cluster)
    np.savez(save_folder + save_clust, spik_clust)   

def plot_hist_dist_ex2intra(save_folder, ext, read_clust = 'extra_clusters.npz'):
    mcr = micro()
    bin_no = 10
    clust = np.load(save_folder + read_clust)
    cluster = clust['arr_0']

    #dist_spwspike, min_dist_spwspike, fs, max_dists_spwspike  = reader.read_dist_SpwfromSpike(save_folder)
    distances, min_distances, fs, max_dist = reader.read_dist_fromSpike(save_folder)
    
    dist_array_all = []

    for electr in range(len(distances)):
        dist_array = []
        for trace in range(len(distances[electr])):

            v, remove = np.where([np.isinf(distances[electr][trace])])#np.where(dist_spwspike[electr][trace] != float('Inf'))
            
            if len(remove) > 0:
                old = 0
                for r in range(len(remove)):
                    del distances[electr][trace][r-old]
                    todel = r + len(dist_array) - old
                    old = old + 1

                    cluster[electr] = np.delete(cluster[electr],todel)

            dist_array = dist_array + distances[electr][trace]
        dist_array_all.append(dist_array) 
    dist_spwspike = ba.idx2recarr(distances)
    

    for electr in range(len(cluster)):
      
        dist_electrode = dist_spwspike[dist_spwspike['electrode']==electr]
        
        min_time, max_time = 0, 100
        n_bins = 300
        bins = np.linspace(min_time, max_time, n_bins)
        
        n_all, _ = np.histogram(dist_electrode['time'], bins)
        #colors = assign_colors(cluster[electr])
        fig = plt.figure()
        plt.title('Distance of extracellular spikes to any intracellular spike')
        n_clu = len(np.unique(cluster[electr]))
        plt.title('Distribution of spikes in each group')

        for clust in range(len(np.unique(cluster[electr]))):
            plt.subplot(n_clu+1, 1, clust+1)

            n1, _ = np.histogram(dist_electrode[cluster[electr]==clust]['time'], bins)
        
            frac1 = n1*1./n_all
            frac1[np.isnan(frac1)]=0
            left = bins[:-1]
            w = bins[1]-bins[0]
            plt.bar(left, frac1, width=w)
            plt.xlim([0,10])
            plt.ylabel('Fraction')        
        if len(range(len(np.unique(cluster[electr])))) > 0:
            plt.subplot(n_clu+1, 1, n_clu+1)
            left = bins[:-1]
            plt.bar(left, n_all, w)
            plt.xlim([0,10])
            plt.title('Distribution of all Spikes')
            plt.ylabel('Number')
            plt.xlabel('Distance in time [ms]')
            nom = save_folder + 'wyniki/clust_distexSpikes' + str(electr) + ext
            fig.savefig(nom,dpi=600)
            plt.close()
            

def plot_hist_dist_ex2SPW(save_folder, ext, read_clust = 'extra_clusters.npz'):
    mcr = micro()
    bin_no = 10
    clust = np.load(save_folder + read_clust)
    cluster = clust['arr_0']

    #dist_spwspike, min_dist_spwspike, fs, max_dists_spwspike  = reader.read_dist_SpwfromSpike(save_folder)
    distances, min_distances, fs, max_dist  = reader.read_dist_fromSPW(save_folder)
    
    dist_array_all = []

    for electr in range(len(distances)):
        dist_array = []
        for trace in range(len(distances[electr])):

            v, remove = np.where([np.isinf(distances[electr][trace])])#np.where(dist_spwspike[electr][trace] != float('Inf'))
            
            if len(remove) > 0:
                old = 0
                for r in range(len(remove)):
                    distances[electr][trace] = distances[electr][trace][r-old + 1:]
                    todel = r + len(dist_array) - old
                    old = old + 1

                    cluster[electr] = np.delete(cluster[electr],todel)

            dist_array = dist_array + distances[electr][trace]
        dist_array_all.append(dist_array) 
    dist_spwspike = ba.idx2recarr(distances)
    

    for electr in range(len(cluster)):
      
        dist_electrode = dist_spwspike[dist_spwspike['electrode']==electr]
        
        min_time, max_time = -5, 100
        n_bins = 100
        bins = np.linspace(min_time, max_time, n_bins)
        
        n_all, _ = np.histogram(dist_electrode['time'], bins)
        #colors = assign_colors(cluster[electr])
        fig = plt.figure()
        #plt.title('Distance of extracellular spikes to any intracellular spike')
        n_clu = len(np.unique(cluster[electr]))
        tit = 'spikes in each cluster as distance from beginning of SPW - 5 ms, e:' + str(electr)
        plt.title(tit)

        for clust in range(len(np.unique(cluster[electr]))):
            plt.subplot(n_clu+1, 1, clust+1)
            if clust == 0:
                plt.title(tit)
            if len(dist_electrode) > 0:
                n1, _ = np.histogram(dist_electrode[cluster[electr]==clust]['time'], bins)
            
                frac1 = n1*1./n_all
                frac1[np.isnan(frac1)]=0
                left = bins[:-1]
                w = bins[1]-bins[0]
                plt.bar(left, frac1, width=w)
                plt.xlim([-5,30])
                plt.ylabel('Fraction')        
        if len(range(len(np.unique(cluster[electr])))) > 0:
            plt.subplot(n_clu+1, 1, n_clu+1)
            left = bins[:-1]
            plt.bar(left, n_all, w)
            plt.xlim([-5,30])
            plt.title('Distribution of all Spikes')
            plt.ylabel('Number')
            plt.xlabel('Distance in time [ms]')
            nom = save_folder + 'wyniki/clust_distexSpikes_to_SPW' + str(electr) + ext
            fig.savefig(nom,dpi=600)
            plt.close()



        
def plot_extra_clusters(save_folder, save_name = 'extra_spikes_clust', save_clust = 'extra_clusters.npz', save_ext = '.png', win = [-2, 5]):
    npzfile = np.load(save_folder + save_clust)
    clusters = npzfile['arr_0']
    data_bas, fs = reader.read_databas(save_folder)
    spike_idxs, spikes_ampl, fs_spikes  = reader.read_extraspikes(save_folder)
    fs_div = fs/fs_spikes
    #if fs != fs_spikes:
    #    values, fs = updater.update_fs(fs_spikes, fs, values = [spike_idxs])
    #spike_idxs = values[0]
    #print np.shape(spike_idxs)
    spike_idxs = ba.idx2recarr(spike_idxs)
    spike_idxs['time'] = spike_idxs['time']*fs_div

    win_pts = (ispw.ms2pts(win[0], fs), ispw.ms2pts(win[1], fs))
    waves, idx = ba.extract_spws(data_bas, spike_idxs, fs, win_pts) 
    t = dat.get_timeline(waves, fs, 'ms')
    ytit = 'Voltage [' + micro() + 'V]'
    for electr in range(len(clusters)):
        fig = plt.figure()
        tit = 'Spikes and their mean for each cluster (baseline remain to be removed), electr: ' + str(electr)
        plt.title(tit)
        for clust in range(len(np.unique(clusters[electr]))):
            #calculate and plot mean spw in each electrode
            plt.subplot(len(np.unique(clusters[electr])), 1, clust)
            spws_clust = waves[:, idx['electrode']==electr][:,clusters[electr]==np.unique(clusters[electr])[clust]] - waves[:, idx['electrode']==electr][:,clusters[electr]==np.unique(clusters[electr])[clust]].mean(0)
            plt.plot(t, spws_clust, 'b', alpha=0.3)
            plt.plot(t, spws_clust.mean(1), 'r', lw=2)
            plt.ylabel(ytit)
        plt.xlabel('Time [ms]')

        nom = save_folder + 'wyniki/' + save_name + str(electr) + save_ext
        fig.savefig(nom)
        plt.close()
    
    

def plot_extra_clusters_vs_distance(save_folder, save_distSPW = 'extra_spikes_clustSPW',  save_distSPIKE = 'extra_spikes_clustSPIKE', save_clust = 'extra_clusters.npz', save_ext = '.png'):
    npzfile = np.load(save_folder + save_clust)
    clusters = npzfile['arr_0']
    
        

def cluster_corrs(save_folder, read_corr = 'corr_electr.npz', save_clust = 'clust_electr.npz', n_clusters = 5):
    
    corrs = np.load(save_folder + read_corr)
    all_corr = corrs['arr_0']

    s_all = []
    cluster_all = []
    medioids_all = []
    for electr in range(len(all_corr)):
        if np.size( all_corr[electr]) > 1:
            cluster, medoids = clust.kmedoids(n_clusters, all_corr[electr], nsteps=1000)
        elif np.size(all_corr[electr]) == 1:
            cluster = [0]
            medoids = [0]
        else:
            cluster = []
            medoids = []  
        
        s_new = clust.reorder_matrix(cluster,all_corr[electr])
        s_all.append(s_new)
        cluster_all.append(cluster)
        medioids_all.append(medoids)
    np.savez(save_folder + save_clust, cluster_all, medioids_all, s_all)   
    
def plot_data_spikes_spws(save_folder, ext):
    data_bas, fs_bas = reader.read_databas(save_folder)
    spw_idxs, spw_maxs, starts_spw, ends_spw, lengths_spw, fs = reader.read_SPWs(save_folder)
    #spike_idxs, spikes_ampl, fs_new = reader.read_extraspikes(save_folder)
    data_intra, fs_intra = reader.read_datafile_intra(save_folder)

    scale = 'ms'
    
    t = dat.get_timeline(data_intra[0][0], fs, scale)
    
    fig = plt.figure()
    plt.subplot(len(data_bas)+1, 1, 1)
    plt.title('All the electrode data (detected SPWs in green)')
    for trace in range(len(data_bas[0])):
        plt.plot(t, data_intra[0][trace])
        
        plt.ylabel('Voltage [mV]')
    plt.xlim([0, max(t)])
    m = micro()
    ytit = 'voltage [' + m + 'V]'
    
    t = dat.get_timeline(data_bas[0][0], fs_bas, scale)
    for electr in range(len(data_bas)):
        plt.subplot(len(data_bas)+1, 1, electr + 2)
        for trace in range(len(data_bas[electr])):
            plt.plot(t, data_bas[electr][trace])
            plt.plot(t[starts_spw[electr][trace]], data_bas[electr][trace][starts_spw[electr][trace]], 'go')
        plt.ylabel(ytit)
        plt.xlim([0, max(t)])
    plt.xlabel('time [ms]')
    nom = save_folder + 'wyniki/normalne_dane' + ext
    fig.savefig(nom,dpi=600)
    plt.close()
    
def assign_colors(groups):
    no_colors = len(np.unique(groups))
    cols = define_colors(no_colors)
    colors = []
    for i in range(len(np.unique(groups))):
        colors.append(cols[i])
        
    return colors

def group_mean(save_folder, ext, read_clust = 'clust_electr.npz', read_idx = 'corr_electr.npz'):
    path = '/home/maja/PhDProject/SPWs/SPWs/saved_data/cell6/gap_free'
    data, fs, idx = ba.load_data(save_folder)
    
    win = (-10, 30)
    win_pts = (ispw.ms2pts(win[0], fs), ispw.ms2pts(win[1], fs))
    corrs = np.load(save_folder + read_idx)

    idx_aligned_all = corrs['arr_1']

    spws, idx = ba.extract_spws(data, idx, fs, win_pts)

    clust= np.load(save_folder + read_clust)
    cluster = clust['arr_0']
    
    #print spws
    ytit = 'Voltage [' + micro() + 'V]'
    for electr in range(len(cluster)):
        fig = plt.figure()
        
        for clust in range(len(np.unique(cluster[electr]))):
            #calculate and plot mean spw in each electrode
            plt.subplot(len(np.unique(cluster[electr])), 1, clust)
            spws_clust = spws[:, idx['electrode']==electr][:,cluster[electr]==np.unique(cluster[electr])[clust]] 
            plt.plot(spws_clust, 'b', alpha=0.3)
            plt.plot(spws_clust.mean(1), 'r', lw=2)
            plt.ylabel(ytit)
        plt.xlabel('Time [ms]')
        nom = save_folder + 'wyniki/mean_clusters' + str(electr) + ext
        fig.savefig(nom)
        plt.close()
        

def plot_corr_with_all(save_folder, ext, read_clust = 'clust_electr.npz'):
    mcr = micro()
    bin_no = 10
    clust = np.load(save_folder + read_clust)
    cluster = clust['arr_0']

    dist_spwspike, min_dist_spwspike, fs, max_dists_spwspike  = reader.read_dist_SpwfromSpike(save_folder)
    dist_array_all = []



    for electr in range(len(dist_spwspike)):
        dist_array = []
        for trace in range(len(dist_spwspike[electr])):

            v, remove = np.where([np.isinf(dist_spwspike[electr][trace])])#np.where(dist_spwspike[electr][trace] != float('Inf'))
            
            if len(remove) > 0:
                old = 0
                for r in range(len(remove)):
                    del dist_spwspike[electr][trace][r-old]
                    todel = r + len(dist_array) - old
                    old = old + 1

                    cluster[electr] = np.delete(cluster[electr],todel)

            dist_array = dist_array + dist_spwspike[electr][trace]
        dist_array_all.append(dist_array) 
    dist_spwspike = ba.idx2recarr(dist_spwspike)
    
    #print dist_spwspike[dist_spwspike['electrode']==0].shape
    #print cluster[0].shape
    #print cluster
    plt.figure()
    for electr in range(len(cluster)):
        #        plt.figure()
        #plt.plot(spw_maxs[electr], cluster , 'o')
#        title = 'amplitudes ' + str(electr)
#        title_x = 'spw amplitudes'
#        title_y = 'correlation clusters'
#        bin_no = 30
#        plot_hist(spw_maxs[electr], cluster, title_x, title_y, 'voltage [' + mcr + 'V]', 'number', bins = bin_no)
#        #plt.hist(lengths_spw[electr])
#        #plt.plot(lengths_spw[electr], cluster , 'o')
#        title = 'widths'
#        title_x = 'spw widths'
#        title_y = 'correlation clusters'        
#        plot_hist(lengths_spw[electr], cluster, title_x, title_y, 'time [ms]', 'number', bins = bin_no)
        #plt.title('widths')
        
        #fig = plt.figure()
        
        #nom = save_folder + 'wyniki/big_fig' + str(electr) + '.png'
        #fig.savefig(nom,dpi=600)
        
        #print cluster[electr]
        dist_electrode = dist_spwspike[dist_spwspike['electrode']==electr]

        #plt.scatter(dist_electrode['time'],  cluster[electr],  c=colors)
        
        min_time, max_time = 0, 100
        n_bins = 100
        bins = np.linspace(min_time, max_time, n_bins)
        
        n_all, _ = np.histogram(dist_electrode['time'], bins)
        colors = assign_colors(cluster[electr])
        fig = plt.figure()
        n_clu = len(np.unique(cluster[electr]))
        plt.title('Distribution of Spws in each group')
        
        
        
        

        
        for clust in range(len(np.unique(cluster[electr]))):
            plt.subplot(n_clu+1, 1, clust+1)

            n1, _ = np.histogram(dist_electrode[cluster[electr]==clust]['time'], bins)
        
            frac1 = n1*1./n_all
            frac1[np.isnan(frac1)]=0
            left = bins[:-1]
            w = bins[1]-bins[0]
            plt.bar(left, frac1, width=w)
            plt.xlim([0,50])
            plt.ylabel('Fraction')        
        if len(range(len(np.unique(cluster[electr])))) > 0:
            plt.subplot(n_clu+1, 1, n_clu+1)
            left = bins[:-1]
            plt.bar(left, n_all, w)
            plt.xlim([0,50])
            plt.title('Distribution of all SPWs')
            plt.ylabel('Number')
            plt.xlabel('Distance in time [ms]')
            nom = save_folder + 'wyniki/clust_distSPW' + str(electr) + ext
            fig.savefig(nom,dpi=600)
            plt.close()
            
            
            
        
#        fig = plt.figure()
#        plt.title ='distance from spike'
#        title_x = 'spw distances from spikes'
#        title_y = 'correlation clusters'
#        plot_hist(dist_spwspike[dist_spwspike['electrode']==electr]['time'], cluster[electr], title_x, title_y, 'time [ms]', 'number', bins = bin_no)
#        nom = save_folder + 'wyniki/big_fig' + str(electr) + '.png'
#        fig.savefig(nom)
#        
#        plot_group_hist(dist_spwspike[dist_spwspike['electrode']==electr]['time'], cluster[electr])
        
        
        
        
        
        
#        title ='distance to widths'
#        title_x = 'spw distances from spikes'
#        title_y = 'spw widths'       
#        plot_hist(dist_spwspike[electr], lengths_spw[electr], title_x, title_y, 'time [ms]', 'time [ms]', bins = bin_no)
#        
#        title ='amplitude and widths'
#        title_x = 'spw amplitudes'
#        title_y = 'spw widths'       
#        plot_hist(spw_maxs[electr], lengths_spw[electr], title_x, title_y, 'voltage [' + mcr + 'V]', 'time [ms]', bins = bin_no)   
    
    
    
    #data_bas, fs_new = reader.read_databas(save_folder)
    #starting_electrode, mean_corr = reader.read_corr_spws_electrwise(save_folder)
    
    #scale = 'ms'
    #t = dat.get_timeline(data_bas[0], fs, scale)
    
    #print np.size(dist_spwspike[0])
    #print len(starting_electrode)
    
    #colors = define_colors(len(min_dist_spwspike))
    #plt.figure()
    #for spw in range(len(starting_electrode)):
    #    electrode_used = starting_electrode[spw]
    #    cl = colors[electrode_used]
        
        #print dist_spwspike[electrode_used][spw]
    #    plt.plot(dist_spwspike[electrode_used][spw], mean_corr[spw] , 'o', color = cl)
    
    #plt.legend(['1', '2', '3', '4', '5', '6', '7'])
    #plt.legend() 


def plot_group_hist(distances, clusters):
    un = np.unique(clusters)
    
    for group in range(len(un)):
        plt.subplot(len(un), 1, group)
        plt.plot(distances[clusters == group])
        #print distances[clusters == group]
    #plt.show()


def plot_hist(x, y, title_x, title_y, xlab = 'Time [ms]', ylab = 'Voltage [mV]', fig = 0, bins = 10):
    y_lim = [0, 5]
    #    import numpy as np
    #import matplotlib.pyplot as plt
    #from matplotlib.ticker import NullFormatter
    
    # the random data
    #x = np.random.randn(1000)
    #y = np.random.randn(1000)
    
    nullfmt   = NullFormatter()         # no labels

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left+width+0.02

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    # start with a rectangular Figure
    if fig == 0:
        fig = plt.figure(figsize=(8,8))
    else: 
        fig = plt.figure(fig, figsize=(8,8)) 
    axScatter = plt.axes(rect_scatter)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)
    #axHistx.set_ylabel(title_x)
    axHistx.set_title(title_x)
    axHisty.set_title(title_y)
    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    
    axHisty.yaxis.set_major_formatter(nullfmt)

    # the scatter plot:
    axScatter.scatter(x, y, alpha = 0.03)
    axScatter.set_xlabel(xlab)
    axScatter.set_ylabel(ylab)
    axScatter.set_xlim([min(x), max(x)])
    #axScatter.set_ylim([max(y), max(y)])
    axScatter.set_ylim(y_lim)
    

    # now determine nice limits by hand:
    #binwidth = 0.25
    xymax = np.max( [np.max(np.fabs(x)), np.max(np.fabs(y))] )
    #lim = ( int(xymax/binwidth) + 1) * binwidth

    #axScatter.set_xlim( (-lim, lim) )
    #axScatter.set_ylim( (-lim, lim) )

    #bins = np.arange(-lim, lim + binwidth, binwidth)
    axHistx.hist(x, bins)
    axHisty.hist(y, bins, orientation='horizontal')
    #axHistx.hist(x, bins=bins)
    #axHisty.hist(y, bins=bins, orientation='horizontal')

    axHistx.set_xlim([0.1, max(x)])# axScatter.get_xlim() )
    axHisty.set_ylim(y_lim)
    #axHisty.set_ylim(y_lim)# axScatter.get_ylim() )
    return fig

#        for h in range(n_clusters):
#            plt.subplot(n_clusters, 1, h)
#            
#            right_spws = [dist_spwspike[electr][i] for i in range(len(cluster)) if cluster[i] == h]
#            right_spws = [right_spws[i] for i in range(len(right_spws)) if right_spws[i]>-1]
#            
#            if h == 0:
#                group0 = group0 + right_spws
#            elif h == 1:
#                group1 = group1 + right_spws
#            elif h == 2:
#                group2 = group2 + right_spws                

    #for electr in range(len(all_corr)):
        #for h in range(n_clusters):
#    max_t = max(group1+group0+group2)
#    bins = np.linspace(0, max_t, 50)   
    


    



#def correlate_spw(save_folder, plot_it = 1, save_file = 'corr_SPWtimewise'):
#    """ the spw will be correlated, and possibly results will be plotted"""
#    use_before = 30 # ms
#    use_after = 40 # ms
#    
#    # read spws and data
#    #data_fast, freq, fs_new = read_filtered_fast(save_folder)
#    data_bas, fs_new = reader.read_databas(save_folder)
#    #data_spw, freq, fs_new = read_filtered_spw(save_folder)
#    spw_idxs, spw_maxs, starts_spw, ends_spw, lengths_spw, fs = reader.read_SPWs(save_folder)
#
#    
#    
#    all_corr = []
#    all_move = []
#    all_idxs_new = []
#    t = dat.get_timeline(data_bas[0][0], fs_new,'ms')
#    plt.plot(t, data_bas[0][0])
#
#    for electr in range(np.size(spw_idxs,0)): 
#        for trace in range(np.size(spw_idxs[electr])):
#            if plot_it == 1:
#                plt.figure()
#            # for every electrode
#
#            electr_corr, electr_move = corr_two(data_bas[electr][trace], fs, spw_idxs[electr][trace], spw_idxs[electr][trace], use_before, use_after)
#    
#            if plot_it == 1:
#                plt.imshow(electr_corr,interpolation='nearest',origin='lower')
#                plt.colorbar()
#        all_corr.append(electr_corr)
#        all_move.append(electr_move)
#        
#    np.savez(save_folder + save_file, all_corr, all_move) 
#    return all_corr, all_move




def plot_data_all(save_folder, electrode = 'all', save_name = 'all_data', save_ext = '.png'):
    """ plots the whole data"""
    data_all, fs = reader.read_datafile(save_folder, save_file = "data.npz")
    plot_data(data_all, fs, save_folder, electrode = 'all', save_name = save_name, save_ext =save_ext)
    
def plot_data_bas(save_folder, electrode = 'all', save_name = 'data_bas', save_ext = '.png'):
    """ plots the baselined data"""
    data_bas, fs = reader.read_databas(save_folder)
    plot_data(data_bas, fs, save_folder, electrode = 'all', save_name = save_name, save_ext = save_ext)

def plot_spw_data(save_folder, electrode = 'all', save_name = 'data_bas', save_ext = '.png'):
    """ plots the spw data"""
    data_spw, freq, fs = reader.read_filtered_spw(save_folder, save_file = "spw_data.npz")
    plot_data(data_spw, fs, save_folder, electrode = 'all', save_name = save_name, save_ext = save_ext)

def plot_ripple_data(save_folder, electrode = 'all', save_name = 'data_bas', save_ext = '.png'):
    """ plots the ripple data"""
    data_ripple, freq, fs =  reader.read_filtered_ripple(save_folder, save_file = "ripple_data.npz")
    plot_data(data_ripple, fs, save_folder, electrode = 'all', save_name = save_name, save_ext = save_ext)

def plot_fast_data(save_folder, electrode = 'all', save_name = 'fast_data', save_ext = '.png'):
    """ plots the fast_data data"""
    data_fast, freq, fs = reader.read_filtered_fast(save_folder, save_file = "fast_data.npz")
    plot_data(data_fast, fs, save_folder, electrode = 'all', save_name = save_name, save_ext = save_ext)









def plot_data(data, fs, save_folder, electrode = 'all', save_name = 'data', save_ext = '.png'):
    scale = 'ms'
    add = 'wyniki/'
    subplots = 0 # if subplots == 0, then each plot is on different figure
    save_names = []
    if electrode == 'all':
        # plot all the electrodes
        electrode = range(len(data))
    else:
        electrode = electrode - 1
    fig = plt.figure()
    t = dat.get_timeline(data[0][0], fs, scale)
    
    for electr in electrode:
        # plot all or chosen electrodes
        #for trace in range(len(data[electr])):
        if subplots == 1:
            plt.subplot(len(electrode), 1, electr+1)
        else:
            if electr != electrode[0]:
                save_nom = save_folder + add + save_name + str(electr) + save_ext
                fig.savefig(save_nom,dpi=600)
                plt.close()
                save_names.append(save_name)
            fig = plt.figure()
        
        plt.plot(t, np.transpose(data[electr]))
        tit = save_name + ', electrode: ' + str(electr + 1)
        plt.title(tit)
        plt.xlim([0, max(t)])
        #plt.ylim([min(data[electr]) - 20, max(data[electr]) + 20])
        plt.xlabel('time [ms]')
        vl = 'voltage [' + micro() + 'V]'
        plt.ylabel(vl)
            
            
    if subplots == 1:
        save_nom = save_folder + save_name + str(electr) + save_ext
        fig.savefig(save_nom,dpi=600)
        plt.close()
        save_names.append(save_nom)
    else:
        save_nom = save_folder + add + save_name + str(len(electrode) + 1) + save_ext
        fig.savefig(save_nom,dpi=600)
        plt.close()
        save_names.append(save_name)
    return save_names



#
#def work_on_swp_spikes(save_folder):
#    dist_fromSpw, min_dist_fromSpw, fs, max_dist_fromSpw = reader.read_dist_fromSPW(save_folder)
#    spike_idxs, spikes_ampl, fs_new  = reader.read_extraspikes(save_folder)
#    dist_fromSpike, min_dist, fs, max_dist = reader.read_dist_fromSpike(save_folder)
#    dist_spwspike, min_dist_spwspike, fs, max_dists_spwspike  = reader.read_dist_SpwfromSpike(save_folder)
#    spw_idxs, spw_maxs, starts_spw, ends_spw, lengths_spw, fs  = reader.read_SPWs(save_folder)
#    spike_idxs, all_valley_to_peak, all_half_valley_width, all_half_peak_width, fs_new = reader.read_ex_sparamas(save_folder)
#    sp_idx_first_intra, sp_idx_all_intra, fs = reader.read_intra_spikes(save_folder)
#    data_bas, fs_new = reader.read_databas(save_folder)
#
#    
#
#    max_dist_pts = ispw.ms2pts(max_dist, fs)
#    min_dist_pts = ispw.ms2pts(0.0, fs)
#    scale = 'ms'
#    t = dat.get_timeline(data_bas[0], fs_new, scale)
#    add_up = 0
#    #print dist_fromSpike
##    for electr in range(len(spike_idxs_extra)):
##        inter = [i for i in range(len(dist_fromSpike[electr])) if dist_fromSpike[electr][i] < max_dist_pts and dist_fromSpike[electr][i] > min_dist_pts]
##        
##        plt.figure(1)
##        #plt.plot(t)
##        add_up = add_up - 100
##        plt.plot(t, data_bas[electr]+add_up)
##        plt.plot(t[starts_spw[electr]], data_bas[electr][starts_spw[electr]]+add_up, 'o')
##        plt.xlim([0,10000])
#    
#        # correlate between the electrodes
#    before = 40
#    after = 60
#    before_pts = ispw.ms2pts(before, fs)
#    after_pts = ispw.ms2pts(after, fs)
#    
#    electr_corr = [] #np.ones((len(one), len(two)))
#    electr_move = [] #np.ones((len(one), len(two)))
#    # first chcek which electrodes begins everything and then calculate correlation from this electrode
#    chosen = 0
#    a = 0
#    starting_electrode = []
#    mean_corr = []
#    trace = 0
#    data_file = 'spw_corr'
#    
#    for spws in range(len(starts_spw[0])):
#        #print spws
#        #print starts_spw
#        #print np.shape(starts_spw)
#        temp_start = starts_spw[chosen][0][spws]
#        for check in range(2):
#            #plt.figure() 
#            electr_corr = []
#            electr_move = []
#            for spw1 in range(len(spike_idxs_extra)):
#                #print chosen
#                #print a
#                #print starts_spw[chosen][0]
#                s1 = data_bas[chosen][0][temp_start-before_pts:temp_start+after_pts]
#                s2 = data_bas[spw1][0][temp_start-before_pts:temp_start+after_pts]
#                s1 = s1 - np.mean(s1)
#                s2 = s2 - np.mean(s2)
#                s12conv = np.convolve(s1, s2[::-1])
#                s12conv = s12conv/np.std(s1) /np.std(s2)
#                corr_best = max(s12conv)
#                corr_idx = np.where(s12conv == corr_best)[0][0]
#                corr_best = corr_best /len(s1)
#    
#                #corr_idx = s12conv[corr_best]
#                electr_corr.append(corr_best)
#                electr_move.append(len(s1) - corr_idx - 1)
##                if spw1 == chosen:
##                    plt.plot(t[starts_spw[chosen][0]-before_pts:starts_spw[chosen][0]+after_pts], data_bas[spw1][starts_spw[chosen][0]-before_pts:starts_spw[chosen][0]+after_pts])
##                    plt.plot(t[starts_spw[chosen][0]], data_bas[spw1][starts_spw[chosen][0]], 'ro')
##                    plt.plot(t[a], data_bas[spw1][a], 'go')                       
#            #print electr_corr    
#            #print electr_move
#            if check == 1:
#                starting_electrode.append(chosen)
#                mean_corr.append(np.mean(electr_corr))     
#            
#            chosen = electr_move.index(min(electr_move))
#            temp_start = starts_spw[chosen][0][spws] - electr_move[chosen]
#
#    np.savez(save_folder + data_file, starting_electrode, mean_corr) 



          
            
    #return electr_corr, electr_move    
#        plt.figure()
#        plt.plot(all_valley_to_peak[electr], all_half_peak_width[electr], 'bo', alpha = 0.1)
#        plt.plot(all_valley_to_peak[electr][inter], all_half_peak_width[electr][inter], 'ro', alpha = 0.5)
#        plt.title(electr)
#        plt.xlabel('all_valley_to_peak [ms]')
#        plt.ylabel('all_half_peak_width')
        
        

        #plt.plot(dist_fromSpike[electr], dist_fromSpw[electr])
        
#        plt.figure()
#        plt.plot(starts_spw[electr])
#        plt.plot(spike_idxs_extra[electr])
#        plt.plot(sp_idx_all_intra)
  
        
        
        #print np.shape(distances_fromSpw[electr])
        #print np.shape(min_distances_fromSpw[electr])  
        #print np.shape(spike_idxs[electr])
        #print np.shape(distances[electr])
        #print np.shape(dist_spwspike[electr])





def plot_extraspikes(save_folder):
    on_top = 0 # if on_top == 1 -> all the traces are plotted on the top of each other, otherwise
    # they are plotted in the row, with horizontal lines in between
    data_bas, fs = reader.read_databas(save_folder, save_file = "data_bas.npz")
    spike_idxs, spikes_ampl, fs_new  = reader.read_extraspikes(save_folder, save_file = "ex_spikes.npz")
    spike_idxs = updater.update_fs(fs_new, fs, values = [spike_idxs])
    spike_idxs = spike_idxs[0]
    
    scale = 'ms'
    t = dat.get_timeline(data_bas[0][0], fs, scale)

   
    for electr in range(len(data_bas)):
        plt.figure()
        h_line = []
        data = []
        spiks = []
        for trace in range(len(data_bas[electr])):
            if on_top == 0:
                #print data_bas[electr][trace].tolist()
                data = data + [data_bas[electr][trace][i] for i in range(len(data_bas[electr][trace]))]
                spiks = spiks + spike_idxs[electr][trace] + len(data)
                h_line.append(len(data))
            else:
                plt.plot(t, data_bas[electr][trace])
                plt.plot(t[spike_idxs[electr][trace]], data_bas[electr][trace][spike_idxs[electr][trace]], 'go')
        if on_top == 0:
            t = dat.get_timeline(data[0][0], fs, scale)
            plt.plot(t, data)
            plt.plot(t[spiks], data[spiks],'go')
            plt.vlines(h_line, 300, -100, color = 'r')


    wyniki = save_folder + 'wyniki/'
    plt.savefig('simple_plot',dpi=600)
    plt.close()



def asa():
    # check if fs - spikes and fs - data is the same
    #plt.plot(t[spike_idxs], data[spike_idxs])
    #ispw.update_extraspikes(data, fs, save_folder)
    
    #ispw.update_expikes_params(data, fs, save_folder, spike_idxs)
    #for electr in electrods:
        #half_widths, teils =fes.find_halfampl(data[electr], fs, spike_idxs)
    data_all, fs = reader.read_datafile(save_folder)
    
    #data_uspl_extra, fs_new_up = read_upsample(save_folder)
    #data_uspl_extra, fs_new_up = ispw.update_upsample(data_all, fs, save_folder, uspl = 2)
    
    spike_idxs, spikes_ampl, whole_spikes, fs_new = ispw.update_extraspikes(data_all, fs, save_folder)
    #spike_idxs, spikes_ampl, whole_spikes, fs_new  = read_extraspikes(save_folder)

    
    
    scale = 'ms'
    t_up = dat.get_timeline(data_all[0], fs, scale)
    t = dat.get_timeline(data_all[0], fs, scale) 
    
    for electr in range(len(data_all)):
        #half_widths, teils =fes.find_halfampl(data[electr], fs, spike_idxs) 
        plt.figure()
        plt.plot(t, data_all[electr])
        #plt.plot(t_up, data_uspl_extra[electr])

        #for o in spike_idxs[electr]:
        #    plt.plot(t_up[o], data_uspl_extra[electr][o], 'o')
        
        


       





   



        



def work_on_spikes(save_folder):
    npzfile = np.load(save_folder + "data.npz")
    #npzfile2 = np.load(save_folder + "fast_data.npz")
    #npzfile3 = np.load(save_folder + "spw_data.npz")
    data = npzfile['arr_0']
    fs = npzfile['arr_1']

    
    #fast_data = npzfile2['arr_0']
    #spw_data = npzfile3['arr_0']
##
    #ispw.update_extraspikes(data, fs, save_folder)
    #ispw.update_events(data, fast_data, spw_data, fs, save_folder)
    #------- 
    # load the extracellular spikes and intracellular spikes
    ex_spikes = "ex_spikes.npz"
    in_spikes = "events.npz"

    # load extracellular spike details
    # spikes, spikes_ampl, half_widths_all, teils_all, fs
    npzfile = np.load(save_folder + ex_spikes)
    ext_spikes = npzfile['arr_0']
    ext_spikeshalf = npzfile['arr_2']
    #fs = npzfile['arr_4']

    #sp_idx_first, sp_idx_all, all_maxs, events_ampl, all_spws, fs
    # load intracellular spike details
    # sp_idx_first, sp_idx_all, all_maxs, events_ampl, fs
    npzfile = np.load(save_folder + in_spikes)
    #int_spikes1 = npzfile['arr_0']
    int_spikes = npzfile['arr_1']
    all_spws = npzfile['arr_4']
    spws_idx = npzfile['arr_2']
    #find the distances from each event to previous intra-spike
    plt.figure(4)
    scale = 'sec'
    t = dat.get_timeline(data[0], fs, scale)
    plt.subplot(np.size(ext_spikes,0)+1, 1, 1)
    plt.plot(t, data[0])
    if len(int_spikes) > 0:
        plt.plot(t[int_spikes], data[0][int_spikes], 'ro')
    plt.xlim([0, t[-1]])
    plt.ylabel('intracellular')
    rang = [1,2] # ms
    rang1 = ispw.ms2pts(rang[0],fs)
    rang2 = ispw.ms2pts(rang[1],fs)
    for electr in range(np.size(ext_spikes,0)):
        plt.figure(1)
        #np.where(distance <= 0)
        distance = ispw.dist_fromSpike(int_spikes, ext_spikes[electr])
        max_dist = 3 # ms
        max_distpts = ispw.ms2pts(max_dist, fs)
        
        #print np.where(distance <= 0)[0]
        #distance[distance[:] < 0] = max_distpts+1
        
        inter = np.where(distance<=max_distpts)[0]
        #print 'here'
        #print np.shape(inter)
        #print inter
        plt.subplot(np.size(ext_spikes,0), 1, electr)

        
        for spik in inter:
            idx1 = ext_spikes[electr][spik] - rang1
            idx2 = ext_spikes[electr][spik] + rang2
            plt.plot(data[electr+1][idx1:idx2])
        
        if len(inter) > 0:
            plt.figure(2)
            plt.subplot(np.size(ext_spikes,0), 1, electr)   
            ext_spikeshalf[electr] = ispw.pts2ms(ext_spikeshalf[electr], fs)
            plt.hist([ext_spikeshalf[electr][i] for i in inter])
            plt.ylabel(electr+1)
        
        plt.figure(3)
        plt.subplot(np.size(ext_spikes,0), 1, electr)  
        plt.hist(ext_spikeshalf[electr])
        plt.ylabel(electr+1)
        
        plt.figure(4)
        plt.subplot(np.size(ext_spikes,0)+1, 1, electr+2)  
        
        plt.plot(t, data[electr+1])        
        plt.plot(t[ext_spikes[electr]], data[electr+1][ext_spikes[electr]], 'go')
#        print 'second round'
#        print electr
#        print len(ext_spikes[electr])
        plt.xlim([0, t[-1]])
        plt.title(electr+1)
        plt.ylabel(len(inter))
        if len(inter) > 0:
            [plt.plot(t[ext_spikes[electr][i]],data[electr+1][ext_spikes[electr][i]], 'ro') for i in inter]
    plt.xlabel('times '+ scale)   
    #distance = ispw.pts2ms(distance, fs2)
    #distances.append(distance)


    
def group_corr(save_folder):
    all_corr, all_move = read_corr_SPWtimewise(save_folder)
    dist_spwspike, min_dist_spwspike, fs, max_dists_spwspike  = read_dist_SpwfromSpike(save_folder)
    spw_idxs, spw_maxs, starts_spw, ends_spw, lengths_spw, fs = read_SPWs(save_folder)   
    
    # data_bas, fs_new = read_databas(save_folder)
    #starting_electrode, mean_corr = read_corr_spws_electrwise(save_folder)
    mcr = micro()
    n_clusters = 3
    cols = define_colors(n_clusters)

    np.savez('korelacje', all_corr)
    group0 = []
    group1 = []
    group2 = []
    for electr in range(len(all_corr)):
        print 'clustering'
        plt.imshow(all_corr[electr],interpolation='nearest',origin='lower')
        plt.colorbar()
        plt.title('old')
        
        cluster, medoids = clust.kmedoids(n_clusters, all_corr[electr], nsteps=1000)
        s_new = clust.reorder_matrix(cluster,all_corr[electr])
#        plt.figure()
#        plt.imshow(s_new,interpolation='nearest',origin='lower')
#        plt.colorbar()
#        plt.title('new')
        for h in range(n_clusters):
            plt.subplot(n_clusters, 1, h)
            
            right_spws = [dist_spwspike[electr][i] for i in range(len(cluster)) if cluster[i] == h]
            right_spws = [right_spws[i] for i in range(len(right_spws)) if right_spws[i]>-1]
            
            if h == 0:
                group0 = group0 + right_spws
            elif h == 1:
                group1 = group1 + right_spws
            elif h == 2:
                group2 = group2 + right_spws                

    #for electr in range(len(all_corr)):
        #for h in range(n_clusters):
    max_t = max(group1+group0+group2)
    bins = np.linspace(0, max_t, 50)   


def group_spws(save_folder):
    spw_idxs, spw_maxs, starts_spw, ends_spw, lengths_spw, fs = reader.read_SPWs(save_folder)
    for electr in range(np.size(spw_idxs,0)): 
        lengths_spw_ms = [ispw.pts2ms(lengths_spw[electr][i], fs) for i in range(len(lengths_spw[electr]))]
        plt.figure()
        plt.plot(spw_maxs[electr], lengths_spw_ms, 'o')
        plt.xlabel('max amplitude [mV]')
        plt.ylabel('length at mean [ms]')







        
def find_nearest(array,value):
    idx=(np.abs(array-value)).argmin()
    return idx, array[idx]





def max_between_spikes(save_folder,  load_chosen  = 'chosen_idx_val.npz', save_name = 'spike_maxs.npz', save_ext = '.png'):
    # read data_all
    # read used spikes
    data_all, fs = reader.read_datafile(save_folder)
    
    npzfile = np.load(save_folder + load_chosen)
    idxs = npzfile['sp_idx']
    #a_s =npzfile['a_s']
    #b_s =npzfile['b_s']
    #ampls = npzfile['ampls']
    #max_left = npzfile['max_left']
    #max_right = npzfile['max_right']
    #norm_fact = npzfile['norm_fact']
    
    spike_idxs = ba.idx2recarr(np.array(idxs), time_dtype=np.float32)
    #print idxs
    max_ipsp = 5 # ms
    max_ipsp = ispw.ms2pts(10, fs)
    
    #a_s = ba.idx2recarr(np.array(a_s), time_dtype=np.float32)
    #b_s = ba.idx2recarr(np.array(b_s), time_dtype=np.float32)
    
    every_max = []
    for trace in range(len(data_all[0])):
        only_this_trace = spike_idxs[(spike_idxs['trace'] == trace)]
        sort_idx = np.argsort(only_this_trace['time'])
        posortowane = only_this_trace[:,:,sort_idx]
        old_spike = -1
        all_maxe = []
        #np.array()
        for (i, spike) in enumerate(posortowane):    

#            if i < len(posortowane)-1:
#                print posortowane[i + 1]
            if old_spike == spike:
                print 'double'
                
            else:
                old_spike = spike
                if i == 0:

                    # if it's the first spike
                    
                    start = int(spike['time']) - max_ipsp
                    if start < 0:
                        start = 0
                    end = int(spike['time'])
    
                elif i == len(posortowane)-1:
                    # if it's the last spike
                    start = int(spike['time'])
                    
                    end = int(spike['time']) + max_ipsp #int(spike['time'])
                    
                    if end > len(data_all[0][0]):
                        end = len(data_all[0][0])
                else:

                    start = int(spike['time'])

                    next_spike = int(posortowane[i + 1]['time'])
                    add1 = max_ipsp #int(spike['time'])
                    add2 = np.abs(int(spike['time']) - next_spike)
                    
                    if add2 ==0:
                        
                        next_spike = int(posortowane[i + 2]['time'])
                        add2 = np.abs(int(spike['time']) - next_spike)
                        
                    
                    end = int(spike['time']) + min(add1, add2)
                    if end > len(data_all[0][0]):
                        end = len(data_all[0][0])
                maxe = []
                for electr in range(len(data_all)):
                    maxe.append(np.max(data_all[electr][trace][start:end]))
                if np.size(all_maxe == 0) or min(add1,add2) == max_ipsp:
                    all_maxe.append(maxe)
                else:
                    all_maxe.append(maxe - all_maxe[-1])
                
        every_max.append(all_maxe)  
    np.savez(save_folder + save_name, every_max)     
    
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! print spike_idxs[(spike_idxs['electrode'] == 0)  & (spike_idxs['trace'] == 0)]

def plot_extra_params_novel(save_folder, save_name = 'chosen_spikes.npz', ext = '.png'):
    """ finds in which of the electrodes spike has the largest amplitude and possibly plot the traces"""
    
    # load the parameters of extracellular spikes
    plot_it     = 0
    allow_shift = 0.15 # ms
    win         = (-5, 5)
    spike_idxs, a_s, b_s, ampls, fs_spikes  =   reader.read_ex_sparamas(save_folder)
    data_all,fs = reader.read_datafile(save_folder)
    
    # transform all the ms data to pts
    allow_shift = ispw.ms2pts(allow_shift, fs)
    before      = ispw.ms2pts(win[0],fs)
    after       = ispw.ms2pts(win[1],fs)
    
    # check if the samplint rantes are compatible - if not write the message
    if fs != fs_spikes:
        print 'warning: the two sampling rates are not compatible!!'
    
    # transform chosen variables to rec array
    spike_idxs  = ba.idx2recarr(spike_idxs)
   
    # initiate variables
    winning_el, winning_idx, winning_trace, nearest_idxs = [], [], [], []
    
    # use every trace 
    for trace in range(len(data_all[0])):
        spike_old = -1
        
        # sort the data for this trace
        only_this_trace = spike_idxs[(spike_idxs['trace'] == trace)]
        sort_idx        = np.argsort(only_this_trace['time'])
        posortowane     = only_this_trace[:,:,sort_idx]
        
        # use old to not use the same spike twice
        old     = -1 # initiate old to be less then idx
        
        # loop trough all the detected spikes sorted timewise
        for (idx, spike) in enumerate(posortowane):
            
            # check if index is not repeating
            if spike['time'] != spike_old:
                spike_old   = spike['time']
                nearest_pts, nearest_el, nearest_idxs = [], [], []
                
                # check if old is less than idx
                if idx >= old:
                    el, tr, sp  = spike['electrode'], spike['trace'], spike['time']
                    add_it      = 0
                
                    plt.clf()
                    # check if spike is not going outside given data
                    if sp > (-before) and sp < (len(data_all[el][tr]) - after):
                        # loop through all the electrodes
                        for electr in range(len(data_all)):
                            
                            # is the trace the same?
                            if tr != trace:
                                print 'error' 
                                
                            # plot the spikes (multiple electrodes                         
                            if plot_it == 1:
                                # calculate timeline for this data
                                t = dat.get_timeline(data_all[electr][tr][before + 1000: 1000+ after], fs, 'ms') + sp - before 
                                data_plot = data_all[electr][tr][sp + before:sp+after]
                                plt.plot(t, data_plot + add_it)
                            
                            # find the spikes for this electrode and this trace
                            other_spikes = spike_idxs[(spike_idxs['electrode'] == electr)  & (spike_idxs['trace'] == trace)]['time']
                            
                            # are we working on this electrode?

                            if electr == el:
                                
                                # plot the spike
                                if plot_it == 1:
                                    plt.plot(t[before], data_plot[before]+ add_it, 'r*')
                                
                                idx_nearest = np.where(other_spikes==sp)[0]
                                #repeated = repeated + len(idx_nearest) - 1
                                idx_nearest = idx_nearest[0]
    
                                nearest_el.append(electr)
                                nearest_pts.append(sp)
                                nearest_idxs.append(idx_nearest)
                            else:
                                # if there were any spikes detected in this electrode, in this trace
                                if len(other_spikes) > 0:
                                    # find the nearest spike to the spike detected in the given electrode
                                    idx_nearest, nearest = find_nearest(other_spikes, sp)
                                    #idx_nearest = np.where(other_spikes==sp)[0]

                                    if (nearest - sp) <= allow_shift and (nearest - sp) >= 0:
                                        
                                        if plot_it == 1:
                                            plt.plot(t[before+nearest-sp], data_plot[before+nearest-sp]+ add_it, 'r*')
                                        #repeated = repeated + len(np.where(other_spikes == nearest)) - 1
                                        nearest_idxs.append(idx_nearest)
                                        nearest_pts.append(nearest)
                                        nearest_el.append(electr)
                                         
                            add_it = add_it + 30
    
                        if len(nearest_pts) > 1:
                            # if more spikes are found close to the given one
                            max_ampls = np.zeros((len(nearest_pts)))
                            # check all the amplitudes for every spike found
                            for maxs in range(len(nearest_pts)):
                                max_ampls[maxs] = ampls[nearest_el[maxs]][trace][nearest_idxs[maxs]] 
                                win = np.argmax(max_ampls)
                            old = idx + len(nearest_pts) #+ repeated
                            
                            if max_ampls[win] > 0: 
                                winning_el.append(nearest_el[win])
                                winning_idx.append(nearest_idxs[win])
                                winning_trace.append(trace) 
    
                            if plot_it == 1:
                                data_plot = data_all[nearest_el[win]][tr][sp + before:sp+after]
                                plt.plot(t[before+nearest_pts[win]-sp], data_plot[before+nearest_pts[win]-sp]+ (30 * nearest_el[win]), 'bo')
                                plt.show()
                                plt.figure()                          
                           
     
                        elif ampls[nearest_el[0]][trace][nearest_idxs[0]] > 0.0:
                            #append this one spike with its data only
                            winning_el.append(nearest_el[0])
                            winning_idx.append(nearest_idxs[0])
                            winning_trace.append(trace)

                            if plot_it == 1:
                                data_plot = data_all[nearest_el[0]][tr][sp + before:sp+after]
                                plt.plot(t[before+nearest_pts[0]-sp], data_plot[before+nearest_pts[0]-sp]+ (30 * nearest_el[0]), 'bo')
                                #plt.show()
                                plt.figure()                          

                        
    # transform all the data so everything is sorted
    np.savez(save_folder + save_name, electrode = winning_el, trace = winning_trace, idx = winning_idx)
    idxs = transform(save_folder, save_name)
    np.savez(save_folder + save_name, idx = idxs)   

#def transform(data_all, electrode, trace, idx):
def transform(save_folder, save_name = 'chosen_spikes.npz'):
#transform all the data so everything is sorted
    all_spikes = []

    npzfile = np.load(save_folder + save_name)
    electrode = npzfile['electrode']
    trace = npzfile['trace']
    idx = npzfile['idx']
    data_all, fs = reader.read_datafile(save_folder)

    for electr in range(len(data_all)): #range(max(np.unique(electrode))+1):

        all_electr = np.where(electrode == electr)[0]
        all_tracespikes = []
        
        #print len(all_electr)
        #print len(data_all[electr])
        if len(all_electr) > 0:
            trace_el = trace[all_electr]
            idx_el = idx[all_electr]

        for tr in range(len(data_all[electr])): #range(max(np.unique(trace_el))):
            if len(all_electr) > 0:
                all_trace = np.where(trace_el == tr)[0]
                if len(all_trace) > 0:
                    idx_el_tr = idx_el[all_trace]
                    all_tracespikes.append(idx_el_tr.tolist())
                else:
                    all_tracespikes.append([])
            else:
                all_tracespikes.append([])


        all_spikes.append(all_tracespikes)
    return all_spikes


def plot_clust(save_folder, load_clust = 'exclust_sep.npz', load_chosen  = 'chosen_idx_val.npz', save_ext = '.png'):
    #np.savez(save_folder + save_clust, spik_clust)   
    no_groups = 2
    npzfile = np.load(save_folder + load_clust)
    clusters = npzfile['clusters']
    npzfile = np.load(save_folder + load_chosen)
    idxs = npzfile['sp_idx']
    a_s =npzfile['a_s']
    b_s =npzfile['b_s']
    ampls = npzfile['ampls']
    max_left = npzfile['max_left']
    max_right = npzfile['max_right']
    norm_fact = npzfile['norm_fact']
    
    groups = np.unique(clusters)
    data_all, fs = reader.read_datafile(save_folder)
    
    before = -2
    after = 2
    before = ispw.ms2pts(before, fs)
    after = ispw.ms2pts(after, fs)

    t = dat.get_timeline(data_all[0][0][0:-before+after], fs, 'ms')
    for electr in range(len(data_all)):
        plt.figure()
        
        plt.title(str(electr))
        for trace in range(len(data_all[electr])):
            
            for spike in range(len(idxs[electr][trace])):
                spike_ten = idxs[electr][trace][spike]
                max_left_ten = max_left[electr][trace][spike]
                baseline = data_all[electr][trace][spike_ten + max_left_ten]
                norm_ten = norm_fact[electr][trace][spike]
                #print spike_ten
                spike_shape = data_all[electr][trace][spike_ten+before:spike_ten+after] - baseline
                
                norm_ten = norm_fact[electr][trace][spike]
                #print 
                #print norm_ten
                #print spike_shape[-before]
                if spike_shape[-before] < 0:
                    
                    if len(spike_shape) >= -before+after:

                        plt.subplot(no_groups,1,clusters[electr][trace][spike] + 1)
                        spike_shape = (spike_shape) /np.abs(spike_shape[-before])
                        plt.plot(t, spike_shape)
                        plt.ylim([-2, 1])
                        #plt.plot(t, spike_shape, 'r')
                        #plt.show()

                        
                    #plt.plot(t[-before], data_all[electr][trace][spike_ten], 'ro')
            
                #plt.plot(t[-before + max_left_ten], data_all[electr][trace][spike_ten + max_left_ten], 'ro')
                #plt.show()

def save_chosen(save_folder, save_name = 'chosen_spikes.npz', save_chosen = 'chosen_idx_val.npz', save_ext = '.png'):

    # load corrected indexes
    plot_it = 0
    npzfile = np.load(save_folder + save_name)
    idxs = npzfile['idx']
    print 'clustering spikes'
    
    # load all indexes
    spike_idxs, a_s, b_s, ampls, fs_spikes = reader.read_ex_sparamas(save_folder)
    npzfile = np.load(save_folder + 'params_' + "ex_sparamas.npz")
    As = npzfile['As']
    Bs = npzfile['Bs']
    norm_fact = npzfile['norm_factor']

    #data_all, fs = reader.read_datafile(save_folder)
    idx = ba.idx2recarr(idxs)
    scale = 'ms'
    a_s = ba.idx2recarr(a_s, time_dtype=np.float32)
    b_s = ba.idx2recarr(b_s, time_dtype=np.float32)
    #amp_rel = ba.idx2recarr(amp_rel, time_dtype=np.float32)
    ampls = ba.idx2recarr(ampls, time_dtype=np.float32)
    As = ba.idx2recarr(As, time_dtype=np.float32)
    Bs  = ba.idx2recarr(Bs, time_dtype=np.float32)
    norm_fact = ba.idx2recarr(norm_fact, time_dtype=np.float32)
    idxs_choise = ba.idx2recarr(spike_idxs, time_dtype=np.float32)
    
    #take only the chosen spikes:
    as_electr, bs_electr, amp_electr, idx_electr =  [], [], [], []
    As_electr, Bs_electr, norm_fact_electr = [], [], []
    for electr in range(len(spike_idxs)):
        as_trace, bs_trace, amp_trace, idx_trace =  [], [], [], []
        As_trace, Bs_trace, norm_fact_trace = [], [], []
        for trace in range(len(spike_idxs[electr])):
            idx_used = idx[(idx['electrode'] == electr) & (idx['trace'] == trace)]['time'].tolist()
            a_s_chosen = a_s[(a_s['electrode'] == electr) & (a_s['trace'] == trace)][idx_used]['time'].tolist()
            b_s_chosen =  b_s[(b_s['electrode'] == electr) & (b_s['trace'] == trace)][idx_used]['time'].tolist()
            ampls_chosen =  ampls[(ampls['electrode'] == electr) & (ampls['trace'] == trace)][idx_used]['time'].tolist()
            idx_chosen = idxs_choise[(idxs_choise['electrode'] == electr) & (idxs_choise['trace'] == trace)][idx_used]['time'].tolist()
            As_chosen = As[(As['electrode'] == electr) & (As['trace'] == trace)][idx_used]['time'].tolist()
            Bs_chosen = Bs[(Bs['electrode'] == electr) & (Bs['trace'] == trace)][idx_used]['time'].tolist()
            norm_fact_chosen = norm_fact[(norm_fact['electrode'] == electr) & (norm_fact['trace'] == trace)][idx_used]['time'].tolist()
            
            as_trace.append(a_s_chosen)
            bs_trace.append(b_s_chosen)
            #amp_rel_trace.append(amp_rel_chosen)
            amp_trace.append(ampls_chosen)
            idx_trace.append(idx_chosen)
            As_trace.append(As_chosen)
            Bs_trace.append(Bs_chosen)
            norm_fact_trace.append(norm_fact_chosen)
            
        as_electr.append(as_trace)
        bs_electr.append(bs_trace)
        #amp_rel_electr.append(amp_rel_trace)
        amp_electr.append(amp_trace)
        idx_electr.append(idx_trace)
        As_electr.append(As_trace)
        Bs_electr.append(Bs_trace)
        norm_fact_electr.append(norm_fact_trace)        
    
    np.savez(save_folder + save_chosen, sp_idx = idx_electr, a_s = as_electr, b_s = bs_electr, ampls = amp_electr, max_left = As_electr, max_right = Bs_electr, norm_fact = norm_fact_electr)  
    
    
    
def exspikes_clust_separate(save_folder, n_clusters = 2, load_chosen  = 'chosen_idx_val.npz', save_clust = 'exclust_sep', save_ext = '.png'):

    npzfile = np.load(save_folder + load_chosen)
    idxs = npzfile['sp_idx']
    a_s =npzfile['a_s']
    b_s =npzfile['b_s']
    ampls = npzfile['ampls']
    max_left = npzfile['max_left']
    max_right = npzfile['max_right']
    norm_fact = npzfile['norm_fact']
    
    plot_it = 0
     
    a_s = ba.idx2recarr(np.array(a_s), time_dtype=np.float32)
    b_s = ba.idx2recarr(np.array(b_s), time_dtype=np.float32)
    #amp_rel = ba.idx2recarr(np.array(amp_rel_electr), time_dtype=np.float32)
    ampls = ba.idx2recarr(np.array(ampls), time_dtype=np.float32)
    
    spik_clust = []
    clusters_electr = []
    for electr in range(len(idxs)):
        a = a_s[a_s['electrode']==electr]['time']
        b = b_s[b_s['electrode']==electr]['time']
        #ampl_rel = amp_rel[amp_rel['electrode']==electr]['time']
        amp = ampls[ampls['electrode']==electr]['time']
        
        # prepare values used for clustering
        X = np.vstack((a, b)).T

        X = X/np.std(X,0);
        X = X-X.mean(0)
        X = X/np.std(X,0)

        cluster = clust.kmeans(n_clusters, X)
        spik_clust.append(cluster)
        cluster_trace = []
        length = 0
        for trace in range(len(idxs[electr])):
            length_new = length + len(idxs[electr][trace])
            cluster_trace.append(cluster[length:length_new].tolist())
            length = length_new
            
        #print len(spike_idxs[electr][trace])
        #print len(cluster_trace)
            
        clusters_electr.append(cluster_trace)
        
        if plot_it == 1:
            symbols = ['>', '<', 'o']
            cols = ['b', 'g', 'k']
            fig = plt.figure()
            plt.subplot(2,1,1)
            for i in range(n_clusters):
                plt.plot(X[cluster==i,0], X[cluster==i,1], symbols[i], color=cols[i], alpha = 0.03)

            tit = 'clustering extracellular spikes, electr: ' + str(electr)
    
            plt.title(tit)
            plt.xlabel('b') #'all valley to peak (a)')
            plt.ylabel('ampls') #'half peaks (c)')

            nom = save_folder + 'wyniki/' + save_clust + str(electr) + save_ext
            fig.savefig(nom,dpi=600)

            plt.subplot(2,1,2)

            x1 = a[np.where(cluster == 0)[0]]
            x2 = a[np.where(cluster == 1)[0]]

        if plot_it == 1:
            plt.hist([x1, x2],  normed=1)
            plot_hist(a, amp, 'a', 'ampls', xlab = 'Time [ms]', ylab = 'Voltage [mV]', fig = 0, bins = 100)
    
    #spik_clust = spik_clust.reshape(np.size(idxs,0), np.size(idxs,1), np.size(idxs,2))
    np.savez(save_folder + save_clust, clusters = clusters_electr)   
    #plt.show()

        
def detect_cell4spikes(save_folder, save_name = 'chosen_spikes.npz'):
    # load corrected indexes
    npzfile = np.load(save_folder + save_name)
    idx = npzfile['idx']
    
    # load all indexes
    spike_idxs, a_s, b_s, ampls, fs_spikes = reader.read_ex_sparamas(save_folder)
    
    data_all, fs = reader.read_datafile(save_folder)
    idx = ba.idx2recarr(idx)
    scale = 'ms'
    data_part_to_plot = [0, 10000]
    t = dat.get_timeline(range(data_part_to_plot[0],data_part_to_plot[1]), fs, scale)
    
    add_it = -30
    # plot all the data
    plt.figure()
    for electr in range(len(data_all)):
        for trace in [0]: #range(len(data_all[electr])):
            add_it = add_it + 300
            data_part = data_all[electr][trace][data_part_to_plot[0]:data_part_to_plot[1]]
            idx_used = idx[(idx['electrode'] == electr) & (idx['trace'] == trace)]['time']
            idx_used = idx_used.tolist()
            #idx_used = idx[(idx['trace'] == trace)]['time']
            #print idx_used
            spikes_used = []
            
            plt.plot(t,  data_part + add_it)
            for spidx in range(len(spike_idxs[electr][trace])):
                sp_temp = spike_idxs[electr][trace][spidx]
                amp_temp = b_s[electr][trace][spidx]
                if sp_temp < data_part_to_plot[1] and sp_temp > data_part_to_plot[0]:
                    plt.plot(t[sp_temp-data_part_to_plot[0]], data_part[sp_temp - data_part_to_plot[0]] + add_it, 'go')
     
                    tex = "%10.1f" % amp_temp #repr(amp_temp)
                    plt.text(t[sp_temp-data_part_to_plot[0]], data_part[sp_temp - data_part_to_plot[0]] + add_it, tex, fontsize=11, ha='center', va='top')
#            for spidx in range(len(spike_idxs[electr][trace])):
#                sp_temp = spike_idxs[electr][trace][spidx]
#                amp_temp = ampls[electr][trace][spidx]
#                if sp_temp < data_part_to_plot[1] and sp_temp > data_part_to_plot[0]:
#                    plt.plot(t[sp_temp-data_part_to_plot[0]], data_part[sp_temp - data_part_to_plot[0]] + add_it, 'go')
#                    t = str(amp_temp)
#                    plt.text(5, 10, t, fontsize=18, ha='center', va='top')
            #plt.text(5, 10, t, fontsize=18, ha='center', va='top')

            
            
            
            for sp in range(len(idx_used)):
                spike = spike_idxs[electr][trace][idx_used[sp]]
                spike_idxs[electr][trace][idx_used[sp]]

                if spike < data_part_to_plot[1] and spike > data_part_to_plot[0]:
                    spikes_used.append(spike)                  
                    plt.plot(t[spike-data_part_to_plot[0]], data_part[spike - data_part_to_plot[0]] + add_it, 'ro')
    plt.show()


def show_mean_root_differance(save_folder, save_file, all = 1, save_root_mean = 'root_mean', save_ext = '.png'):
    """ if all is 1 it means that indexes were calculated separately for each single electrode and then their equivalent in each single 
    electrode; if it's 0 there are all indexes put into one """
    
    add = 'wyniki/'
    npzfile = np.load(save_folder + 'init_' + save_file)
    init_root_mean = npzfile['root_mean']
    
    npzfile = np.load(save_folder + 'spont_' + save_file)
    spont_root_mean = npzfile['root_mean']
    
    fig = plt.figure(figsize=(20, 20))
    
    
    #import pdb; pdb.set_trace()
    if all == 1:
        for electr in range(len(spont_root_mean)):
            plt.subplot(np.ceil(np.sqrt(len(spont_root_mean))), np.ceil(np.sqrt(len(spont_root_mean))), electr+1)
            tit = 'detected in el: ' + str(electr)
            plt.title(tit)
            # all electrodes, all
            index = np.arange(len(init_root_mean))
    
            width = 0.40       # the width of the bars
            if len(init_root_mean[electr]) > 0:
                p1 = plt.bar(index - width, init_root_mean[electr], width, color='yellow', label ='initiated')[0]
            if len(spont_root_mean[electr]) > 0:
                p2 = plt.bar(index, spont_root_mean[electr], width, color='blue', label = 'spontaneous')[0]
    
            #plt.title('root mean square difference between SPWs')
            plt.legend(loc=4)
            plt.ylabel('root mean square')
            plt.xlabel('electrode')
            #save_nom = save_folder + add + save_root_mean + str(el) + save_ext
    else:
        
        index1 = np.arange(len(init_root_mean))
        index2 = np.arange(len(spont_root_mean))

        width = 0.40       # the width of the bars
        if len(init_root_mean) > 0:
            p1 = plt.bar(index1 - width, init_root_mean, width, color='yellow', label ='initiated')[0]
            
        if len(spont_root_mean) > 0:
            p2 = plt.bar(index2, spont_root_mean, width, color='blue', label = 'spontaneous')[0]

        plt.title('root mean square difference between SPWs')
        plt.legend(loc=4)
        plt.ylabel('root mean square')
        plt.xlabel('electrode')
        
    save_nom = save_folder + add + save_root_mean + save_ext
    fig.savefig(save_nom) #,dpi=600)
        
    plt.clf()



def spw_spike_plotter(save_folder,   save_name, save_ext, load_chosen  = 'chosen_idx_val.npz', load_maxs = 'spike_maxs.npz'):
    data_all, fs = reader.read_datafile(save_folder)
    spw_idxs, spw_maxs, starts_spw, ends_spw, lengths_spw, fs_spw = reader.read_SPWs(save_folder)

    distances, min_distances, fs, max_dist = reader.read_dist_SpwfromSpike(save_folder)

    max_spw2intraspike = 8 # ms
    max_spw2intraspike = ispw.ms2pts(max_spw2intraspike, fs)
    
    #fs_div = fs/fs_spw
    if fs != fs_spw:
        values, fs = updater.update_fs(fs_spw, fs, values = [starts_spw])
        starts_spw = values[0]
    
    # separate spikes to two different groups - those which are closer and those which are further than max_spw.. 
    distances = ba.idx2recarr(distances, time_dtype=np.float32)
    starts_spw = ba.idx2recarr(np.array(starts_spw))
    
    initiated = ((distances['time'] <= max_spw2intraspike) & (distances['time'] > 0))
    spontaneous = ((distances['time'] > max_spw2intraspike)) # & (distances['time'] !=float('inf'))]

    distances_in = distances[initiated]
    distances_sp = distances[spontaneous]
    group1 = starts_spw[initiated]
    group2 = starts_spw[spontaneous]
    
    win = (-3, 30)
    
    print 'ploting avg spw'
    if len(group1) > 0:
        plot_mean_wave_elecrwise(data_all, fs, group1, save_folder, save_name = 'avgSPWs_init', save_ext = save_ext, title = 'mean SPWs for each electrode initiated', win = win)
    if len(group2) > 0:
        plot_mean_wave_elecrwise(data_all, fs, group2, save_folder, save_name = 'avgSPWs_spont', save_ext = save_ext, title = 'mean SPWs for each electrode spontaneous', win = win)
    
    print 'ploting all SPWs and avg'
    if len(group1) > 0:
        plot_all_wave_elecrwise(data_all, fs, group1, save_folder, save_name = 'allSPWs_init', save_ext = save_ext, title = 'SPWs for each electrode initiated', win = win, save_root = 'init_mean_root.npz')
    if len(group2) > 0:
        plot_all_wave_elecrwise(data_all, fs, group2, save_folder, save_name = 'allSPWs_spont', save_ext = save_ext, title = 'SPWs for each electrode spontaneous', win =win, save_root = 'spont_mean_root.npz')
    
    print 'ploting all SPWs electrowise'
    if len(group1) > 0:
        plot_all_wave_elecrwise_all(data_all, fs, group1, save_folder, save_name = 'allSPWs_all_electrinit', save_ext = save_ext, title = 'SPWs for each electrode initiated', win = win,  save_root_all = 'init_mean_root_all.npz')
    if len(group2) > 0:
        plot_all_wave_elecrwise_all(data_all, fs, group2, save_folder, save_name = 'allSPWs_all_electrspont', save_ext = save_ext, title = 'SPWs for each electrode spontaneous', win = win,  save_root_all = 'spont_mean_root_all.npz')
    
    npzfile = np.load(save_folder + load_chosen)
    sp_idxs = npzfile['sp_idx']

    sp_idxs = ba.idx2recarr(np.array(sp_idxs))

    print 'ploting all SPWs and extra spikes'
    if len(group1) > 0:
        plot_each_wave_and_exspikes(data_all, fs, group1, sp_idxs,  distances_in, save_folder, save_name = 'each_SPWs_init', save_ext = save_ext, title = 'spikes SPWs for each electrode induced', win = win, save_file = 'dane_init.npz')
    if len(group2) > 0:
        plot_each_wave_and_exspikes(data_all, fs, group2, sp_idxs,  distances_sp, save_folder, save_name = 'each_SPWs_spont', save_ext = save_ext, title = 'spikes SPWs for each electrode spontaneous', win = win, save_file = 'dane_spont.npz')
#    
    if len(group1) > 0 and len(group2) > 0:
        show_mean_root_differance(save_folder, 'mean_root_all.npz', 1, 'root_mean_all', save_ext)
    if len(group1) > 0 and len(group2) > 0:  
        show_mean_root_differance(save_folder, 'mean_root.npz', 0, 'root_mean', save_ext)
    
    print 'loading the data'
    if len(group1) > 0:    
         plot_temp_hist(data_all, save_folder,  add_folder = 'all_spws/', save_name = 'init', save_ext = save_ext, title = 'mean SPWs for each electrode', win =win, save_file = 'dane_init.npz')
    if len(group2) > 0:    
         plot_temp_hist(data_all, save_folder,  add_folder = 'all_spws/', save_name = 'spont', save_ext = save_ext, title = 'mean SPWs for each electrode', win =win, save_file = 'dane_spont.npz')
    if len(group1) > 0:
        calc_waves_used(data_all, fs, group1, win, save_folder, add_name = 'induced')
        a = calc_cum_root_mean(data_all, fs, group1, win, save_folder, add_name  = 'induced.npz')
    if len(group2) > 0:
        calc_waves_used(data_all, fs, group2, win, save_folder, add_name = 'spont')
        b = calc_cum_root_mean(data_all, fs, group2, win, save_folder, add_name  = 'spont.npz')
   
    plt.plot(np.transpose(a), linestyle='--')
    plt.plot(np.transpose(b))
    plt.show()
    plt.legend()


def plot_all_wave_elecrwise_all(data, fs, starts_spw, save_folder, save_name = 'SPWs', save_ext = '.png', title = 'mean SPWs for each electrode', win = (-10, 30), save_root_all = 'mean_root_all.npz'):    
    add = 'wyniki/'
    scale = 'ms'
    add_it = np.linspace(0, 1000, len(data))
    
    #colors = ['b', 'g', 'r', 'c', 'm', 'y', ]
    colors = define_colors(no_colors = len(data))
    win_pts = (ispw.ms2pts(win[0], fs), ispw.ms2pts(win[1], fs))
     
    t = []
    root_mean_all = []
    for el in range(len(data)):
        st_spw = starts_spw[starts_spw['electrode']==el]
        l = len(st_spw)
        fig = plt.figure()
        tity = 'RMS: '
        root_mean = 0
        root_mean_electr = []
        
        for electr in range(len(data)):
            st_spw['electrode'] = electr
            spws, idx = ba.extract_spws(data, st_spw, fs, win_pts)   
            spw_mean = np.array([spws[:, idx['electrode']==i].mean(1) + add_it[i] 
                         for i in np.unique(idx['electrode'])]).T
            t = dat.get_timeline(spw_mean, fs, scale)
            #import pdb; pdb.set_trace()
            #if spws[:, idx['electrode']==electr].
            if electr in np.unique(idx['electrode']):
                s_mean = spws[:, idx['electrode']==electr].mean(1)
                substr_mean = spws - s_mean[:, np.newaxis]
                sqrt_substr = substr_mean**2
                sum_sqrt = np.sum(sqrt_substr, 0)
                mean_sum = sum_sqrt.mean(0)
                root_mean = np.sqrt(mean_sum)     
                root_mean_electr.append(root_mean)       
                tity = tity + ',el:' + str(electr) + ':' +str(str("{0:.0f}".format(root_mean)))

            if np.size(spws) > 0:
                plt.plot(t, spws[:, idx['electrode']==electr]+add_it[electr], color = colors[electr], alpha = 0.1)
                
         
            lines = plt.plot(t, spw_mean, linewidth=3.0)
            labels = []
        root_mean_all.append(root_mean_electr)
        #plt.xlabel(tit)
        #root_mean_electr.append(root_mean)
        #for i in np.unique(range(len(data))):
        
        if len(t) > 0:    
            lab = str(l) + ' SPWs' #+ ' in e.' + str(i)
            labels.append(lab)
            
            #labels = ['e. %d' % i for i in np.unique(idx['electrode'])]
            plt.legend(lines, labels)
            tit = title + ' found in el: ' + str(el)
            plt.title(tit)
            plt.xlim([0, t[-1]])
            xlab = 'time (' + scale + ')'
            ylab = 'voltage ('+ micro() + 'V)' 
            plt.ylabel(ylab)
            plt.xlabel(tity)
            plt.ylim([-100, 1200])
            save_nom = save_folder + add + save_name + str(el) + save_ext
        
            fig.savefig(save_nom,dpi=600)
            plt.clf()
    
    plt.close()     
    root_mean_all
    np.savez(save_folder + save_root_all, root_mean = root_mean_all)   


def find_ipsps(save_folder, load_file = 'IPSPs.npz', save_name = 'isps_clean', load_chosen  = 'chosen_idx_val.npz', title = 'mean IPSPs for each electrode', save_ext = '.png', win = [-1, 4]):
    add = 'wyniki/'
    scale = 'ms'
    start_range = 3 # ms

    print_all = 1
    
    data, fs = reader.read_datafile(save_folder)
    ipsp_idxs, ipsp_maxs, starts_ipsp, ends_ipsp, lengths_ipsp, fs_ipsp = reader.read_SPWs(save_folder, save_file = load_file)
    #import pdb; pdb.set_trace()
    npzfile = np.load(save_folder + load_chosen)
    sp_idxs = npzfile['sp_idx']
    sp_idxs = ba.idx2recarr(np.array(sp_idxs))
    
    #fs_div = fs/fs_spw
    if fs != fs_ipsp:
        values, fs = updater.update_fs(fs_ipsp, fs, values = [starts_ipsp])
        starts_ipsp = values[0]
        
    colors = define_colors(no_colors = len(data))
    start_range = ispw.ms2pts(start_range, fs)        
    after = ispw.ms2pts(win[1], fs)
    before = ispw.ms2pts(win[0], fs)
    

    add_it = np.linspace(0, 1000, len(data))
    
    starts_ipsp = ba.idx2recarr(np.array(starts_ipsp))
    win_pts = (ispw.ms2pts(win[0], fs), ispw.ms2pts(win[1], fs))
    ipsps, idx = ba.extract_spws(data, starts_ipsp, fs, win_pts) 
    
    g = 0
    plt.figure()
    for trace in range(len(data[0])):
        only_this_trace = starts_ipsp[(starts_ipsp['trace'] == trace)]
        sort_idx = np.argsort(only_this_trace['time'])
        posortowane = only_this_trace[sort_idx]
        
        for (i, spw) in enumerate(posortowane):
            last_end = 0
            sp_idxs_trace  = sp_idxs[(sp_idxs['trace'] == trace)]
            sp_idxs_spw    = sp_idxs_trace[(sp_idxs_trace['time'] < spw['time']+after) & 
                                           (sp_idxs_trace['time'] > spw['time']+before)]    
            
            if last_end < spw['time'] and len(sp_idxs_spw) < 3:

                old_g = 0
                g = g+1
                last_end = spw['time'] + start_range
                
                #plt.clf()
                
                for electr in range(len(data)):
                    #t = dat.get_timeline(spws[:, idx['electrode']==electr]+add_it[electr], fs, scale) +spw['time'] + before 
                    
                    t = dat.get_timeline(ipsps[:, idx['electrode']==electr], fs, scale) + win[0] 
                    t2 = dat.get_timeline(ipsps[:, idx['electrode']==electr]+add_it[electr], fs, scale) + win[0]
                    if len(t) == len(data[electr][trace][spw['time']+before:spw['time']+after]):
                        if print_all == 1:
                            fig = plt.figure(1)
                            plt.plot(t, data[electr][trace][spw['time']+before:spw['time']+after] + add_it[electr], color = colors[electr])
                        
                        
                        
                        sp_idxs_electr = sp_idxs_spw[(sp_idxs_spw['electrode'] == electr)]
                        
                        for sp_idx in sp_idxs_electr['time']:  #np.where(sp_idxs[electr][trace] < spw['time']+after and sp_idxs[electr][trace] > spw['time']+before):
                            if print_all == 1:
                                plt.figure(1)
                                plt.plot(t[sp_idx-spw['time']-before],data[electr][trace][sp_idx] +  add_it[electr], 'go' , color = colors[electr], alpha = 0.03)
                                #plt.plot(t[sp_idx-spw['time']-before],-100, 'g|' , color = colors[electr])
                            
                            plt.figure(2)
                            #import pdb; pdb.set_trace()
                            plt.plot(t2[sp_idx-spw['time']-before],g, 'g|' , color = colors[electr])
                            #spike_times.append(ispw.pts2ms(sp_idx-spw['time']-before, fs))
                            if print_all == 1:
                                plt.figure(1)
                                tit = title + 'at time: ' + str(ispw.pts2ms(spw['time'], fs)) + 'ms'
                                plt.title(tit)
                                plt.xlabel('time (ms)')  
        plt.show()              
        

def calc_waves_used(data, fs, starts_spw, win, save_folder, add_name = 'induced'):
    
    win_pts = (ispw.ms2pts(win[0], fs), ispw.ms2pts(win[1], fs))
    spws, idx = ba.extract_spws(data, starts_spw, fs, win_pts)  
    start_range = 5 # ms
    start_range = ispw.ms2pts(start_range, fs)
    starts_used = []
    electr_used = []
    trace_used = []
    
    # extract the data - later put into the separate function
    for trace in range(len(data[0])):
        last_end = 0
        only_this_trace = starts_spw[(starts_spw['trace'] == trace)]
        sort_idx = np.argsort(only_this_trace['time'])
        posortowane = only_this_trace[sort_idx]
        sp_times = []
        for (i, spw) in enumerate(posortowane):
            if last_end < spw['time']:  
                last_end = spw['time'] + start_range
                starts_used.append(spw['time'])
                electr_used.append(spw['electrode'])
                trace_used.append(spw['trace'])
                
    #ba.idx2recarr(idx, time_dtype)             
    wave_used = np.rec.fromarrays([electr_used, trace_used, starts_used], names='electrode, trace, time')
    #import pdb; pdb.set_trace()
    save_file = 'waves_used_' + add_name
    np.savez(save_folder + save_file, waves_used = wave_used) 
    
    
def calc_cum_root_mean(data, fs, starts_spw, win, save_folder, add_name = 'induced.npz'):
    # calculate culumative root mean 
    
    win_pts = (ispw.ms2pts(win[0], fs), ispw.ms2pts(win[1], fs))
    save_file = 'waves_used_' + add_name
    npzfile = np.load(save_folder + save_file)
    waves_used = npzfile['waves_used']
    
    no_of_pts = 50
    rms_all = []
    t = dat.get_timeline(np.linspace(win[0], win[1], no_of_pts), fs, 'ms')
    
    #for idx_spw in range(len(waves_used)):
    for electr in range(len(data)):
        root_mean_all = []
        for next_pt in np.linspace(win_pts[0] + 1, win_pts[1], no_of_pts):
            
            next_pt = int(next_pt)
            used_pts = [win_pts[0], next_pt]
            spws, idx = ba.extract_spws(data, waves_used, fs, used_pts) 
            
            s_mean = spws[:, idx['electrode']==electr].mean(1)
            substr_mean = spws[:, idx['electrode']==electr] - s_mean[:, np.newaxis]
            sqrt_substr = substr_mean**2
            sum_sqrt = np.sum(sqrt_substr, 0)
            mean_sum = sum_sqrt.mean(0)
            root_mean = np.sqrt(mean_sum)  
            root_mean_all.append(root_mean)
        rms_all.append(root_mean_all)
    return rms_all
    
        
            
            
            
            

        #mean_given = 
     
    
def plot_each_wave_and_exspikes(data, fs, starts_spw, sp_idxs,  distances, save_folder, save_name = 'SPWs', save_ext = '.png', title = 'mean SPWs for each electrode', win = (-10, 20), save_file = 'dane.npz'):
    add = 'wyniki/'
    scale = 'ms'
    print_all = 1
    
    #if len(data[0]) == 1:
    #only one trace
    after = ispw.ms2pts(win[1], fs)
    before = ispw.ms2pts(win[0], fs)
    start_range = 5 # ms
    start_range = ispw.ms2pts(start_range, fs)
    add_folder = 'all_spws/'
    
    add_it = np.linspace(0, 1000, len(data))
    colors = define_colors(no_colors = len(data))
    #colors = ['b', 'g', 'r', 'c', 'm', 'y', (0.5, 0.5, 0.5)]
    #colors = define_colors(no_colors = len(data))
    
    
    win_pts = (ispw.ms2pts(win[0], fs), ispw.ms2pts(win[1], fs))
    spws, idx = ba.extract_spws(data, starts_spw, fs, win_pts)    
    g = 0
    spike_times = []
    temp_sp_times = []
    plt.figure()
    spike_contour = []
    for trace in range(len(data[0])):
        last_end = 0
        only_this_trace = starts_spw[(starts_spw['trace'] == trace)]
        sort_idx = np.argsort(only_this_trace['time'])
        posortowane = only_this_trace[sort_idx]
        sp_times = []

        for (i, spw) in enumerate(posortowane):

            sp_idxs_trace  = sp_idxs[(sp_idxs['trace'] == trace)]
            sp_idxs_spw    = sp_idxs_trace[(sp_idxs_trace['time'] < spw['time']+after) & 
                                           (sp_idxs_trace['time'] > spw['time']+before)]
            if last_end < spw['time']:
                old_g = 0
                g = g+1
                last_end = spw['time'] + start_range
                
                #plt.clf()
                
                for electr in range(len(data)):
                    #t = dat.get_timeline(spws[:, idx['electrode']==electr]+add_it[electr], fs, scale) +spw['time'] + before 
                    
                    t = dat.get_timeline(spws[:, idx['electrode']==electr], fs, scale) + win[0] 
                    t2 = dat.get_timeline(spws[:, idx['electrode']==electr]+add_it[electr], fs, scale) + win[0]
                    if len(t) == len(data[electr][trace][spw['time']+before:spw['time']+after]):
                        if print_all == 1:
                            plt.figure(1)
                            plt.plot(t, data[electr][trace][spw['time']+before:spw['time']+after] + add_it[electr], color = colors[electr])
                        
                        sp_idxs_electr = sp_idxs_spw[(sp_idxs_spw['electrode'] == electr)]
                        
                        for sp_idx in sp_idxs_electr['time']:  #np.where(sp_idxs[electr][trace] < spw['time']+after and sp_idxs[electr][trace] > spw['time']+before):
                            if print_all == 1:
                                plt.figure(1)
                                plt.plot(t[sp_idx-spw['time']-before],data[electr][trace][sp_idx] +  add_it[electr], 'go' , color = colors[electr])
                                plt.plot(t[sp_idx-spw['time']-before],-100, 'g|' , color = colors[electr])
                            
                            plt.figure(2)
                            #import pdb; pdb.set_trace()
                            plt.plot(t2[sp_idx-spw['time']-before],g, 'g|' , color = colors[electr])
                            #spike_times.append(ispw.pts2ms(sp_idx-spw['time']-before, fs))
                            if print_all == 1:
                                plt.figure(1)
                                tit = title + 'at time: ' + str(ispw.pts2ms(spw['time'], fs)) + 'ms'
                                plt.title(tit)
                                plt.xlabel('time (ms)')
                            sp_contour = np.rec.fromarrays([ispw.pts2ms(sp_idx-spw['time']-before, fs), electr], names='time,electrode')
                            spike_contour.append(sp_contour)
                            
                plt.figure(1)
                save_nom = save_folder + add_folder + save_name + str(ispw.pts2ms(spw['time'], fs)) + save_ext
                plt.savefig(save_nom,dpi=600)
                plt.clf()
                
                sp_t = ispw.pts2ms(sp_idxs_spw['time']-spw['time']-before, fs)
                sp_data = np.rec.fromarrays([sp_t, sp_idxs_spw['electrode']], names='time,electrode')
                spike_times.append(sp_data)
                
                
#                            if old_g == g:
#                                temp_sp_times.append(ispw.pts2ms(sp_idx-spw['time']-before, fs))
#                            else:
#                                old_g = g
#                                if len(temp_sp_times) > 0:
#                                    sp_times.append(temp_sp_times)
#                                    temp_sp_times = []
#                                    print np.size(sp_times, 1)
                            
  
                    
    plt.figure(2)
    plt.ylim([0,g])
    plt.title(title)
    plt.xlim([win[0], win[1]])
    plt.xlabel('time (ms)')
    save_nom = save_folder + add_folder + 'raster_'+ save_name  + save_ext
    plt.savefig(save_nom,dpi=600)
    
    
    np.savez(save_folder + save_file, spike_contour = spike_contour, spike_times = spike_times)


def plot_distSPWfromSpike(save_folder, add_folder = 'all_spws/', save_ext = '.png'):
    bins = 5000
    #rang = np.linspace(0, 101, 500)
    data_all, fs = reader.read_datafile(save_folder)
    spw_idxs, spw_maxs, starts_spw, ends_spw, lengths_spw, fs_spw = reader.read_SPWs(save_folder)

    distances, min_distances, fs, max_dist = reader.read_dist_SpwfromSpike(save_folder)

    max_spw2intraspike = 5 # ms
    max_spw2intraspike = ispw.ms2pts(max_spw2intraspike, fs)
    
    #fs_div = fs/fs_spw
    if fs != fs_spw:
        values, fs = updater.update_fs(fs_spw, fs, values = [starts_spw])
        starts_spw = values[0]
    
    # separate spikes to two different groups - those which are closer and those which are further than max_spw.. 
    distances = ba.idx2recarr(distances, time_dtype=np.float32)
#    print distances
#    import pdb; pdb.set_trace()
    distances = distances[(distances['time'] >= 0) & (distances['time']!= float('Inf'))]['time']
    plt.hist(distances, 200, range = [0,101], normed = True)
    plt.xlim([0, 100])
    title = 'distanceSPWspike'
    plt.title(title)
    save_nom = save_folder + add_folder + title  + save_ext
    plt.savefig(save_nom,dpi=600)  
    plt.clf()       
    
    #starts_spw = ba.idx2recarr(np.array(starts_spw))
    
    #initiated = ((distances['time'] <= max_spw2intraspike) & (distances['time'] > 0))
    #spontaneous = ((distances['time'] > max_spw2intraspike)) # & (distances['time'] !=float('inf'))]

    #distances_in = distances[initiated]
    #distances_sp = distances[spontaneous]
    
    
    
    
def plot_temp_hist(data, save_folder,  add_folder, save_name = 'SPWs', save_ext = '.png', title = 'mean SPWs for each electrode', win = [-10, 30], save_file = 'dane_init.npz'):
    # must improve this function
    
    
    npzfile = np.load(save_folder + save_file)
    spike_contour = npzfile['spike_contour']
    spike_times = npzfile['spike_times']
    #times = npzfile['times']
    #import pdb; pdb.set_trace()
    #no_electr = np.unique(spike_contour['electrode'])
    no_electr = len(data)
    
    colors = define_colors(no_colors = no_electr)

    plt.figure()
    spt = np.concatenate([t['time'] for t in spike_times])
    plt.hist(spt, 100)
    #plt.ylim([0, 400])
    plt.title(title)
    save_nom = save_folder + add_folder + 'hist_' + save_name  + save_ext
    plt.savefig(save_nom,dpi=600)


    plt.figure()
    for electr in range(len(data)):
        plt.subplot(len(data), 1, electr)
        spike_electr = [int(spike_contour[i]['time']) for i in range(len(spike_contour)) if spike_contour[i]['electrode'] == electr]
        plt.hist(spike_electr, 100, win, color = colors[electr])
    plt.xlim([win[0], win[1]])
    plt.title(title)
    save_nom = save_folder + add_folder + 'hist_precise_'+ save_name  + save_ext
    plt.savefig(save_nom,dpi=600)  
    plt.clf()          
    
    #import pdb; pdb.set_trace()
    plt.figure()
    # figura bartkowa
    spike_count = [len(i) for i in spike_times]
    spw_order = np.argsort(spike_count)
    
    colors = np.array(colors)
    
    for idx, row in enumerate(spw_order):
        el = spike_times[row]['electrode']
        #import pdb; pdb.set_trace()
        times = spike_times[row]['time'] + win[0]
        #for spikes in range(len(times)):
        plt.scatter(times, np.ones(len(times))*idx, marker='|' , c = colors[el])
        
    plt.xlim([win[0], win[1]])
    plt.ylim([0, len(spw_order)])
    plt.title(title)
    save_nom = save_folder + add_folder + 'ordered_'+ save_name  + save_ext
    plt.savefig(save_nom,dpi=600)  
    plt.clf()  
    


    
    plt.figure()
    for electr in range(no_electr):
        plt.subplot(no_electr, 1, electr)
        spike_electr = [int(spike_contour[i]['time'] + win[0]) for i in range(len(spike_contour)) if spike_contour[i]['electrode'] == electr]
        plt.hist(spike_electr, 100, win, color = colors[electr])
    plt.xlim([win[0], win[1]])
    plt.title(title)
    save_nom = save_folder + add_folder + 'hist_precise_'+ save_name  + save_ext
    plt.savefig(save_nom,dpi=600)  
    plt.clf() 
    
        #import pdb; pdb.set_trace()
    plt.figure()
    # counterplot
    spike_electros = []
    for electr in range(no_electr):
        spike_electr = [float(spike_contour[i]['time'] + win[0]) for i in range(len(spike_contour)) if spike_contour[i]['electrode'] == electr]
        (hist), bin_edges = np.histogram(spike_electr, np.linspace(win[0], win[1], 100))
        spike_electros.append(hist)
    plt.imshow(spike_electros) 
    save_nom = save_folder + add_folder + 'colourbar_'+ save_name  + save_ext
    plt.title(title)
    plt.savefig(save_nom,dpi=600)  
    plt.clf()        
    
    

def plot_all_wave_elecrwise(data, fs, starts_spw, save_folder, save_name = 'SPWs', save_ext = '.png', title = 'mean SPWs for each electrode', win = (-10, 30), save_root = 'mean_root.npz'):    
    add = 'wyniki/'
    scale = 'ms'
    add_it = np.linspace(0, 1000, len(data))
    
    colors = define_colors(no_colors = len(data))
    
    win_pts = (ispw.ms2pts(win[0], fs), ispw.ms2pts(win[1], fs))
    spws, idx = ba.extract_spws(data, starts_spw, fs, win_pts)    
    spw_mean = np.array([spws[:, idx['electrode']==i].mean(1) + add_it[i] 
                         for i in np.unique(idx['electrode'])]).T
    t = dat.get_timeline(spw_mean, fs, scale)
    fig = plt.figure()
    tity = 'RMS: '
    root_mean_all = []
    
    for electr in range(len(data)):

        if np.size(spws[:, idx['electrode']==electr]+add_it[electr]) > 0:
            plt.plot(t, spws[:, idx['electrode']==electr]+add_it[electr], color = colors[electr], alpha = 0.1)
            #import pdb; pdb.set_trace()
            if electr in np.unique(idx['electrode']):
                s_mean = spws[:, idx['electrode']==electr].mean(1)
                substr_mean = spws[:, idx['electrode']==electr] - s_mean[:, np.newaxis]
                sqrt_substr = substr_mean**2
                sum_sqrt = np.sum(sqrt_substr, 0)
                mean_sum = sum_sqrt.mean(0)
                root_mean = np.sqrt(mean_sum)  
                root_mean_all.append(root_mean)
                tity = tity + ',el:' + str(electr) + ':' + str("{0:.0f}".format(root_mean))
    
    lines = plt.plot(t, spw_mean, linewidth=3.0)
    labels = []
    for i in np.unique(idx['electrode']):
        l = np.size(spws[:, idx['electrode']==i], 1)
        lab = str(l) + ' in e.' + str(i)
        labels.append(lab)
    
    #labels = ['e. %d' % i for i in np.unique(idx['electrode'])]
    plt.legend(lines, labels)
    plt.title(title)
    plt.xlim([0, t[-1]])
    xlab = 'time (' + scale + ')'
    ylab = 'voltage ('+ micro() + 'V)' 
    plt.xlabel(tity)
    plt.ylabel(ylab)
    #plt.ylabel()
    plt.ylim([-100, 1200])
    save_nom = save_folder + add + save_name + save_ext
    fig.savefig(save_nom,dpi=600)
    
    
    plt.close()     
    np.savez(save_folder + save_root, root_mean = root_mean_all)   


    
def plot_mean_wave_elecrwise(data, fs, starts_spw, save_folder, save_name = 'avgSPWs', save_ext = '.png', title = 'mean SPWs for each electrode', win = (-10, 30)):    
    add = 'wyniki/'
    scale = 'ms'
    
    win_pts = (ispw.ms2pts(win[0], fs), ispw.ms2pts(win[1], fs))
    if len(starts_spw) > 0:
        spws, idx = ba.extract_spws(data, starts_spw, fs, win_pts)    
        spw_mean = np.array([spws[:, idx['electrode']==i].mean(1) 
                             for i in np.unique(idx['electrode'])]).T
        
        #import pdb; pdb.set_trace()
        t = dat.get_timeline(spw_mean, fs, scale)
        fig = plt.figure()
        lines = plt.plot(t, spw_mean)
        labels = []
        for i in np.unique(idx['electrode']):
            l = np.size(spws[:, idx['electrode']==i], 1)
            lab = str(l) + ' in e.' + str(i)
            labels.append(lab)
        
        #labels = ['e. %d' % i for i in np.unique(idx['electrode'])]
        plt.legend(lines, labels)
        plt.title(title)
        plt.xlim([0, t[-1]])
        xlab = 'time (' + scale + ')'
        ylab = 'voltage ('+ micro() + 'V)' 
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        save_nom = save_folder + add + save_name + save_ext
        fig.savefig(save_nom,dpi=600)
        
        plt.close()     


def plot_extra_params(save_folder, ext = '.png', intra = 0):
    if intra == 1:
        distances, min_distances, fs, max_dist = reader.read_dist_fromSpike(save_folder)
    spike_idxs, all_valley_to_peak, all_half_valley_width, all_half_peak_width, ampls, fs_new = reader.read_ex_sparamas(save_folder)

    data_all, fs = reader.read_datafile(save_folder)
    spw_idxs, spw_maxs, starts_spw, ends_spw, lengths_spw, fs_spw = reader.read_SPWs(save_folder)
    np.shape(spike_idxs)
    #if fs != fs_SPWs
    spike_idxs = ba.idx2recarr(spike_idxs) 
    
    starts_spw = ba.idx2recarr(starts_spw)
    if intra == 1:
        distances = ba.idx2recarr(distances, time_dtype=np.float32)
    all_valley_to_peak = ba.idx2recarr(all_valley_to_peak, time_dtype=np.float32)
    all_half_peak_width = ba.idx2recarr(all_half_peak_width, time_dtype=np.float32)

    scale = 'ms'     

    t = dat.get_timeline(data_all[0], fs, scale)
    for electr in range(len(spike_idxs)): 

        if intra == 1:
            dist_ok = []
            all_val_to_peak_ok_dist = []
            all_half_peak_ok_dist = []
            for eliminate in range(len(distances[distances['electrode']==electr]['time'])):
                if ~np.isinf(distances[distances['electrode']==electr]['time'][eliminate]):
                    #print distances[distances['electrode']==electr]['time'][eliminate]
                    dist_ok.append(distances[distances['electrode']==electr]['time'][eliminate])
                    all_val_to_peak_ok_dist.append(all_valley_to_peak[all_valley_to_peak['electrode']==electr]['time'][eliminate])
                    all_half_peak_ok_dist.append(all_half_peak_width[all_half_peak_width['electrode']==electr]['time'][eliminate])
        all_val_to_peak_ok = all_valley_to_peak[all_valley_to_peak['electrode']==electr]['time']
        all_half_peak_ok = all_half_peak_width[all_half_peak_width['electrode']==electr]['time']

        if  len(all_val_to_peak_ok) >0 and len(all_half_peak_ok) > 0:
            fig = plot_hist(all_val_to_peak_ok, all_half_peak_ok, 'all val to peak (a)', 'all half peak (c) (ms)', xlab = 'Time [ms]', ylab = 'Voltage [mV]', fig = 0, bins = 20)
            nom = save_folder + 'wyniki/ext_spikes_a_vs_c' + str(electr) + ext
            fig.savefig(nom,dpi=600)
            plt.close()
        
            fig = plt.figure()
            hist,xedges,yedges = np.histogram2d(all_val_to_peak_ok, all_half_peak_ok,bins=40,range=[[0.01,max(all_val_to_peak_ok)],[0.01,max(all_half_peak_ok)]])
            extent = [xedges[0], xedges[-1], yedges[0], yedges[-1] ]
            plt.imshow(hist.T,extent=extent,interpolation='nearest',origin='lower')
            plt.colorbar()
            plt.ylabel('half peaks [ms]')
            plt.xlabel('all valley to peak [ms]')
            tit = 'ext spikes, electr: ' + str(electr)
            plt.title(tit)
            nom = save_folder + 'wyniki/ext_spikes_a_vs_cv2_' + str(electr) + ext
            fig.savefig(nom,dpi=600)
            plt.close()        
        
        
        if intra == 1 and len(dist_ok) > 0:
            
            plt.figure()
            hist,xedges,yedges = np.histogram2d(dist_ok, all_val_to_peak_ok_dist,bins=40,range=[[0.01,max(dist_ok)],[0.01,max(all_val_to_peak_ok_dist)]])
            extent = [xedges[0], xedges[-1], yedges[0], yedges[-1] ]
            plt.imshow(hist.T,extent=extent,interpolation='nearest',origin='lower')
            tit = 'ext spikes, electr: ' + str(electr)
            plt.xlabel('distances from any intracellular spike [ms]')
            plt.ylabel('all valley to peak [ms]')
            plt.title(tit)
            plt.colorbar()   
            plot_hist(dist_ok, all_val_to_peak_ok_dist, 'distances from any intracellular spike [ms]', 'all half peak (c) (ms)', xlab = 'Time [ms]', ylab = 'Voltage [mV]', fig = 0, bins = 50)
            nom = save_folder + 'wyniki/ext_spikes_a_vs_distToSpike_' + str(electr) + ext
            fig.savefig(nom,dpi=600)
            plt.close()            
    