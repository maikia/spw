                                                                     
                                                                     
                                                                     
                                             
import numpy as np
import data_mang as dat
import pylab as plt
import update_data as updater
import induc_SPW as ispw
import scipy.signal as signal
import b_analyser as ba
from matplotlib.ticker import NullFormatter
from mpl_toolkits.mplot3d import Axes3D
import folder_manager as fold_mng
from matplotlib.widgets import Slider, Button
from matplotlib.transforms import Bbox
from matplotlib.path import Path
import logging 
import scipy.cluster as cluster
import gc #garbage collector

def plot_dist_spw2spike(save_folder, plot_folder, save_plots, dist_file, ext):
    """ plots the histogram of the distribution of spws from the spike"""
    npzfile         = np.load(save_folder + dist_file)
    dist = npzfile['dist_spwspike']
    
    distance = dist['distance']
    npzfile.close()
    
    xlim = [-3, 70]
    no_bins = 200
    distance = distance[(distance > xlim[0]) & (distance < xlim[1])]
    #import pdb; pdb.set_trace() 
    fig = plt.figure()                
    plt.hist(distance, no_bins, normed=1, color = 'k')
    plt.title('Distribution of SPWs from the spike')
    plt.xlim(xlim)
    plt.ylabel('Fraction of SPWs')
    plt.xlabel('Distance from spike (ms)')
    fig.savefig(save_folder + plot_folder + save_plots + ext,dpi=600)         
    #plt.show() 
    #plt.close()  
    
    
def plot_noIpsps2distance(save_folder, plot_folder, save_plots, spw_file, dist_file, ext):
    """ it takes previously calculated distances spw to spike and calculates number of Ipsp for each spike
    and plots them against each other"""
    
    npzfile         = np.load(save_folder + spw_file)
    
    ipsps      = npzfile['spw_ipsps'] # spikes_all
    npzfile.close()    
    
    npzfile         = np.load(save_folder + dist_file)
    distances      = npzfile['dist_spwspike'] # spikes_all
    npzfile.close()    
    
    all_spws = np.unique(distances['spw_no'])
    numb_ipsps = np.zeros([len(all_spws)])
    dists = np.zeros([len(all_spws)])
    
    fig = plt.figure()
    for idx, spw_no in enumerate(all_spws):
        # count how many IPSPs are detected in each SPW
        ipsps_used = ipsps[ipsps['spw_no'] == spw_no]['group']
        uniq = np.unique(ipsps_used)
        no_ipsps = len(uniq)
        
        # check how far this spw is from the spike
        dis = np.unique(distances[distances['spw_no'] == spw_no]['distance'])
        
        numb_ipsps[idx] = no_ipsps
        dists[idx] = dis
        
    plt.plot(dists, numb_ipsps, 'ro', alpha = 0.2)
    plt.title('Relation of distance and number of IPSPs in SPW')
    plt.xlabel('Distance from spike [ms]')
    plt.ylabel('Number of IPSPs')
    
    save_fold = save_folder + plot_folder
    fold_mng.create_folder(save_fold)
    #plt.show()
    fig.savefig(save_fold + save_plots + ext,dpi=600)      
    plt.close()
        

def display_data(save_folder, plot_folder, save_plots, data_file, trace = 0, part = [0, 10000], ext = '.pdf'):
    """ it plots given data trace"""
    npzfile        = np.load(save_folder + data_file)
    data = npzfile['data']
    fs = npzfile['fs']
    npzfile.close() 
    
    fig = plt.figure()
    data_temp = data[:,trace, part[0]: part[1]]
    t = dat.get_timeline(data_temp[0,:], fs, 'ms')
    #import pdb; pdb.set_trace() 
    add_it = 100
    for electr in range(np.size(data_temp, 0)):
        plt.plot(t, data_temp[electr, :] + electr * add_it)
    plt.title('data')
    plt.ylabel('Voltage (' + micro() + 'V)')
    plt.xlabel('Time (ms)')
        

    save_fold = save_folder + plot_folder
    fold_mng.create_folder(save_fold)
    fig.savefig(save_fold + save_plots + ext,dpi=600)     
    fig.savefig(save_fold + save_plots + '.eps',dpi=600)    
    #plt.show() 
    plt.close()   

def display_data_no_electrode(save_folder, plot_folder, save_plots, data_file, part = [0, 10000], ext = '.pdf'):
    """ it plots given data trace"""
    npzfile        = np.load(save_folder + data_file)
    data = npzfile['data']
    fs = npzfile['fs']
    npzfile.close() 
    
    if part == 'all':
        part = [0, np.size(data,1)]
    fig = plt.figure()
    import pdb; pdb.set_trace()
    for trace in range(len(data)):
        plt.subplot(2,2,trace)
        data_temp = data[trace, part[0]: part[1]]
        t = dat.get_timeline(data_temp[:], fs, 'ms')
    #import pdb; pdb.set_trace() 
    
        plt.plot(t, data_temp[:])
        plt.ylabel('Voltage (' + micro() + 'V)')
        plt.xlabel('Time (ms)')
        plt.xlim([t[0], t[-1]])
        
    save_fold = save_folder + plot_folder
    fold_mng.create_folder(save_fold)
    fig.savefig(save_fold + save_plots + ext,dpi=600)       
    plt.show() 
    plt.close()   

def plot_data_interactive(save_folder, load_datafile, load_spw_ipsps, load_spikefile, load_spikesall, load_ipspsOld, spw_base, load_dataintrafile, load_intraSpikes):
    """ Plots interactively all the data and detected SPWs, spikes and IPSPs. does not save anything"""
    

    npzfile        = np.load(save_folder + load_datafile)
    data = npzfile['data']
    fs = npzfile['fs']
    npzfile.close() 
    
    npzfile        = np.load(save_folder + load_spw_ipsps)
    ipsps = npzfile['spw_ipsps']
    npzfile.close()
    #print np.unique(ipsps['electrode'])
    #import pdb; pdb.set_trace()
    npzfile        = np.load(save_folder + load_spikefile)
    #import pdb; pdb.set_trace()
    spikes = npzfile['spike_idx'] 
    #spikes = npzfile['chosen_spikes']
    npzfile.close()  
    
    npzfile        = np.load(save_folder + load_spikesall)
    spikes_alll = npzfile['spike_idx'] 
    npzfile.close()  
    
    npzfile = np.load(save_folder + load_ipspsOld)
    ipsps_old = npzfile['spw_ipsps'] 
    npzfile.close()
    
    #import pdb; pdb.set_trace()

    npzfile = np.load(save_folder + spw_base)
    primitive_starts = npzfile['spw_details']
    #primitive_starts = npzfile['starts']
    npzfile.close()
    
    #import pdb; pdb.set_trace()
    npzfile = np.load(save_folder + load_dataintrafile)
    intra_data = npzfile['data']
    npzfile.close()
    
    
    npzfile = np.load(save_folder + load_intraSpikes)
    intra_spikes = npzfile['spikes_first']
    npzfile.close()    
    

    ax = plt.subplot(111)
    plt.subplots_adjust(bottom=0.2)
     
    
    class Update_Plot():
        # initiate all the variable used throughout the class
        win = [-20, 80]
        ipsps_used = ipsps
        used_spikes = spikes
        used_primStarts = primitive_starts
        used_ipsps_old = ipsps_old
        spikes_all = spikes_alll
        all_data = data
        all_intradata = intra_data
        intraspikes_used =intra_spikes

        used_trace =  min(ipsps_used['trace'])
        used_spw_no = min(ipsps_used[ipsps_used['trace'] == used_trace]['spw_no'])
        full = False
        used_spw_idx = 0
        add_it = 150
        
        data_used = all_data[:, used_trace, :]
        t = dat.get_timeline(data_used[0,:], fs, 'ms')
        ylims = [-add_it, add_it * (np.size(data_used,0) + 3)]
        xs1 = [0, win[1] - win[0]]
        xs2 = [0, max(t)]
        xlims = [0, win[1] - win[0]]  
        
        actual_ipsps = ipsps_used[ipsps_used['trace'] == used_trace]
        actual_spikes = used_spikes[used_spikes['trace'] == used_trace]
        actual_intraspikes = intraspikes_used[intraspikes_used['trace'] == used_trace]
        actual_spikes_all = spikes_all[spikes_all['trace'] == used_trace]
        actual_ipsps_old = used_ipsps_old[used_ipsps_old['trace'] == used_trace]
        actual_prim_starts = used_primStarts[used_primStarts['trace'] == used_trace]
        
        group_colors = np.array([(1.0,0., 0.), 
                                 (0., 1., 0.), 
                                 (0., 0., 1.), 
                                 (0.6,0., 0.6),
                                 (0., 0.6, 0.6),
                                 (0.6, 0.6, 0)]) 
        ipsp_group_all = actual_ipsps['group']
        ipsps_group_colors = group_colors[ipsp_group_all % len(group_colors)]
        
        #import pdb; pdb.set_trace()
        used_spws = ipsps_used[ipsps_used['trace'] == used_trace]
        used_spw_start = used_spws['spw_no' == used_spw_no]['spw_start'] 
        no_spw_in_trace = np.unique(used_spws['spw_no'])
        no_spw_all = np.unique(ipsps_used['spw_no'])

        # purly plot variables
        spike_sign = '^'
        ipsps_sign = 'o'
        spike_color = 'g'
        ipsps_color = 'r'
      
        plt.xlim([max(0,used_spw_start + win[0]), min(used_spw_start + win[1], np.size(data_used,1))])   
        plt.title('trace' + str(used_trace))        
                   
 
 

        
        def full_part(self, event):
            if self.full:
                self.full = False
                self.xlims = self.xs1
            else:
                self.full = True
                self.xlims = self.xs2
            self.redraw_data()                
        
        def next_spw(self):
            #import pdb; pdb.set_trace()
            self.used_spw_no = self.no_spw_all[self.used_spw_idx]
            if len(self.no_spw_all) != 0:
                self.used_spw_idx = min(self.used_spw_idx + 1, len(self.no_spw_all)-1)
                spw_no = self.no_spw_all[self.used_spw_idx]
                self.used_spw_start = self.ipsps_used[self.ipsps_used['spw_no'] == spw_no]['spw_start'][0]
                temp_trace = self.ipsps_used[self.ipsps_used['spw_no'] == spw_no]['trace'][0]
            return temp_trace
        
        def previous_spw(self):
            self.used_spw_no = self.no_spw_all[self.used_spw_idx]
            if len(self.no_spw_all) != 0:
                self.used_spw_idx = max(self.used_spw_idx -1, 0)
                spw_no = self.no_spw_all[self.used_spw_idx]
                self.used_spw_start = self.ipsps_used[self.ipsps_used['spw_no'] == spw_no]['spw_start'][0]
                temp_trace = self.ipsps_used[self.ipsps_used['spw_no'] == spw_no]['trace'][0]
            return temp_trace
        
        def next(self, event): #self, event):         
            temp_trace = self.next_spw()
            if temp_trace == self.used_trace:
                # trace is unchanged
                central = self.used_spw_start
                self.set_lim(central)
            else:
                self.used_trace = temp_trace
                self.redraw_data()
            plt.draw()
        
        def set_lim(self, center):
            if self.full:
                ax.set_xlim([self.xlims[0], self.xlims[1]]) 
            else:
                ax.set_xlim([max(0,center + self.win[0]), min(center + self.win[1], np.size(self.data_used,1))]) 
        
        def prev(self, event): #self, event):
            temp_trace = self.previous_spw()
            if temp_trace == self.used_trace:
                # trace is unchanged
                central = self.used_spw_start
                self.set_lim(central)
            else:
                self.used_trace = temp_trace
                self.redraw_data()
            plt.draw()
            
        def next_trace(self, event):
            temp_trace = min(self.used_trace +1, np.size(self.all_data,1)-1)
            if temp_trace != self.used_trace:
                self.used_trace = temp_trace
                self.redraw_data()
            # update data used, 
            
        def prev_trace(self, event):
            #import pdb; pdb.set_trace()
            temp_trace = max(self.used_trace -1, 0)
            if temp_trace != self.used_trace:
                self.used_trace = temp_trace
                self.redraw_data()
            #used_trace = max(used_trace - 1, 0)
            #sed_spw_no = min(ipsps['trace' == used_trace]['spw_no'])
            
       
        def redraw_data(self):  
            #import pdb; pdb.set_trace()
            all_spws = ipsps[ipsps['trace'] == self.used_trace]['spw_no']
            ax.cla()
            ax.set_ylim(self.ylims[0], self.ylims[1])     
            
            ax.set_title('trace: ' + str(self.used_trace))
            self.data_used = self.all_data[:, self.used_trace, :]
            self.actual_ipsps = self.ipsps_used[self.ipsps_used['trace'] == self.used_trace]
            self.actual_spikes = self.used_spikes[self.used_spikes['trace'] == self.used_trace]
            self.actual_intraspikes = self.intraspikes_used[self.intraspikes_used['trace'] == self.used_trace]
            self.actual_spikes_all = self.spikes_all[self.spikes_all['trace'] == self.used_trace]
            self.actual_ipsps_old = self.used_ipsps_old[self.used_ipsps_old['trace'] == self.used_trace]
            self.actual_prim_starts = self.used_primStarts[self.used_primStarts['trace'] == self.used_trace]
            self.ipsp_group_all = self.actual_ipsps['group']
            self.ipsps_group_colors = self.group_colors[self.ipsp_group_all % len(self.group_colors)]
            self.intradata_used = self.all_intradata[:, self.used_trace, :]

            if len(all_spws) == 0:
                self.no_spw_in_trace = 0
                ax.set_xlim([self.xlims[0], self.xlims[1]]) 
            else:
                
                self.used_spw_no = min(ipsps[ipsps['trace'] == self.used_trace]['spw_no'])  
                
                self.used_spws = ipsps[ipsps['trace'] == self.used_trace]
                self.no_spw_in_trace = np.unique(self.used_spws['spw_no'])
                
                #self.used_spw_idx = 0
                self.used_spw_start = self.used_spws[self.used_spws['spw_no']==self.used_spw_no]['spw_start'][0]

                self.set_lim(self.used_spw_start)
                
                for sp in range(len(self.used_spws['spw_start'])):
                    rect = plt.Rectangle((self.used_spws['spw_start'][sp], self.ylims[0]), 1, self.ylims[1] + self.add_it, facecolor='y', lw=0)
                    #plt.gca().add_patch(rect)
                    ax.add_patch(rect)
            
            for electr in range(np.size(self.data_used,0)):
                actual_ipsps_electr = self.actual_ipsps[self.actual_ipsps['electrode'] == electr]['ipsp_start']
                actual_ipsps_pts = ispw.ms2pts(actual_ipsps_electr, fs).astype('i4')
                
                ax.plot(self.t, self.data_used[electr,:] + electr * self.add_it, 'k')
                #ax.plot(self.t[actual_ipsps_pts], self.data_used[electr,actual_ipsps_pts] + electr * self.add_it, self.ipsps_color + self.ipsps_sign)       
                if len(self.actual_ipsps) > 0:
                    actual_ipsps_selected = self.actual_ipsps[self.actual_ipsps['electrode'] == electr]

                    actual_ipsps_electr = actual_ipsps_selected['ipsp_start']
                    ipsp_color = self.ipsps_group_colors[self.actual_ipsps['electrode'] == electr]

                    actual_ipsps_pts = ispw.ms2pts(actual_ipsps_electr, fs).astype('i4')
                    ax.scatter(self.t[actual_ipsps_pts], self.data_used[electr,actual_ipsps_pts] + electr * self.add_it, 
                           c=ipsp_color, marker=self.ipsps_sign, s=70,zorder=10)
                      
                    
                actual_spikes_allelectr = self.actual_spikes_all[self.actual_spikes_all['electrode'] == electr]['time']
                actual_spikes_allpts = ispw.ms2pts(actual_spikes_allelectr, fs).astype('i4')
                ax.plot(self.t[actual_spikes_allpts], self.data_used[electr,actual_spikes_allpts] + electr * self.add_it, 
                        self.spike_color + self.spike_sign, mfc='none', ms=7, linewidth = 4)
                
                actual_spikes_electr = self.actual_spikes[self.actual_spikes['electrode'] == electr]['time']
                actual_spikes_pts = ispw.ms2pts(actual_spikes_electr, fs).astype('i4')
                ax.plot(self.t[actual_spikes_pts], self.data_used[electr,actual_spikes_pts] + electr * self.add_it, self.spike_color + self.spike_sign)
                
                actual_ipsps_old_electr = self.actual_ipsps_old[self.actual_ipsps_old['electrode'] == electr]['ipsp_start']
                actual_ipsps_old_pts = ispw.ms2pts(actual_ipsps_old_electr, fs).astype('i4')
                ax.plot(self.t[actual_ipsps_old_pts], self.data_used[electr,actual_ipsps_old_pts] + electr * self.add_it, 
                        self.ipsps_color + self.ipsps_sign,  mfc='none', ms=7, linewidth = 4)
                
                actual_prim_starts_electr = self.actual_prim_starts[self.actual_prim_starts['electrode'] == electr]['time']
                actual_prim_starts_pts = ispw.ms2pts(actual_prim_starts_electr, fs).astype('i4')
                ax.plot(self.t[actual_prim_starts_pts], self.data_used[electr,actual_prim_starts_pts] + electr * self.add_it, 'r*', 
                        ms=7, linewidth = 4)  
               
             
            ax.plot(self.t, self.intradata_used[0,:] + (electr+3) * self.add_it, 'k')     
            
            actual_intraspikes_electr = self.actual_intraspikes['time']
            actual_intraspikes_pts = ispw.ms2pts(actual_intraspikes_electr, fs).astype('i4')
            
            ax.plot(self.t[actual_intraspikes_pts], self.intradata_used[0,actual_intraspikes_pts] + (electr+3) * self.add_it, 'r*') #self.spike_color + self.spike_sign)
            plt.draw()
            
            #import pdb; pdb.set_trace()
                
        
    
    callback = Update_Plot()
    callback.redraw_data()
    
    axprev = plt.axes([0.6, 0.05, 0.15, 0.075])
    axnext = plt.axes([0.75, 0.05, 0.15, 0.075])
    axprev_trace = plt.axes([0.05, 0.05, 0.15, 0.075])
    axnext_trace = plt.axes([0.20, 0.05, 0.15, 0.075])
    axfull_trace = plt.axes([0.40, 0.05, 0.15, 0.075])
    
    
    bnext = Button(axnext, 'Next')
    bnext.on_clicked(callback.next)
    bprev = Button(axprev, 'Previous')
    bprev.on_clicked(callback.prev)
    
    bfull = Button(axfull_trace, 'All/part')
    bfull.on_clicked(callback.full_part)
    
    bnext_trace = Button(axnext_trace, 'Next trace')
    bnext_trace.on_clicked(callback.next_trace)
    bprev_trace = Button(axprev_trace, 'Previous trace')
    bprev_trace.on_clicked(callback.prev_trace)
    

    plt.show()

#def plot_by_ipsp_amplitude(save_folder, plot_folder, spw_data, ext):
#
#    npzfile        = np.load(save_folder + spw_data)
#    spontaneous = npzfile['spontaneous']
#    initiated = npzfile['initiated']
#    npzfile.close()  
#    #import pdb; pdb.set_trace()      
#    
#    for spw_no in initiated['spw_no']:
#        import pdb; pdb.set_trace()     
#        for ipsp_no in initiated[initiated['spw_no'] == spw_no]['group']


def find_max_corr(x, y):
    
    nomean_x = x - np.mean(x)
    nomean_y = y - np.mean(y)

    cxy = np.correlate(nomean_x,nomean_y,mode='full')
    i = cxy.argmax()
    lag = i - len(cxy)/2
    return lag
 
def PCA(data,ncomps=2, start_comp=0):
    """Perfrom a principle component analysis.
    
    Parameters
    ----------
    data : array
    (n_vars, n_obs) array where `n_vars` is the number of
    variables (vector dimensions) and `n_obs` the number of
    observations
    
    Returns
    -------
    evals : array
    sorted eigenvalues
    evecs : array
    sorted eigenvectors
    score : array
    projection of the data on `ncomps` components
    """

    #norm=data/np.std(data,1)[:,np.newaxis]
    #norm[np.isnan(norm)]=0
    #norm = data
    data = data.astype(np.float64)
    K=np.cov(data)
    evals,evecs=np.linalg.eig(K)
    order=np.argsort(evals)[::-1]
    evecs=np.real(evecs[:,order])
    evals=np.abs(evals[order])
    score= np.dot(evecs[:,start_comp:start_comp+ncomps].T,data)
    score = score/np.sqrt(evals[:ncomps, np.newaxis])
    return evals,evecs,score

def calculate_PCA(new_starts_pts, traces, data, fs, electr_to_use, ncomps, window, window_base):
    # give new beginning of SPWs and data
 
    remove_baseline = True
    pc = []
    s_p_temp1 = ispw.ms2pts(window_base[0], fs).astype('i4') + new_starts_pts 
    s_p_temp2 = ispw.ms2pts(window_base[1], fs).astype('i4') + new_starts_pts
    window_pts = [s_p_temp1, s_p_temp2]
    
    for electr in electr_to_use:
        #import pdb; pdb.set_trace() 
        data_used = get_spws(data, traces, electr, window, fs, new_starts_pts, remove_baseline, window_pts)
        #import pdb; pdb.set_trace() 

        a1, b1, principal_components = PCA(data_used,ncomps=ncomps)
        
        #plt.figure()
        #plt.plot(b1[:,:4]) # plot first 4 comps
        #plt.figure()

        pc.append(principal_components)
    return np.vstack(pc).T

def get_spws(data, traces, electr, window, fs, start_pts, remove_baseline = False, window_pts = [0, 0]):
    data_used = []
    win_pts = [ispw.ms2pts(window[0], fs), ispw.ms2pts(window[1], fs)]
    for spw in range(len(traces)):
        t = start_pts[spw]
        
        data_temp = data[electr, traces[spw], t + win_pts[0]:t + win_pts[1]]

        data_used.append(data_temp)
    data_used = np.transpose(np.array(data_used))   
    
    if remove_baseline:
        baselines = np.zeros(len(start_pts))
        for idx, trace in enumerate(traces):
            #import pdb; pdb.set_trace() 
            baseline = data[electr, trace, window_pts[0][idx] : window_pts[1][idx]]
            baseline = np.mean(baseline)
            baselines[idx] = baseline
        
        #import pdb; pdb.set_trace() 
        data_used = data_used - baselines
    return data_used

def compare_spws(data, traces, electrodes, window, fs, start_pts, labels, window_base = [0,0]): 
    plt.figure()    
    n_electr = len(electrodes)
    colors = ['r','b','g']
    s_p_temp1 = ispw.ms2pts(window_base[0], fs).astype('i4') + start_pts 
    s_p_temp2 = ispw.ms2pts(window_base[1], fs).astype('i4') + start_pts
    window_pts = [s_p_temp1, s_p_temp2]
    
    remove_baseline = True
    for i,el in enumerate(electrodes):
        plt.subplot(n_electr,1,i+1)
        spw_traces = get_spws(data, traces, el, window, fs, start_pts, remove_baseline, window_pts)
        for l in np.unique(labels):
            plt.plot(spw_traces[:, labels==l], colors[l],  alpha=0.2)

    
def compare_clusters(sc, labels, names=None):
    n_fets = sc.shape[0]
    plt.figure()
    colors = ['r', 'b', 'g']
    ulabs = np.unique(labels)
    for i in range(n_fets):
        for j in range(n_fets):
            f1 = sc[i,:]
            f2 = sc[j,:]
            plt.subplot(n_fets, n_fets, i+j*n_fets+1)
            if i==j:
                for l in ulabs:
                    plt.hist(f1[labels==l], 30, histtype='stepfilled', 
                             normed=True,alpha=0.4, color=colors[l])
            else:
                for l in ulabs:                                
                    plt.plot(f1[labels==l],f2[labels==l], colors[l]+'.')
            if i==0 and names is not None:
                plt.ylabel(names[j])
            if j==0 and names is not None:
                plt.title(names[i])
                            

def remove_with_less_ipsps(save_folder, save_file, spw_data ,min_ipsps_group):
    npzfile        = np.load(save_folder + spw_data)
    spws = npzfile['spw_ipsps']
    npzfile.close()          
   
    new_spw = []
    for spw_no in np.unique(spws['spw_no']):
        spw_used = spws[spws['spw_no'] == spw_no]
        groups = np.unique(spw_used['group'])
        number_groups = len(groups)
        if len(min_ipsps_group) == 2: 
            if number_groups >= min_ipsps_group[0] and number_groups <= min_ipsps_group[1]:
                new_spw.append(spw_used)
        else:
            if number_groups >= min_ipsps_group[0]:
                new_spw.append(spw_used)
    new_spw = np.concatenate(new_spw)
    np.savez(save_folder + save_file, spw_ipsps = new_spw) 

def plot_dendograms(save_folder, plot_folder, plot_file, data_file, spw_groups, spw_details,
                    spike_data , ext, win):
    #import scipy as sc
    #import pylab
    import scipy.cluster.hierarchy as sch
    #import clust.hierarchy as sch
    #import numpy as np
    npzfile        = np.load(save_folder + spw_groups)
    matrix = npzfile['group1']
    npzfile.close()  
    D = matrix['group']
    #import pdb; pdb.set_trace()
    #try:
    #    spws = [npzfile['spw_ipsps']]
    #    types = ['all']
    #    double = False
    #except:
    #    spont = npzfile['spontaneous']
    #    init = npzfile['initiated']
    #    types = ['spontaneous', 'initiated']
    
    #D = np.genfromtxt('LtoR.txt', dtype=None)
    def llf(id):
        return str(id)
    fig = plt.figure(figsize=(10,10))
    Y = sch.linkage(D, method='single')
    Z1 = sch.dendrogram(Y,leaf_label_func=llf,leaf_rotation=90)
    fig.show()
    fig.savefig('dendrogram.png')


def plot_spikes(save_folder, save_name, distances, ext = '.png'):
    
    # plot histogram
    fig = plt.figure()
    plt.hist(distances, normed = True)
    plt.title('normal, blue = spont, green = init')
    fig.savefig(save_folder + 'init_distance' + ext, dpi=600)  
    
    import pdb; pdb.set_trace()
    
    # errorbars
    fig, ax = plt.subplots()
    fig.canvas.draw()
    plt.boxplot(np.transpose(np.array(distances)))
    labels = [item.get_text() for item in ax.get_xticklabels()]
    labels[0] = 'spontaneous'
    labels[1] = 'initiated'
    ax.set_xticklabels(labels)
    plt.ylabel('Distance of from SPW initiation site* from intracellular electrode')
    
    fig.savefig(save_folder + 'init_distance_box' + ext, dpi=900)  
    import pdb; pdb.set_trace()
    #plt.show()

def plot_all_cum_change_var(plot_folder, plot_file, 
                            all_var_spont, all_var_init, timeline, fs, 
                            ext = '.png'):
    plot_it = True

    # save the data in the txt files
    all_var_spont_array = [np.array(all_var_spont[var]) for var in range(len(all_var_spont))]
    all_var_init_array = [np.array(all_var_init[var]) for var in range(len(all_var_init))]
    np.savetxt(plot_folder + 'cum_spontaneous.txt', np.array(all_var_spont_array),delimiter='\t')
    np.savetxt(plot_folder + 'cum_initiated.txt', np.array(all_var_init_array),delimiter='\t')
    
    
#    import csv
#    
#    #cum_spont= csv.writer(open(plot_folder + "MYFILE.csv", "wb"))
#    #spamWriter.writerow([1,2,3])
#    
#    file = open(plot_folder + "MYFILE3.csv", "wb")
#    #fileWriter = csv.writer(file , delimiter='\n',quotechar='|', quoting=csv.QUOTE_MINIMAL)
#    #fileWriter.writerow([1,2,3])
#    #spamWriter = csv.writer(file , delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
#    spamWriter = csv.writer(file)
#    for row in range(len(all_var_spont)):
#    #spamWriter = csv.writer(file , delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
#        import pdb; pdb.set_trace()
#        spamWriter.writerow(all_var_spont[row].tolist())
#    file.close()
    
    #temp = np.loadtxt(plot_folder + 'cum_spontaneous.txt')

    from scipy.stats import nanmean
    #plt.plot(np.transpose(all_var_spont), 'b', label = 'spontaneous', alpha = 0.3, ms = 12)
    #plt.plot(np.transpose(all_var_init), 'g', label = 'initiated', alpha = 0.3, ms = 12)
    all_vars = [all_var_spont, all_var_init]
    names = ['spontaneous', 'initiated']
    colors = ['b', 'g']
    line_widt = 4
    
    for idx, var in enumerate(all_vars):
        # every result separately
        #plt.plot(timeline, np.transpose(var), colors[idx], alpha =0.2, lw = 3)
        # all the mean
        plt.plot(timeline, np.transpose(np.mean(var, 0)), colors[idx], label = names[idx], lw = line_widt)
        box_pos = [10, 20, 30, 40, 50, 60]
        vs = [var[:,int(ispw.ms2pts(mils, fs))] for mils in box_pos]
        plt.boxplot(vs, positions=box_pos) 
    plt.legend(loc=2)
    plt.xlim([min(timeline), max(timeline)])
    plt.xlabel('Time (ms)')
    plt.ylabel('Cumulative change of variance')
    #import pdb; pdb.set_trace()
    
    #plt.box
    #plt.boxplot(vs, positions=box_pos) 
    #plt.er
    #plt.plot(t, nanmean(all_cums[:, typ, :], 0), colors[typ], label = types[typ])
    #if plot_it: 
    #    plt.show()
    plt.savefig(plot_folder +plot_file + ext, dpi=600) 
    if plot_it:
        plt.show()
     

    plt.clf()
    gc.collect() 

    #import pdb; pdb.set_trace()


    
 
def create_scatter_synch(ampl, synch, group, name, save_file, ext = '.png'):
    # plots the scatter plot
    
    # define variables
    
    #groups_for_colors = np.array([1, 1, 2, 2, 3, 3, 4, 4]) # less groups, no 1 IPSP
    #ticks_labels = ['2-3 IPSPs','4-5 IPSPs', '6-7 IPSPs', '8 or more'] # less groups, no 1 IPSP
    
    #groups_for_colors = np.array([1, 2, 2, 3, 3, 4, 4, 5, 5]) # less groups, with 1 IPSP
    #ticks_labels = ['1 IPSPs','2-3 IPSPs','4-5 IPSPs', '6-7 IPSPs', '8 or more'] # less groups, with 1 IPSP
    
    #groups_for_colors = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8]) # more groups, no 1 IPSP
    #ticks_labels = ['2 IPSPs','3 IPSPs', '4 IPSPs', '5 IPSPs','6 IPSPs', '7 IPSPs', '8 or more']
    
    groups_for_colors = np.array([1, 2, 3, 4, 5, 6, 7, 8]) # more groups, with 1 IPSP
    ticks_labels = ['1 IPSP', '2 IPSPs','3 IPSPs', '4 IPSPs', '5 IPSPs','6 IPSPs', '7 IPSPs', '8 or more']
    
    max_color = max(groups_for_colors) * 1.0
    groups_for_colors = groups_for_colors / max_color
    
    # define colorbar ticks
    #ticks = [0,1,2, 3, 4, 5]
    ticks = np.linspace(min(groups_for_colors), max(groups_for_colors), len(np.unique(groups_for_colors)))
    
    #ticks_labels = ['1 IPSPs','2-3 IPSPs','4-5 IPSPs', '6-7 IPSPs', '8 or more']
   
    #ticks_labels = ['2 IPSPs','3 IPSPs', '4 IPSPs', '5 IPSPs','6 IPSPs', '7 IPSPs', '8 or more']
    #ticks_labels = ['1 IPSP', '2 IPSPs','3 IPSPs', '4 IPSPs', '5 IPSPs','6 IPSPs', '7 IPSPs', '8 or more']
    #ticks_labels = ticks_labels[::-1]
    marker_size = 60
    import matplotlib as mpl
    #import pdb; pdb.set_trace() 
    
    fig = plt.figure()

    group[group > len(groups_for_colors) -1] = len(groups_for_colors) - 1
    colors_group = groups_for_colors[group-1] 
    #plt.scatter(ampl, synch, c = colors_group, s = marker_size, cmap=mpl.cm.gray)
    plt.scatter(ampl, synch, c = colors_group, s = marker_size) #, alpha = 0.2)
    cbar = plt.colorbar()
    cbar.set_ticks(ticks)
    cbar.set_ticklabels(ticks_labels)

    plt.ylabel('synchrony [(no_ipsps all together)/(no_of_electrodes * no_group_ipsp)]')
    plt.xlabel('amplitude [micro V]')
    plt.title(name + ' ,no of SPWs: ' + str(len(ampl)))
    #import pdb; pdb.set_trace()      
    #fig.colorbar(im)       
    plt.ylim([0,1.05])
    plt.xlim([0,1200])
    plt.savefig(save_file + name + ext, dpi=600) 
    
     
    
def plot_amplitude_vs_synchrony_all(plot_folder, plot_file, cells, 
                                                 amplitudes, synchronise, 
                                                 group_nos, names, ext = '.pdf'):
    plot_it = True
    #import pdb; pdb.set_trace()
    all_cells = np.unique(cells)
   
    #import pdb; pdb.set_trace()
    cells = np.array(cells)
    
    #amplitudes = np.array(amplitudes)
    #synchronise = np.array(synchronise)
    spw_nos_group = 0
    linestyl = '.'
    all_ampl = []
    all_synch = []
    all_group = []

    for group in range(len(names)):
        
        ampl_group = amplitudes[group]
        ampl_group = np.concatenate(ampl_group)
        all_ampl.append(ampl_group)
        
        synch_group = synchronise[group]
        synch_group = np.concatenate(synch_group)
        all_synch.append(synch_group)
        
        group_group = group_nos[group]
        group_group = np.concatenate(group_group)
        #import pdb; pdb.set_trace()
        all_group.append(group_group) 
        create_scatter_synch(ampl_group, synch_group, group_group, names[group], plot_folder +plot_file, ext)
        #spw_nos_group = spw_nos_group + len(ampl)
    np.savetxt(plot_folder + 'all_synch_transposed.txt', np.transpose(all_synch),delimiter='\t')  
    np.savetxt(plot_folder + 'all_ampl_transposed.txt', np.transpose(all_ampl),delimiter='\t')   
    import pdb; pdb.set_trace()
    all_ampl = np.concatenate(all_ampl)
    all_synch = np.concatenate(all_synch)
    all_group = np.concatenate(all_group)
    create_scatter_synch(all_ampl, all_synch, all_group, 'spont_init_', plot_folder +plot_file, ext)
    plt.close()
    
    plt.figure() 
        
    # fit the data points
    fitfunc = lambda p, x: p[0]* np.sqrt(x) +p[1] # Target function
    all_ampl = all_ampl[all_ampl >= 0 ]
    all_synch = all_synch[all_ampl >= 0]
    
    plot_fit(fitfunc, all_ampl, all_synch, guess_values = [0.1, 0.2]) 

    # linear function only
    fitfunc = lambda p, x: p[0]*x +p[1] # Target function
    plot_fit(fitfunc, all_ampl, all_synch, guess_values = [0.1, 0.2]) 

    fitfunc = lambda p, x: x * p[0] # Target function
    plot_fit(fitfunc, all_ampl, all_synch, guess_values = [0.1]) 
        
    if plot_it: 
        plt.show()
    import pdb; pdb.set_trace()
    plt.clf()
    gc.collect()    
    
def plot_fit(fitfunc, x_data, y_data, guess_values):
    
    from scipy import optimize
    plt.figure()
    errfunc = lambda p, x, y: fitfunc(p, x) - y
 
    idxs = range(0,len(x_data))
    
    sort_idx = np.argsort(x_data)
    x_data = x_data[sort_idx]
    y_data = y_data[sort_idx]

    p1, success = optimize.leastsq(errfunc, guess_values, args=(x_data, y_data))
    #import pdb; pdb.set_trace()
    plt.plot(x_data, y_data, "ro")
    plt.plot(x_data, fitfunc(p1, x_data), 'r-')
    plt.ylim([0, 1.1])
    #plt.show()      
        

def plot_amplitude_vs_synchrony(save_folder, save_file, plot_folder,plot_file, data_file, spw_groups,spw_details, ext):
    
    npzfile        = np.load(save_folder + data_file)
    data = npzfile['data']
    fs = npzfile['fs']
    npzfile.close() 
    
    plot_it = False
    
    npzfile = np.load(save_folder + spw_details)
    try:
        
        spws = npzfile['spw_ipsps']
        types = ['all']
        npzfile.close()    
        
        names = [0]    
        groups = [0]  
         
    except:
        spont = npzfile['spontaneous']
        init = npzfile['initiated']
        types = ['spontaneous', 'initiated']
        spws = [spont, init]     
        npzfile.close()

        npzfile = np.load(save_folder + spw_groups)
        group1 = npzfile['group1']
        group2 = npzfile['group2']
        groups = [group1, group2]
        names = npzfile['names']
    npzfile.close()       
    # go through every type possible
    remove_baseline = True
    win_base = [-3, -1]
    win_base_pts = [ispw.ms2pts(win_base[0], fs),ispw.ms2pts(win_base[1], fs)]
    save_fold = save_folder + plot_folder
    fold_mng.create_folder(save_fold)
    save_base = save_fold + plot_file
    
    window_to_plot = [0, 5] #win
    win_pts = [ispw.ms2pts(window_to_plot[0], fs), ispw.ms2pts(window_to_plot[1], fs)]
    size_win_pts = win_pts[1] - win_pts[0] + 1
     
     
    #no_of_colors = 5
    #groups_for_colors = np.array([0, 0, 1, 1, 2, 2, 3, 3, 4, 4])
    #cols = np.array(define_colors(no_colors = no_of_colors + 1, type = 'grey'))
    #cols = cols[:-1]
    #cols = cols[::-1]
    #cols = np.array(cols)
    
    #import pdb; pdb.set_trace() 
    all_ampls = []
    all_syncs = []
    all_groups = []
    for typ in range(len(names)):
        spw_group_typ = groups[typ]
        #import pdb; pdb.set_trace() 
        if spw_group_typ != 0:
            group_nos = np.unique(spw_group_typ['group'])
            spw_type = spws[typ]
        else:
            group_nos = [0]
            spw_type = spws
        temp_ampls = []
        temp_syncs = []
        temp_groups = []
        # go through every group detected
        for group_no in group_nos:    

            if spw_group_typ != 0:
                spw_nos_used = spw_group_typ[spw_group_typ['group'] == group_no]['spw_no']
            else:
                spw_nos_used = np.unique(spws['spw_no'])
            hist_electr_all = []
            all_spws = np.zeros([len(spw_nos_used), np.size(data,0), size_win_pts])
            
            no_bins = 150
            electro_bins = np.zeros([np.size(data,0), no_bins])
            all_sync = []
            all_ampl = []
            all_group = []
            # go through every spw used in this group
            for idx, spw_no in enumerate(spw_nos_used):
                spw_used = spw_type[spw_type['spw_no'] == spw_no]
                spw_start = spw_used['spw_start'][0]
                last_ipsp = np.max(spw_used['ipsp_start'])
                trace = spw_used['trace'][0]
                data_used = data[:,trace,:]

                
                start_trace = spw_start + window_to_plot[0]
                end_trace = last_ipsp + window_to_plot[1]

                 

                start_trace_pts = ispw.ms2pts(start_trace, fs).astype('i4')
                end_trace_pts = ispw.ms2pts(end_trace, fs).astype('i4')

                data_spw = data_used[:, start_trace_pts: end_trace_pts]
                #import pdb; pdb.set_trace() 
                if remove_baseline:
                    base = data_used[:, start_trace_pts + win_base_pts[0]: start_trace_pts + win_base_pts[1]]
                    base = np.mean(base, axis = 1)
                    data_spw = np.transpose(data_spw) - base
                    data_spw = np.transpose(data_spw)
                    
                #import pdb; pdb.set_trace()
                # calculate amplitude
                ampls = max(np.max(data_spw, 1))
                if ampls < 0:
                    print 1
                # calculate synchrony
                # (no_ipsps all together)/(no_of_electrodes * no_group_ipsp)
                no_ipsps = len(spw_used)
                all_electr = np.size(data, 0) * 1.0
                no_group_ipsp = len(np.unique(spw_used['group']))
                #import pdb; pdb.set_trace()
                sync = (no_ipsps * 1.0)/(all_electr * no_group_ipsp * 1.0)
                
                
                
                all_ampl.append(ampls)
                all_sync.append(sync)
                all_group.append(no_group_ipsp)
                #import pdb; pdb.set_trace() 
            
            #import pdb; pdb.set_trace()
            #group_group = np.copy(all_group)
            #group_group[group_group > len(groups_for_colors) -1] = len(groups_for_colors) -1
            #colors_group = groups_for_colors[group_group]
            #all_group.append(colors_group)
            #import pdb; pdb.set_trace()
            name  = '_group_' + str(group_no) + '_' + types[typ]
            #import pdb; pdb.set_trace()
            create_scatter_synch(all_ampl, all_sync, np.array(all_group), name, save_base, ext = '.png')
            
            
            #plt.scatter(all_ampl, all_sync, color = cols[colors_group,:], s = 3)
            #import pdb; pdb.set_trace() 
            #im = ax1.scatter(ampl_group, synch_group, color = cols[colors_group,:], s = marker_size)
            #plt.plot(all_ampl, all_sync, '.', alpha = 0.4)
            #plt.ylabel('synchrony [(no_ipsps all together)/(no_of_electrodes * no_group_ipsp)]')
            #plt.xlabel('amplitude [micro V]')
            #plt.title('Group: ' + str(group_no) + ', ' + types[typ] + ',no of SPWs: ' + str(len(spw_nos_used)))
            #plt.xlim([0, 1600])
            #plt.ylim([0, 1])
            temp_ampls.append(all_ampl)
            temp_syncs.append(all_sync)
            temp_groups.append(all_group)
        #import pdb; pdb.set_trace() 
        all_syncs.append(np.concatenate(np.array(temp_syncs)))
        all_ampls.append(np.concatenate(np.array(temp_ampls)))
        all_groups.append(np.concatenate(np.array(temp_groups)))
        #import pdb; pdb.set_trace() 
        #fig.savefig(save_base + '_group_' + str(group_no) + '_' + types[typ] + ext, dpi=600)   
        #plt.savefig(save_file + name + ext, dpi=600) 
    if plot_it: 
        plt.show() 
    plt.close()
    #import pdb; pdb.set_trace() 
    np.savez(save_folder + save_file, group = types, all_ampls = all_ampls, all_syncs = all_syncs, groups_ipsp = all_groups) 
    gc.collect()    
    
    
def cum_distribution_funct(save_folder, save_file, plot_folder, plot_file, data_file, spw_details, ext, win):
    """ it calculate mean root mean square of each of the SPWs (separately for spontaneous and induced,
    and from this create comulative distribution function for each of the groups"""
    npzfile        = np.load(save_folder + data_file)
    data = npzfile['data']
    fs = npzfile['fs']
    npzfile.close() 
    
    npzfile = np.load(save_folder + spw_details)
    spont = npzfile['spontaneous']
    init = npzfile['initiated']
    types = ['spontaneous', 'initiated']
    spws = [spont, init]     
    npzfile.close()    
    
    # make sure that there is equal number of spontaneous and initiated SPWs
    assert len(np.unique(init['spw_no'])) == len(np.unique(spont['spw_no']))
    
    remove_mean = True
    remove_baseline = True
    
    win_base = [-10, -5]
    win_base_pts = [ispw.ms2pts(win_base[0], fs),ispw.ms2pts(win_base[1], fs)]
    save_fold = save_folder + plot_folder
    fold_mng.create_folder(save_fold)
    save_base = save_fold + plot_file
    
    window_to_plot = [0, 70] #win
    win_pts = [ispw.ms2pts(window_to_plot[0], fs), ispw.ms2pts(window_to_plot[1], fs)]
    size_win_pts = win_pts[1] - win_pts[0]
    
    no_of_spws_in_group =  len(np.unique(spont['spw_no']))
    all_cums = np.zeros([len(data), len(types), size_win_pts])

    # go through every type possible
    for typ in range(len(types)):
        spws_used = spws[typ]   
        # go through every group detected  
        all_spws = np.ones([no_of_spws_in_group, np.size(data,0), size_win_pts]) * np.nan
        spw_nos_used = np.unique(spws_used['spw_no'])
        #no_bins = 150
        #electro_bins = np.zeros([np.size(data,0), no_bins])
        
        # go through every spw used in this group
        for idx, spw_no in enumerate(spw_nos_used):
            spw_used = spws_used[spws_used['spw_no'] == spw_no]
            spw_start = spw_used['spw_start'][0]
            trace = spw_used['trace'][0]
            data_used = data[:,trace,:]
            
            start_trace = spw_start + window_to_plot[0]
            end_trace = spw_start + window_to_plot[1]
            
            # plot every trace of this spw
            start_trace_pts = ispw.ms2pts(start_trace, fs).astype('i4')
            end_trace_pts = ispw.ms2pts(end_trace, fs).astype('i4')
            
            data_spw = data_used[:, start_trace_pts: end_trace_pts]
            #import pdb; pdb.set_trace() 
            if remove_baseline:
                base = data_used[:, max(0, start_trace_pts + win_base_pts[0]): max(20, start_trace_pts + win_base_pts[1])]
                base = np.mean(base, axis = 1)
                data_spw = np.transpose(data_spw) - base
                data_spw = np.transpose(data_spw)

            all_spws[idx, :, 0:np.size(data_spw,1)] = data_spw[:, 0:size_win_pts]              
            
        all_root_means = np.zeros([len(data_spw), np.size(all_spws, 2)])   
        all_root_meaned = np.zeros([len(data_spw), np.size(all_spws, 2)]) 
        # calculate variablility - cumulative sum - of all spws in this type
        from scipy.stats import nanmean
        for electr in range(len(data_spw)):
            # normal equation
            #import pdb; pdb.set_trace() 
            if not remove_mean:
                s_mean_across = nanmean(all_spws[:, electr, :], 0) # mean across all the spq in this electrode
                variance = all_spws[:, electr, :] - s_mean_across[None, :] # subtracts mean from each SPW
                squared = variance ** 2 # power of every point
                meaned = nanmean(squared, 0) # calculate mean from the powers
                sqruted = np.sqrt(meaned) # sqrt of 
                all_root_means[electr, :] = np.cumsum(sqruted) 
            else:
               
                # normalised by amplitude of mean
                s_mean_across = nanmean(all_spws[:, electr, :], 0) # mean across all the spq in this electrode
                max_mean = max(s_mean_across)
                normalize_by_mean = all_spws[:, electr, :] / max_mean
                #import pdb; pdb.set_trace() 
                s_mean_across = s_mean_across/max_mean
                #import pdb; pdb.set_trace() 
                variance = normalize_by_mean - s_mean_across[None, :] # subtracts mean from each SPW
                squared = variance ** 2 # power of every point
                meaned = nanmean(squared, 0) # calculate mean from the powers
                sqruted = np.sqrt(meaned) # sqrt of 
                all_root_means[electr, :] = np.cumsum(sqruted)             

            #all_root_means[electr, :] = all_root_means[electr, :]/ all_root_means[electr, -1]           
        # save all the root mean squares for this function    
        all_cums[:, typ, :] = all_root_means
    colors = ['b', 'g']
    
    #from scipy.stats import ks_2samp
    t = dat.get_timeline(all_cums[0, 0, :], fs, 'ms') + window_to_plot[0]
    
    for electr in range(len(data_spw)):
        fig = plt.figure()
        #plt.subplot(len(data_spw), 1, electr + 1 )
        for typ in range(len(types)):
            plt.plot(t, all_cums[electr, typ, :], colors[typ], label = types[typ])
        plt.legend()
        plt.xlabel('Time (ms)')
        plt.ylabel('Cumulative change of variance')
        
        if remove_mean:
            plt.title('no of SPWs: ' + str(len(spw_nos_used)) + 'corrected by mean , electrode: ' + str(electr))
            fig.savefig(save_base + 'corrected_electr_' + str(electr) + ext, dpi=600)
            
        else:
            plt.title('no of SPWs: ' + str(len(spw_nos_used)) + ', electrode: ' + str(electr))
            fig.savefig(save_base + '_electr_' + str(electr) + ext, dpi=600)

        plt.close()
        
    fig = plt.figure()
    for typ in range(len(types)):
        plt.plot(t, nanmean(all_cums[:, typ, :], 0), colors[typ], label = types[typ])
    plt.xlabel('Time (ms)')
    plt.legend()    
    plt.ylabel('Cumulative change of variance')  
    if remove_mean:
        plt.title('Corrected, no of SPWs: ' + str(len(spw_nos_used)) + ', mean of all electrodes')
        fig.savefig(save_base + '_corrected_all_'+ ext, dpi=600) 
        
    else:
        plt.title('No of SPWs: ' + str(len(spw_nos_used)) + ', mean of all electrodes')
        fig.savefig(save_base + '_all_'+ ext, dpi=600) 
    
    np.savez(save_folder + save_file, cum_change_spont = nanmean(all_cums[:, 0, :]), cum_change_init = nanmean(all_cums[:, 1, :]), timeline = t, fs = fs) 
    
    #plt.show()
    #import pdb; pdb.set_trace() 
    plt.close()    
    del all_root_meaned, all_root_means
    gc.collect()   

def plot_groups_w_fr(save_folder, plot_folder, plot_file, data_file, spw_groups, spw_details, spike_data, ext, win):
    """ makes the plot of every given group and finds the firing rate for it"""
    npzfile        = np.load(save_folder + data_file)
    data = npzfile['data']
    fs = npzfile['fs']
    npzfile.close() 
    
    npzfile        = np.load(save_folder + spike_data)
    spike_idx = npzfile['spike_idx']
    npzfile.close()
    
    
    npzfile = np.load(save_folder + spw_details)
    try:
        spws = [npzfile['spw_ipsps']]
        types = ['all']
        npzfile.close()     
        
        npzfile = np.load(save_folder + spw_groups)
        groups = [npzfile['group']]
        names = npzfile['names']       
         
    except:
        spont = npzfile['spontaneous']
        init = npzfile['initiated']
        types = ['spontaneous', 'initiated']
        spws = [spont, init]     
        npzfile.close()

        npzfile = np.load(save_folder + spw_groups)
        group1 = npzfile['group1']
        group2 = npzfile['group2']
        groups = [group1, group2]
        names = npzfile['names']
    npzfile.close()     
     
    remove_baseline = True
    win_base = [-10, -5]
    win_base_pts = [ispw.ms2pts(win_base[0], fs),ispw.ms2pts(win_base[1], fs)]
    save_fold = save_folder + plot_folder
    fold_mng.create_folder(save_fold)
    save_base = save_fold + plot_file
    
    window_to_plot = [-10, 70] #win
    win_pts = [ispw.ms2pts(window_to_plot[0], fs), ispw.ms2pts(window_to_plot[1], fs)]
    size_win_pts = win_pts[1] - win_pts[0] + 1
    
    # go through every type possible
    for typ in range(len(names)):
        spw_group_typ = groups[typ]
        group_nos = np.unique(spw_group_typ['group'])
        spw_type = spws[typ]
        
        # go through every group detected
        for group_no in group_nos:
            
            fig = plt.figure()
            spw_nos_used = spw_group_typ[spw_group_typ['group'] == group_no]['spw_no']
            all_spws = np.zeros([len(spw_nos_used), np.size(data,0), size_win_pts])
            
            no_bins = 150
            electro_bins = np.zeros([np.size(data,0), no_bins])
            # go through every spw used in this group
            for idx, spw_no in enumerate(spw_nos_used):
                spw_used = spw_type[spw_type['spw_no'] == spw_no]
                spw_start = spw_used['spw_start'][0]
                trace = spw_used['trace'][0]
                data_used = data[:,trace,:]
                spikes_used = spike_idx[spike_idx['trace'] == trace]
                
                start_trace = spw_start + window_to_plot[0]
                end_trace = spw_start + window_to_plot[1]

                bins = np.linspace(start_trace, end_trace, no_bins + 1)
                
                # plot every trace of this spw
                #import pdb; pdb.set_trace() 
                for electr in range(np.size(data,0)):
                    
                    # find spikes in this electrode in this trace
                    spikes_electr = spikes_used[spikes_used['electrode'] == electr]
                    spikes_spw = spikes_electr[(spikes_electr['time']>=start_trace) & (spikes_electr['time']<=end_trace)]['time']
                    
                    #if len(spikes_spw) > 0:
                    #    import pdb; pdb.set_trace() 
                    in_bin, _ = np.histogram(spikes_spw, bins)
                    #import pdb; pdb.set_trace() 
                    electro_bins[electr, :] = electro_bins[electr, :] + in_bin

                start_trace_pts = ispw.ms2pts(start_trace, fs).astype('i4')
                end_trace_pts = ispw.ms2pts(end_trace, fs).astype('i4')
                add_it = 300
                

                data_spw = data_used[:, start_trace_pts: end_trace_pts]
                #import pdb; pdb.set_trace() 
                if remove_baseline:
                    base = data_used[:, max(0, start_trace_pts + win_base_pts[0]): max(20, start_trace_pts + win_base_pts[1])]
                    base = np.mean(base, axis = 1)
                    data_spw = np.transpose(data_spw) - base
                    data_spw = np.transpose(data_spw)
                    
                t = dat.get_timeline(data_spw[0,:], fs, 'ms')
                # plot data trace 
                #import pdb; pdb.set_trace() 
                for electr in range(np.size(data, 0)):
                    plt.plot(t, data_spw[electr,:] + add_it * electr, color = 'k', alpha = max(0.05, 1.0/len(spw_nos_used))) #, lw = 4)
                

                all_spws[idx, :, 0:np.size(data_spw,1)] = data_spw

            #import pdb; pdb.set_trace()  
            #all_spws[all_spws == 0] = 0.01   
            from scipy.stats import nanmean
            mean_spw = nanmean(all_spws, 0)
            
            #mean_spw = np.mean(all_spws, 0)
            #import pdb; pdb.set_trace()  
            t = dat.get_timeline(mean_spw[0,:], fs, 'ms')
            for electr in range(np.size(data, 0)):
                #import pdb; pdb.set_trace() 
                plt.plot(t, mean_spw[electr,:] + add_it * electr, color = 'r', lw = 2)

               

                    #plt.boxplot(all_spws[:,electr,ispw.ms2pts(mils, fs)], positions=[mils])
                #import pdb; pdb.set_trace() 
                
            #import pdb; pdb.set_trace() 
            #plt.plot(t, np.mean(data_spw[electr,:]) + add_it * electr, color = 'r')
            
            bar_lin = np.linspace(0, t[-1], no_bins + 1)
            bar_width = bar_lin[1]-bar_lin[0]
            # for this group plot the histogram of the spikes
            for electr in range(np.size(data,0)):
                plt.bar(bar_lin[:-1], (electro_bins[electr,:]/len(spw_nos_used))*50, bottom = add_it * electr, alpha = 0.8, width = bar_width)  
             


                
                
                 
            spike_distribution = np.sum(electro_bins, 0)
            #import pdb; pdb.set_trace() 
            plt.bar(bar_lin[:-1], (spike_distribution/len(spw_nos_used))*50, bottom = add_it * (-1), width = bar_width) #, alpha = 0.7)  
            plt.xlabel('time (ms)')
            plt.title('Group: ' + str(group_no) + ', ' + types[typ] + ',no of SPWs: ' + str(len(spw_nos_used)))
            #import pdb; pdb.set_trace()
            plt.xlim([t[0], t[-1]])
            fig.savefig(save_base + '_group_' + str(group_no) + '_' + types[typ] + ext, dpi=600)    
            
            
            #plt.show()             
            # for drawing box plot
#            add_it = 500
#            plt.figure()
#            for electr in range(np.size(data,0)):   
#                box_pos = [10,20,30, 40, 50, 60, 70]
#                vs = [all_spws[:,electr,ispw.ms2pts(mils, fs)]+add_it*electr for mils in box_pos]
#                plt.boxplot(vs, positions=box_pos)
#            plt.title(types[typ])
#            plt.xlabel('time (ms)')      
    #plt.show()

def calculate_dist_extra_spikes_to_intra_spike(intra_spikes, extra_spikes):    
    distance = []
    
    for trace in np.unique(intra_spikes['trace']):
        intra_spikes_used = intra_spikes[intra_spikes['trace'] == trace]
        spikes_used = extra_spikes[extra_spikes['trace'] == trace]
        
        if len(intra_spikes_used['time']) > 1:
            # only one SPW in this trace
            for spik in range(len(intra_spikes_used['time'])):
                if spik ==0:
                    # it's the first one
                    current_spik = intra_spikes_used[spik]['time']
                    next_spik = intra_spikes_used[spik + 1]['time']
                    spiki_used = spikes_used[spikes_used['time'] < next_spik]
                    
                    
                elif spik != len(intra_spikes_used['time']) - 1:

                    current_spik = intra_spikes_used[spik]['time']
                    next_spik = intra_spikes_used[spik + 1]['time']
                    
                    spiki_used = spikes_used[(spikes_used['time'] >= current_spik) & (spikes_used['time'] < next_spik)]

                else:
                    # it's the last spike
                    current_spik = intra_spikes_used[spik]['time']
                    spiki_used = spikes_used[spikes_used['time'] >= current_spik]
                dist = spiki_used['time']  -current_spik
                distance.append(dist)
            #import pdb; pdb.set_trace()
        else:       
            dist = spikes_used['time'] - intra_spikes_used['time']
            distance.append(dist)
        
    distance = np.concatenate(distance)
    return distance

def plot_fr_after_spike_and_distances_after_spike(save_folder, plot_folder, 
                                      plot_file, intra_spikes, dist_file,
                                      spike_data , ext):
    """ plots the histogram of the distribution of spws from the spike
    and the histogram of firing rate after spike"""
    npzfile         = np.load(save_folder + dist_file)
    dist = npzfile['dist_spwspike']
    distance = dist['distance']
    npzfile.close()

    npzfile        = np.load(save_folder + intra_spikes)
    intra_spikes = npzfile['spikes_first']
    npzfile.close()
    
    npzfile        = np.load(save_folder + spike_data)
    spikes = npzfile['spike_idx']
    npzfile.close()
    
    xlim = [-3, 100]
    no_bins = 200
    distance = distance[(distance > xlim[0]) & (distance < xlim[1])]
    
    distance_spike = calculate_dist_extra_spikes_to_intra_spike(intra_spikes, spikes)
    distance_spike = distance_spike[(distance_spike > xlim[0]) & (distance_spike < xlim[1])]
    #import pdb; pdb.set_trace() 

    fig = plt.figure()   
    plt.hist(distance_spike, no_bins, xlim, normed=1, color = 'b', alpha = 0.7, label = 'firing rate')             
    plt.hist(distance, no_bins, xlim, normed=1,  color = 'k', alpha = 0.6, label = 'SPW beginnings')
    plt.legend()
    
    plt.title('Distribution of SPWs from the spike')
    plt.xlim(xlim)
    plt.ylabel('Fraction of SPWs')
    plt.xlabel('Distance from spike (ms)')
    fig.savefig(save_folder + plot_folder + plot_file + ext,dpi=600)         
    #plt.show() 
    plt.close()  
   


def plot_fr_after_spike(save_folder, plot_folder, 
                                      plot_file, intra_spikes,
                                      spike_data , ext, win):
    npzfile        = np.load(save_folder + intra_spikes)
    intra_spikes = npzfile['spikes_first']
    npzfile.close()
    
    npzfile        = np.load(save_folder + spike_data)
    spikes = npzfile['spike_idx']
    npzfile.close()
    
    xlim = [-3, 70]
    no_bins = 200
    plt.rcParams.update({'font.size': 22})


    distance = calculate_dist_extra_spikes_to_intra_spike(intra_spikes, spikes)
    
    distance = distance[(distance > xlim[0]) & (distance < xlim[1])]
    #import pdb; pdb.set_trace() 
    fig = plt.figure()                
    plt.hist(distance, no_bins, normed=1, color = 'k')
    plt.title('Distribution of extracelluar spikes from the intracellular spike')
    plt.xlim(xlim)
    plt.ylabel('Fraction of SPWs')
    plt.xlabel('Distance from spike (ms)')
    fig.savefig(save_folder + plot_folder + plot_file + ext,dpi=600)         
    plt.show() 
    plt.close()     
    #import pdb; pdb.set_trace()
     



def plot_spw_ipsps_no_groups_all(save_folder, save_file, data_file, spw_data, ext):
    """ similar to plot_spw_ipsps_no_groups but does not divide first group into
    origin of the first ipsp"""
    
    npzfile        = np.load(save_folder + spw_data)
    try:
        spws = [npzfile['spw_ipsps']]
        types = ['all']
        double = False
    except:
        spont = npzfile['spontaneous']
        init = npzfile['initiated']
        types = ['spontaneous', 'initiated']
        spws = [spont, init]    
        double = True 
    npzfile.close()  
    
    npzfile        = np.load(save_folder + data_file)
    data = npzfile['data']
    fs = npzfile['fs']
    npzfile.close()      
    all_group_spws = []
    window_base = [-10, -5]
    for typ in range(len(spws)):
        # if there spws are divided into spontaneous and initiated, it will
        # analyze both groups separately, otherwise it will do everything as the whole
        type = spws[typ] #spws[len(np.unique(spws['group'])) == 1]
        spw_nos = []
        groups = []
        
        spw_nos = np.unique(type['spw_no']) #np.array(spw_nos)
        groups = np.zeros(len(spw_nos))
        all_spw = []
        
        #import pdb; pdb.set_trace() 
        group_idcs, = np.where(groups == 0)
        spw_used = spw_nos[group_idcs]
        
        # don't event bother if the group is has only 1 spw
        if np.size(spw_used,0) != 1:
            amplitudes = []
            answer = False
            subgroups = np.zeros(len(spw_used))
            spw_used_subgroup = spw_used.copy()
#            for idx, spw_no in enumerate(spw_used_subgroup):
#                spw_temp = type[type['spw_no'] == spw_no]
#                min_group = spw_temp[spw_temp['ipsp_start'] == min(spw_temp['ipsp_start'])]['group'][0]
#                ampls = spw_temp[spw_temp['group']==min_group]['amplitude']
#                amplitudes.append(ampls)
            
            #import pdb; pdb.set_trace() 
            sub = 0
            already_clustered = subgroups > 0
            while not answer:
                window = [-0.5, 7]
                print 'analysing subgroup: ' + str(sub)
                print subgroups
                group_name = str(0) + '.' +  str(sub) + ' ' + types[typ]
                #import pdb; pdb.set_trace() 
                ampls_used, new_starts_pts, traces, output = display_group_data(type, spw_used[subgroups == sub], data, fs, tit = group_name, window = window, window_base = window_base)
                #import pdb; pdb.set_trace() 
                #import pdb; pdb.set_trace()
                if output ==True:
                    # group is alright
                    already_clustered[subgroups == sub] = True
                    sub = sub + 1
                    
                elif output == False:
                    # group has to be further divided
                    
                    # variables set for pca calculations
                    electr_to_use = [1 ,4]
                    n_comps = 3
                    pcs = calculate_PCA(new_starts_pts, traces, data, fs, electr_to_use, n_comps, window = window, window_base = window_base)
                    
                    #characteristics = np.concatenate([pcs, ampls_used[:,[1, 4, 7]]], axis = 1)
                    _, actual_clusters = cluster.vq.kmeans2(pcs, 2,minit='points')
                    
                    names = ["E%d:P%d" % (e, p) for e in electr_to_use for p in range(n_comps) ]
                    compare_clusters(pcs.T, actual_clusters, names)
                    compare_spws(data, traces, electr_to_use, window, fs, new_starts_pts, actual_clusters, window_base = window_base)

                    #import pdb; pdb.set_trace()
                    assert len(actual_clusters)==len(traces),  "Cluster analysis failed-wrong feature dimensions"
                    actual_clusters[actual_clusters==1] =  np.max(subgroups) + 1
                    actual_clusters[actual_clusters==0] =  sub
                    
                    subgroups[subgroups == sub] = actual_clusters
                                    
                    plt.show()
                if sub > np.max(subgroups):
                    answer = True  
                #import pdb; pdb.set_trace()

                while len(spw_used[subgroups == sub]) == 1:
                    sub = sub+1
                    if sub > np.max(subgroups):
                        answer = True  
                        break;
        else:
            subgroups = np.zeros(len(spw_used))            
        # add the different spws to their groups and subgroups
        #import pdb; pdb.set_trace()  
        spw_no_temp = spw_used.astype('i4')
        subgroups_temp = (subgroups).astype('f8')
        
        new_spw_groups =  np.rec.fromarrays([spw_no_temp, subgroups_temp], names='spw_no, group')
        
        all_spw.append(new_spw_groups)
        #import pdb; pdb.set_trace()
        all_spw = np.concatenate(all_spw)  
        all_group_spws.append(all_spw)
       
    #all_group_spws = np.concatenate(all_group_spws)
    if double:
        group1 = all_group_spws[0]
        group2 = all_group_spws[1]
        np.savez(save_folder + save_file, group1 = group1, group2 = group2, names = types)  
        #np.savez(save_folder + save_file, group = all_spw, names = types)  
    else:
        np.savez(save_folder + save_file, group = all_spw, names = types)  

       

def plot_spw_ipsps_no_groups(save_folder, save_file, data_file, spw_data, ext):
    
    npzfile        = np.load(save_folder + spw_data)
    spws = npzfile['spw_ipsps']
    #initiated = npzfile['initiated']
    npzfile.close()           
    
    npzfile        = np.load(save_folder + data_file)
    data = npzfile['data']
    fs = npzfile['fs']
    npzfile.close()      
    
    types = ['spontaneous', 'initiated']
    all_ampls = []
    
    #type = spontaneous.copy()
    type = spws #spws[len(np.unique(spws['group'])) == 1]
    all_starts = []
    spw_nos = []
    #import pdb; pdb.set_trace() 
    max_group = -1
    groups = []
    # divide into groups looking at which electrode was the ipsp detected
    for spw_no in np.unique(type['spw_no']):
        spw_used = type[type['spw_no'] == spw_no]

        min_ipsps_group = spw_used[spw_used['ipsp_start'] == min(spw_used['ipsp_start'])]['group'][0]
        
        starting_electrodes = spw_used[spw_used['group'] == min_ipsps_group]['electrode']

        trace = spw_used['trace'][0]
        spw_start = spw_used['spw_start'][0] - 10
        spw_end = spw_start + 60
        spw_start_pts = ispw.ms2pts(spw_start, fs)
        spw_end_pts = ispw.ms2pts(spw_end,fs)

        # convert all to the different groups (using binary system coding)
        #if np.sum(2**(starting_electrodes + 1)) == 2:
        #    import pdb; pdb.set_trace()
            
        groups.append(np.sum(2**(starting_electrodes + 1)))
        spw_nos.append(spw_no)
        all_starts.append(starting_electrodes)
    
    # divide groups into subgroups by the amplitude of the first ipsp in different
    # electrodes (kmeans)
    
    
    groups = np.array(groups)
    spw_nos = np.array(spw_nos)
    
    
    all_spw = []
    
    for group in np.unique(groups):
        group_idcs, = np.where(groups == group)
        spw_used = spw_nos[group_idcs]
        if np.size(spw_used,0) != 1:
            amplitudes = []
            #subgroup = '.'
            answer = False
            subgroups = np.zeros(len(spw_used))
            spw_used_subgroup = spw_used.copy()
            for idx, spw_no in enumerate(spw_used_subgroup):
                spw_temp = spws[spws['spw_no'] == spw_no]
                min_group = spw_temp[spw_temp['ipsp_start'] == min(spw_temp['ipsp_start'])]['group'][0]
                ampls = spw_temp[spw_temp['group']==min_group]['amplitude']
                #if len(ampls) != 3:
                #    import pdb; pdb.set_trace()
                amplitudes.append(ampls)
            print amplitudes
            try:
                ampls = np.vstack(amplitudes)
            except:
                import pdb; pdb.set_trace()
            sub = 0
            already_clustered = subgroups > 0
            print np.size(spw_used,0)
            while not answer:
                print 'analysing subgroup: ' + str(sub)
                print subgroups
                group_name = str(group) + '.' +  str(sub)
                ampls_used, output = display_group_data(spws, spw_used[subgroups == sub], data, fs, tit = group_name)
                
                if output ==True:
                    # group is alright
                    already_clustered[subgroups == sub] = True
                    sub = sub + 1
                    
                elif output == False:
                    # group has to be further divided
                    ampls_used2 = ampls[subgroups == sub]
                    #try:
                    klastry = cluster.vq.kmeans2(ampls_used, 2,minit='points')
                    #except:
                    #import pdb; pdb.set_trace()
                    actual_clusters = klastry[1]
                    print actual_clusters
                    actual_clusters[actual_clusters==1] =  np.max(subgroups) + 1
                    actual_clusters[actual_clusters==0] =  sub
                    
                    print actual_clusters
                    subgroups[subgroups == sub] = actual_clusters
                
                if sub > np.max(subgroups):
                    answer = True  
                #import pdb; pdb.set_trace()

                while len(spw_used[subgroups == sub]) == 1:
                    sub = sub+1
                    if sub > np.max(subgroups):
                        answer = True  
                        break;
        else:
            subgroups = np.zeros(len(spw_used))            
        # add the different spws to their groups and subgroups
        #import pdb; pdb.set_trace()  
        spw_no_temp = spw_used
        subgroups_temp = group + subgroups/10. 
        
        new_spw_groups =  np.rec.fromarrays([spw_no_temp, subgroups_temp], names='spw_no, group')
        
        all_spw.append(new_spw_groups)
    all_spw = np.concatenate(all_spw)       
    np.savez(save_folder + save_file, group = all_spw)  



def display_group_data(spws, spw_to_use, data, fs, tit, window, window_base = [-10, -5]):    
   
    
    version = 1 # version 1 - keep original amplitudes and original data for plotting
        # version 2 - aligns everything to the peak of the first IPSP and calculates
        # its' amplitude
        
    align_on_max = False
    remove_baseline = True
    #if version == 1:
    #    window = [-10, 50]
    #else:
    #window = [-5, 5]
    window_pts0 = ispw.ms2pts(window[0], fs)   
    window_pts1 = ispw.ms2pts(window[1], fs) 
    add_it = 150
    
    #win0 = ispw.ms2pts(window[0], fs)
    window_plot = [-5, 30]
    win_pts0 = ispw.ms2pts(window_plot[0], fs)
    win_pts1 = ispw.ms2pts(window_plot[1], fs)
    baselin_length = 5 # ms
    baselin_len_pts = ispw.ms2pts(baselin_length, fs)
        
    #for group in np.unique(groups):
        #import pdb; pdb.set_trace() 
    #group_idcs, = np.where(groups == group)
    #spw_used = spw_nos[group_idcs]
    colors = define_colors(len(spw_to_use))
    ampls_used = []
    #fig = plt.figure(figsize=(20,10))
    fig = plt.figure(figsize = (15, 8))
    #fig.set_size_inches(18.5,10.5)
    ax = plt.subplot(111)
    plt.subplots_adjust(bottom=0.2)
    new_starts_pts = []
    traces = []
    
    # plot every spw given by spw_used
    for idx, spw_no in enumerate(spw_to_use):    
        #import pdb; pdb.set_trace() 
        spw_used = spws[spws['spw_no'] == spw_no]
        spw_used = np.sort(spw_used, order = 'electrode')
        trace = spw_used['trace'][0]

        if remove_baseline:
            s_p_temp1 = ispw.ms2pts(spw_used['spw_start'][0]+window_base[0], fs).astype('i4') 
            s_p_temp2 = ispw.ms2pts(spw_used['spw_start'][0]+window_base[1], fs).astype('i4')
            baseline = data[:, trace, s_p_temp1 : s_p_temp2]
            baseline = np.mean(baseline, axis = 1)
            #import pdb; pdb.set_trace()

        spw_start = spw_used['spw_start'][0]
        #spw_end = spw_used['spw_start'][0] + window[1]
        
        spw_start_pts = ispw.ms2pts(spw_start, fs).astype('i4')
        #spw_end_pts = ispw.ms2pts(spw_end, fs).astype('i4')
        
        
        data_used = data[:,trace,spw_start_pts + window_pts0:
                         spw_start_pts + window_pts1].copy()
        if remove_baseline:
            data_used[:,:] = data_used[:,:] - baseline[:,np.newaxis]
                
        electr_max = np.argmax(np.max(data_used, axis = 1))
#        if version == 1:
#            min_group = spw_used[spw_used['ipsp_no'] == min(spw_used['ipsp_no'])]['group'][0]
#            ampls = spw_used[spw_used['group'] == min_group]['amplitude']
#            ampls_used.append(ampls)
#            maxs = ampls
#        else:
        maxs = np.max(data_used, axis = 1)
        ampls_used.append(maxs)
        #import pdb; pdb.set_trace() 
        if version == 1:
            peak = spw_start_pts
        else:
            #import pdb; pdb.set_trace() 
            peak = spw_start_pts + np.argmax(data_used[electr_max, :]) +window_pts0
        new_starts_pts.append(peak)
        #new_starts_pts.append(ispw.ms2pts(spw_used['spw_start'][0], fs))
        traces.append(trace)
        #import pdb; pdb.set_trace() 
    
        data_to_plot = data[:,trace,peak + win_pts0:peak + win_pts1].copy()
        
        t = dat.get_timeline(data_to_plot[0,:], fs, 'ms')
        for electr in range(np.size(data,0)):
            if remove_baseline:
                #import pdb; pdb.set_trace()
                baseline = np.mean(data_to_plot[electr,:baselin_len_pts])
                #import pdb; pdb.set_trace()
                data_to_plot[electr,:] = data_to_plot[electr,:] - baseline

            if align_on_max:
                ax.plot(t, data_to_plot[electr, :] - maxs[electr] + add_it * electr, color = colors[idx], alpha = 0.2)

                
            else:
                ax.plot(t, data_to_plot[electr, :] + add_it * electr, color = colors[idx], alpha = 0.2)
            #data_used
        #import pdb; pdb.set_trace()    
   
    #class Update_Plot():
    # initiate all the variable used throughout the class
    #global answer
    def yes_button(event):
        #ax.Destroy()
        plt.close()
    #sys.exit(0)
        global answer
        answer = True
                      
    def no_button(event):
        plt.close()
        global answer 
        answer = False 
    plt.hlines(-add_it, t[window_pts0 - win_pts0], t[window_pts1 - win_pts0], 'r') #, width = 2.0)       
    plt.title('Group: ' + tit)
    axnext = plt.axes([0.75, 0.05, 0.15, 0.075])    
    bnext = Button(axnext, 'Update')
    bnext.on_clicked(no_button)
    axprev = plt.axes([0.6, 0.05, 0.15, 0.075])  
    bprev = Button(axprev, 'Keep')
    bprev.on_clicked(yes_button)
    
    plt.show()
    
    return np.vstack(ampls_used), new_starts_pts, traces, answer


def plot_spw_amplitude(save_folder, plot_folder, save_plots, data_file, spw_data, ext):
    npzfile        = np.load(save_folder + spw_data)
    spontaneous = npzfile['spontaneous']
    initiated = npzfile['initiated']
    npzfile.close()           
    
    npzfile        = np.load(save_folder + data_file)
    data = npzfile['data']
    fs = npzfile['fs']
    npzfile.close()           
    
    print 'plotting spw amplitude'
    types = ['spontaneous', 'initiated']
    all_ampls = []
    type = spontaneous[np.unique(spontaneous['group']) > 2]
    spw_numbs = []
    for spw_no in np.unique(type['spw_no']):
        spw_used = type[type['spw_no'] == spw_no]
        
        no_ipsps = len(np.unique(spw_used['group']))
        spw_start = spw_used['spw_start'][0] - 10
        spw_end = spw_start + 70#max(spw_used['ipsp_start'])
        if spw_start == spw_end:
            spw_end = spw_start + 10
        spw_start_pts = ispw.ms2pts(spw_start, fs)
        spw_end_pts = ispw.ms2pts(spw_end, fs)
        trace = spw_used['trace'][0]
        
        data_used = data[:, trace, spw_start_pts:spw_end_pts]
        ampls = [max(data_used[electr, :]) for electr in range(len(data_used))]
        ampls = np.argsort(ampls)
        all_ampls.append(ampls)
        spw_numbs.append(spw_no)
        electr_no =len(data_used)
        np.zeros([electr_no, electr_no])
        
        all_shifts = []
        all_similarities = []
        for electr1 in range(electr_no - 1):
            electr2 = electr1+1
            #import pdb; pdb.set_trace()
            #print electr1
            #print electr2
            #corr = signal.correlate(data_used[electr1, :],data_used[electr2, :])
            dat1 =  data_used[electr1, :]# - data_used[electr1, :].mean()
            #dat1/=dat1.std()
            dat2 =  data_used[electr1, :] #- data_used[electr2, :].mean()
            #dat2/=dat2.std()
            
            corr = find_max_corr(dat1,dat2)
            #plt.plot(corr)
            #plt.show()
            
            nsamples = dat1.size
            shift = corr.argmax()
            all_similarities.append(np.max(corr))
            #dt = np.arange(nsamples) #np.arange(1-nsamples, nsamples) - len(dat1)
            #recovered_time_shift = dt[corr.argmax()]

            recovered_time_shift = shift - nsamples/2.
            all_shifts.append(recovered_time_shift)
            plt.plot(dat1 + 50 * electr1)
            #import pdb; pdb.set_trace()
            
        print
        print all_shifts
        print all_similarities
        t = dat.get_timeline(data_used[0,:], fs, 'ms')
        #import pdb; pdb.set_trace()
        #for electr in range(len(data_used)):
        #    plt.plot(data_used[electr, :] + 50 * electr)
        plt.show()
                        

## Load datasets, taking mean of 100 values in each table row
#A = numpy.loadtxt("vb-sync-XReport.txt")[:,1:].mean(axis=1)
#B = numpy.loadtxt("vb-sync-YReport.txt")[:,1:].mean(axis=1)
#
#nsamples = A.size
#
## Put in an artificial time shift between the two datasets
#time_shift = 20
#A = numpy.roll(A, time_shift)
#
## Find cross-correlation
#xcorr = correlate(A, B)
#
## delta time array to match xcorr
#dt = numpy.arange(1-nsamples, nsamples)
#
#recovered_time_shift = dt[xcorr.argmax()]
#
#print "Added time shift: %d" % (time_shift)
#print "Recovered time shift: %d" % (recovered_time_shift)
#
## SAMPLE OUTPUT:
## Added time shift: 20
## Recovered time shift: 20


            
    #import pdb; pdb.set_trace()
    corr_matrix = np.array([[np.sum(x==y) for x in all_ampls] for y in all_ampls])
    #all_ampls_redistributed = np.vstack(all_ampls)
    #corr_matrix = np.dot(all_ampls_redistributed, all_ampls_redistributed.T)
    #max_value = np.diag(corr_matrix)[0]
    idx_i, idx_j = np.where(np.triu(corr_matrix,1)> 6)
    for value in np.unique(np.append(idx_i,idx_j)):
        the_same = np.concatenate([idx_j[idx_i==value], idx_i[idx_j==value]])
        plt.figure()
        for spw in the_same:
            spw_no = spw_numbs[spw]
            spw_used =  type[type['spw_no'] == spw_no]
            trace = spw_used['trace'][0]
            spw_start = spw_used['spw_start'][0] - 10
            spw_end = spw_start + 70
            spw_start_pts = ispw.ms2pts(spw_start, fs)
            pw_end_pts = ispw.ms2pts(spw_end, fs)
            data_used = data[:, trace, spw_start_pts:pw_end_pts]
            t = dat.get_timeline(data_used[0,:], fs, 'ms')
            #import pdb; pdb.set_trace()
            for electr in range(len(data_used)):
                plt.plot(t, data_used[electr, :] + 150 * electr)
        plt.show()
        #amls = [max(data_used[])]
            
def plot_imshow_origin(dane, window, save_name, title, electrodes, vrange = [0, 0]):
    # plots and saves imshow form the given parameters
    fig = plt.figure()  
    if vrange[0] == vrange[1] == 0: 
        plt.imshow(dane, aspect = 'auto', interpolation='nearest', origin='lower', extent=[window[0],window[1],0.5,len(electrodes)+0.5]) #, vmin=0, vmax=0.3) #, interpolation='bilinear', aspect = 'auto') #interpolation='nearest', aspect='auto')
    else: 
        plt.imshow(dane, aspect = 'auto', interpolation='nearest', origin='lower', extent=[window[0],window[1],0.5,len(electrodes)+0.5], vmin=vrange[0], vmax=vrange[1]) #, interpolation='bilinear', aspect = 'auto') #interpolation='nearest', aspect='auto')

    plt.colorbar()

    plt.title(title)
    plt.ylabel('electrode number')
    plt.xlabel('time from beginning of the detected start of SPW (ms)')
        
    fig.savefig(save_name, dpi=600)     
        
    

def plot_spike(save_folder, plot_folder, save_plots, save_file, save_name_max_electr, spike_data = 'spike.npz', spw_data = 'spw.npz', ext = '.pdf', win = [-20, 20]):
    """ counts how many spikes are in each electrode during all the spws and do the image shows """
    #import pdb; pdb.set_trace()
    npzfile        = np.load(save_folder + spike_data)
    spikes = npzfile['spike_idx']
    npzfile.close()    
    
    npzfile        = np.load(save_folder + spw_data)
    spontaneous = npzfile['spontaneous']
    initiated = npzfile['initiated']
    npzfile.close()           
    types = ['spontaneous', 'initiated']
    
    n_bins = 6
    all_p_dist = []
    #window = [-.5, 1.]
    window = [-1, 5.]
    all_dists_hist = []
    for_imshow = []
    bins = np.linspace(window[0], window[1], n_bins + 1)
    considered_bin = np.where(bins <= 0)[0][-1]
    
    save_fold = save_folder + plot_folder
    fold_mng.create_folder(save_fold)
    save_base = save_fold + save_plots
    
    all_spikes = []
    # calculate all the spikes which are within given window
    for idx_type, typ in enumerate([spontaneous, initiated]):
        
        spikes_list = [[] for x in range(0,len(np.unique(typ['electrode'])))]
        
        # go through every spw in this group
        for spw_start in np.unique(typ['spw_start']):
            if spw_start > -window[0]:
                trace = typ[typ['spw_start'] == spw_start]['trace'][0]
                spikes_used = spikes[spikes['trace'] == trace]
                spikes_spw = spikes_used[(spikes_used['time'] > spw_start + window[0]) & 
                                         (spikes_used['time'] < spw_start + window[1])]
                
                # look for spikes in each electrode
                for electr in np.unique(spikes_spw['electrode']):
                    spikes_times = (spikes_spw[spikes_spw['electrode'] == electr]['time'] - spw_start).tolist()
                    spikes_list[electr] = spikes_list[electr] + spikes_times
        
        for_histogram = []
        dist_electrode = []
        numb_spikes = 0
        
        
        # check how many spikes are in every bin
        for electr in np.unique(typ['electrode']):
            n_all, _ = np.histogram(spikes_list[electr], bins) #, normed = True)
            numb_spikes = (numb_spikes + len(spikes_list[electr])) * 1.0
            dist_electrode.append(n_all)
            for_histogram.append(n_all[considered_bin])
        
        # dist_electrode - just simple number of spikes in each electrode
        dist_electrode = [dist_electrode[i]/numb_spikes for i in range(len(dist_electrode))]
        
        # numb_spikes - number of spikes all together
        numb_spikes = numb_spikes * 1.0
        all_spikes.append(numb_spikes)
        # for_histogram - only those dist_electrode from considered_bin
        for_histogram = np.array(for_histogram)
        
        # same as for_histogram, but will not be normalized - used for chi-sqare later
        for_p_distribution = for_histogram.copy()
        
        # divide by all number of spikes in this histogram
        for_histogram = for_histogram / (for_histogram.sum() * 1.0)
        
        for_imshow.append(dist_electrode)
        all_dists_hist.append(for_histogram)
        all_p_dist.append(for_p_distribution)
    
    # plot imshow - not normalized
    for idx_type, typ in enumerate([spontaneous, initiated]):  
        title = types[idx_type] + ', found spws: ' + str(len(np.unique(typ['spw_no']))) + ', found spikes: ' + str(all_spikes[idx_type])
        save_name =  save_base + types[idx_type] + ext
        electrs = np.unique(typ['electrode'])
        plot_imshow_origin(for_imshow[idx_type], window,save_name, title, electrs) 
    
    # calculate chi_spare between the two distributions
    import scipy.stats.mstats as mst
    fig = plt.figure()
    width = 0.5
    left = np.arange(0.5, len(for_histogram) + 0.5)
    find_zeros = np.sum(all_p_dist,0)
    x_squared, p_value = mst.chisquare(np.array(all_p_dist[0][find_zeros != 0]), np.array(all_p_dist[1][find_zeros != 0]))

    #import pdb; pdb.set_trace()
    # plot bar plot for the distributions
    rects1 = plt.bar(left , all_dists_hist[0], width, color='r')
    rects2 = plt.bar(left + width, all_dists_hist[1], width, color = 'b')
    plt.title('p: ' + str(p_value))
    plt.legend( (rects1[0], rects2[0]), ('Spontaneous', 'Induced') )
    plt.title(str(p_value))
    plt.xlim([left[0], 1 + left[-1]])     
    fig.savefig(save_fold + save_plots + types[idx_type] + '_hist' + ext, dpi=600)     
    
    
    # normalize by sum in each pixel
    #import pdb; pdb.set_trace()
    sum_in_each = (np.sum(for_imshow,0)) * 1.0
    for idx_type, typ in enumerate([spontaneous, initiated]):  
        title = types[idx_type] + ', norm_pixel, found spws: ' + str(len(np.unique(typ['spw_no']))) + ', found spikes: ' + str(all_spikes[idx_type])
        save_name =  save_base + types[idx_type]+'norm_pixel' + ext
        electrs = np.unique(typ['electrode'])
        
        dane = np.array(for_imshow[idx_type])/sum_in_each
        dane[sum_in_each == 0.] = 0.
        dane[dane == np.nan] = 0
        plot_imshow_origin(dane, window,save_name, title, electrs) 
    
    # plot bar for normalized data        
    plt.figure()
    column_used = np.where(bins <= 0)[0][-1]
    rects1 = plt.bar(left , all_dists_hist[0]/sum_in_each[:,column_used] , width, color='r')
    rects2 = plt.bar(left + width, all_dists_hist[1]/sum_in_each[:,column_used], width, color = 'b')
    plt.legend( (rects1[0], rects2[0]), ('Spontaneous', 'Induced'))
    plt.title(str(p_value))
    plt.xlim([left[0], 1 + left[-1]])
    fig.savefig(save_fold + save_plots + types[idx_type] + '_hist_norm_pixel' + ext, dpi=600)   
    
    #import pdb; pdb.set_trace()
    # normalize by sum in each electrode
    sum_in_electr = np.sum(for_imshow,2)
    sum_between = np.sum(sum_in_electr, 0)
    for idx_type, typ in enumerate([spontaneous, initiated]):  
        title = types[idx_type] + ', norm_electr, found spws: ' + str(len(np.unique(typ['spw_no']))) + ', found spikes: ' + str(all_spikes[idx_type])
        save_name =  save_base + types[idx_type]+'norm_electr' + ext
        electrs = np.unique(typ['electrode'])
        
        dane = []
        for electr in range(len(for_imshow[idx_type])):
            if sum_between[electr] != 0:
                dane.append(np.array(for_imshow[idx_type][electr])/sum_between[electr])
            else:
                dane.append(np.array(for_imshow[idx_type][electr] - for_imshow[idx_type][electr]))
        
        plot_imshow_origin(np.array(dane), window,save_name, title, electrs) 

    # plot bar by substracting not normalized data  
    fig = plt.figure()
    column_used = np.where(bins <= 0)[0][-1]
    difference_between = all_dists_hist[0]-all_dists_hist[1]
    spont_part = difference_between.copy()
    spont_part[spont_part < 0] = 0
    induc_part = difference_between.copy()
    induc_part[induc_part > 0] = 0    
    
    rects1 = plt.bar(left + width/2 , spont_part, width, color='r')
    rects2 = plt.bar(left + width/2, induc_part, width, color = 'b')
    plt.legend( (rects1[0], rects2[0]), ('Spontaneous', 'Induced'))
    plt.title(str(p_value))
    plt.xlim([left[0], 1 + left[-1]])
    fig.savefig(save_fold + save_plots + types[idx_type] + '_hist_difference' + ext, dpi=600)          
    #import pdb; pdb.set_trace()
    
    electr_spont_diff = np.argmax(spont_part) + 1
    electr_init_diff = np.argmin(induc_part) + 1
    electr_spont =  np.argmax(all_dists_hist[0]) + 1
    electr_init = np.argmax(all_dists_hist[1]) + 1
    
    #import pdb; pdb.set_trace()
    np.savez(save_folder + save_name_max_electr, spont_diff = electr_spont_diff, init_diff = electr_init_diff, spont = electr_spont, init = electr_init)
    # plot imshow - difference between the two
    #for idx_type, typ in enumerate([spontaneous, initiated]):  
    title = 'found spws: ' + str(len(np.unique(typ['spw_no']))) + ', found spikes: ' + str(numb_spikes)
    save_name =  save_base + types[idx_type] + '_difference_' + ext
    electrs = np.unique(typ['electrode'])
    #import pdb; pdb.set_trace()
    
    difference = [for_imshow[0][electr] - for_imshow[1][electr] for electr in range(len(for_imshow[0]))]
    max_value = max(abs(difference_between))
    plot_imshow_origin(difference, window,save_name, title, electrs, vrange = [-max_value, max_value])   
    plt.close('all')
    #import pdb; pdb.set_trace()
    #plt.show()

def plot_spikes4spw(save_folder, plot_folder, save_plots = 'saved', data_file = 'data.npz', spike_data = 'spikes.npz', spw_data = 'spw.npz', spikes_filter = [], ext = '.pdf', win = [-20, 20], filt = 600.0):
    """ plots every spw separately (in all electrodes)"""

    npzfile        = np.load(save_folder + data_file)
    data = npzfile['data']
    fs = npzfile['fs']
    npzfile.close()     
    colors = define_colors(no_colors = np.size(data,0))
    
    npzfile        = np.load(save_folder + spw_data)
    spontaneous = npzfile['spontaneous']
    initiated = npzfile['initiated']
    npzfile.close()         
    
    npzfile        = np.load(save_folder + spike_data)
    #import pdb; pdb.set_trace()
    spikes = npzfile['spike_idx']
    npzfile.close()     
    #import pdb; pdb.set_trace()
    types = ['spontaneous', 'initiated']
    #add_before_after = 20 # ms
    before_pts = ispw.ms2pts(win[0], fs)
    after_pts = ispw.ms2pts(win[1], fs)
    
    use_filter = [200, 500]
    plot_only_last = False
    
    print 'plotting the data and their spikes'
    add_it = 500
    for idx, typ in enumerate([spontaneous, initiated]):
        #data[:, ]
        #t = 
        spw_noms = np.unique(typ['spw_no'])
        for spw_used in spw_noms:
            spw = typ[typ['spw_no'] == spw_used][0]
            trace = spw['trace']
            spw_start = ispw.ms2pts(spw['spw_start'], fs).astype('i4')
            spw_end = ispw.ms2pts(spw['spw_end'], fs).astype('i4')
            data_used = data[:, trace, spw_start+before_pts:spw_end + after_pts]
            
            t = dat.get_timeline(data_used[0, :], fs, 'ms') + win[0]
            #import pdb; pdb.set_trace()
            spikes_used = spikes[spikes['trace'] == trace]
            spikes_used['time'] = ispw.ms2pts(spikes_used['time'], fs)
            spikes_used = spikes_used[(spikes_used['time'] > spw_start+before_pts) & (spikes_used['time'] <spw_end + after_pts)] 
            
            fig = plt.figure()
            
            
            
            mult_factor = 10
            #import pdb; pdb.set_trace()
            #
            
            #import pdb; pdb.set_trace()
            if plot_only_last and spikes_filter != []:
                mult_factor = 10
                npzfile = np.load(save_folder + spikes_filter + str(filt) + '_0_' + str(trace) + '.npz')
                #import pdb; pdb.set_trace()
                filt_dat = npzfile['data']
                npzfile.close() 
                data_fast = filt_dat[spw_start+before_pts:spw_end + after_pts]
                plt.plot(t, data_fast * mult_factor + add_it * -1.5,  'k')
                #filt_dat = filt_dat[spw_start+before_pts:spw_end + after_pts]
                #plt.plot(t, data_fast[0, :]*mult_factor + add_it * -1.5,  'k')
                #import pdb; pdb.set_trace()
            #else:
            #    filt_dat, fs = ispw.load_create(save_folder + spikes_filter, str(use_filter) + '_0_' + str(trace) + '.npz', use_filter, fs, data[:, trace, :], 100)
            #    data_fast = filt_dat[:,spw_start+before_pts:spw_end + after_pts]
            #    import pdb; pdb.set_trace()
            for electr in range(np.size(data_used, 0)):
                filt_dat, fs = ispw.load_create(save_folder + spikes_filter, str(use_filter) + '_' + str(electr) + '_' + str(trace) + '.npz', use_filter, fs, data[electr, trace, :], 100)
                data_fast = filt_dat[spw_start+before_pts:spw_end + after_pts]
                
                if not plot_only_last:
                    #filt_dat = filt_dat[spw_start+before_pts:spw_end + after_pts]
                    plt.plot(t, data_fast*mult_factor + add_it* electr + (0.5 * max(data_used[electr, :])) ,  'k', alpha = 0.3)
                    
                plt.plot(t, data_used[electr, :] + add_it * electr,  'k') #color = colors[electr],
                spikes_electr = spikes_used[spikes_used['electrode'] == electr]['time'].astype('i4')- spw_start
                #plt.plot(t[spikes_electr], data_used[electr, spikes_electr], )
                #spikes_plot = spikes_used[spikes_used['electrode'] == electr]['time']
                #spikes_idx = (-ispw.ms2pts(spw['spw_start'],fs).astype('i4') +
                #              ispw.ms2pts(spikes_plot, fs).astype('i4')-
                #              before_pts).astype(int)
                                #(spw['spw_start'] - spikes_used)#spw_start - ispw.ms2pts(spikes_used[spikes_used['electrode'] == electr]['time'],fs) - ispw.ms2pts(win[0], fs)
                

        
    
                
                #import pdb; pdb.set_trace()    
                #spikes_idx = spikes_idx[spikes_idx < np.size(data_used,1)]
                #plt.plot(t[spikes_idx], data_used[electr, spikes_idx] + add_it * electr - 5, 'k*', linewidth = 6.)
                plt.scatter(t[spikes_electr], data_used[electr, spikes_electr] + add_it * electr - 20, s = 25., marker = '*')
                
            
            plt.xlim([t[0], t[-1]])
            #import pdb; pdb.set_trace()
            if plot_only_last and spikes_filter != []:
                plt.ylim([min(data_fast[0, :] * mult_factor + add_it * -1.5)-100, data_used.max() + electr *add_it + 100])
                plt.ylabel('Voltage (' + micro() + 'V); filter multiplied by: '  + str(mult_factor))
            else:
                plt.ylim([data_used.min()-100, data_used.max() + electr *add_it + 100])
                plt.ylabel('Voltage (' + micro() + 'V)')
            plt.title(types[idx] + ', spw no: ' + str(spw_used))
            plt.xlabel('time(ms) from the detected beginning of the SPW')
            save_fold = save_folder + plot_folder
            fold_mng.create_folder(save_fold)
            fig.savefig(save_fold + save_plots + str(spw_used) + types[idx] + ext,dpi=600)     
            fig.savefig(save_fold + save_plots + str(spw_used) + types[idx] + '.eps',dpi=600) 
            #plt.+() 
            plt.close()   

        
def plot_alignedSPW(save_folder, plot_folder, save_plots, data_file, intra_data_file, induc_spont, intra_spikes, ext):
    """ it divides SPWs to two groups - close and far from the spike and alignes them, and plots together"""
    
    from scipy.stats import nanmean
    win = [-5, 80] #ms
    
    npzfile         = np.load(save_folder + induc_spont)
    spontaneous      = npzfile['spontaneous'] # spikes_all
    initiated      = npzfile['initiated'] # spikes_all
    npzfile.close()    

    npzfile        = np.load(save_folder + intra_data_file)
    data_intra = npzfile['data']
    fs = npzfile['fs']
    npzfile.close() 
     
    
    npzfile        = np.load(save_folder + data_file)
    data = npzfile['data']
    fs = npzfile['fs']
    npzfile.close() 
    
    npzfile        = np.load(save_folder + intra_spikes)
    intra_spikes = npzfile['spikes_first']
    npzfile.close() 
      
    win_base = [-10, -5]
    win_base_pts = [ispw.ms2pts(win_base[0], fs),ispw.ms2pts(win_base[1], fs)]  
    
    
    before_pts = ispw.ms2pts(win[0], fs)
    after_pts = ispw.ms2pts(win[1], fs)
    
    #initiated_no = distances[distances['distance'] <= min_dist]['spw_no']
    #spont_no = distances[distances['distance'] > min_dist]['spw_no']
    
    #import pdb; pdb.set_trace() 
    #init_set = np.in1d(ipsps['spw_no'], initiated_no, assume_unique= False)
    #initiated = ipsps[init_set]
    #spont_set = np.in1d(ipsps['spw_no'], spont_no, assume_unique=False)
    #spontaneous = ipsps[spont_set]
    #import pdb; pdb.set_trace() 
    spws_all = [initiated, spontaneous]
    all_data_traces = []
    all_in_spikes = []
    
    #import pdb; pdb.set_trace() 
    add_nans = np.ones([np.size(data, 0), np.size(data,1), after_pts - before_pts]) * np.nan  
    all_data = np.zeros([np.size(data, 0), np.size(data,1), np.size(data,2) + np.size(add_nans, 2)])
    #all_data[:, :, 0:np.size(add_nans, 2)] = add_nans
    all_data[:, :, 0:-np.size(add_nans, 2)] = data
    all_data[:, :, -np.size(add_nans, 2):] = add_nans
    temp_len = np.size(data,2)
    data = all_data
    del all_data
    
    save_fold = save_folder + plot_folder
    fold_mng.create_folder(save_fold)
    
    all_data_intra = np.zeros([np.size(data_intra, 0), np.size(data_intra,1), np.size(data_intra,2) + np.size(add_nans, 2)])
    all_data_intra[:, :, 0:-np.size(add_nans, 2)] = data_intra
    all_data_intra[:, :, -np.size(add_nans, 2):] = add_nans[0,:,:]
    data_intra = all_data_intra
    del all_data_intra
    for spws in spws_all:
        # take out those spws which are too early or too late (don't fit in the data size)
        #used = ispw.ms2pts(spws['spw_start'], fs).astype(int)
        #spws_used = spws[(used > -before_pts) & (used + after_pts < np.size(data,2))]
        
        spws_used = spws
        spw = np.unique(spws_used['spw_no'])
        spw_traces = np.zeros([np.size(data,0) + 1, len(spw), after_pts - before_pts])
        in_spikes = []

        for spw_idx, spw_n in enumerate(spw):
            spw_start = spws_used[spws_used['spw_no'] == spw_n]['spw_start'][0]
            spw_start_pts = ispw.ms2pts(spw_start, fs).astype(int)
            trace = spws_used[spws_used['spw_no'] == spw_n]['trace'][0]
            
            base = data[:, trace, spw_start_pts + win_base_pts[0]: spw_start_pts + win_base_pts[1]]
            base = np.mean(base, axis = 1)
            
            #data_spw = np.transpose(data_spw) - base
            #data_spw = np.transpose(data_spw)

            
            data_temp = data[:, trace, spw_start_pts + before_pts: spw_start_pts + after_pts]
            
            data_temp = np.transpose(data_temp) - base
            data_temp = np.transpose(data_temp)
            try:
                spw_traces[1:, spw_idx, :] = data_temp
            except:
                import pdb; pdb.set_trace() 
            data_temp_intra = data_intra[:, trace, spw_start_pts + before_pts: spw_start_pts + after_pts]
            try:
                spw_traces[0, spw_idx, :] = data_temp_intra
            except:
                import pdb; pdb.set_trace() 
            
            spikes_detected = intra_spikes[(intra_spikes['time'] < spw_start + win[1]) & (intra_spikes['time'] > spw_start + win[0])]['time']
            
            spikes_detected = (ispw.ms2pts(spikes_detected, fs) - spw_start_pts - before_pts).astype(int)
            in_spikes.append(spikes_detected)
        all_in_spikes.append(in_spikes)    
        all_data_traces.append(spw_traces)
    
            
    titles = ['Induced', 'Spontaneous']
    #import pdb; pdb.set_trace() 
    add_it = 400
    t = dat.get_timeline(data_temp[0], fs, 'ms') + win[0]
    #in_spikes = np.concatenate(in_spikes).astype(int)
    for idx, data_spw in enumerate(all_data_traces): 
        #print 'tak'
        fig = plt.figure()
        in_spikes = all_in_spikes[idx]
        for electr in range(len(data_spw)):
            
            data_used = data_spw[electr,:]
            for s in range(len(data_used)):

                plt.plot(t, data_used[s] + electr * add_it, 'b', alpha=0.2)
                
                if electr == 0:
                    print s
                    #import pdb; pdb.set_trace()
                    plt.plot(t[in_spikes[s]], data_used[s, in_spikes[s]] + electr * add_it, 'r.')
                
                
            if electr != 0:

                plt.plot(t, nanmean(data_used, 0) + electr * (add_it), 'r')
                
            if electr == 1:
                plt.xlim([t[0], t[-1]])
                fig_fname = save_fold + save_plots + titles[idx] + 'electr_0' + ext
                fig.savefig(fig_fname,dpi=600)   
        #import pdb; pdb.set_trace()
        
        plt.title(titles[idx] + ', spws: ' + str(len(data_used)))     

        fig_fname = save_fold + save_plots + titles[idx] + str(win[0]) + '_' + str(win[1]) + ext
        logging.info("saving figure %s" % fig_fname)
        fig.savefig(fig_fname,dpi=600)    
        #logging.info("saving figure %s" % fig_fname)       
        plt.close() 
        
    
#---------------------------old -----------------------------
def define_colors(no_colors = 8, type = 'color'):
    import colorsys
    if type == 'color':
        RGB_tuples = []
        #import pdb; pdb.set_trace()
        for i in range(no_colors):
            i = i * 1.0
            RGB_tuples.append(colorsys.hsv_to_rgb(i/no_colors, 0.5, 0.5))
        #RGB_tuples = [(x*1.0/no_colors, x*0.5/no_colors,x*0.5/no_colors) for x in range(no_colors)]
        RGB_tuples = np.array(RGB_tuples)
        #RGB_tuples = [(1, 0, 0), (1, 0.5, 0), (1, 1, 0), (0.5, 1, 0), (0, 1, 0.5), (0, 1, 1), (0, 0.5, 1), (0, 0, 1)]
        #RGB_tuples = [(0.5, 1, 0), (0, 1, 0.5), (0, 1, 1), (0, 0.5, 1), (0, 0, 1), (0.5, 0, 1), (1, 0, 1), (1, 0, 0.5)]
    
    elif type == 'grey':
        #import pdb; pdb.set_trace()
        use_vals = np.linspace(0, 1.0, no_colors + 1)
        RGB_tuples = []
        #import pdb; pdb.set_trace()
        for i in use_vals[1:]:
            #i = i*1.0 + 1
            RGB_tuples.append((i, i, i))
    
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
    

