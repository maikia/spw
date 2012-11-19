import numpy as np

#datanpz = np.load('data_bas.npz')
#data = datanpz['data']
#fs = datanpz['fs']
def rec2array(save_folder, file_name_old, file_name_new):
    datanpz = np.load(save_folder+ file_name_old)
    data = datanpz['data']
    fs = datanpz['fs']
    electrodes = np.unique(data['electrode'])
    traces = np.unique(data['trace'])
    
    n_el = electrodes.max()+1
    n_traces = traces.max()+1
    n_pts = np.sum((data['electrode']==1) & (data['trace']==0))
    data3d = np.zeros((n_el, n_traces, n_pts))
    
    for i in electrodes:
        for j in traces:
            data3d[i,j,:] = data['time'][(data['electrode']==i) & (data['trace']==j)]
    
    np.savez(save_folder + file_name_new, data=data3d, fs=fs)