import numpy as np

datanpz = np.load('data_bas.npz')
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

np.savez('data_bas_3d.npz', data=data3d, fs=fs)