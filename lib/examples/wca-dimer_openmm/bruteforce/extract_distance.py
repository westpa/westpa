import netCDF4 as netcdf
import numpy as np

f = netcdf.Dataset('data/md-solvent-langevin.nc', 'r')

dis = f.variables['distance']

chunksize = 50000
data = []
maxstep = dis.shape[0]
i = list(range(0, maxstep + chunksize, chunksize))

for k in range(len(i)-1):
    print(i[k], i[k+1])
    data.append(dis[i[k]:i[k+1]])

d = np.hstack(data)
np.save('data/md-solvent-langevin-distance.npy', d)


