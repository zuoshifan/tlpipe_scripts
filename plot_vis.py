import os
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


in_dir1 = './output_ps_subtract/r2t/'
in_dir2 = './output_ps_subtract/data/ps_cal/'
out_dir = './plot/vis/'
if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

with h5py.File(in_dir1 + 'file.hdf5', 'r') as f:
    vis1 = f['vis'][:]
    vis_mask1 = f['vis_mask'][:]
    bl_order = f['blorder'][:]
    lhour = f['local_hour'][:]

with h5py.File(in_dir2 + 'file.hdf5', 'r') as f:
    vis2 = f['vis'][:]
    vis_mask2 = f['vis_mask'][:]

nt, nf, npol, nbl = vis1.shape

bi = 1000
v1 = np.where(vis_mask1[:, 0, 0, bi], complex(np.nan, np.nan), vis1[:, 0, 0, bi])
v1 = np.where((lhour>8.0)&(lhour<19.5), complex(np.nan, np.nan), v1)

v2 = np.where(vis_mask2[:, 0, 0, bi], complex(np.nan, np.nan), vis2[:, 0, 0, bi])
v2 = np.where((lhour>8.0)&(lhour<19.5), complex(np.nan, np.nan), v2)

plt.figure()
# plt.plot(lhour, np.abs(v1))
plt.plot(lhour, np.abs(v2))
plt.savefig('vis_sub_%d_%d.png' % tuple(bl_order[bi]))
# plt.savefig('vis_%d_%d_before.png' % tuple(bl_order[bi]))
plt.close()
