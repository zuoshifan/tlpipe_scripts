import os
import itertools
import argparse
from collections import defaultdict
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description='Check visibilities of redundant baselines.')
parser.add_argument('--vis_file', type=str, default='./vis.hdf5', help='Visibility files')
parser.add_argument('--gain_file', type=str, default='./gain.hdf5', help='Gain files')
parser.add_argument('--cylinder', type=int, choices=[1, 2, 3], default=1, help='Which cylinder to plot.')
parser.add_argument('-o', '--out_dir', type=str, default='./', help='Output directory where the plotted figures should be put in.')
parser.add_argument('-t', '--time', type=int, default=0,  help='Time indices to be plot.')
parser.add_argument('-f', '--freq', type=int, default=0,  help='Frequency indices to be plot.')

args = parser.parse_args()


out_dir = args.out_dir if args.out_dir.endswith('/') else ('%s/' % args.out_dir)
if not os.path.isdir(out_dir):
    os.makedirs(out_dir)



vis_file = args.vis_file
gain_file = args.gain_file

# read in vis
with h5py.File(vis_file, 'r') as f:
    sky_vis = f['sky_vis'][:]
    src_vis = f['src_vis'][:]
    nt, nf, npol, nfeed, _ = sky_vis.shape
    sky_vis_xx = sky_vis[args.time, args.freq, 0]
    sky_vis_yy = sky_vis[args.time, args.freq, 1]
    src_vis_xx = src_vis[args.time, args.freq, 0]
    src_vis_yy = src_vis[args.time, args.freq, 1]

# read in gain
with h5py.File(gain_file, 'r') as f:
    gain_xx = f['gain'][args.freq, 0]
    gain_yy = f['gain'][args.freq, 1]

feed_ranges = [(0, 31), (31, 31+32), (31+32, 31+32+33)]
feed_range = feed_ranges[args.cylinder - 1]

colors = ['b', 'g', 'r', 'm', 'y', 'c', 'k']
markers = ['o', '^', 's', '*', 'd']
cm = list(itertools.product(colors, markers))
# print len(list(cm))

# check visibilities observed by a single cylinder
# for xx
for vis_xx, vis_xx_str in zip([sky_vis_xx, src_vis_xx], ['sky_vis_xx', 'src_vis_xx']):
    vis_xx1 = defaultdict(list)
    vis_xx1_cal = defaultdict(list)
    for d1 in range(*feed_range):
        for d2 in range(d1, feed_range[1]):
            vis_xx1[abs(d2 - d1)].append(vis_xx[d1, d2])
            g_xx = gain_xx[d1] * np.conj(gain_xx[d2])
            if np.isfinite(g_xx):
                vis_xx1_cal[abs(d2-d1)].append(vis_xx[d1, d2] / g_xx)

    # plt vis
    plt.figure()
    for i, (k, v) in enumerate(vis_xx1.iteritems()):
        vis_real = np.array(v).real
        vis_imag = np.array(v).imag
        c, m = cm[i][0], cm[i][1]
        plt.scatter(vis_real, vis_imag, c=c, marker=m)
    plt.xlim(-1.0e9, 1.0e9)
    plt.ylim(-1.0e9, 1.0e9)
    plt.xlabel('Re(Visibility)')
    plt.ylabel('Im(Visibility)')
    plt.savefig(out_dir+'%s_cyl%d.png' % (vis_xx_str, args.cylinder))
    plt.close()

    plt.figure()
    for i, (k, v) in enumerate(vis_xx1_cal.iteritems()):
        vis_real = np.array(v).real
        vis_imag = np.array(v).imag
        c, m = cm[i][0], cm[i][1]
        plt.scatter(vis_real, vis_imag, c=c, marker=m)
    plt.xlim(-0.8e4, 0.8e4)
    plt.ylim(-0.8e4, 0.8e4)
    plt.xlabel('Re(Visibility)')
    plt.ylabel('Im(Visibility)')
    plt.savefig(out_dir+'%s_cyl%d_cal.png' % (vis_xx_str, args.cylinder))
    plt.close()

# for yy
for vis_yy, vis_yy_str in zip([sky_vis_yy, src_vis_yy], ['sky_vis_yy', 'src_vis_yy']):
    vis_yy1 = defaultdict(list)
    vis_yy1_cal = defaultdict(list)
    for d1 in range(*feed_range):
        for d2 in range(d1, feed_range[1]):
            vis_yy1[abs(d2 - d1)].append(vis_yy[d1, d2])
            g_yy = gain_yy[d1] * np.conj(gain_yy[d2])
            if np.isfinite(g_yy):
                vis_yy1_cal[abs(d2-d1)].append(vis_yy[d1, d2] / g_yy)

    # plt vis
    plt.figure()
    for i, (k, v) in enumerate(vis_yy1.iteritems()):
        vis_real = np.array(v).real
        vis_imag = np.array(v).imag
        c, m = cm[i][0], cm[i][1]
        plt.scatter(vis_real, vis_imag, c=c, marker=m)
    plt.xlim(-1.0e9, 1.0e9)
    plt.ylim(-1.0e9, 1.0e9)
    plt.xlabel('Re(Visibility)')
    plt.ylabel('Im(Visibility)')
    plt.savefig(out_dir+'%s_cyl%d.png' % (vis_yy_str, args.cylinder))
    plt.close()

    plt.figure()
    for i, (k, v) in enumerate(vis_yy1_cal.iteritems()):
        vis_real = np.array(v).real
        vis_imag = np.array(v).imag
        c, m = cm[i][0], cm[i][1]
        plt.scatter(vis_real, vis_imag, c=c, marker=m)
    plt.xlim(-0.8e4, 0.8e4)
    plt.ylim(-0.8e4, 0.8e4)
    plt.xlabel('Re(Visibility)')
    plt.ylabel('Im(Visibility)')
    plt.savefig(out_dir+'%s_cyl%d_cal.png' % (vis_yy_str, args.cylinder))
    plt.close()
