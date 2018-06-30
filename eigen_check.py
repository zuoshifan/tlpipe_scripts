import os
import argparse
import numpy as np
from scipy import linalg as la
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description='Eigen-analysis.')
parser.add_argument('-i', '--in_file', type=str, default='./vis.hdf5', help='Input visibility file.')
parser.add_argument('-o', '--out_dir', type=str, default='./', help='Output directory where the plotted figures should be put in.')
parser.add_argument('-t', '--time', type=int, default=0,  help='Frequency indices to be plot.')
parser.add_argument('-f', '--freq', type=int, default=0,  help='Frequency indices to be plot.')

args = parser.parse_args()


out_dir = args.out_dir if args.out_dir.endswith('/') else ('%s/' % args.out_dir)

if not os.path.isdir(out_dir):
    os.makedirs(out_dir)


with h5py.File(args.in_file, 'r') as f:
    src_vis_xx = f['src_vis'][args.time, args.freq, 0]
    sky_vis_xx = f['sky_vis'][args.time, args.freq, 0]
    otl_vis_xx = f['outlier_vis'][args.time, args.freq, 0]
    src_vis_yy = f['src_vis'][args.time, args.freq, 1]
    sky_vis_yy = f['sky_vis'][args.time, args.freq, 1]
    otl_vis_yy = f['outlier_vis'][args.time, args.freq, 1]


if (~np.isfinite(sky_vis_xx)).all():
    print 'All values are invalid for ti = %d, fi = %d, pol = xx...' % (args.time, args.freq)
else:
    e, U = la.eigh(sky_vis_xx)
    e1, U1 = la.eigh(src_vis_xx)

    # print e
    plt.figure()
    plt.plot(range(1, 97), e[::-1], 'ro', label=r'$V$')
    plt.plot(range(1, 97), e1[::-1], 'go', label=r'$V_0$')
    plt.xlim(0, 97)
    yl, yh = plt.ylim()
    plt.ylim(0, yh)
    plt.xlabel('Feed number')
    plt.legend()
    plt.savefig(out_dir+'e_%d_%d_xx.png' % (args.time, args.freq))
    plt.close()

    plt.figure()
    plt.plot(range(1, 97), np.abs(U[:, -1]), 'r', label=r'$V$')
    # plt.plot(range(1, 97), np.abs(U[:, -2]), 'r')
    # plt.plot(range(1, 97), np.abs(U[:, -3]), 'r')
    plt.plot(range(1, 97), np.abs(U1[:, -1]), 'g', label=r'$V_0$')
    plt.xlim(0, 97)
    plt.xlabel('Feed number')
    plt.legend()
    plt.savefig(out_dir+'u_%d_%d_xx.png' % (args.time, args.freq))
    plt.close()

    plt.figure()
    plt.plot(range(1, 97), e[-1]**0.5 * np.abs(U[:, -1]), 'r', label=r'$V$')
    plt.plot(range(1, 97), e1[-1]**0.5 * np.abs(U1[:, -1]), 'g', label=r'$V_0$')
    plt.xlim(0, 97)
    plt.xlabel('Feed number')
    plt.legend()
    plt.savefig(out_dir+'eu_%d_%d_xx.png' % (args.time, args.freq))
    plt.close()

if (~np.isfinite(sky_vis_yy)).all():
    print 'All values are invalid for ti = %d, fi = %d, pol = yy...' % (args.time, args.freq)
else:
    e, U = la.eigh(sky_vis_yy)
    e1, U1 = la.eigh(src_vis_yy)

    # print e
    plt.figure()
    plt.plot(range(1, 97), e[::-1], 'ro', label=r'$V$')
    plt.plot(range(1, 97), e1[::-1], 'go', label=r'$V_0$')
    plt.xlim(0, 97)
    yl, yh = plt.ylim()
    plt.ylim(0, yh)
    plt.xlabel('Feed number')
    plt.legend()
    plt.savefig(out_dir+'e_%d_%d_yy.png' % (args.time, args.freq))
    plt.close()

    plt.figure()
    plt.plot(range(1, 97), np.abs(U[:, -1]), 'r', label=r'$V$')
    # plt.plot(range(1, 97), np.abs(U[:, -2]), 'r')
    # plt.plot(range(1, 97), np.abs(U[:, -3]), 'r')
    plt.plot(range(1, 97), np.abs(U1[:, -1]), 'g', label=r'$V_0$')
    plt.xlim(0, 97)
    plt.xlabel('Feed number')
    plt.legend()
    plt.savefig(out_dir+'u_%d_%d_yy.png' % (args.time, args.freq))
    plt.close()

    plt.figure()
    plt.plot(range(1, 97), e[-1]**0.5 * np.abs(U[:, -1]), 'r', label=r'$V$')
    plt.plot(range(1, 97), e1[-1]**0.5 * np.abs(U1[:, -1]), 'g', label=r'$V_0$')
    plt.xlim(0, 97)
    plt.xlabel('Feed number')
    plt.legend()
    plt.savefig(out_dir+'eu_%d_%d_yy.png' % (args.time, args.freq))
    plt.close()