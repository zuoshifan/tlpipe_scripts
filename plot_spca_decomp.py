import os
import argparse
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description='Plot visibility matrix components separated by stable PCA algorithm.')
parser.add_argument('-i', '--in_file', type=str, default='src_vis.hdf5', help='Input file containts the separated visibilities.')
parser.add_argument('-o', '--out_dir', type=str, default='./', help='Output directory where the plotted figures should be put in.')
parser.add_argument('-t', '--time', type=int, default=None, nargs='+',  help='Time indices to be plot.')
parser.add_argument('-f', '--freq', type=int, default=None, nargs='+',  help='Frequency indices to be plot.')
parser.add_argument('--min', type=float, default=None,  help='Min color scale.')
parser.add_argument('--max', type=float, default=None,  help='Max color scale.')

args = parser.parse_args()


out_dir = args.out_dir if args.out_dir.endswith('/') else ('%s/' % args.out_dir)
if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

with h5py.File(args.in_file, 'r') as f:
    time = f.attrs['time']
    freq = f.attrs['freq']
    pols = f.attrs['pol']
    feed = f.attrs['feed'].tolist()
    sky_vis = f['sky_vis'][:]
    src_vis = f['src_vis'][:]
    otl_vis = f['outlier_vis'][:]

nt, nf = time.shape[0], freq.shape[0]
print 'nt = %d, nf = %d...' % (nt, nf)
t_inds = args.time if not args.time is None else range(nt)
f_inds = args.freq if not args.freq is None else range(nf)

for ti in t_inds:
    for fi in f_inds:
        for pi, pol in enumerate(pols):
            Vmat = sky_vis[ti, fi, pi]
            V0 = src_vis[ti, fi, pi]
            S = otl_vis[ti, fi, pi]
            N = Vmat - V0 - S
            if (~np.isfinite(Vmat)).all():
                print 'All invalid values for ti = %d, fi = %d, pol = %s, skip...' % (ti, fi, pol)
                continue

            # plot Vmat
            plt.figure(figsize=(13, 5))
            plt.subplot(121)
            plt.imshow(Vmat.real, aspect='equal', origin='lower', interpolation='nearest', vmin=args.min, vmax=args.max)
            plt.colorbar(shrink=1.0)
            plt.subplot(122)
            plt.imshow(Vmat.imag, aspect='equal', origin='lower', interpolation='nearest', vmin=args.min, vmax=args.max)
            plt.colorbar(shrink=1.0)
            plt.savefig(out_dir + '%s_%d_%d_%s.png' % ('V', ti, fi, pol))
            plt.close()

            # plot V0
            plt.figure(figsize=(13, 5))
            plt.subplot(121)
            plt.imshow(V0.real, aspect='equal', origin='lower', interpolation='nearest', vmin=args.min, vmax=args.max)
            plt.colorbar(shrink=1.0)
            plt.subplot(122)
            plt.imshow(V0.imag, aspect='equal', origin='lower', interpolation='nearest', vmin=args.min, vmax=args.max)
            plt.colorbar(shrink=1.0)
            plt.savefig(out_dir + '%s_%d_%d_%s.png' % ('V0', ti, fi, pol))
            plt.close()

            # plot S
            plt.figure(figsize=(13, 5))
            plt.subplot(121)
            plt.imshow(S.real, aspect='equal', origin='lower', interpolation='nearest', vmin=args.min, vmax=args.max)
            plt.colorbar(shrink=1.0)
            plt.subplot(122)
            plt.imshow(S.imag, aspect='equal', origin='lower', interpolation='nearest', vmin=args.min, vmax=args.max)
            plt.colorbar(shrink=1.0)
            plt.savefig(out_dir + '%s_%d_%d_%s.png' % ('S', ti, fi, pol))
            plt.close()

            # plot N
            plt.figure(figsize=(13, 5))
            plt.subplot(121)
            plt.imshow(N.real, aspect='equal', origin='lower', interpolation='nearest', vmin=-1.1e8, vmax=1.1e8)
            plt.colorbar(shrink=1.0)
            plt.subplot(122)
            plt.imshow(N.imag, aspect='equal', origin='lower', interpolation='nearest', vmin=-1.1e8, vmax=1.1e8)
            plt.colorbar(shrink=1.0)
            plt.savefig(out_dir + '%s_%d_%d_%s.png' % ('N', ti, fi, pol))
            plt.close()
