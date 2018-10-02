import os
import argparse
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description='Plot Gain.')
parser.add_argument('-i', '--in_file', type=str, default='gain.hdf5', help='Input file containts the gain data.')
parser.add_argument('-o', '--out_dir', type=str, default='./', help='Output directory where the plotted figures should be put in.')
parser.add_argument('-t', '--time', type=int, default=None, nargs='+',  help='Time indices to be plot.')
parser.add_argument('-f', '--freq', type=int, default=None, nargs='+',  help='Frequency indices to be plot.')

args = parser.parse_args()


out_dir = args.out_dir if args.out_dir.endswith('/') else ('%s/' % args.out_dir)
if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

with h5py.File(args.in_file, 'r') as f:
    Gain = f['Gain'][:]
    time = f['Gain'].attrs['time']
    freq = f['Gain'].attrs['freq']
    pols = f['Gain'].attrs['pol']
    feed = f['Gain'].attrs['feed']

nt, nf = time.shape[0], freq.shape[0]
print 'nt = %d, nf = %d...' % (nt, nf)
t_inds = args.time if not args.time is None else range(nt)
f_inds = args.freq if not args.freq is None else range(nf)

for ti in t_inds:
    for fi in f_inds:
        for pi, pol in enumerate(pols):
            g = Gain[ti, fi, pi]

            if (~np.isfinite(g)).all():
                print 'All invalid values for ti = %d, fi = %d, pol = %s, skip...' % (ti, fi, pol)
                continue

            plt.figure()
            plt.plot(feed, g.real, 'b-', label='real')
            plt.plot(feed, g.real, 'bo')
            plt.plot(feed, g.imag, 'g-', label='imag')
            plt.plot(feed, g.imag, 'go')
            plt.plot(feed, np.abs(g), 'r-', label='abs')
            plt.plot(feed, np.abs(g), 'ro')
            plt.xlim(np.min(feed)-1, np.max(feed)+1)
            yl, yh = plt.ylim()
            plt.ylim(yl, yh+(yh-yl)/5)
            plt.xlabel('Feed number')
            plt.ylabel(r'$\mathbf{G}$ / arbitrary units')
            plt.legend()
            plt.savefig(out_dir + 'Gain_%d_%d_%s.png' % (ti, fi, pol))
            plt.close()
