import os
import argparse
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description='Plot gain.')
parser.add_argument('-i', '--in_file', type=str, default='gain.hdf5', help='Input file containts the gain data.')
parser.add_argument('-o', '--out_dir', type=str, default='./', help='Output directory where the plotted figures should be put in.')
parser.add_argument('-f', '--freq', type=int, default=None, nargs='+',  help='Frequency indices to be plot.')

args = parser.parse_args()


out_dir = args.out_dir if args.out_dir.endswith('/') else ('%s/' % args.out_dir)
if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

with h5py.File(args.in_file, 'r') as f:
    gain = f['gain'][:]
    freq = f['Gain'].attrs['freq']
    pols = f['Gain'].attrs['pol']
    feed = f['Gain'].attrs['feed']

nf, nfd = freq.shape[0], feed.shape[0]
print 'nf = %d...' % nf
f_inds = args.freq if not args.freq is None else range(nf)

for fi in f_inds:
    for pi, pol in enumerate(pols):
        if np.isfinite(gain[fi, pi, 0]):
            g = gain[fi, pi] * np.sign(gain[fi, pi, 0].real)
        else:
            print 'gain[0] = ', gain[fi, pi, 0], ' is invalid...'
            g = gain[fi, pi]

        if (~np.isfinite(g)).all():
            print 'All invalid values for fi = %d, pol = %s, skip...' % (fi, pol)
            continue

        # plot real, imag, abs
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
        plt.legend()
        plt.savefig(out_dir + 'gain_%d_%s.png' % (fi, pol))
        plt.close()

        # plot phase
        plt.figure()
        plt.plot(feed, np.angle(g), 'b-')
        plt.plot(feed, np.angle(g), 'bo')
        plt.xlabel('Feed number')
        plt.savefig(out_dir + 'gain_phase_%d_%s.png' % (fi, pol))
        plt.close()

        # plot scatter
        gg = []
        for di in range(nfd):
            # for dj in range(nfd):
            for dj in range(di, nfd):
                gg.append(g[di] * np.conj(g[dj]))

        gg = np.array(gg)
        plt.figure()
        plt.scatter(gg.real, gg.imag)
        plt.savefig(out_dir + 'gain_scatter_%d_%s.png' % (fi, pol))
        plt.close()
