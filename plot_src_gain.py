import os
import os
import argparse
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description='Plot gain solved by different point source calibrators.')
parser.add_argument('-i', '--in_dir', type=str, default='./', help='Input directory where the visibility files are in.')
parser.add_argument('-o', '--out_dir', type=str, default='./', help='Output directory where the plotted figures should be put in.')
parser.add_argument('-s', '--srcs', type=str, default=['cyg'], nargs='+',  help='Point source names to be plotted.') # -s cyg cas crab

args = parser.parse_args()


in_dir = args.in_dir if args.in_dir.endswith('/') else ('%s/' % args.in_dir)
out_dir = args.out_dir if args.out_dir.endswith('/') else ('%s/' % args.out_dir)

if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

gr = []
gi = []
ga = []
gp = []
for si, src in enumerate(args.srcs):
    in_fl = '%s_gain.hdf5' % src
    with h5py.File(in_dir + in_fl, 'r') as f:
        gain = f['gain'][:]
        if si == 0:
            freq = f['gain'].attrs['freq']
            pols = f['gain'].attrs['pol']
            feed = f['gain'].attrs['feed']

    nf, npol, nfd = gain.shape

    for fi in range(nf):
        gr.append([])
        gi.append([])
        ga.append([])
        gp.append([])
        for pi, pol in enumerate(pols):
            gr[fi].append([])
            gi[fi].append([])
            ga[fi].append([])
            gp[fi].append([])

            # g = gain[fi, pi]
            g = gain[fi, pi] * np.sign(gain[fi, pi, 0].real)

            gr[fi][pi].append(g.real)
            gi[fi][pi].append(g.imag)
            ga[fi][pi].append(np.abs(g))
            phs = np.angle(g)
            # print phs[0], g[0]
            gp[fi][pi].append(phs)


# colors to use
cs = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
nc = len(cs)

for fi in range(nf):
    for pi, pol in enumerate(pols):
        plt.figure()
        f, axarr = plt.subplots(3, sharex=True)
        for si, src in enumerate(args.srcs):
            axarr[0].plot(feed, gr[fi][pi][si], cs[si%nc], label=src)
            axarr[0].plot(feed, gr[fi][pi][si], cs[si%nc]+'o')
            axarr[1].plot(feed, gi[fi][pi][si], cs[si%nc], label=src)
            axarr[1].plot(feed, gi[fi][pi][si], cs[si%nc]+'o')
            axarr[2].plot(feed, ga[fi][pi][si], cs[si%nc], label=src)
            axarr[2].plot(feed, ga[fi][pi][si], cs[si%nc]+'o')
        axarr[0].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=4, mode="expand", borderaxespad=0.)
        axarr[2].set_xlim([np.min(feed)-1, np.max(feed)+1])
        plt.savefig(out_dir + 'gain_%d_%s.png' % (fi, pol))
        plt.close()

        plt.figure()
        for si, src in enumerate(args.srcs):
            plt.plot(feed, gp[fi][pi][si], cs[si%nc], label=src)
            plt.plot(feed, gp[fi][pi][si], cs[si%nc]+'o')
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=4, mode="expand", borderaxespad=0.)
        plt.xlim([np.min(feed)-1, np.max(feed)+1])
        plt.ylabel('phase / radian')
        plt.savefig(out_dir + 'phase_%d_%s.png' % (fi, pol))
        plt.close()
