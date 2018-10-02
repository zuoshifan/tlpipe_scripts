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
parser.add_argument('-f', '--freq', type=int, default=None, nargs='+',  help='Frequency indices to be plot.')

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
    ifreqs = range(nf) if args.freq is None else args.freq

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

for fi in ifreqs:
    for pi, pol in enumerate(pols):
        # plot real, imag, abs
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
        axarr[2].set_xlabel('Feed number')
        plt.savefig(out_dir + 'gain_%d_%s.png' % (fi, pol))
        plt.close()

        # plot phs
        plt.figure()
        for si, src in enumerate(args.srcs):
            plt.plot(feed, gp[fi][pi][si], cs[si%nc], label=src)
            plt.plot(feed, gp[fi][pi][si], cs[si%nc]+'o')
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=4, mode="expand", borderaxespad=0.)
        plt.xlim([np.min(feed)-1, np.max(feed)+1])
        plt.xlabel('Feed number')
        plt.ylim([-4, 4])
        plt.gca().set_yticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
        plt.gca().set_yticklabels([r'$-\pi$', r'$-\frac{\pi}{2}$', r'$0$', r'$\frac{\pi}{2}$', r'$\pi$'])
        plt.ylabel('phase / radian')
        plt.savefig(out_dir + 'phase_%d_%s.png' % (fi, pol))
        plt.close()

        # plot abs and phs
        plt.figure()
        f, axarr = plt.subplots(2, sharex=True)
        for si, src in enumerate(args.srcs):
            axarr[0].plot(feed, ga[fi][pi][si], cs[si%nc], label=src)
            axarr[0].plot(feed, ga[fi][pi][si], cs[si%nc]+'o')
            axarr[1].plot(feed, gp[fi][pi][si], cs[si%nc], label=src)
            axarr[1].plot(feed, gp[fi][pi][si], cs[si%nc]+'o')
        axarr[0].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=4, mode="expand", borderaxespad=0.)
        axarr[0].set_ylabel('Amplitude')
        axarr[1].set_xlabel('Feed number')
        axarr[1].set_ylabel('Phase / radian')
        axarr[1].set_xlim([np.min(feed)-1, np.max(feed)+1])
        axarr[1].set_ylim([-4, 4])
        axarr[1].set_yticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
        axarr[1].set_yticklabels([r'$-\pi$', r'$-\frac{\pi}{2}$', r'$0$', r'$\frac{\pi}{2}$', r'$\pi$'])
        plt.subplots_adjust(hspace=0.1)
        plt.savefig(out_dir + 'gain_amp+phs_%d_%s.png' % (fi, pol))
        plt.close()

        # plot abs and phs and phs diff
        plt.figure()
        f, axarr = plt.subplots(3, sharex=True)
        for si, src in enumerate(args.srcs):
            axarr[0].plot(feed, ga[fi][pi][si], cs[si%nc], label=src)
            axarr[0].plot(feed, ga[fi][pi][si], cs[si%nc]+'o')
            axarr[1].plot(feed, gp[fi][pi][si], cs[si%nc], label=src)
            axarr[1].plot(feed, gp[fi][pi][si], cs[si%nc]+'o')
            pdiff = gp[fi][pi][si] - gp[fi][pi][0]
            pdiff = np.where(pdiff>np.pi, pdiff-2*np.pi, pdiff)
            pdiff = np.where(pdiff<-np.pi, pdiff+2*np.pi, pdiff)
            axarr[2].plot(feed, pdiff, cs[si%nc], label=src)
            axarr[2].plot(feed, pdiff, cs[si%nc]+'o')
        axarr[0].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=4, mode="expand", borderaxespad=0.)
        axarr[0].set_ylabel('Amplitude')
        # axarr[1].set_xlabel('Feed number')
        axarr[1].set_ylabel('Phase / radian')
        axarr[1].set_xlim([np.min(feed)-1, np.max(feed)+1])
        axarr[1].set_ylim([-4, 4])
        axarr[1].set_yticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
        axarr[1].set_yticklabels([r'$-\pi$', r'$-\frac{\pi}{2}$', r'$0$', r'$\frac{\pi}{2}$', r'$\pi$'])
        axarr[2].set_xlabel('Feed number')
        axarr[2].set_ylabel('Phase diff / degree')
        # plt.subplots_adjust(hspace=0.1)
        plt.savefig(out_dir + 'gain_amp+phs+diff_%d_%s.png' % (fi, pol))
        plt.close()

        # compare RMS between pairs of calibrators
        fl = open(out_dir+'RMS.txt', 'w')
        for s1 in range(len(args.srcs)):
            for s2 in range(s1+1, len(args.srcs)):
                pdiff = gp[fi][pi][s1] - gp[fi][pi][s2]
                pdiff = np.where(pdiff>np.pi, pdiff-2*np.pi, pdiff)
                pdiff = np.where(pdiff<-np.pi, pdiff+2*np.pi, pdiff)
                pdiff1 = pdiff[np.isfinite(pdiff)]
                rms = np.sqrt(np.sum(pdiff1**2) / len(pdiff1))
                # rms = np.sqrt(np.sum(pdiff1**2) / len(feed))
                msg = 'For pol %s RMS between %s and %s: %f radian, %f degree' % (pol, args.srcs[s1], args.srcs[s2], rms, np.degrees(rms))
                print msg
                fl.write('%s\n' % msg)
        fl.close()
