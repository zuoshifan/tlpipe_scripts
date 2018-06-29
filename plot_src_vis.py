import os
from datetime import datetime
import argparse
import numpy as np
import h5py
import ephem
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import MaxNLocator, AutoMinorLocator


parser = argparse.ArgumentParser(description='Plot visibility components separated by stable PCA algorithm.')
parser.add_argument('-i', '--in_dir', type=str, default='./', help='Input directory where the visibility files are in.')
parser.add_argument('-o', '--out_dir', type=str, default='./', help='Output directory where the plotted figures should be put in.')
parser.add_argument('-f', '--freq', type=int, default=None, nargs='+',  help='Frequency indices to be plot.')
parser.add_argument('-s', '--srcs', type=str, default=['cyg'], nargs='+',  help='Point source names to be plotted.') # -s cyg cas crab
parser.add_argument('--feed1', type=int, default=[1], nargs='+',  help='Feed No. of the first element of baseline to be plot.') # --feed1 1 3 15 31 32 45 62 78 96
parser.add_argument('--feed2', type=int, default=[1], nargs='+',  help='Feed No. of the second element of baseline to be plot.') # --feed2 1 2 10 31 40 63 80 95

args = parser.parse_args()


in_dir = args.in_dir if args.in_dir.endswith('/') else ('%s/' % args.in_dir)
out_dir = args.out_dir if args.out_dir.endswith('/') else ('%s/' % args.out_dir)


def juldate2ephem(num):
    """Convert Julian date to ephem date, measured from noon, Dec. 31, 1899."""
    return ephem.date(num - 2415020.)


for src in args.srcs:
    in_fl = '%s_vis.hdf5' % src
    src_out_dir = '%s%s/' % (out_dir, src)
    if not os.path.isdir(src_out_dir):
        os.makedirs(src_out_dir)
    with h5py.File(in_dir + in_fl, 'r') as f:
        time = f.attrs['time']
        freq = f.attrs['freq']
        pol = f.attrs['pol']
        feed = f.attrs['feed'].tolist()
        sky_vis = f['sky_vis'][:]
        src_vis = f['src_vis'][:]
        otl_vis = f['outlier_vis'][:]

    dt = np.array([ ephem.Date(juldate2ephem(t) + 8*ephem.hour).datetime() for t in time ]) # python datetime
    xlabel = '%s' % dt[0].date()
    xval = mdates.date2num(dt)

    ifreqs = range(len(freq)) if args.freq is None else args.freq
    for fi in ifreqs:
        for pi, pol in enumerate(pol):
            for fd1 in args.feed1:
                fd1_idx = feed.index(fd1)
                for fd2 in args.feed2:
                    fd2_idx = feed.index(fd2)
                    sky_vis_plot = sky_vis[:, fi, pi, fd1_idx, fd2_idx]
                    src_vis_plot = src_vis[:, fi, pi, fd1_idx, fd2_idx]
                    otl_vis_plot = otl_vis[:, fi, pi, fd1_idx, fd2_idx]

                    plt.figure()
                    f, axarr = plt.subplots(3, sharex=True)
                    axarr[0].plot(xval, sky_vis_plot.real, label='total')
                    axarr[0].plot(xval, src_vis_plot.real, label=src)
                    axarr[0].plot(xval, otl_vis_plot.real, label='outlier')
                    axarr[0].plot(xval, (sky_vis_plot - src_vis_plot - otl_vis_plot).real, label='noise')
                    axarr[0].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=4, mode="expand", borderaxespad=0.)
                    axarr[1].plot(xval, sky_vis_plot.imag, label='total')
                    axarr[1].plot(xval, src_vis_plot.imag, label=src)
                    axarr[1].plot(xval, otl_vis_plot.imag, label='outlier')
                    axarr[1].plot(xval, (sky_vis_plot - src_vis_plot - otl_vis_plot).imag, label='noise')
                    axarr[2].plot(xval, np.abs(sky_vis_plot), label='total')
                    axarr[2].plot(xval, np.abs(src_vis_plot), label=src)
                    axarr[2].plot(xval, np.abs(otl_vis_plot), label='outlier')
                    axarr[2].plot(xval, np.abs(sky_vis_plot - src_vis_plot - otl_vis_plot), label='noise')
                    yl, yh = axarr[2].get_ylim()
                    ylim = (yl-0.05*(yh-yl), yh)
                    axarr[2].set_ylim(ylim)
                    axarr[2].set_xlim([xval[0], xval[-1]])
                    axarr[2].xaxis_date()
                    date_format = mdates.DateFormatter('%H:%M')
                    axarr[2].xaxis.set_major_formatter(date_format)
                    locator = MaxNLocator(nbins=6)
                    axarr[2].xaxis.set_major_locator(locator)
                    axarr[2].xaxis.set_minor_locator(AutoMinorLocator(2))
                    axarr[2].set_xlabel(xlabel)
                    plt.savefig(src_out_dir + 'vis_%d_%d_%d_%s.png' % (fd1, fd2, fi, pol))
                    plt.close()