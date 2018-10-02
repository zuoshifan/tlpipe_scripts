import os
from datetime import datetime
import argparse
import numpy as np
from scipy import linalg as la
import h5py
import ephem
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import MaxNLocator, AutoMinorLocator


parser = argparse.ArgumentParser(description='Eigen-analysis.')
parser.add_argument('-i', '--in_file', type=str, default='./vis.hdf5', help='Input visibility file.')
parser.add_argument('-o', '--out_dir', type=str, default='./', help='Output directory where the plotted figures should be put in.')
parser.add_argument('-f', '--freq', type=int, default=0,  help='Frequency indices to be plot.')
parser.add_argument('--utc', type=str, default='8h',  help='Time zone of time axis, default 8h for Beijing time, else 0h.') # --utc 0h

args = parser.parse_args()


out_dir = args.out_dir if args.out_dir.endswith('/') else ('%s/' % args.out_dir)

if not os.path.isdir(out_dir):
    os.makedirs(out_dir)


def juldate2ephem(num):
    """Convert Julian date to ephem date, measured from noon, Dec. 31, 1899."""
    return ephem.date(num - 2415020.)


with h5py.File(args.in_file, 'r') as f:
    src_vis_xx = f['src_vis'][:, args.freq, 0]
    sky_vis_xx = f['sky_vis'][:, args.freq, 0]
    otl_vis_xx = f['outlier_vis'][:, args.freq, 0]
    src_vis_yy = f['src_vis'][:, args.freq, 1]
    sky_vis_yy = f['sky_vis'][:, args.freq, 1]
    otl_vis_yy = f['outlier_vis'][:, args.freq, 1]
    time = f.attrs['time']

nt = src_vis_xx.shape[0]

if args.utc == '0h':
    dt = np.array([ ephem.Date(juldate2ephem(t)).datetime() for t in time ]) # python datetime
elif args.utc == '8h':
    dt = np.array([ ephem.Date(juldate2ephem(t) + 8*ephem.hour).datetime() for t in time ]) # python datetime
else:
    raise ValueError('Unsupported utc: %s' % args.utc)
xlabel = '%s' % dt[0].date()
if args.utc == '0h':
    xlabel = 'UT+0h: ' + xlabel
xval = mdates.date2num(dt)


es = []
es2 = []
e1s = []
e2s = []
for ti in range(nt):
    if (~np.isfinite(sky_vis_xx[ti])).all():
        print 'All values are invalid for ti = %d, fi = %d, pol = xx...' % (ti, args.freq)
        es.append(np.nan)
        es2.append(np.nan)
        e1s.append(np.nan)
        e2s.append(np.nan)
    else:
        e, U = la.eigh(sky_vis_xx[ti])
        e1, U1 = la.eigh(src_vis_xx[ti])
        e2, U2 = la.eigh(sky_vis_xx[ti] - src_vis_xx[ti])

        es.append(e[-1])
        es2.append(e[-2])
        e1s.append(e1[-1])
        e2s.append(e2[-1])

# rescale es
es = np.array(es) * 1.0e-10
es2 = np.array(es2) * 1.0e-10
e1s = np.array(e1s) * 1.0e-10
e2s = np.array(e2s) * 1.0e-10

# plot es
plt.figure()
# plt.plot(range(nt), es, 'b', label=r'$V \, 1$')
# plt.plot(range(nt), es2, 'b--', label=r'$V \, 2$')
# plt.plot(range(nt), e1s, 'g', label=r'$V_0$')
# plt.plot(range(nt), e2s, 'r', label=r'$V-V_0$')
plt.plot(xval, es, 'b', label=r'$V \, 1$')
plt.plot(xval, es2, 'b--', label=r'$V \, 2$')
plt.plot(xval, e1s, 'g', label=r'$V_0$')
plt.plot(xval, e2s, 'r', label=r'$V-V_0$')
plt.legend()
plt.xlim(xval[0], xval[-1])
ax = plt.gca()
ax.xaxis_date()
date_format = mdates.DateFormatter('%H:%M')
ax.xaxis.set_major_formatter(date_format)
locator = MaxNLocator(nbins=6)
ax.xaxis.set_major_locator(locator)
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.set_xlabel(xlabel)
# plt.xlabel('Integral time / 4 second')
ax.set_ylabel('Variation of the largest eigenvalue')
plt.savefig(out_dir+'es_%d_xx.png' % args.freq)
plt.close()


es = []
es2 = []
e1s = []
e2s = []
for ti in range(nt):
    if (~np.isfinite(sky_vis_yy[ti])).all():
        print 'All values are invalid for ti = %d, fi = %d, pol = yy...' % (ti, args.freq)
        es.append(np.nan)
        es2.append(np.nan)
        e1s.append(np.nan)
        e2s.append(np.nan)
    else:
        e, U = la.eigh(sky_vis_yy[ti])
        e1, U1 = la.eigh(src_vis_yy[ti])
        e2, U2 = la.eigh(sky_vis_yy[ti] - src_vis_yy[ti])

        es.append(e[-1])
        es2.append(e[-2])
        e1s.append(e1[-1])
        e2s.append(e2[-1])

# rescale es
es = np.array(es) * 1.0e-10
es2 = np.array(es2) * 1.0e-10
e1s = np.array(e1s) * 1.0e-10
e2s = np.array(e2s) * 1.0e-10

# plot es
plt.figure()
# plt.plot(range(nt), es, 'b', label=r'$V \, 1$')
# plt.plot(range(nt), es2, 'b--', label=r'$V \, 2$')
# plt.plot(range(nt), e1s, 'g', label=r'$V_0$')
# plt.plot(range(nt), e2s, 'r', label=r'$V-V_0$')
plt.plot(xval, es, 'b', label=r'$V \, 1$')
plt.plot(xval, es2, 'b--', label=r'$V \, 2$')
plt.plot(xval, e1s, 'g', label=r'$V_0$')
plt.plot(xval, e2s, 'r', label=r'$V-V_0$')
plt.legend()
plt.xlim(xval[0], xval[-1])
ax = plt.gca()
ax.xaxis_date()
date_format = mdates.DateFormatter('%H:%M')
ax.xaxis.set_major_formatter(date_format)
locator = MaxNLocator(nbins=6)
ax.xaxis.set_major_locator(locator)
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.set_xlabel(xlabel)
# plt.xlabel('Integral time / 4 second')
ax.set_ylabel('Variation of the largest eigenvalue')
plt.savefig(out_dir+'es_%d_yy.png' % args.freq)
plt.close()