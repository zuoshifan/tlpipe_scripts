import numpy as np
from scipy import optimize
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


plot = True
# plot = False
plot_gain = True
plot_An = True
plot_Gain = True

# file_dir = '/tianlai/zuoshifan/cyl_20160927/output_cal_cyg1/gain/'
# file_dir = '/tianlai/zuoshifan/cyl_20160927/output_cal_cyg1/gain_new_rpca_1200_link/'
# file_dir = '/tianlai/zuoshifan/cyl_20160927/output_map/gain/'
# file_dir = '/tianlai/zuoshifan/cyl_20160927/output_cal_spca/gain/'
# file_dir = '/tianlai/zuoshifan/cyl_20160927/output_plot_phase/gain/'
# file_dir = '/tianlai/zuoshifan/cyl_20170903/output_cal/gain/'
# file_dir = '/tianlai/zuoshifan/cyl_20160927/output_cal_timing/gain/'
file_dir = '/tianlai/zuoshifan/cyl_20171209/output_check/gain/'
file_name = 'gain.hdf5'
with h5py.File(file_dir+file_name, 'r') as f:
    Gain = f['Gain'][:]
    gain = f['gain'][:]
    time = f['Gain'].attrs['time']
    freq = f['Gain'].attrs['freq']
    pols = f['Gain'].attrs['pol']
    feed = f['Gain'].attrs['feed']


nt, nf, npol, nfeed = Gain.shape
extent = [0, nt, 1, nfeed+1]
An = np.abs(Gain/gain[np.newaxis, :, :, :])

# one sidereal day in second
sday = 86164.0905

# plot gain
if plot and plot_gain:
    for pi, pol in enumerate(pols):
        for fi in range(nf):
        # for fi in range(300, 303):
            plt.figure()
            plt.plot(feed, gain[fi, pi].real, 'b-', label='real')
            plt.plot(feed, gain[fi, pi].real, 'bo')
            plt.plot(feed, gain[fi, pi].imag, 'g-', label='imag')
            plt.plot(feed, gain[fi, pi].imag, 'go')
            plt.plot(feed, np.abs(gain[fi, pi]), 'r-', label='abs')
            plt.plot(feed, np.abs(gain[fi, pi]), 'ro')
            plt.xlim(feed[0]-1, feed[-1]+1)
            yl, yh = plt.ylim()
            plt.ylim(yl, yh+(yh-yl)/5)
            plt.xlabel('Feed number')
            plt.legend()
            fig_name = 'gain_%d_%s.png' % (fi, pol)
            plt.savefig(file_dir+fig_name)
            plt.close()


# plot An
if plot and plot_An:
    for pi, pol in enumerate(pols):
        for fi in range(nf):
        # for fi in range(300, 303):
            plt.figure()
            # plt.imshow(np.abs(An[:, fi, pi, :]).T, extent=extent, origin='lower', aspect='auto', vmax=1.0)
            plt.imshow(np.abs(An[:, fi, pi, :]).T, extent=extent, origin='lower', aspect='auto')
            plt.axhline(31+1, color='k')
            plt.axhline(31+32+1, color='k')
            plt.colorbar()
            plt.xlim(0, nt)
            plt.ylim(1, nfeed+1)
            plt.xlabel('Integral time / 4 second')
            plt.ylabel('Feed number')
            fig_name = 'An_wf_%d_%s.png' % (fi, pol)
            plt.savefig(file_dir+fig_name)
            plt.close()


# plot Gain
if plot and plot_Gain:
    for pi, pol in enumerate(pols):
        for fi in range(nf):
        # for fi in range(300, 303):
            plt.figure()
            # plt.imshow(np.abs(Gain[:, fi, pi, :]).T, extent=extent, origin='lower', aspect='auto', vmax=450)
            plt.imshow(np.abs(Gain[:, fi, pi, :]).T, extent=extent, origin='lower', aspect='auto')
            plt.axhline(31+1, color='k')
            plt.axhline(31+32+1, color='k')
            plt.colorbar()
            plt.xlim(0, nt)
            plt.ylim(1, nfeed+1)
            plt.xlabel('Integral time / 4 second')
            plt.ylabel('Feed number')
            fig_name = 'Gain_wf_%d_%s.png' % (fi, pol)
            plt.savefig(file_dir+fig_name)
            plt.close()


# Equation for Gaussian
def fg(x, a, b, c, d):
    return a * np.exp(-(x - b)**2.0 / (2 * c**2)) + d

def fc(x, a, b, c, d):
    return a * np.sinc(c * (x - b)) + d

def shift(a, val=1, pad=0):
    if val == 0:
        return a

    b = np.roll(a, val)
    if val > 0:
        b[:val] = pad
    if val < 0:
        b[val:] = pad

    return b

# choose data slice near the transit time
c = nt/2 # center ind
# li = max(0, c - 100)
# hi = min(nt, c + 100 + 1)
li = max(0, c - 150)
hi = min(nt, c + 150 + 1)
# li = max(0, c - 90)
# hi = min(nt, c + 90 + 1)
x = np.arange(li, hi)
# save normalized Gain and its center
Gain_norm = np.zeros_like(Gain)
center = np.zeros((2, nf, nfeed))
center[:] = np.nan
FWHM_xx = []
FWHM_yy = []
fl = open(file_dir+'fwhm.txt', 'w')
for pi, pol in enumerate(pols):
    for fi in range(nf):
    # for fi in range(300, 303):
        # flag exceptional value
        data = np.abs(Gain[:, fi, pi, :]).T
        median = np.ma.median(data, axis=0)
        abs_diff = np.ma.abs(data - median[np.newaxis, :])
        mad = np.ma.median(abs_diff, axis=0) / 0.6745
        data = np.where(abs_diff>3.0*mad[np.newaxis, :], np.nan, data)

        # plot
        if plot:
            plt.figure()
            plt.imshow(data, extent=extent, origin='lower', aspect='auto', interpolation='nearest')
            plt.axhline(31+1, color='k')
            plt.axhline(31+32+1, color='k')
            plt.colorbar()
            plt.xlabel('Integral time / 4 second')
            plt.ylabel('Feed number')
            fig_name = 'Gain_flag_wf_%d_%s.png' % (fi, pol)
            plt.savefig(file_dir+fig_name)
            plt.close()

        # gaussian/sinc fit
        for idx, fd in enumerate(feed):
            plot_fit = False
            y = data[idx, li:hi]
            inds = np.where(np.isfinite(y))[0]
            if len(inds) > 0.75 * len(y):
                # get the best estimate of the central val
                cval = y[inds[np.argmin(np.abs(inds-c))]]
                try:
                    # gaussian fit
                    # popt, pcov = optimize.curve_fit(fg, x[inds], y[inds], p0=(cval, c, 90, 0))
                    # sinc function seems fit better
                    popt, pcov = optimize.curve_fit(fc, x[inds], y[inds], p0=(cval, c, 1.0e-2, 0))
                    print popt
                    center[pi, fi, idx] = popt[1]
                    Gain_norm[:, fi, pi, idx] = shift(Gain[:, fi, pi, idx] / fc(popt[1], *popt), int(np.around(c - popt[1])))
                    plot_fit = True
                except RuntimeError:
                    print 'curve_fit failed for fi = %d, pol = %s, feed = %d' % (fi, pol, feed[idx])
                    pass

            # # plot fit
            # plt.figure()
            # plt.plot(x, y)
            # if plot_fit:
            #     plt.plot(x, fc(x, *popt))
            # fig_name = 'Gain_sinc_fit_%d_%d_%s.png' % (fd, fi, pol)
            # plt.savefig(file_dir+fig_name)
            # plt.close()

        # plot
        if plot:
            plt.figure()
            plt.imshow(data, extent=extent, origin='lower', aspect='auto', interpolation='nearest')
            plt.axhline(31+1, color='k')
            plt.axhline(31+32+1, color='k')
            plt.colorbar()
            plt.scatter(center[pi, fi], range(1, 97))
            plt.xlim(0, nt)
            plt.ylim(1, nfeed+1)
            plt.xlabel('Integral time / 4 second')
            plt.ylabel('Feed number')
            fig_name = 'Gain_flag_wf_center_%d_%s.png' % (fi, pol)
            plt.savefig(file_dir+fig_name)
            plt.close()

        # flag Gin_norm
        data_norm = np.abs(Gain_norm[:, fi, pi, :]).T
        median = np.ma.median(data_norm, axis=0)
        abs_diff = np.ma.abs(data_norm - median[np.newaxis, :])
        mad = np.ma.median(abs_diff, axis=0) / 0.6745
        data_norm = np.where(abs_diff>3.0*mad[np.newaxis, :], np.nan, data_norm)

        # plot Gain_norm
        if plot:
            plt.figure()
            plt.imshow(data_norm, extent=extent, origin='lower', aspect='auto', interpolation='nearest')
            plt.axhline(31+1, color='k')
            plt.axhline(31+32+1, color='k')
            plt.colorbar()
            plt.xlim(0, nt)
            plt.ylim(1, nfeed+1)
            plt.xlabel('Integral time / 4 second')
            plt.ylabel('Feed number')
            fig_name = 'Gain_norm_wf_%d_%s.png' % (fi, pol)
            plt.savefig(file_dir+fig_name)
            plt.close()

        # plot sum of Gain_norm
        data_mask = np.ma.masked_invalid(data_norm)
        data_sum = data_mask.sum(axis=0) / data_mask.count(axis=0)
        if plot:
            plt.figure()
            plt.plot(data_sum)
            plt.xlabel('Integral time / 4 second')
            fig_name = 'Gain_norm_sum_%d_%s.png' % (fi, pol)
            plt.savefig(file_dir+fig_name)
            plt.close()

        # combined beam profile fit
        xs = []
        xss = []
        ys = []
        yss = []
        for idx, fd in enumerate(feed):
            inds = np.where(np.isfinite(data_norm[idx, li:hi]))[0]
            xs = xs + x[inds].tolist()
            ys = ys + data_norm[idx, li:hi][:, inds].tolist()
            inds1 = np.where(np.isfinite(data_norm[idx]))[0]
            xss = xss + inds1.tolist()
            yss = yss + data_norm[idx, inds1].tolist()
        # sinc fit
        poptc, pcovc = optimize.curve_fit(fc, xs, ys, p0=(1.0, c, 1.0e-2, 0))
        # gauss fit
        poptg, pcovg = optimize.curve_fit(fg, xs, ys, p0=(1.0, c, 90, 0))

        # print fc(x, *poptc) - fg(x, *poptg)

        # plot combined fit
        if plot:
            plt.figure()
            # plt.scatter(xs, ys)
            # plt.plot(x, fc(x, *poptc), 'r-', linewidth=2)
            # plt.plot(x, fg(x, *poptg), 'g-', linewidth=2)
            x1 = np.arange(nt)
            plt.scatter(xss, yss)
            slable = r'sinc:  $a \, \sinc{(c \, (x - b))} + d$' + '\n' + r'$a = %.2f, \, b = %.2f, \, c = %.2f, \, d = %.2f$' % tuple(poptc)
            glable = r'Gaussian:  $a \, \exp{(- (x - \mu)^2 / {2 \sigma^2})} + d$' + '\n' + r'$a = %.2f, \, \mu = %.2f, \, \sigma = %.2f, \, d = %.2f$' % tuple(poptg)
            plt.plot(x1, fc(x1, *poptc), 'r-', linewidth=2, label=slable)
            plt.plot(x1, fg(x1, *poptg), 'g-', linewidth=2, label=glable)
            plt.axhline(0, color='k')
            plt.ylim(-0.5, 1.05)
            lg = plt.legend(loc=8) # lower center
            lg.draw_frame(False)
            plt.xlabel('Integral time / 4 second')
            fig_name = 'combined_fit_%d_%s.png' % (fi, pol)
            plt.savefig(file_dir+fig_name)
            plt.close()

        # compute FWHM of Gaussian fit
        FWHM = 2 * (2 * np.log(2))**0.5 * poptg[2]
        FWHM_degree = FWHM * 4.0 * 360.0 / sday
        if pol == 'xx':
            FWHM_xx.append(FWHM_degree)
        elif pol == 'yy':
            FWHM_yy.append(FWHM_degree)
        msg = 'FWHM of the Gaussian fit for freq = %g, pol = %s: %g degree' % (freq[fi], pol, FWHM_degree)
        print msg
        fl.write('%s\n' % msg)


fl.write('\n\n\n')
fl.write('%s\n\n' % freq)
fl.write('%s\n\n' % FWHM_xx)
fl.write('%s\n\n' % FWHM_yy)
fl.close()

# plot FWHM vs freq
plt.figure()
plt.plot(freq, FWHM_xx, label='XX-pol')
plt.plot(freq, FWHM_yy, label='YY-pol')
plt.xlabel(r'$\nu$ / MHz')
plt.ylabel('FWHM of beam / degree')
plt.legend()
plt.savefig(file_dir+'fwhm.png')
plt.close()

# print center.shape
# print np.ma.masked_invalid(center.flatten()).mean()
# print np.ma.masked_invalid(center[0].flatten()).mean()
# print np.ma.masked_invalid(center[1].flatten()).mean()
# print np.ma.masked_invalid(center[0, :, :31].flatten()).mean()
# print np.ma.masked_invalid(center[0, :, 31:63].flatten()).mean()
# print np.ma.masked_invalid(center[0, :, 63:].flatten()).mean()
# print np.ma.masked_invalid(center[1, :, :31].flatten()).mean()
# print np.ma.masked_invalid(center[1, :, 31:63].flatten()).mean()
# print np.ma.masked_invalid(center[1, :, 63:].flatten()).mean()
