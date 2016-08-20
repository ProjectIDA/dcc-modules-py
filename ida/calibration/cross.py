#######################################################################################################################
# Copyright (C) 2016  Regents of the University of California
#
# This is free software: you can redistribute it and/or modify it under the terms of the
# GNU General Public License (GNU GPL) as published by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# A copy of the GNU General Public License can be found in LICENSE.TXT in the root of the source code repository.
# Additionally, it can be found at http://www.gnu.org/licenses/.
#
# NOTES: Per GNU GPLv3 terms:
#   * This notice must be kept in this source file
#   * Changes to the source must be clearly noted with date & time of change
#
# If you use this software in a product, an explicit acknowledgment in the product documentation of the contribution
# by Project IDA, Institute of Geophysics and Planetary Physics, UCSD would be appreciated but is not required.
#######################################################################################################################

from math import floor, atan2, sqrt
import logging
from numpy import ndarray, full, pi, sqrt, array, zeros, float64, concatenate, unwrap, sign, add, subtract, \
    rad2deg, max, min
from numpy.fft import fft

"""Python port of subst of cross.f Fortran code tailored with IDA-specific
parameter values.
"""

def cross_correlate(sampling_rate, ts1, ts2, smoothing_factor=2.0):
    """
    Compute coherence of and transfer function between two time series

    :param sampling_rate: Digitizing sampling rate
    :type sampling_rate: float
    :param ts1: First time series
    :type ts1: numpy.ndarry
    :param ts2: Second time series
    :type ts2: numpy.ndarray
    :return: Tuple of
        freqs: ndarray of frequencies,
        gain: ndarray of transfer function gain values at freqs frequencies
        phase: ndarray of transfer function phase values (in degrees) at freqs frequencies
        coh: ndarray of square of coherence between the two time series at freqs frequencies
    :rtype: (ndarray, ndarray, ndarray, ndarray)
    """

    if (type(ts1) != ndarray) or (type(ts2) != ndarray):
        msg = 'ERROR: Timeseries need to be of type list or numpy.ndarray.'
        raise TypeError(msg)

    logging.debug('Making copies of time series...')
    ts1_data = ts1.copy()
    ts2_data = ts2.copy()

    logging.debug('De-mean timeseries and capture variance...')
    # remove mean then calc variance
    ts1_mean = ts1_data.mean()
    ts1_data.__isub__(ts1_mean)
    ts1_var = ts1_data.var()
    ts2_mean = ts2_data.mean()
    ts2_data.__isub__(ts2_mean)
    # ts2_var = ts2_data.var()

    smoothing_factor = 2.0
    # calculate number of tapers
    # NOTE: inner floor probably not ideal
    taper_cnt = floor(floor((3.0 + 0.3 * sqrt(ts1_data.size))) * sqrt(smoothing_factor))
    logging.debug('Cross taper cnt: ' + str(taper_cnt))

    # with open('py-cross-pre-fft.txt', 'wt') as ofl:
    #     for r in range(ts_info['ts1']['data'].size):
    #         ofl.write('{} {}\n'.format(ts_info['ts1']['data'][r], ts_info['ts2']['data'][r]))

    logging.debug('calling spcmat...')

    sxy, fft_usable_len = spcmat(ts1_data, ts2_data, taper_cnt)

    logging.debug('calling spcmat... complete (fft-len: ' + str(fft_usable_len) + ')')

    power = 0.5 * (sxy[0, 0] + sxy[-1, 0]) + sxy[1:fft_usable_len - 1, 0].sum()
    # for ndx in range(1,ts_info['fft_usable_len']-1):
    #     power += sxy[ndx, 0]

    const = ts1_var / (power * (sampling_rate * 0.5 / (fft_usable_len - 1)))

    sxy *= const

    # not currently supporting 'adaptive tapers ( see Fortran adapt2 subroutine). so setting kopt to it's constant val
    kopt = full(fft_usable_len, taper_cnt, dtype=int)

    freq_bin_size = (sampling_rate / 2.0) / (fft_usable_len - 1)

    freqs = array([freq_bin_size * ndx for ndx in range(fft_usable_len)], dtype=float64)
    gain = zeros(fft_usable_len, dtype=float64)
    coh = zeros(fft_usable_len, dtype=float64)
    phase = zeros(fft_usable_len, dtype=float64)

    for freqndx in range(fft_usable_len):
        coh[freqndx] = (sxy[freqndx, 2] ** 2 + sxy[freqndx, 3] ** 2) / (sxy[freqndx, 0] * sxy[freqndx, 1])
        gain[freqndx] = sqrt(coh[freqndx] * sxy[freqndx, 1] / sxy[freqndx, 0])
        phase[freqndx] = atan2(sxy[freqndx, 3], sxy[freqndx, 2])  # phase in radians

    # # kmin=nf
    # #     kmax=0
    # #     kbar=0
    # #     do 2000 j=1, nf
    # #       gamsq=(sxy(j,3)**2 + sxy(j,4)**2)/(sxy(j,1)*sxy(j,2))
    # #       gain=sqrt(gamsq*sxy(j,2)/sxy(j,1))
    # #       phase=deg* atan2( sxy(j,4),  sxy(j,3))
    # #       if (named.ge. 0)
    # #    $    write(1,'(1p,e11.5,2e11.4,3e12.4,0p,f8.5,f8.2,i5)')
    # #    $    (j-1)*df,(sxy(j,k),k=1,4),gain,gamsq,phase,kopt(j)
    # #       kmin=min(kopt(j), kmin)
    # #       kmax=max(kopt(j), kmax)
    # #       kbar=kbar + kopt(j)
    # #       call lag(j, phase, gamsq)
    # # SKIPPING LAG SINCE WE KNOW TIME SERIES ARE IN SYNC wrt TIME
    # # 2000 continue
    # #     kbar=kbar/nf

    print('Phase: 1st, min, max:', phase[0], min(phase), max(phase))
    phase = unwrap(phase, discont=pi*1.05)  # phase in degrees
    print('Unwrapped Phase: 1st, min, max:', phase[0], min(phase), max(phase))
    # if max(phase) > pi:
    #     phase  = subtract(phase, pi)
    # elif min(phase) < -pi:
    #     phase = add(phase, pi)
    # print(min(phase), max(phase))
    phase = rad2deg(phase)  # phase in degrees

    del ts1_data
    del ts2_data

    return freqs, gain, phase, coh, sxy[:0], sxy[:1], sxy[:2], sxy[:3], kopt

# def lag(ts_info, fndx, phase, gamsq):
#
#     if fndx == 0:
#         ts_info['sw'] = 0.0
#         ts_info['swx'] = 0.0
#     else:
#         x = phase - ts_info['pho']
#         if abs(x) <= 200.0:
#             ts_info['sw'] += gamsq + ts_info['gmo']
#             ts_info['swx'] += (gamsq + ts_info['gmo']) * x
#     ts_info['pho'] = phase
#     ts_info['gmo'] = gamsq

#      subroutine lag(j, phase, gamsq)
# c$$$$ calls nothing
# c  If j=0 return mean and std deviation in phase, gamsq
# c  if j=1 initialize
# c  if j>1 find various sums related to the sum of delta-phi
# c     weighted by coherence
# c
#      common /hold/ pho,gmo,sw,swx
# c
# c  Estimates mean and its variance
#      if (j .eq. 0) then
#        phase = swx/sw
# c  Initializes sums and other parameters
#      elseif (j .eq. 1) then
#        pho = phase
#        gmo = gamsq
#        sw = 0.0
#        swx = 0.0
#      else
#        x = phase - pho
# c  Exclude phase wrap
#        if (abs(x) .le. 200.0) then
#          w = gamsq + gmo
#          sw = sw + w
#          swx = swx + w*x
#        endif
#        pho = phase
#        gmo = gamsq
#      endif
#      return
#      end
#


def spcmat(ts1, ts2, taper_cnt):
    """Python implementation of the cross.f spcmat() routine with IDA fixed parameters"""

    opt_len = int((ts1.size // 2) * 2)
    # opt_len = 139968  # just for testing aginst cross.f
    pad_len = 2 * opt_len

    # ts_info['fft_usable_len'] = opt_len // 2
    fft_usable_len = opt_len // 2

    logging.debug('cross.spcmat() fft_usable_len: ' + str(fft_usable_len))

    padded = concatenate([ts1[0:opt_len], zeros(opt_len)])
    ts1_fft = fft(padded).conjugate()
    padded = concatenate([ts2[0:opt_len], zeros(opt_len)])
    ts2_fft = fft(padded).conjugate()

    sxy = zeros([fft_usable_len, 4], dtype=float64)
    klim = taper_cnt
    ck = 1.0 / klim ** 2
    wt = 6.0 * klim / (4 * klim ** 2 + 3 * klim - 1)

    for freqndx in range(fft_usable_len):
        freqndx2 = freqndx * 2

        for taperndx in range(1, klim + 1):
            j1 = (freqndx2 + pad_len - taperndx) % pad_len
            j2 = (freqndx2 + taperndx) % pad_len
            z1 = ts1_fft[j1] - ts1_fft[j2]
            z2 = ts2_fft[j1] - ts2_fft[j2]

            totwt = wt * (1.0 - ck * (taperndx - 1) ** 2)

            sxy[freqndx, 0] += (totwt * (z1.real ** 2 + z1.imag ** 2))
            sxy[freqndx, 1] += (totwt * (z2.real ** 2 + z2.imag ** 2))
            sxy[freqndx, 2] += (totwt * (z1.real * z2.real + z1.imag * z2.imag))
            sxy[freqndx, 3] += (totwt * (z2.real * z1.imag - z1.real * z2.imag))

    # ojzfl.close()
    # c  Loop over frequency[
    #       d] 1500 m[0, nf]+1
    #         j=m+1
    #         klim=kopt(j)
    #         ck=1.0/klim**2
    #         m2=2*m
    # c
    #         do 1300 i=1, 4
    #           sxy(j,i)=0.0
    #  1300   continue
    # c
    # c  Average over the tapers with parabolic weight wk
    #         do 1400 k=1, klim
    #           j1=mod(m2+n2-k, n2)
    #           j2=mod(m2+k, n2)
    #           z1=y(j1+1,1) - y(j2+1,1)
    #           z2=y(j1+1,2) - y(j2+1,2)
    #           wk=1.0 - ck*(k-1)**2
    #           sxy(j,1)=sxy(j,1)+wk*(real(z1)**2 + aimag(z1)**2)
    #           sxy(j,2)=sxy(j,2)+wk*(real(z2)**2 + aimag(z2)**2)
    #           sxy(j,3)=sxy(j,3)+wk*(real(z1)*real(z2)+aimag(z1)*aimag(z2))
    #           sxy(j,4)=sxy(j,4)+wk*(real(z2)*aimag(z1)-real(z1)*aimag(z2))
    #  1400   continue
    # *       wt=1.0/real(klim)

    # c  Exact normalization for parabolic weighting
    #         wt=6.0*klim/real(4*klim**2+3*klim-1)
    #         do 1450 i=1, 4
    #           sxy(j,i)=sxy(j,i)*wt
    #  1450   continue
    # c
    #  1500 continue

    return sxy, fft_usable_len, kopt


# def prepare_for_cross(seis_model, input_strm, output_strm, cal_log, paz):
#     MODELS_SUPPORTED = ['STS2.5']
#
#     if seis_model not in MODELS_SUPPORTED:
#         raise ValueError('Unsupported seismometer model: {}, must be one of {}'.format(
#             seis_model,
#             MODELS_SUPPORTED))
#
#     if len(output_strm) != 1:
#         raise ValueError('There should be a single Trace in the output Stream object.')
#
#     if output_strm[0].stats.channel[2] != 'Z':
#         raise ValueError('Output Trace must be a vertical channel')
#
#     if len(input_strm) != 1:
#         raise ValueError('There should be a single Trace in the input Stream object.')
#
#     # invert qcal signals if needed
#     if qcal.invert_signal(input_strm[0].stats.channel, seis_model):
#         input_strm[0].data *= -1.0
#     if qcal.invert_signal(output_strm[0].stats.channel, seis_model):
#         output_strm[0].data *= -1.0
#
#     # trim settling and trailing times from traces
#     op.trim_stream(input_strm, left=cal_log['settling_time'], right=cal_log['trailing_time'])
#     op.trim_stream(output_strm, left=cal_log['settling_time'], right=cal_log['trailing_time'])
#     samp_count = input_strm[0].stats.npts
#
#     samprate = input_strm[0].stats.sampling_rate
#     taper_size = 0.1
#     taper = signal.tukey(samp_count, alpha=(taper_size * 2))
#
#     rb_tapered = np.multiply(input_strm[0].data, taper)
#
#     # read normative PAZ response and create FAP freq resp
#     resp, resp_freqs = ida.response.utils.pz2fap(input_strm[0].data.size,
#                                                  paz.poles('acc', 'rad'),
#                                                  paz.zeros('acc', 'rad'),
#                                                  1.0,
#                                                  samprate)
#     resp_freqs /= 2.0 * np.pi
#     resp_norm, scale, ndx = ida.response.utils.normalize(resp, resp_freqs, 0.05)
#
#     # convolve input RB with sensor normative response, inv transform and divide out taper
#     rb_fft = np.fft.rfft(rb_tapered)
#     rb_fft_conv = np.multiply(rb_fft, resp_norm)
#     rb_with_resp = np.fft.irfft(rb_fft_conv, samp_count)
#     rb_processed = np.divide(rb_with_resp, taper)
#
#     # remove taper bin_count
#     taper_bin_count = taper_size * samp_count
#     data_start = taper_bin_count
#     data_end = rb_processed.size - taper_bin_count
#
#     rb_cross_trimmed = rb_processed[data_start:data_end]
#
#     # DEBUG
#     # with open(os.path.join('', 'py-resp.txt'), 'w') as ofl:
#     #     abs(resp).tofile(ofl, sep="\n")
#     #     ofl.write('\n')
#     #
#     # with open(os.path.join('', 'py-resp-norm.txt'), 'w') as ofl:
#     #     resp_norm.tofile(ofl, sep="\n")
#     #     ofl.write('\n')
#     #
#     # with open(os.path.join('', 'py-resp-norm-mag.txt'), 'w') as ofl:
#     #     abs(resp_norm).tofile(ofl, sep="\n")
#     #     ofl.write('\n')
#     #
#     # with open(os.path.join('', 'py-rbin-FFT.txt'), 'w') as ofl:
#     #     rb_fft.tofile(ofl, sep="\n")
#     #     ofl.write('\n')
#     #
#     # with open(os.path.join('', 'py-rbin-FFT-with-resp.txt'), 'w') as ofl:
#     #     rb_fft_conv.tofile(ofl, sep="\n")
#     #     ofl.write('\n')
#     #
#     # with open(os.path.join('', 'py-rbin-with-resp.txt'), 'w') as ofl:
#     #     rb_with_resp.tofile(ofl, sep="\n")
#     #     ofl.write('\n')
#     #
#     # rbin w/resp, detaperred
#     # with open(os.path.join('', 'py-rbin-processed.txt'), 'w') as ofl:
#     #     rb_processed.tofile(ofl, sep="\n")
#     #     ofl.write('\n')
#
#     # with open(os.path.join('', 'py-rbin-processed-trimmed.txt'), 'w') as ofl:
#     #     rb_cross_trimmed.tofile(ofl, sep="\n")
#     #     ofl.write('\n')
#
#     with open(os.path.join('', 'py-rbout.txt'), 'w') as ofl:
#         output_strm[0].data.tofile(ofl, sep="\n")
#         ofl.write('\n')
#
#     rb_cross_trimmed.__itruediv__(rb_cross_trimmed.std())
#     # rb_cross_trimmed.__isub__(rb_cross_trimmed.mean())
#
#     outdata = output_strm[0].data.copy()
#     # change to after taper sections removed, per email 2016-04-06 from PD
#     # outdata.__itruediv__(outdata.std())
#     # outdata.__isub__(outdata.mean())
#     outdata = outdata[data_start:data_end]
#     outdata.__itruediv__(outdata.std())
#     outdata.__isub__(outdata.mean())
#
#     with open(os.path.join('', 'py-rbout-normed.txt'), 'w') as ofl:
#         outdata.tofile(ofl, sep="\n", format="%15.5f")
#         ofl.write('\n')
#
#     # DEBUG
#     with open(os.path.join('', 'py-cross-data.txt'), 'w') as ofl:
#         for ndx, val in enumerate(rb_cross_trimmed):
#             ofl.write("{:15.10f} {:15.10f}\n".format(val, outdata[ndx]))
#
#     snr = 1.0 / (outdata - rb_cross_trimmed).std()
#
#     plt.figure()
#     plt.grid(True)
#     stats = output_strm[0].stats
#     plt.title('{}.{}.{}.{}\nSNR: {:3.3f}'.format(stats['network'],
#                                                  stats['station'],
#                                                  stats['location'],
#                                                  stats['channel'], snr))
#     plt.plot(outdata)
#     plt.plot(outdata - rb_cross_trimmed, '-r')
#     plt.show()
#
#     return rb_cross_trimmed, outdata
#
#
