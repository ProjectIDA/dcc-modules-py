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

from numpy import angle, pi, abs

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import matplotlib.gridspec as gridspec

from numpy import linspace, ceil

from ida.instruments import CALTYPE_RBLF, CALTYPE_RBHF

"""Convenience methods for plotting response and calibration results."""

def save_comp_response_comparison_plot(sta, chan, loc, resp1_fn, resp2_fn, seis_model, timestamp,
                                   operating_sample_rate, num_freqs, norm_freq, resp1,
                                   resp2, adev, pdev):
    """Generate plots of measured response deviations between two respoonses for
     a single component channel from a nominal response.

    :param sta: Station code
    :type sta: str
    :param chancodes: Channle codes for all three components
    :type chancodes: ComponentsTpl
    :param loc: Location code
    :type loc: str
    :param amp_fn: Fully qualified pathname for file name in which to save Amplitude plots
    :type amp_fn: str
    :param pha_fn: Fully qualified pathname for file name in which to save Phase plots
    :type pha_fn: str
    :param seis_model: Seismometer model code. Must be one of instruments.SEISMOMETER_MODELS
    :type seis_model: str
    :param timestamp: Timestamp to use on plot title. Typically time of data acquisition
    :type timestamp: datetime
    :param operating_sample_rate: Sample rates for channels being plotted. Assumed to all be the same rate.
    :type operating_sample_rate: float
    :param num_freqs: How many frequencies to plot. Frequency bins will be linearly spaced from 0 to nyquist.
    :type num_freqs: int
    :param resp1: Nominal/baseline frequency response being copmared to
    :type resp1: ndarray
    :param resp2: North/south channel measured frequency response
    :type resp2: ndarray
    :param adev: North/south channel response amplitude deviation from nominal response
    :type adev: ndarray
    :param pdev: North/south channel response phase deviation from nominal response
    :type pdev: ndarray
    """

    freqs = linspace(0, operating_sample_rate/2, num_freqs)  # must start with 0hz
    nyquist = operating_sample_rate / 2.0
    nyq_90pct_freqndx = int(ceil(0.9 * len(freqs)))

    datestr = timestamp.strftime('%Y-%m-%d %H:%M UTC')

    f100 = plt.figure(100, figsize=(14, 18))
    gspec = gridspec.GridSpec(2,1, hspace=.3)

    ax = plt.subplot(gspec[0])
    plt.tick_params(labelsize=13)
    title_fmt = 'Amplitude Response Comparison\n' +  \
                '{:<}-{:<}-{:<} ({:<} @ {}hz)\n' + \
                '{} / {}\n(Calibration Date: {})'
    plt.title(title_fmt.format(sta, loc, chan,
                               seis_model, int(round(operating_sample_rate, 0)),
                               resp1_fn, resp2_fn, datestr),
              fontsize=14,
              fontweight='bold')

    plt.ylabel('Amplitude (Normalized @ {} hz, V/m/s)'.format(norm_freq), fontsize=12, fontweight='bold')
    plt.xlabel('Frequency (hz)', fontsize=12, fontweight='bold')
    ax.grid(which='both')
    line2, = plt.loglog(freqs, abs(resp2), 'r:', linewidth=1.5)
    line1, = plt.loglog(freqs, abs(resp1), 'k', linewidth=0.5)
    plt.legend((line2, line1), (resp2_fn, resp1_fn), loc=0, handlelength=5, fontsize=12)
    #
    axlim = plt.axis()
    plt.axis([1e-3, nyquist, axlim[2], axlim[3]])
    #
    ax = plt.subplot(gspec[1])
    plt.tick_params(labelsize=13)
    plt.ylabel('Amplitude Deviation (%)\n(up to 90% of Nyquist)', fontsize=12, fontweight='bold')
    plt.xlabel('Frequency (hz)', fontsize=12, fontweight='bold')
    ax.grid(which='both')
    plt.semilogx(freqs[1:nyq_90pct_freqndx], adev[:nyq_90pct_freqndx-1], 'r', linewidth=0.75)
    amp_toler_pcnt = 5.0
    plt.axis([1e-3, nyquist, -amp_toler_pcnt * 2, amp_toler_pcnt * 2])
    axlim = plt.axis()
    within_tolerance_verts = [(axlim[0], -amp_toler_pcnt),
                              (axlim[0], amp_toler_pcnt),
                              (axlim[1], amp_toler_pcnt),
                              (axlim[1], -amp_toler_pcnt)]
    poly = Polygon(within_tolerance_verts, facecolor='#E6E6E6', edgecolor='0.9', label='Acceptable Tolerance Band')
    ax.add_patch(poly)
    # plt.clf()

    # Now plot for Phase
    f101 = plt.figure(101, figsize=(14, 18))
    ax = plt.subplot(211)
    plt.tick_params(labelsize=13)
    title_fmt = 'Phase Response Comparison\n' +  \
                '{:<}-{:<}-{:<} ({:<} @ {}hz)\n' + \
                '{} / {}\n(Calibration Date: {})'
    plt.title(title_fmt.format(sta, loc, chan,
                               seis_model, int(round(operating_sample_rate, 0)),
                               resp1_fn, resp2_fn, datestr),
              fontsize=14,
              fontweight='bold')
    plt.ylabel('Phase (deg)', fontsize=12, fontweight='bold')
    plt.xlabel('Frequency (hz)', fontsize=12, fontweight='bold')
    ax.grid(which='both')
    line2, = plt.semilogx(freqs, angle(resp2) * 180 / pi, 'r:', linewidth=1.5)
    line1, = plt.semilogx(freqs, angle(resp1) * 180 / pi, 'k', linewidth=0.5)
    plt.legend((line2, line1), (resp2_fn, resp1_fn), loc=0, handlelength=5, fontsize=12)
    axlim = plt.axis()
    plt.axis([1e-3, nyquist, axlim[2], axlim[3]])

    ax = plt.subplot(212)
    plt.tick_params(labelsize=13)
    plt.ylabel('Phase Deviation (deg)\n(up to 90% of Nyquist)', fontsize=12, fontweight='bold')
    plt.xlabel('Frequency (hz)', fontsize=12, fontweight='bold')
    ax.grid(which='both')
    line2, = plt.semilogx(freqs[1:nyq_90pct_freqndx], pdev[:nyq_90pct_freqndx-1], 'r', linewidth=0.75)
    pha_toler_degs = 5.0
    plt.axis([1e-3, nyquist, -pha_toler_degs * 2, pha_toler_degs * 2])
    poly = Polygon(within_tolerance_verts, facecolor='#E6E6E6', edgecolor='0.9', label='Acceptable Tolerance Band')
    ax.add_patch(poly)

    return f100, f101

def save_response_comparison_plots(sta, chancodes, loc, amp_fn, pha_fn, seis_model, timestamp,
                                   operating_sample_rate, num_freqs, norm_freq, nom_resp,
                                   n_resp, n_adev, n_pdev,
                                   e_resp, e_adev, e_pdev,
                                   v_resp, v_adev, v_pdev):
    """Generate plots of measured response deviations of three component channels from a nominal response.

    :param sta: Station code
    :type sta: str
    :param chancodes: Channle codes for all three components
    :type chancodes: ComponentsTpl
    :param loc: Location code
    :type loc: str
    :param amp_fn: Fully qualified pathname for file name in which to save Amplitude plots
    :type amp_fn: str
    :param pha_fn: Fully qualified pathname for file name in which to save Phase plots
    :type pha_fn: str
    :param seis_model: Seismometer model code. Must be one of instruments.SEISMOMETER_MODELS
    :type seis_model: str
    :param timestamp: Timestamp to use on plot title. Typically time of data acquisition
    :type timestamp: datetime
    :param operating_sample_rate: Sample rates for channels being plotted. Assumed to all be the same rate.
    :type operating_sample_rate: float
    :param num_freqs: How many frequencies to plot. Frequency bins will be linearly spaced from 0 to nyquist.
    :type num_freqs: int
    :param nom_resp: Nominal/baseline frequency response being copmared to
    :type nom_resp: ndarray
    :param n_resp: North/south channel measured frequency response
    :type n_resp: ndarray
    :param n_adev: North/south channel response amplitude deviation from nominal response
    :type n_adev: ndarray
    :param n_pdev: North/south channel response phase deviation from nominal response
    :type n_pdev: ndarray
    :param e_resp: North/south channel measured frequency response
    :type e_resp: ndarray
    :param e_adev: East/west channel response amplitude deviation from nominal response
    :type e_adev: ndarray
    :param e_pdev: East/West channel response phase deviation from nominal response
    :type e_pdev: ndarray
    :param v_resp: Vertical channel measured frequency response
    :type v_resp: ndarray
    :param v_adev: Vertical channel response amplitude deviation from nominal response
    :type v_adev: ndarray
    :param v_pdev: Vertical channel response phase deviation from nominal response
    :type v_pdev: ndarray
    """

    freqs = linspace(0, operating_sample_rate/2, num_freqs)  # must start with 0hz
    nyquist = operating_sample_rate / 2.0
    nyq_90pct_freqndx = int(ceil(0.9 * len(freqs)))

    datestr = timestamp.strftime('%Y-%m-%d %H:%M UTC')

    f100 = plt.figure(100, figsize=(15, 20))
    ax = plt.subplot(211)
    plt.tick_params(labelsize=13)
    if loc.strip() != '':
        plt.title('Amplitude Responses - ({:<} - {:<} - {:<} - {:<})\n(sampling rate: {} hz)'.format(sta,
                                                                                                     loc,
                                                                                                     seis_model,
                                                                                                     datestr,
                                                                                                     int(round(operating_sample_rate, 0))),
                  fontsize=14,
                  fontweight='bold')
    else:
        plt.title('Amplitude Responses - ({:<} - {:<} - {:<})\n(sampling rate: {} hz)'.format(sta,
                                                                                             seis_model,
                                                                                             datestr,
                                                                                             int(round(operating_sample_rate, 0))),
                  fontsize=14,
                  fontweight='bold')
    plt.ylabel('Amplitude (Normalized @ {} hz, V/m/s)'.format(norm_freq), fontsize=12, fontweight='bold')
    plt.xlabel('Frequency (hz)', fontsize=12, fontweight='bold')
    ax.grid(which='both')
    bh1, = plt.loglog(freqs, abs(n_resp), 'g', linewidth=0.5)
    bh2, = plt.loglog(freqs, abs(e_resp), 'y', linewidth=0.75)
    bhz, = plt.loglog(freqs, abs(v_resp), 'r', linewidth=0.5)
    nom, = plt.loglog(freqs, abs(nom_resp), 'k', linewidth=0.5)
    plt.legend((bh1, bh2, bhz, nom), (chancodes.north, chancodes.east, chancodes.vertical, 'Nominal'), loc=0, handlelength=5, fontsize=12)
    #
    axlim = plt.axis()
    plt.axis([1e-3, nyquist, axlim[2], axlim[3]])
    #
    ax = plt.subplot(212)
    plt.tick_params(labelsize=13)
    if loc.strip() != '':
        plt.title('Amplitude Deviations from Nominal - ({:<} - {:<} - {:<} - {:<})'.format(sta, loc, seis_model, datestr), fontsize=14, fontweight='bold')
    else:
        plt.title('Amplitude Deviations from Nominal - ({:<} - {:<} - {:<})'.format(sta, seis_model, datestr), fontsize=14, fontweight='bold')
    plt.ylabel('Amplitude Deviation (%)\n(up to 90% of Nyquist)', fontsize=12, fontweight='bold')
    plt.xlabel('Frequency (hz)', fontsize=12, fontweight='bold')
    ax.grid(which='both')
    bh1, = plt.semilogx(freqs[1:nyq_90pct_freqndx], n_adev[:nyq_90pct_freqndx-1], 'g', linewidth=0.75)
    bh2, = plt.semilogx(freqs[1:nyq_90pct_freqndx], e_adev[:nyq_90pct_freqndx-1], 'y', linewidth=0.75)
    bhz, = plt.semilogx(freqs[1:nyq_90pct_freqndx], v_adev[:nyq_90pct_freqndx-1], 'r', linewidth=0.75)
    plt.legend((bh1, bh2, bhz), (chancodes.north, chancodes.east, chancodes.vertical), loc=0, handlelength=5, fontsize=12)
    amp_toler_pcnt = 5.0
    plt.axis([1e-3, nyquist, -amp_toler_pcnt * 2, amp_toler_pcnt * 2])
    axlim = plt.axis()
    within_tolerance_verts = [(axlim[0], -amp_toler_pcnt),
                              (axlim[0], amp_toler_pcnt),
                              (axlim[1], amp_toler_pcnt),
                              (axlim[1], -amp_toler_pcnt)]
    poly = Polygon(within_tolerance_verts, facecolor='#D8FFD8', edgecolor='0.9', label='Acceptable Tolerance Band')
    ax.add_patch(poly)
    f100.savefig(amp_fn, dpi=400)
    plt.clf()

    f101 = plt.figure(101, figsize=(15, 20))
    ax = plt.subplot(211)
    plt.tick_params(labelsize=13)
    if loc.strip() != '':
        plt.title('Phase Responses - ({:<} - {:<} - {:<} - {:<})\n(sampling rate: {} hz)'.format(sta,
                                                                                                 loc,
                                                                                                 seis_model,
                                                                                                 datestr,
                                                                                                 int(round(operating_sample_rate, 0))),
                  fontsize=14,
                  fontweight='bold')
    else:
        plt.title('Phase Responses - ({:<} - {:<} - {:<})\n(sampling rate: {} hz)'.format(sta,
                                                                                         seis_model,
                                                                                         datestr,
                                                                                         int(round(operating_sample_rate, 0))),
                  fontsize=14,
                  fontweight='bold')

    plt.ylabel('Phase (deg)', fontsize=12, fontweight='bold')
    plt.xlabel('Frequency (hz)', fontsize=12, fontweight='bold')
    ax.grid(which='both')
    bh1, = plt.semilogx(freqs, angle(n_resp) * 180 / pi, 'g', linewidth=0.5)
    bh2, = plt.semilogx(freqs, angle(e_resp) * 180 / pi, 'y', linewidth=0.75)
    bhz, = plt.semilogx(freqs, angle(v_resp) * 180 / pi, 'r', linewidth=0.5)
    nom, = plt.semilogx(freqs, angle(nom_resp) * 180 / pi, 'k', linewidth=0.5)
    plt.legend((bh1, bh2, bhz, nom), (chancodes.north, chancodes.east, chancodes.vertical, 'Nominal'), loc=0, handlelength=5, fontsize=12)
    axlim = plt.axis()
    plt.axis([1e-3, nyquist, axlim[2], axlim[3]])

    ax = plt.subplot(212)
    plt.tick_params(labelsize=13)
    if loc.strip() != '':
        plt.title('Phase Deviations from Nominal - ({:<} - {:<} - {:<} - {:<})'.format(sta, loc, seis_model, datestr), fontsize=14, fontweight='bold')
    else:
        plt.title('Phase Deviations from Nominal - ({:<} - {:<} - {:<})'.format(sta, seis_model, datestr), fontsize=14, fontweight='bold')
    plt.ylabel('Phase Deviation (deg)\n(up to 90% of Nyquist)', fontsize=12, fontweight='bold')
    plt.xlabel('Frequency (hz)', fontsize=12, fontweight='bold')
    ax.grid(which='both')
    bh1, = plt.semilogx(freqs[1:nyq_90pct_freqndx], n_pdev[:nyq_90pct_freqndx-1], 'g', linewidth=0.75)
    bh2, = plt.semilogx(freqs[1:nyq_90pct_freqndx], e_pdev[:nyq_90pct_freqndx-1], 'y', linewidth=0.75)
    bhz, = plt.semilogx(freqs[1:nyq_90pct_freqndx], v_pdev[:nyq_90pct_freqndx-1], 'r', linewidth=0.75)
    plt.legend((bh1, bh2, bhz), (chancodes.north, chancodes.east, chancodes.vertical), loc=0, handlelength=5, fontsize=12)
    pha_toler_degs = 5.0
    plt.axis([1e-3, nyquist, -pha_toler_degs * 2, pha_toler_degs * 2])
    poly = Polygon(within_tolerance_verts, facecolor='#D8FFD8', edgecolor='0.9', label='Acceptable Tolerance Band')
    ax.add_patch(poly)
    #
    f101.savefig(pha_fn, dpi=400)
    plt.clf()

def apc_plot(sampling_freq, freqs, amp, pha, coh):

    """Generate simple plot that mirrors Matlab go_parker.m plots"""

    freq_max = 0.4 * sampling_freq
    plt.grid(which='both')

    plt.subplot(311)
    plt.semilogx(freqs, amp)
    plt.xlim(1e-3, freq_max)
    plt.ylim()
    plt.grid(which='both')

    plt.subplot(312)
    plt.semilogx(freqs, pha)
    plt.xlim(1e-3, freq_max)
    plt.grid(which='both')

    plt.subplot(313)
    plt.semilogx(freqs, coh)
    plt.xlim(1e-3, freq_max)
    plt.ylim(0.95, 1.05)
    plt.show()

def cross_tf_plot(sta: object, loc: object, chn: object, sensor: object, ondate: object, cal_type: object,
                  samp_rate: object, freqs: object, cr_amp: object, cr_pha: object, cr_coh: object, green_tol_lims: object = None, grey_tol_lims: object = None) -> object:
    '''Python port of go_parker.m plots of cross.f output'''

    band_limit = 0.9  # plot to 70% of nyquist
    nyq = samp_rate * 0.5
    freq_limit = nyq * band_limit

    freq_plt = [f for f in freqs if f <= freq_limit]
    amp_plt = cr_amp[:len(freq_plt)]
    pha_plt = cr_pha[:len(freq_plt)]
    coh_plt = cr_coh[:len(freq_plt)]

    fig = plt.figure(figsize=(8.5,11))

    gspec = gridspec.GridSpec(3,1, hspace=.5)
    title_substr = 'LOW' if cal_type == CALTYPE_RBLF else 'HIGH'
    # subp = plt.subplot(3,1,1)
    subp = plt.subplot(gspec[0])
    plt.title('{} TF on {}\n{} {}-{} ({})'.format(
        title_substr + ' Freq', ondate, sta.upper(), loc, chn.upper(), sensor.upper()
    ))
    plt.xlabel('Frequency')
    plt.ylabel('TF Amp')
    plt.grid(which='both')
    plt.semilogx(freq_plt, amp_plt)
    plt.xlim(freq_plt[0], freq_limit)
    plt.axis([samp_rate / 1e4, nyq, .94, 1.06])
    ax = plt.axis()
    if grey_tol_lims:
        tol_verts = [(ax[0], 1.0-(grey_tol_lims[0]/100.0)),
                     (ax[0], 1.0+(grey_tol_lims[0]/100.0)),
                     (ax[1], 1.0+(grey_tol_lims[0]/100.0)),
                     (ax[1], 1.0-(grey_tol_lims[0]/100.0))]
        poly = Polygon(tol_verts, facecolor='#E6E6E6', edgecolor='0.9', label='Grey Acceptable Tolerance Band')
        subp.add_patch(poly)
    if green_tol_lims:
        tol_verts = [(ax[0], 1.0-(green_tol_lims[0]/100.0)),
                     (ax[0], 1.0+(green_tol_lims[0]/100.0)),
                     (ax[1], 1.0+(green_tol_lims[0]/100.0)),
                     (ax[1], 1.0-(green_tol_lims[0]/100.0))]
        poly = Polygon(tol_verts, facecolor='#D8FFD8', edgecolor='0.9', label='Green Acceptable Tolerance Band')
        subp.add_patch(poly)


    subp = plt.subplot(gspec[1])
    plt.title('{} TF on {}\n{} {}-{} ({})'.format(
        title_substr + ' Freq', ondate, sta.upper(), loc, chn.upper(), sensor.upper()
    ))
    plt.xlabel('Frequency')
    plt.ylabel('TF Pha')
    plt.grid(which='both')
    plt.semilogx(freq_plt, pha_plt)
    plt.xlim(freq_plt[0], freq_limit)
    ax = plt.axis()
    plt.axis([samp_rate / 1e4, nyq, min(ax[2] - 1, -green_tol_lims[1] - 1), max(ax[3] + 1, green_tol_lims[1] + 1)])
    ax = plt.axis()
    if grey_tol_lims:
        tol_verts = [(ax[0], -grey_tol_lims[1]),
                     (ax[0], grey_tol_lims[1]),
                     (ax[1], grey_tol_lims[1]),
                     (ax[1], -grey_tol_lims[1])]
        poly = Polygon(tol_verts, facecolor='#E6E6E6', edgecolor='0.9', label='Acceptable Tolerance Band')
        subp.add_patch(poly)
    if green_tol_lims:
        tol_verts = [(ax[0], -green_tol_lims[1]),
                     (ax[0], green_tol_lims[1]),
                     (ax[1], green_tol_lims[1]),
                     (ax[1], -green_tol_lims[1])]
        poly = Polygon(tol_verts, facecolor='#D8FFD8', edgecolor='0.9', label='Acceptable Tolerance Band')
        subp.add_patch(poly)

    subp = plt.subplot(gspec[2])
    plt.title('{} TF on {}\n{} {}-{} ({})'.format(
        title_substr + ' Freq', ondate, sta.upper(), loc, chn.upper(), sensor.upper()
    ))
    plt.xlabel('Frequency')
    plt.ylabel('TF COH')
    plt.grid(which='both')
    plt.semilogx(freq_plt, coh_plt)
    plt.xlim(freq_plt[0], freq_limit)
    plt.axis([samp_rate / 1e4, nyq, 0.95, 1.0])

    return fig
