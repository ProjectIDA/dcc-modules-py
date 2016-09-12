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
from datetime import datetime
import os.path
import yaml
import collections

from numpy import float32, logical_and, less_equal, greater_equal, greater, \
    polyfit, polyval, subtract, log10, ceil, floor
from numpy.fft import rfft, irfft
import matplotlib.pyplot as plt

from fabulous.color import red, bold
from obspy.core import read, UTCDateTime
from obspy.signal.invsim import evalresp

from ida.calibration.cross import cross_correlate


# for cleanup of (mostly) past crud.
bhz_chan_list = ['UN1', 'BHZ', 'HHZ', 'HNZ']
bhn_chan_list = ['UN2', 'BHN', 'HHN', 'HHY', 'HNY', 'BH1']
bhe_chan_list = ['UN3', 'BHE', 'HHE', 'HHX', 'HNX', 'BH2']

def rename_chan(inchan):
    if inchan in bhz_chan_list:
        return 'BHZ'
    elif inchan in bhn_chan_list:
        return 'BH1'
    elif inchan in bhe_chan_list:
        return 'BH2'
    else:
        return inchan


def read_intimes(fn):

    intimes_data = collections.OrderedDict()
    with open(fn, 'rt') as ifl:
        for line in ifl:
            cline = line.upper().strip()
            if cline.startswith('BH'):
                recitems = cline.split()
                # print(recitems)
                # print(cline)
                if (len(recitems) != 6) and (len(recitems) != 4):
                    print(red(bold('Error reading record: ' + line)))
                    print(red('from file ' + fn + '. '))
                    print(red('Skipping record.'))
                    continue

                chan = rename_chan(recitems[0][:3])
                if chan not in (bhz_chan_list + bhn_chan_list + bhe_chan_list):
                    print(red(bold('Unknown channel code: {}. Skipping record.'.format(chan))))
                    continue
                start_t =  UTCDateTime(datetime.strptime(recitems[2][:17], '%Y:%j-%H:%M:%S'))
                end_t = UTCDateTime(datetime.strptime(recitems[3][:17], '%Y:%j-%H:%M:%S'))
                intimes_data[chan] = {
                    'loc': recitems[0][3:5],
                    'ref_chan': rename_chan(recitems[1][:3]),
                    'starttime': start_t,
                    'endtime':  end_t,
                    'station': 'TRI',   # may be overwritten below
                    'network': 'II',    # ditto
                }
                if (len(recitems) == 6):
                    intimes_data[chan]['station'] = recitems[5]
                    intimes_data[chan]['network']= recitems[4]

                # print(chan, ':', intimes_data[chan])
    return intimes_data


def prepare_traces(meta, msfn):

    try:
        strm = read(msfn)
    except:
        print(red(bold('Error reading miniseed file: ' + msfn)))
        return {}

    # cleanup chan codes in miniseed
    for tr in strm:
        if tr.stats.channel in bhz_chan_list:
            tr.stats.channel = 'BHZ'
        elif tr.stats.channel in bhn_chan_list:
            tr.stats.channel = 'BH1'
        elif tr.stats.channel in bhe_chan_list:
            tr.stats.channel = 'BH2'
        else:
            raise ValueError('Unknown channel code encoutnered: ' + tr.stats.channel)

    # loop through metadata record with trace start/end time info
    # and drop and channels without waveform 'wf' ot 'wf_ref' data
    good_meta = collections.OrderedDict()
    for chan, vals in meta.items():
        wf = strm.select(channel=chan).copy().trim(starttime=vals['starttime'], endtime=vals['endtime'])
        wf_ref = strm.select(channel=vals['ref_chan']).copy().trim(starttime=vals['starttime'], endtime=vals['endtime'])
        # print('CHAN:', chan, wf)
        # print('REF:', vals['ref_chan'], wf_ref)
        if wf and wf_ref:
            # print('GOOD FOR ', chan)
            good_meta[chan] = vals
            good_meta[chan]['wf'] = wf[0]
            good_meta[chan]['wf_ref'] = wf_ref[0]
            good_meta[chan]['wf'].stats.network = vals['network']
            good_meta[chan]['wf'].stats.loc = vals['loc']
            good_meta[chan]['wf'].stats.station = vals['station']
            good_meta[chan]['wf_ref'].stats.network = vals['network']
            good_meta[chan]['wf_ref'].stats.loc = vals['loc']
            good_meta[chan]['wf_ref'].stats.station = vals['station']
        else:
            # print('BAD FOR ', chan)
            if not wf:
                print(red('No trace in {} for channel: {} during {} - {}'.format(msfn, chan,
                                                                                 vals['starttime'], vals['endtime'])))
            if not wf_ref:
                print(red('No trace in {} for reference channel: {} during {} - {}'.format(msfn, vals['ref_chan'],
                                                                                           vals['starttime'], vals['endtime'])))
            print(red('Skipping {} record.'.format(chan)))

    return good_meta

def save_chan_traces(chancode, fn, strm):

    try:
        strm.write(fn, format='MSEED')
    except:
        return False
    else:
        return True

def correlate_channel_traces(chan_trace, ref_trace, sample_rate, shake_m_per_volt, digi_sens_cnts_per_volt, **kwargs):

    cross_results = {}
    npts = ref_trace.stats.npts
    # comp = chan_trace.stats.channel[2]
    # ref_comp = ref_trace.stats.channel[2]
    # thedate = chan_trace.stats.starttime
    if not os.environ.get('SEEDRESP'):
        print(red('Error: Can not find RESP files without SEEDRESP env var being set.'))
        return False, cross_results
    else:
        resp_dir = os.environ.get('SEEDRESP')

    # construct RESP file filename for ref_trace
    resp_file = 'RESP.{}.{}.{}.{}'.format(
        ref_trace.stats.network,
        ref_trace.stats.station,
        ref_trace.stats.loc,
        ref_trace.stats.channel
    )
    resp_filepath = os.path.join(resp_dir, resp_file)

    fresp, f = evalresp(1/sample_rate,
                        ref_trace.stats.npts,
                        resp_filepath,
                        ref_trace.stats.starttime,
                        station=ref_trace.stats.station,
                        channel=ref_trace.stats.channel,
                        network=ref_trace.stats.network,
                        locid=ref_trace.stats.loc,
                        units='DIS', freq=True )

    # Convolving ref data with nominal response...
    refdata    = ref_trace.data.astype(float32)
    mean       = refdata.mean()
    refdata   -= mean
    refdata_fft    = rfft(refdata)
    refdata_fft   *= fresp
    ref_wth_resp  = irfft(refdata_fft, npts)
    ref_wth_resp *= (shake_m_per_volt/digi_sens_cnts_per_volt)


    # trim 20 smaples off both ends.
    ref_wth_resp = ref_wth_resp[20:-20]
    outdata = chan_trace.data[20:-20].astype(float32)

    if 'smoothing_factor' in kwargs:
        sf = kwargs['smoothing_factor']
    else:
        sf = 0.5
    # noinspection PyTupleAssignmentBalance
    freqs, amp, pha, coh, psd1, psd2, _, _, _ = cross_correlate(sample_rate,
                                                                outdata,
                                                                ref_wth_resp, smoothing_factor=sf)

    cross_results = {
        'freqs': freqs,
        'amp': amp,
        'pha': pha,
        'coh': coh,
        'psd1': psd1,
        'psd2': psd2,
        'cospect': [],
        'quadspect': [],
    }

    return cross_results

def shake_table_chan_plots(datadir, comp, cross_res_dict, coh_min=0.98, freq_band=(0.1, 10)):

    freqs = cross_res_dict['freqs']
    amp = cross_res_dict['amp']
    pha = cross_res_dict['pha']
    coh = cross_res_dict['coh']
    psd1 = cross_res_dict['psd1']
    psd2 = cross_res_dict['psd2']

    # get ndxs of good coh in freq_band
    use_freqs = logical_and(less_equal(freqs, freq_band[1]), greater_equal(freqs, freq_band[0]))

    fr = freqs[use_freqs]
    co = coh[use_freqs]
    ps1 = psd1[use_freqs]
    ps2 = psd2[use_freqs]
    am = amp[use_freqs]
    ph = pha[use_freqs]

    # comp = 'Z'
    fig_ndx = 0
    if comp == 'BHZ':
        fig_ndx = 0
    elif comp == 'BH1':
        fig_ndx = 10
    elif comp == 'BH2':
        fig_ndx = 30

    # plot psd for both time series
    fig1 = plt.figure(fig_ndx + 1, figsize=(8.5, 11))
    plt.subplot(311)
    plt.title('Fig 1: {} : {}'.format(datadir, comp))
    plt.suptitle('Big Title')
    plt.xlim(1e-1, 10)
    # plt.ylim(75, 110)
    plt.ylabel('PSD (dB)')
    plt.xlabel('Freq (Hz)')
    plt.grid(which='both')
    plt.semilogx(fr, 10 * log10(ps1))
    plt.semilogx(fr, 10 * log10(ps2))

    plt.subplot(312)
    plt.xlim(1e-1, 10)
    plt.ylim(0.97, 1.01)
    plt.ylabel('Coh**2')
    plt.xlabel('Freq (Hz)')
    plt.grid(which='both')
    plt.semilogx(fr, co)

    fig2 = plt.figure(fig_ndx + 2, figsize=(8.5, 11))
    plt.subplot(311)
    plt.title('Fig 2: {} : {}'.format(datadir, comp))
    plt.xlim(1e-1, 10)
    # plt.ylim(0.99, 1.01)
    plt.ylabel('TF Gain')
    plt.xlabel('Freq (Hz)')
    plt.grid(which='both')
    plt.semilogx(fr, am)

    plt.subplot(312)
    plt.xlim(1e-1, 10)
    plt.ylim(-10, 5)
    plt.ylabel('TF Phase')
    plt.xlabel('Freq (Hz)')
    plt.grid(which='both')
    plt.semilogx(fr, ph)

    plt.subplot(313)
    plt.xlim(1e-1, 10)
    plt.ylim(0.97, 1.01)
    plt.xlabel('Freq (Hz)')
    plt.ylabel('Coh**2')
    plt.grid(which='both')
    plt.semilogx(fr, co)

    # now just coh >= coh_min
    good_coh = greater(co, coh_min)
    fr = fr[good_coh]
    am = am[good_coh]
    ph = ph[good_coh]
    co = co[good_coh]

    # and remove trend/time shift from phase
    coeffs = polyfit(fr, ph, 1)
    # set y-intercept to 0
    coeffs = (coeffs[0], 0)
    correction = polyval(coeffs, fr)
    ph = subtract(ph, correction)

    fig3 = plt.figure(fig_ndx + 3, figsize=(8.5, 11))
    plt.subplot(311)
    plt.title('Fig 3: {} : {}'.format(datadir, comp))
    plt.xlim(1e-1, 10)
    # plt.ylim(0.99, 1.01)
    plt.ylabel('TF Gain')
    plt.xlabel('Freq (Hz)')
    plt.grid(which='both')
    plt.semilogx(fr, am)

    plt.subplot(312)
    plt.xlim(1e-1, 10)
    # plt.ylim(-10, 5)
    plt.ylabel('TF Phase')
    plt.xlabel('Freq (Hz)')
    plt.grid(which='both')
    plt.semilogx(fr, ph)

    plt.subplot(313)
    plt.xlim(1e-1, 10)
    plt.ylim(0.97, 1.01)
    plt.xlabel('Freq (Hz)')
    plt.ylabel('Coh**2')
    plt.grid(which='both')
    plt.semilogx(fr, co)

    return fig1, fig2, fig3

def shake_table_psd_plot(fignum, datadir, comp, freqs, psd1, psd2, coh, add_title_text=''):

    # plot psd for both time series
    fig = plt.figure(fignum, figsize=(8.5, 11))
    plt.subplot(311)
    plt.title('Shake Table PSDs: {} : {}'.format(datadir, comp) + add_title_text)
    plt.xlim(1e-1, 10)
    # plt.ylim(75, 110)

    plt.ylabel('PSD (dB)')
    plt.xlabel('Freq (Hz)')
    plt.grid(which='both')
    plt.semilogx(freqs, 10 * log10(psd1))
    plt.semilogx(freqs, 10 * log10(psd2))

    plt.subplot(312)
    plt.xlim(1e-1, 10)
    plt.ylim(0.9, 1.01)
    plt.ylabel('Coh**2')
    plt.xlabel('Freq (Hz)')
    plt.grid(which='both')
    plt.semilogx(freqs, coh)

    return fig

def shake_table_tf_plot(fignum, datadir, comp, freqs, amp, pha, coh, add_title_text=''):

    fig = plt.figure(fignum, figsize=(8.5, 11))
    plt.subplot(311)
    plt.title('Shake Table TF: {} : {}'.format(datadir, comp) + add_title_text)
    plt.xlim(1e-1, 10)
    plt.ylim(floor(amp.min() * 100) / 100, ceil(amp.max() * 100) / 100)
    axlim = plt.axis()

    plt.ylabel('TF Gain')
    plt.xlabel('Freq (Hz)')
    plt.grid(which='both')
    plt.semilogx(freqs, amp)

    plt.subplot(312)
    plt.xlim(1e-1, 10)
    # plt.ylim(-10, 5)
    plt.ylabel('TF Phase')
    plt.xlabel('Freq (Hz)')
    plt.grid(which='both')
    plt.semilogx(freqs, pha)

    plt.subplot(313)
    plt.xlim(1e-1, 10)
    plt.ylim(0.9, 1.01)
    plt.xlabel('Freq (Hz)')
    plt.ylabel('Coh**2')
    plt.grid(which='both')
    plt.semilogx(freqs, coh)

    return fig

class ShakeConfig(object):

    def __init__(self, fn):
        self._config = yaml.load('')
        with open(fn, 'rt') as cfl:
            config_txt = cfl.read()

        try:
            self._config = yaml.load(config_txt)
        except:
            print(red(bold('Error parsing config file: ' + fn)))
            print(config_txt)

    def digi_cnts_per_volt(self, digi_sn=None):
        return self._config['digi_cnts_per_volt']

    def shake_table_meters_per_volt(self, comp, thedate, shaketable_sn=None):

        if comp.upper() not in ['Z', '1', '2']:
            raise ValueError('Invalid component: ', comp)

        # noinspection PyUnresolvedReferences
        if not isinstance(thedate, datetime.datetime):
            raise TypeError('Invalid date: ' + str(thedate))

        datenum = int(thedate.strftime('%Y%j'))

        res = 0
        if comp.upper() == 'Z':
            for per in self._config['shaketable_vert_resp']:
                if (datenum >= per['startdate']) and (datenum <= per['enddate']):
                    res = per['meters_per_volt']
                    break
        else:
            for per in self._config['shaketable_hori_resp']:
                if (datenum >= per['startdate']) and (datenum <= per['enddate']):
                    res = per['meters_per_volt']
                    break

        return res

    def sample_rate(self):
        return self._config['sample_rate']

    def plot_min_freq(self):
        return self._config['plots']['start_freq']

    def plot_max_freq(self):
        return self._config['plots']['end_freq']

    def plot_coh_min(self):
        return self._config['plots']['coh_min']

    def analysis_smoothing_factor(self):
        return self._config['analysis']['smoothing_factor']
