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
import logging

from obspy.core import Trace, Stream  #, UTCDateTime
from obspy.signal.cross_correlation import correlate, xcorr_max
from scipy.signal import freqs
from scipy.signal.ltisys import zpk2tf
from numpy import array, ndarray, isclose, abs, mod, divide, multiply, pi, exp, cos, sin, \
        angle, absolute, complex128, asarray, less
from numpy.fft import rfft

from fabulous.color import red, bold

import ida.signals.paz
# from ida.signals.trace import IDATrace
import ida.calibration.qcal_utils
from ida.instruments import SEIS_INVERT_CAL_CHAN, SEIS_INVERT_NORTH_CHAN, SEIS_INVERT_EAST_CHAN


TAPER_TYPES = [
    'tukey'
]


def time_offset(trace1, trace2, bpfreqmin=0.1, bpfreqmax=1.0, winsize=1000):
    """ Timeseries time offset computed using correlation computation.

        xcorr() does not look at start/end times, just aligns the samples starting at 0th sample.

        Returns offset in seconds that when added to trace2 the actual time of events
        will be the same.

    """
    if trace1.stats.sampling_rate != trace2.stats.sampling_rate:
        print(red(bold('The supplied traces must have the same sampling rate. Cannot compute time_offset.')))
        return 0, 0.0, [], 'The supplied traces must have the same sampling rate. Cannot compute time_offset.'

    if Stream([trace1, trace2]).get_gaps():
        print(red(bold('The supplied traces have gaps. Cannot compute time_offset.')))
        return 0, 0.0, [], 'The supplied traces have gaps. Cannot compute time_offset.'

    # leave originals alone...
    tr1 = trace1.copy()
    tr2 = trace2.copy()

    # xcorr deprecated in obspy 1.3: 
    # index, val, corrfun = xcorr(tr1, tr2, winsize, full_xcorr=True)
    corrfun = correlate(tr1, tr2, winsize)
    index, val = xcorr_max(corrfun)

    offset = index / trace1.stats.sampling_rate

    return offset, val, corrfun, ''

def taper_high_freq_resp(resp, taper_fraction):
    """
    This comes from idaresponse/resp.c where highest 5% of complex frequency response
    is tapered using exponential function.

    resp: complex response to be tapered
    taper_fraction: decimal fraction of high freq to taper

    Tapering is done IN-PLACE and updated resp is type np.complex128.
    """

    if (taper_fraction < 0.0) or (taper_fraction > 1.0):
        raise ValueError('taper_high_freq_resp: taper_fraction must be between 0.0 and 1.0')
    resp = asarray(resp, dtype=complex128)
    if len(resp) == 0:
        raise ValueError('taper_high_freq_resp: resp must not be empty')

    # exponentially taper off highest fraction of freqs before nyquist
    # taken from idaresponse/resp.c
    resplen = len(resp)
    startndx = int( taper_fraction * resplen)
    decay_rng = range(startndx, resplen)
    factors = [exp(5 * (ndx - startndx) / (resplen - startndx)) for ndx in decay_rng]

    pha = angle(resp[decay_rng])
    amp = factors * absolute(resp[decay_rng])
    resp[decay_rng] = complex128(amp*cos(pha) + amp*sin(pha))


def dynlimit_resp_min(resp, dyn_range):
    """
    This comes from idaresponse/resp.c where dynamic range of response is clipped
    by setting minimum value for amplitude response.

    resp: complex response to be tapered
    dyn_range: amplitude dynamic range limit. Lower resp amplitudes clipped to 
        (max amplitude)/dyn_range.

    resp is modified IN-PLACE and updated resp is type np.complex128.
    """

    resp = asarray(resp, dtype=complex128)
    if len(resp) == 0:
        raise ValueError('dynlimit_resp_min: resp must not be empty')

    amp = absolute(resp)
    max_amp = amp.max()
    limit_min = max_amp / dyn_range
    amp_flag = less(amp, limit_min)
    # save phase for freqs below minimum amp limit
    phavals = angle(resp[amp_flag])
    resp[amp_flag] = limit_min * complex128(cos(phavals) + sin(phavals)*1j)


def check_and_fix_polarities(strm, seis_model):

    for tr in strm:
        if invert_signal(tr.stats.channel, seis_model):
            tr.data *= -1.0

def invert_signal(channel, seis_model):

    if len(channel) < 3:
        msg = "Invalid channel code '{}'".format(channel)
        logging.error(msg)
        raise Exception(msg)

    seis_model = seis_model.upper()
    component = channel[2].upper()

    invert = ((component in ['S', 'F']) and (seis_model in SEIS_INVERT_CAL_CHAN)) or \
             ((component in ['1', 'N']) and (seis_model in SEIS_INVERT_NORTH_CHAN)) or \
             ((component in ['2', 'E']) and (seis_model in SEIS_INVERT_EAST_CHAN))

    return invert


# def taper_trace(trace, taper_type, both_sides=True, fraction=0.1, remove=False):
#
#     if not isinstance(trace, obspy.core.trace.Trace):
#         msg = 'Invalid type for trace. Expected an ObsPy Trace object.'
#         logging.error(msg)
#         raise TypeError(msg)
#
#     if not isinstance(taper_type, str):
#         msg = 'Invalid value type for taper_type: taper_type must be a string.'
#         logging.error(msg)
#         raise TypeError(msg)
#
#     if taper_type == 'tukey':
#         taper = scipy.signal.tukey(trace.data.size, alpha=fraction*2, sym=both_sides)
#         taper_bin_cnt = trace.data.size * fraction
#     else:
#         msg = "Invalid Taper Type requested: '{}'. Valid values: {}".format(
#             taper_type,
#             TAPER_TYPES)
#         logging.error(msg)
#         raise ValueError(msg)
#
#     if not remove:
#         trace.data *= taper
#     else:
#         # will introduce NaN into trace data where taper == 0
#         trace.data /= taper
#
#     return taper_bin_cnt
#
#
def trim_stream(src_stream, left=0, right=0):

    for trace in src_stream:
        trace.trim(trace.stats.starttime + left, trace.stats.endtime - right)


def ntrim_stream(traces, left=0, right=0):

    for trace in traces:
        trace.trim(left, right)


def compute_response(freqlist, paz, mode='vel'):

    if isinstance(freqlist, float) or isinstance(freqlist, int):
        freqlist = array([freqlist])

    b, a = zpk2tf(paz.zeros(units='rad', mode=mode),
                                      paz.poles(units='rad', mode=mode),
                                      paz.h0)
    # a has to be a list for the scipy.signal.freqs() call later but zpk2tf()
    # strangely returns it as an integer.
    if not isinstance(a, ndarray) and a == 1.0:
        a = [1.0]

    _, h = freqs(b, a, freqlist * 2 * pi)

    return h


def compute_response_fir(fir_coeffs, fft_len):

    fir_fft = rfft(fir_coeffs, fft_len)

    return fir_fft


def normalize_response(freq_resp, freqlist, norm_freq):

    normed = None
    # find the index in freqs of the first freq >= nom_freq

    ndx = min([freq[0] for freq in enumerate(freqlist) if freq[1] >= norm_freq])
    if not (ndx is None):

        scale = abs(freq_resp[ndx])
        normed = divide(freq_resp, scale)
    else:
        scale = 1.0

    return normed, scale, ndx


def channel_xform(trace_tpl, xfrm):

    #  traces and xfrms assumed to be ENZ (aka 21Z and XYZ for triaxial seis output) order
    chn_2chrs = trace_tpl[0].stats.channel[:2]

    tr_2 = Trace(header=trace_tpl[0].stats.copy(), data=(multiply(trace_tpl[0].data, xfrm[0][0]) +
                                                         multiply(trace_tpl[1].data, xfrm[0][1]) +
                                                         multiply(trace_tpl[2].data, xfrm[0][2]))
                   )
    tr_2.stats.channel = chn_2chrs + '2'

    tr_1 = Trace(header=trace_tpl[0].stats.copy(), data=(multiply(trace_tpl[0].data, xfrm[1][0]) +
                                                         multiply(trace_tpl[1].data, xfrm[1][1]) +
                                                         multiply(trace_tpl[2].data, xfrm[1][2]))
    )
    tr_1.stats.channel = chn_2chrs + '1'

    tr_z = Trace(header=trace_tpl[0].stats.copy(), data=(multiply(trace_tpl[0].data, xfrm[2][0]) +
                                                         multiply(trace_tpl[1].data, xfrm[2][1]) +
                                                         multiply(trace_tpl[2].data, xfrm[2][2]))
    )
    tr_z.stats.channel = chn_2chrs + 'Z'

    # if isinstance(trace_tpl, ida.calibration.qcal_utils.QCalData):
    #     output_trace_tpl = ida.calibration.qcal_utils.QCalData(east=tr_E, north=tr_N, vertical=tr_Z, input=trace_tpl.input.copy())
    # else:
    output_trace_tpl = (tr_2, tr_1, tr_z)

    return output_trace_tpl


def unpack_paz(paz, paz_map):

    nump = len(paz_map[0])
    numz = len(paz_map[1])

    poles = paz._poles[paz_map[0]]
    zeros = paz._zeros[paz_map[1]]

    data = []
    flags = ([],[])

    prev = ''
    for ndx in range(0, nump):
        if isclose(abs(poles[ndx]), 0):  # check for 0+0j
            cur = 'zero'
        elif not isclose(poles[ndx].imag, 0):  # check for complex value
            if prev == 'complex':
                if isclose(poles[ndx].imag, -poles[ndx-1].imag):
                    cur = 'conjugate'
                else:
                    cur = 'complex'
                    data.extend([poles[ndx].real, poles[ndx].imag])
            else:
                cur = 'complex'
                data.extend([poles[ndx].real, poles[ndx].imag])
        else:
            if (prev == 'real') and isclose(poles[ndx].real, poles[ndx-1].real):
                cur = 'real-double'   # if not zer and not complex must be real
            else:
                data.append(poles[ndx].real)
                cur = 'real'

        flags[0].append(cur)
        prev = cur

    prev = ''
    for ndx in range(0, numz):
        if isclose(abs(zeros[ndx]), 0):  # check for 0+0j
            cur = 'zero'
        elif not isclose(zeros[ndx].imag, 0):  # check for complex value
            if prev == 'complex':
                if isclose(zeros[ndx].imag, zeros[ndx-1].imag):
                    cur = 'conjugate'
                else:
                    cur = 'complex'
                    data.extend([zeros[ndx].real, zeros[ndx].imag])
            else:
                cur = 'complex'
                data.extend([zeros[ndx].real, zeros[ndx].imag])
        else:
            if (prev == 'real') and isclose(zeros[ndx].real, zeros[ndx-1].real):
                cur = 'real-double'   # if not zer and not complex must be real
            else:
                data.append(zeros[ndx].real)
                cur = 'real'

        flags[1].append(cur)
        prev = cur

    data.append(1.0)
    return array(data), flags


def pack_paz(data, flags):

    nump = len(flags[0])
    numz = len(flags[1])

    paz_partial = ida.signals.paz.PAZ('vel', 'hz')

    datandx = 0
    for ndx in range(0, nump):
        if flags[0][ndx] == 'zero':
            paz_partial.add_pole(complex(0,0))
        elif flags[0][ndx] == 'complex':
            paz_partial.add_pole(complex(data[datandx], data[datandx+1]))
            datandx += 2
        elif flags[0][ndx] == 'conjugate':
            paz_partial.add_pole(complex(data[datandx-2], -data[datandx-1]))
        elif flags[0][ndx] == 'real':
            paz_partial.add_pole(complex(data[datandx],0))
            datandx += 1
        elif flags[0][ndx] == 'real-double':
            paz_partial.add_pole(complex(data[datandx-1],0))
        else:
            logging.error('Invalid pole type: ' + flags[0][ndx])

    for ndx in range(0, numz):
        if flags[1][ndx] == 'zero':
            paz_partial.add_zero(complex(0,0))
        elif flags[1][ndx] == 'complex':
            paz_partial.add_zero(complex(data[datandx], data[datandx+1]))
            datandx += 2
        elif flags[1][ndx] == 'conjugate':
            paz_partial.add_zero(complex(data[datandx-2], -data[datandx-1]))
        elif flags[1][ndx] == 'real':
            paz_partial.add_zero(complex(data[datandx],0))
            datandx += 1
        elif flags[1][ndx] == 'real-double':
            paz_partial.add_zero(complex(data[datandx-1],0))
        else:
            logging.error('Invalid zero type: ' + flags[1][ndx])

    paz_partial.h0 = data[-1]

#     print('Partial Poles:',paz_partial._poles)
#     print('Partial Zeros:',paz_partial._zeros)
    return paz_partial

