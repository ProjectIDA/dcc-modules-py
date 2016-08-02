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
from functools import reduce

from numpy import pi

from ida.instruments import *

"""Methods and structures supporting construction of CTBTO IMS 2.0 messages"""


# Structure holding individual channel calibration results needed for CALIBRATE_RESULT message type
CTBTChannelResult = namedtuple('CTBTChannelResult', [
    'channel',
    'calib',
    'calper',
    'sample_rate',
    'in_spec',
    'paz',
    'A0'
])


def ims2_calibrate_result_msg_header(cal_timestamp):
    ts_str = cal_timestamp.strftime('%Y/%m/%d %H:%M:%S')

    return "BEGIN IMS2.0\n\n" \
           "MSG_TYPE COMMAND_RESPONSE\n" \
           "MSG_ID <REPLACE WITH VALUE>\n" \
           "REF_ID <REPLACE WITH VALUE>\n" \
           "TIME_STAMP {}\n".format(ts_str)

def ims2_calibrate_result_msg_comp_info(sta, chan, inspec, calib, calper):

    return "\nSTA_LIST {}\n" \
           "CHAN_LIST {}\n" \
           "CALIBRATE_RESULT\n" \
           "IN_SPEC {}\n" \
           "CALIB {:<15.8f}\n" \
           "CALPER {}\n".format(sta,
                                chan,
                                inspec,
                                calib,
                                calper)

def ims2_PAZ2_msg(sta, loc, chan, seis_model, cal_timestamp,
                  opsr, calib, calper, paz, sys_scale_factor):
    # sys_scale_factor is A0 * gnom * gcalib (multiplier
    date_str = cal_timestamp.strftime('%Y/%m/%d')
    time_str = cal_timestamp.strftime('%H:%M')

    msg = "DATA_TYPE RESPONSE\n"

    # CAL2 header record...
    msg = msg + '{:<4} {:<5} {:<3} {:<4} {:<6} {:<15.8e} {:<7.3f} {:<11.5f} {:<10} {:<5}\n'.format(
        'CAL2',
        sta,
        chan,
        loc,
        seis_model[:6],
        calib,
        calper,
        opsr,
        date_str,
        time_str)

    # PAZ2 header records...
    # assumed at 1hz
    # covert to displacement and radians
    msg = msg + '{:<4} {:<2} {:1} {:<15.8e} {:<4} {:<8.3f} {:<3d} {:<3d} {:<25}\n'.format(
        'PAZ2',
        1,
        'V',
        sys_scale_factor * 2 * pi/ 1e9,
        '',
        0,
        paz.num_poles,
        len(paz.zeros(mode='disp', units='hz')),
        '{} {} Disp; Rad'.format(seis_model[:9], chan))

    # PAZ values - in displacement units and radians/sec
    poles = paz.poles(units='rad')
    for pole in poles:
        msg = msg + ' {:15.8e} {:15.8e}\n'.format(pole.real, pole.imag)

    zeros = paz.zeros(mode='disp', units='rad')
    for zero in zeros:
        msg = msg + ' {:15.8e} {:15.8e}\n'.format(zero.real, zero.imag)

    return msg


def calibration_result_msg(sta, loc, seis_model, cal_timestamp, channel_results, dig2_msg_fn=None, fir2_msg_fn=None, msgfn=None):
    """
    High level function used to create IMS2.0 calibration result message. This is a CALIBRATE_RESULT message for
    each channel in the channels_results list with embedded DATA_TYPE RESPONSE information for each channel.

    :param sta: Staiuon code (blank ok for ctbto)
    :type sta: str
    :param loc: Location code
    :type loc: str
    :param seis_model: Seismometer model key
    :type seis_model: str
    :param cal_timestamp: Timestamp of calibration
    :type cal_timestamp: timestamp
    :param channel_results: List of CTBTChannelResult tuples
    :type channel_results: [CTBTOChannelResult...]
    :param dig2_msg_fn: filename of preformatted IMS2.0 digitizer stage 2 response info
    :type dig2_msg_fn: str
    :param fir2_msg_fn: filename of preformatted IMS2.0 digitizer FIR filter stage 3 coefficients
    :type fir2_msg_fn: str
    :param msgfn: Filename of option file to write IMS message text to. (Optional)
    :type msgfn: str
    :return: IMS 2.0 Message text
    :rtype: str
    """

    ts_str = cal_timestamp.strftime('%Y/%m/%d %H:%M:%S')

    dig2_msg = ''
    if dig2_msg_fn:
        with open(dig2_msg_fn, 'rt') as ifl:
            dig2_msg = ifl.read()

    fir2_msg = ''
    if fir2_msg_fn:
        with open(fir2_msg_fn, 'rt') as ifl:
            fir2_msg = ifl.read()

    msg = "BEGIN IMS2.0\n\n" \
          "MSG_TYPE COMMAND_RESPONSE\n" \
          "MSG_ID <REPLACE WITH VALUE>\n" \
          "REF_ID <REPLACE WITH VALUE>\n" \
          "TIME_STAMP {}\n".format(ts_str)

    for ndx, chn_res in enumerate(channel_results):

        assert isinstance(chn_res, CTBTChannelResult)

        msg = msg + \
            "\nSTA_LIST {}\n" \
            "CHAN_LIST {}\n" \
            "CALIBRATE_RESULT\n" \
            "IN_SPEC {}\n" \
            "CALIB {:<15.8f}\n" \
            "CALPER {}\n".format(sta, chn_res.channel, chn_res.in_spec, chn_res.calib, chn_res.calper)

        msg = msg + _sensor_response_msg(sta, loc, seis_model, cal_timestamp, chn_res)

        msg = msg + dig2_msg + fir2_msg

    msg = msg + "\nSTOP\n\n"

    if msgfn:
        with open(msgfn, 'wt') as mfl:
            mfl.write(msg)

    return msg


def _sensor_response_msg(sta, loc, seis_model, cal_timestamp, chn_res):
    """
    Creates DATA_TYPE RESPONSE portion channel response for IMS 2.0 messages.
    Results need to be in radians and diplacement.

    :param sta: Station code
    :type sta: str
    :param loc: Location code
    :type loc: str
    :param seis_model: Seismometer model key
    :type seis_model: str
    :param cal_timestamp: Time of calibration
    :type cal_timestamp: timestamp
    :param chn_res: CTBTOChannelResult tuple with cal results info
    :type chn_res:
    :return: Message text
    :rtype: str
    """

    date_str = cal_timestamp.strftime('%Y/%m/%d')
    time_str = cal_timestamp.strftime('%H:%M')

    msg = "DATA_TYPE RESPONSE\n"

    # CAL2 header record...
    msg = msg + '{:<4} {:<5} {:<3} {:<4} {:<6} {:<15.8e} {:<7.3f} {:<11.5f} {:<10} {:<5}\n'.format(
        'CAL2',
        sta,
        chn_res.channel,
        loc,
        seis_model[:6],
        chn_res.calib,
        chn_res.calper,
        chn_res.sample_rate,
        date_str,
        time_str)

    # PAZ2 header records...
    msg = msg + '{:<4} {:<2} {:1} {:<15.8e} {:<4} {:<8.3f} {:<3d} {:<3d} {:<25}\n'.format(
        'PAZ2',
        1,
        'V',
        chn_res.A0 * INSTRUMENT_NOMINAL_GAINS[seis_model] * 2 * pi/ 1e9,
        '',
        0,
        chn_res.paz.num_poles,
        len(chn_res.paz.zeros(mode='disp', units='hz')),
        '{} {} Disp; Rad'.format(seis_model[:9], chn_res.channel))

    # PAZ values - in displacement units and radians/sec
    poles = chn_res.paz.poles(units='rad')
    for pole in poles:
        msg = msg + ' {:15.8e} {:15.8e}\n'.format(pole.real, pole.imag)

    zeros = chn_res.paz.zeros(mode='disp', units='rad')
    for zero in zeros:
        msg = msg + ' {:15.8e} {:15.8e}\n'.format(zero.real, zero.imag)


    return msg

def ims2_DIG2_msg(stage_seq_num, nom_sens, samp_rate, description):
    txt = 'DIG2 {:<2} {:<15.8e} {:<11.5f} {:<}'.format(stage_seq_num,
                                                       nom_sens,
                                                       samp_rate,
                                                       description[:25])

    return txt.ljust(61, ' ')

def ims2_FIR2_msg(stage_seq_num, gain, decimation, group_correction,
                 symmetry, description, factors):
    fmt = 'FIR2 {:<2} {:<10.2e} {:<4} {:<8.3f} {} {:<4} {:<}'
    txthead = fmt.format(stage_seq_num,
                         gain,
                         decimation,
                         group_correction,
                         symmetry,
                         str(len(factors)),
                         description[:25]).ljust(65, ' ')

    txtdatalist = [txthead]
    row_count = (len(factors) // 5) + 1
    row = 0
    while row < row_count:
        txtrow = reduce(lambda x, y: x + ' {:<15.8e}'.format(y), factors[row*5:row*5+5], '')
        txtrow = txtrow.ljust(80, ' ')
        txtdatalist.append(txtrow)
        row += 1

    return '\n'.join(txtdatalist)

