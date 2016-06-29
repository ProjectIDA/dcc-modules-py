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
import os.path
import sys
from numpy import pi, linspace, abs
import ida.calibration.plots
from ida.calibration.process import analyze_cal_component, \
    compare_component_response, \
    prepare_cal_data, \
    nominal_sys_sens_1hz
import ida.signals.paz
from ida.signals.utils import compute_response
import ida.ctbto.messages
from ida.instruments import ComponentsTpl, SEISMOMETER_RESPONSES

"""Methods for performing CTBTO/Sandia specific calibration analysis and data processing"""

def process_qcal_data(sta, chancodes, loc, data_dir, lf_fnames, hf_fnames, seis_model, full_paz_fn):
    """Main method for analyzing random binary calibration data for Sandia/CTBTO STS2.5 (w/Q330HR digi) sensors.
    It processes both low and high frequency frequency random binary data sets. Supplied timeseries files must
    contain all 3 components plus input cal signal.

    :param sta: Station code
    :type sta: str
    :param chancodes: ComponentsTpl tuple containing chan codes for each component.
    :type chancodes: ida.instruments.ComponentsTpl
    :param loc: Location code, 2 decimal digits or spaces
    :type loc: str
    :param data_dir: Directory where miniseed and qcal log files are located
    :type data_dir: str
    :param lf_fnames: For Low frequency: (miniseed_filename, qcal_log_filemane)
    :type lf_fnames: (str, str)
    :param hf_fnames: For High frequency: (miniseed_filename, qcal_log_filemane)
    :type hf_fnames: (str, str)
    :param seis_model: Seismometer model key
    :type seis_model: str
    :param full_paz_fn: PAZ FULL response filename or Tuple of filenames for each component
    :type full_paz_fn: str or ComponentsTpl
    :return: Results: (IMS2.0_msg_filename, amp_plot_filename, pha_plot_filename)
    :rtype: (str, str, str)
    """

    if not isinstance(full_paz_fn, str) and not isinstance(full_paz_fn, ComponentsTpl):
        msg = 'full_paz_fn must be a response file name or ComponentsTpl with filenames of responses for each component'
        raise TypeError(msg)

    if isinstance(full_paz_fn, str):
        full_paz_fn_tpl = ComponentsTpl(vertical=full_paz_fn, north=full_paz_fn, east=full_paz_fn)
    else:
        full_paz_fn_tpl = full_paz_fn

    logging.debug('Reading FULL response files: {}'.format(str(full_paz_fn_tpl)))

    # read FULL resp paz
    full_paz_v = ida.signals.paz.PAZ('vel', 'hz', pzfilename=full_paz_fn_tpl.vertical, fileformat='ida')
    full_paz_n = ida.signals.paz.PAZ('vel', 'hz', pzfilename=full_paz_fn_tpl.north, fileformat='ida')
    full_paz_e = ida.signals.paz.PAZ('vel', 'hz', pzfilename=full_paz_fn_tpl.east, fileformat='ida')

    full_paz_tpl = ComponentsTpl(vertical=full_paz_v, north=full_paz_n, east=full_paz_e)

    # perturbed paz must not be in both lf anf hf sets
    lf_fit_paz_map = (SEISMOMETER_RESPONSES[seis_model]['fit']['lf_poles'],
                       SEISMOMETER_RESPONSES[seis_model]['fit']['lf_zeros'])

    hf_fit_paz_map = (SEISMOMETER_RESPONSES[seis_model]['fit']['hf_poles'],
                       SEISMOMETER_RESPONSES[seis_model]['fit']['hf_zeros'])

    # make fitting paz response objs and reset gain to 1.0 for fitting
    lf_fit_paz_v = full_paz_tpl.vertical.make_partial(lf_fit_paz_map, 1.0)
    lf_fit_paz_v.h0 = 1.0
    lf_fit_paz_n = full_paz_tpl.north.make_partial(lf_fit_paz_map, 1.0)
    lf_fit_paz_n.h0 = 1.0
    lf_fit_paz_e = full_paz_tpl.east.make_partial(lf_fit_paz_map, 1.0)
    lf_fit_paz_e.h0 = 1.0

    lf_fit_paz_tpl = ComponentsTpl(vertical=lf_fit_paz_v, north=lf_fit_paz_n, east=lf_fit_paz_e)

    #ToDo: refactor to only use one paz tpl since LF and HF come from same response files.
    hf_fit_paz_v = full_paz_tpl.vertical.make_partial(hf_fit_paz_map, 1.0)
    hf_fit_paz_v.h0 = 1.0
    hf_fit_paz_n = full_paz_tpl.north.make_partial(hf_fit_paz_map, 1.0)
    hf_fit_paz_n.h0 = 1.0
    hf_fit_paz_e = full_paz_tpl.east.make_partial(hf_fit_paz_map, 1.0)
    hf_fit_paz_e.h0 = 1.0

    hf_fit_paz_tpl = ComponentsTpl(vertical=hf_fit_paz_v, north=hf_fit_paz_n, east=hf_fit_paz_e)

    logging.debug('Reading responses complete.')

    logging.debug('Preparing cal data for coherence analysis...')
    # read, trim, transpose, invert and convolve input with response

    samp_rate_lf, lf_start_time, lfinput, lfmeas, \
    samp_rate_hf, hf_start_time, hfinput, hfmeas, \
    _, freqs_lf, \
    _, freqs_hf = prepare_cal_data(os.path.abspath(data_dir),
                                   lf_fnames, hf_fnames,
                                   seis_model,
                                   lf_fit_paz_tpl, hf_fit_paz_tpl)

    logging.debug('Preparing cal data for coherence analysis complete.')

    # perturbed paz must not be in both lf anf hf sets
    lf_paz_pert_map = (SEISMOMETER_RESPONSES[seis_model]['perturb']['lf_poles'],
                       SEISMOMETER_RESPONSES[seis_model]['perturb']['lf_zeros'])
    hf_paz_pert_map = (SEISMOMETER_RESPONSES[seis_model]['perturb']['hf_poles'],
                       SEISMOMETER_RESPONSES[seis_model]['perturb']['hf_zeros'])

    logging.debug('Analyzing cal data and calculating new VERTICAL response...')

    operating_sample_rate = 40.0
    new_full_paz_v = analyze_cal_component(full_paz_tpl.vertical,
                                          lf_paz_pert_map, hf_paz_pert_map,
                                          samp_rate_lf, samp_rate_hf, operating_sample_rate,
                                          lfinput.vertical, hfinput.vertical,
                                          lfmeas.vertical, hfmeas.vertical)

    logging.debug('Analyzing cal data and calculating new VERTICAL response complete.')
    logging.debug('Analyzing cal data and calculating new NORTH response...')

    new_full_paz_n = analyze_cal_component(full_paz_tpl.north,
                                          lf_paz_pert_map, hf_paz_pert_map,
                                          samp_rate_lf, samp_rate_hf, operating_sample_rate,
                                          lfinput.north, hfinput.north,
                                          lfmeas.north, hfmeas.north)

    logging.debug('Analyzing cal data and calculating new NORTH response complete.')
    logging.debug('Analyzing cal data and calculating new EAST response...')

    new_full_paz_e = analyze_cal_component(full_paz_tpl.east,
                                          lf_paz_pert_map, hf_paz_pert_map,
                                          samp_rate_lf, samp_rate_hf, operating_sample_rate,
                                          lfinput.east, hfinput.east,
                                          lfmeas.east, hfmeas.east)

    logging.debug('Analyzing cal data and calculating new EAST response complete.')

    # lets find "nice" number for resp length: i.e. a multiple of 20 * sample rate. So 0.05 and 1 hz in freqs exactly.
    resp_len = ((hfinput.vertical.size * 2) // (20 * operating_sample_rate)) * (20 * operating_sample_rate)
    freqs, _ = linspace(0, operating_sample_rate/2, resp_len+1, retstep=True)  # must start with 0hz
    # resp_len = ((hfinput.size * 2) // (20 * 20)) * (20 * 20)
    # freqs = linspace(0, 20, resp_len)  # must start with 0hz

    logging.debug('Comparing VERTICAL response with FULL RESP normed at 0.05...')

    bl_resp, v_resp, v_adev, v_pdev, v_adev_max, v_pdev_max = \
        compare_component_response(freqs, full_paz_tpl.vertical, new_full_paz_v,
                                   norm_freq=0.05, mode='vel', phase_detrend=True)

    logging.debug('Comparing VERTICAL response with FULL RESP complete.')
    logging.debug('Comparing NORTH response with FULL RESP normed at 0.05...')

    bl_resp, n_resp, n_adev, n_pdev, n_adev_max, n_pdev_max = \
        compare_component_response(freqs, full_paz_tpl.north, new_full_paz_n,
                                   norm_freq=0.05, mode='vel', phase_detrend=True)

    logging.debug('Comparing NORTH response with FULL RESP complete.')
    logging.debug('Comparing EAST response with FULL RESP normed at 0.05...')

    bl_resp, e_resp, e_adev, e_pdev, e_adev_max, e_pdev_max = \
        compare_component_response(freqs, full_paz_tpl.east, new_full_paz_e,
                                   norm_freq=0.05, mode='vel', phase_detrend=True)

    logging.debug('Comparing EAST response nom_paz FULL RESP complete.')
    logging.debug('Finding sensitivities @ 1hz...')

    # sens X_resp responses normed at 0.05hz, need to get resp at 1hz
    # find 1hz bin
    #todo: find the better way
    bin_1hz_ndx = min([freq[0] for freq in enumerate(freqs) if freq[1] >= 1.0])
    n_sys_sens = nominal_sys_sens_1hz(n_resp[bin_1hz_ndx], seis_model)
    e_sys_sens = nominal_sys_sens_1hz(e_resp[bin_1hz_ndx], seis_model)
    v_sys_sens = nominal_sys_sens_1hz(v_resp[bin_1hz_ndx], seis_model)

    # convert cts/(m/s) to cts/m to nm/ct
    n_sys_sens = 1 / (n_sys_sens * (2 * pi) / 1e9)
    e_sys_sens = 1 / (e_sys_sens * (2 * pi) / 1e9)
    v_sys_sens = 1 / (v_sys_sens * (2 * pi) / 1e9)

    logging.debug('Finding sensitivities @ 1hz... complete.')

    logging.debug('Finding A0 at 1hz...')
    resp_1hz = compute_response(1.0, new_full_paz_v)
    v_a0 = 1/abs(resp_1hz[0])
    resp_1hz = compute_response(1.0, new_full_paz_n)
    n_a0 = 1/abs(resp_1hz[0])
    resp_1hz = compute_response(1.0, new_full_paz_e)
    e_a0 = 1/abs(resp_1hz[0])
    logging.debug('Finding A0 at 1hz complete.')

    # logging.debug('Reprocessing with new response for snr analysis...')
    # # read, trim, transpose, invert and convolve input with response
    # samp_rate_lf, lf_start_time, lfinput, lfmeas, \
    # samp_rate_hf, hf_start_time, hfinput, hfmeas, \
    # _, freqs_lf, \
    # _, freqs_hf = prepare_cal_data(os.path.abspath(data_dir), lf_fnames, hf_fnames, seis_model, start_paz)
    # logging.debug('Reprocessing with new response for snr analysis complete.')

    logging.debug('Computing IN_SPEC status...')

    if ((n_adev_max <= 5.0) and (n_pdev_max <= 5)):
        north_inspec = 'YES'
    else:
        north_inspec = 'NO'
    if ((e_adev_max <= 5.0) and (e_pdev_max <= 5)):
        east_inspec = 'YES'
    else:
        east_inspec = 'NO'
    if ((v_adev_max <= 5.0) and (v_pdev_max <= 5)):
        vert_inspec = 'YES'
    else:
        vert_inspec = 'NO'

    logging.debug('Computing IN_SPEC status... complete')
    logging.debug('Creating result tuples...')

    vert_chan_result = ida.ctbto.messages.CTBTChannelResult(channel=chancodes.vertical, calib=v_sys_sens, calper=1.0,
                                                        sample_rate=40.0,
                                                        in_spec=vert_inspec, paz=new_full_paz_v, A0=v_a0)
    north_chan_result = ida.ctbto.messages.CTBTChannelResult(channel=chancodes.north, calib=n_sys_sens, calper=1.0,
                                                             sample_rate=40.0,
                                                             in_spec=north_inspec, paz=new_full_paz_n, A0=n_a0)
    east_chan_result = ida.ctbto.messages.CTBTChannelResult(channel=chancodes.east, calib=e_sys_sens, calper=1.0,
                                                        sample_rate=40.0,
                                                        in_spec=east_inspec, paz=new_full_paz_e, A0=e_a0)
    logging.debug('Creating result tuples... complete.')

    report_files_basename = os.path.splitext(lf_fnames[0])[0].replace('rblf-', '')

    logging.info('Generating Response plots...')
    amp_fn = os.path.join(data_dir, report_files_basename + '_AMP_Resp.png')
    pha_fn = os.path.join(data_dir, report_files_basename + '_PHA_Resp.png')

    ida.calibration.plots.save_response_comparison_plots(sta, chancodes, loc,
                                                         amp_fn, pha_fn,
                                                         seis_model,
                                                         lf_start_time, operating_sample_rate,
                                                         resp_len+1, 1.0,
                                                         bl_resp,
                                                         n_resp, n_adev, n_pdev,
                                                         e_resp, e_adev, e_pdev,
                                                         v_resp, v_adev, v_pdev)
    logging.info('Generating Response plots... complete.')

    if getattr(sys, 'frozen', False):
        bundle_dir = sys._MEIPASS
        dig2_msg_fn = os.path.abspath(os.path.join(bundle_dir, 'IDA', 'data', 'Q330_DIG2.txt'))
        fir2_msg_fn = os.path.abspath(os.path.join(bundle_dir, 'IDA', 'data', 'Q330_FIR2.txt'))
    else:
        # bundle_dir = os.path.dirname(os.path.abspath(__file__))
        dig2_msg_fn = os.path.abspath(os.path.join('.', 'data', 'Q330_DIG2.txt'))
        fir2_msg_fn = os.path.abspath(os.path.join('.', 'data', 'Q330_FIR2.txt'))

    ims_calres_txt_fn = os.path.join(data_dir, report_files_basename + '_CALRESULT.msgtxt')

    logging.info('Generating IMS 2.0 messages...')

    calres_msg = ida.ctbto.messages.calibration_result_msg(
        sta, loc, seis_model, hf_start_time,
        [north_chan_result, east_chan_result, vert_chan_result],
        dig2_msg_fn, fir2_msg_fn,
        msgfn=ims_calres_txt_fn)

    logging.debug(calres_msg)
    logging.info('Generating IMS 2.0 messages... complete.')

    return ims_calres_txt_fn, amp_fn, pha_fn



