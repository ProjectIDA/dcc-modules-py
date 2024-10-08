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
from os.path import join

from numpy import complex128, pi, ceil, sin, cos, angle, abs, linspace, multiply, \
    logical_and, less_equal, polyfit, polyval, \
    divide, subtract, concatenate
from numpy.fft import rfft, irfft
from scipy.signal.windows import tukey
import scipy.signal as ss

from scipy.optimize import least_squares
import obspy.signal.filter as osf

import ida.calibration.qcal_utils
from ida.signals.paz import PAZ
import ida.signals.utils
from ida.instruments import CALTYPE_RBLF, CALTYPE_RBHF, TRIAXIAL_SEIS_MODELS, \
        TRIAXIAL_TRANSFORMS, XFRM_TYPE_XYZ2UVW, XFRM_TYPE_UVW2ENZ_ABS

import matplotlib.pyplot as plt

"""utility functions for processing of IDA Random Binary calibration data"""

def compare_component_response(freqs, paz1, paz2, norm_freq=0.05, mode='vel', phase_detrend=False):
    """Compute amp and pha response of paz1 against paz2.

    :param freqs: Frequencies at which to compare responses. phase_detrend assumes freqs go to Nyquist.
    :type freqs: ndarray
    :param paz1: First PAZ repsonse
    :type paz1: PAZ instance
    :param paz2: Second PAZ response
    :type paz2: PAZ instance
    :param norm_freq: Frequency at which to normalize the two responses
    :type norm_freq: float
    :param mode: Units to use when computing responses
    :type mode: str ['disp', 'vel', 'acc']
    :param phase_detrend: Boolean indicating whether to remove linear trend from phase. Assumes nyquist in freqs[-1]
    :type phase_detrend: Bool
    :return: Normalized responses of paz1,
         Normalized responses of paz2,
         Amplitude deviation of paz2 from paz1 in percent,
         Phase deviation of paz2 from paz1 in degrees,
         Absolute max Amplitutde deviation,
         Absolute max Phase deviation,
    :rtype: ndarray, ndarray, ndarray, ndarray, float, float, float
    """

    if (paz1 == None):
        msg = 'You must supply a nominal PAZ to compare.'
        logging.error(msg)
        raise Exception(msg)

    if (paz2 == None):
        msg = 'You must supply a measured component PAZ to compare.'
        logging.error(msg)
        raise Exception(msg)

    resp1 = ida.signals.utils.compute_response(freqs, paz1, mode=mode)
    resp1_norm, _, _ = ida.signals.utils.normalize_response(resp1, freqs, norm_freq)

    resp2 = ida.signals.utils.compute_response(freqs, paz2, mode=mode)
    resp2_norm, _, _ = ida.signals.utils.normalize_response(resp2, freqs, norm_freq)

    if phase_detrend:
        freq_range = less_equal(freqs, 0.9 * freqs[-1])  # 90% of nyquist

        lin_trend_coeff = polyfit(freqs[freq_range], angle(resp1_norm[freq_range]), 1)
        trend_vals = polyval(lin_trend_coeff, freqs)
        amp = abs(resp1_norm)
        pha = angle(resp1_norm) - trend_vals
        resp1_norm = amp * complex128(cos(pha) + 1j*sin(pha))

        lin_trend_coeff = polyfit(freqs[freq_range], angle(resp2_norm[freq_range]), 1)
        trend_vals = polyval(lin_trend_coeff, freqs)
        amp = abs(resp2_norm)
        pha = angle(resp2_norm) - trend_vals
        resp2_norm = amp * complex128(cos(pha) + 1j*sin(pha))


    # calculate percentage deviations
    resp2_a_dev = (divide(abs(resp2_norm[1:]), abs(resp1_norm[1:])) - 1.0) * 100.0
    resp2_p_dev = subtract(angle(resp2_norm[1:])*180/pi, angle(resp1_norm[1:])*180/pi)

    resp2_a_dev_max = abs(resp2_a_dev).max()
    resp2_p_dev_max = abs(resp2_p_dev).max()

    return resp1_norm, resp2_norm, resp2_a_dev, resp2_p_dev, resp2_a_dev_max, resp2_p_dev_max


def analyze_cal_component(fullpaz, lfpertndxs, hfpertndxs, opsr, lftf_f, lf_tf, hftf_f, hf_tf, cal_type):
    """Analyze both high and low frequency calibration component timeseries output with calibration input
    using starting paz fitting_paz.

    Find improved PAZ fit to reduce transfer function based on fitting_paz between
    input/output time series using a least_squares minimization approach.

    """

    def resp_cost(p, paz_partial_flags, freqs, normfreq, tf_target, resp_pert0):
        """

        Args:
            p (list): initial values of poles and zeros in list form
            paz_partial_flags (list): flags identifying which poles and zeros to perturb
            freqs (list): list of frequencies used in frequnecy response calc
            normfreq (float): frequency of normalization
            tf_target (list): list of real followed by imag components of target TF
            resp_pert0 (list(complex)): initial frequency response of paz to be perturbed

        Returns:
            resid (list): array of difference between target and new TF fit with real followed by imag components
        """

        # pack up into PAZ instances
        paz_pert = ida.signals.utils.pack_paz(p, paz_partial_flags)

        # compute perturbed response andnormalize
        resp = ida.signals.utils.compute_response(freqs, paz_pert)
        resp_norm, scale, ndx = ida.signals.utils.normalize_response(resp, freqs, normfreq)

        # calc new TF
        new_tf = divide(resp_norm, resp_pert0)
        # simple dif for residuals array. Assuems tf_target is np.concatenate((tf_target.real, tf_target.imag))
        new_real_imag_tf = concatenate((new_tf.real, new_tf.imag))
        resid = subtract(new_real_imag_tf, tf_target)

        return resid


    # trim freqs and norm freqs
    lflo = 1e-04
    lfhi = 0.45  # 90% of Nyquist of 1hz LF cal signal input
    lf_norm_freq = 0.05

    hflo = 0.45
    hfhi = opsr * 0.45  # 90% of Nyquist of channel operating sample rate
    hf_norm_freq = 1.0

    new_paz = fullpaz.copy()

    if cal_type == CALTYPE_RBLF:
        lf_range = logical_and(lftf_f <= lfhi, lftf_f > lflo)
        lfmeas_f_t = lftf_f[lf_range]
        lfmeas_tf = lf_tf[lf_range]
        lfmeas_tf_norm, _, _ = ida.signals.utils.normalize_response(lfmeas_tf, lfmeas_f_t, lf_norm_freq)

        # Setting paz perturbation map and splitting...
        lf_paz_pert = fullpaz.make_partial(lfpertndxs, lf_norm_freq)

        # computing frequency response of perturbed paz...
        # initial response of paz_pert over freq_band of interest
        lf_resp0 = ida.signals.utils.compute_response(lfmeas_f_t, lf_paz_pert)

        lf_paz_pert_flat, lf_paz_pert_flags = ida.signals.utils.unpack_paz(lf_paz_pert,
                                                                           (list(range(0, lf_paz_pert.num_poles)),
                                                                            list(range(0, lf_paz_pert.num_zeros))))
        # set upper/lower bounds for fitting algorithm
        lf_pazpert_lb = lf_paz_pert_flat - 0.5 * abs(lf_paz_pert_flat)
        lf_pazpert_ub = lf_paz_pert_flat + 0.5 * abs(lf_paz_pert_flat)

        # Fitting new HF response...
        lf_res = least_squares(resp_cost,  # cost function
                                lf_paz_pert_flat,  # initial values
                                bounds=(lf_pazpert_lb,
                                        lf_pazpert_ub),  # lb, ub for each parameter
                                method='trf',
                                jac='3-point',  # I think this matches MATLAB FiniteDifferenceType='central'
                                xtol=1e-6,
                                ftol=1e-4,
                                diff_step=0.001,
                                max_nfev=300,  # max number of function evaluations
                                # list addl args to fittung fun
                                args=(lf_paz_pert_flags,
                                      lfmeas_f_t,
                                      lf_norm_freq,
                                      concatenate((lfmeas_tf_norm.real, lfmeas_tf_norm.imag)),
                                      lf_resp0
                                      ),
                               verbose=0)
        print(lf_res.message)

        # re-combine perturbed poles/zeros into full set of poles and zeros
        new_lf_paz_pert = ida.signals.utils.pack_paz(lf_res.x, lf_paz_pert_flags)
        new_paz.merge_paz_partial(new_lf_paz_pert, lfpertndxs, hf_norm_freq)

    if cal_type == CALTYPE_RBHF:
        hf_range = logical_and(hftf_f <= hfhi, hftf_f > hflo)
        hfmeas_f_t = hftf_f[hf_range]
        hfmeas_tf = hf_tf[hf_range]
        hfmeas_tf_norm, _, _ = ida.signals.utils.normalize_response(hfmeas_tf, hfmeas_f_t, hf_norm_freq)

        # Setting paz perturbation map and splitting...
        hf_paz_pert = fullpaz.make_partial(hfpertndxs, hf_norm_freq)

        # computing frequency response of perturbed paz...
        # initial response of paz_pert over freq_band of interest
        hf_resp0 = ida.signals.utils.compute_response(hfmeas_f_t, hf_paz_pert)

        hf_paz_pert_flat, hf_paz_pert_flags = ida.signals.utils.unpack_paz(hf_paz_pert,
                                                                           (list(range(0, hf_paz_pert.num_poles)),
                                                                            list(range(0, hf_paz_pert.num_zeros))))

        # set upper/lower bounds for fitting algorithm
        hf_pazpert_lb = hf_paz_pert_flat - 0.5 * abs(hf_paz_pert_flat)
        hf_pazpert_ub = hf_paz_pert_flat + 0.5 * abs(hf_paz_pert_flat)

        # Fitting new HF response...
        hf_res = least_squares(resp_cost,
                               hf_paz_pert_flat,
                               bounds=(hf_pazpert_lb,
                                       hf_pazpert_ub),  # lb, ub for each parameter
                               method='trf',
                               jac='3-point',  # I think this matches MATLAB FiniteDifferenceType='central'
                               xtol=1e-6,
                               ftol=1e-4,
                               diff_step=0.001,
                               max_nfev=300,  # max number of function evaluations
                               # list addl args to fittung fun
                               args=(hf_paz_pert_flags,
                                     hfmeas_f_t,
                                     hf_norm_freq,
                                     concatenate((hfmeas_tf_norm.real, hfmeas_tf_norm.imag)),
                                     hf_resp0
                                     ),
                               verbose=0)
        print(hf_res.message)

        # re-combine perturbed poles/zeros into full set of poles and zeros
        new_hf_paz_pert = ida.signals.utils.pack_paz(hf_res.x, hf_paz_pert_flags)
        new_paz.merge_paz_partial(new_hf_paz_pert, hfpertndxs, hf_norm_freq)

    return new_paz


def prepare_cal_data(lfpath, lffile, hfpath, hffile, sensor, comp, fullpaz, opsr):
    """Prepares cal input and measured (output) timeseries for analysis.


    Each components is:
        - trimmed to remove settling and trailing portions of time series
        - corrected for polarity
        - If triaxial seismometer, UVW components are transformed to get horizontal signals
        - Input timeseries are:
            - tapered
            - xformed to frequency space
            - convolved with paz response
            - inverse xformed
            - taper portions trimmed off
            - normed and de-meaned
        - each observed compopnent output time series is then:
            - trimmed just as input timeseries
            - normed and de-meaned

    The algorithm follows the Matlab scripts previously used for calibration processing with the following exceptions:
        1) All 3 components are handled at once
        2) Components for triaxial seismometer are transformed to UVW and then to abs() values for each XYZ component

    Args:
        lfpath (str): path to LF MS and LOG files
        lffile (str): name (stem) of LF MS and LOG files
        hfpath (str): path to LF MS and LOG files
        hffile (str): name (stem) of LF MS and LOG files
        sensor (str): sensor key string (should be one of SEISMOMETER_MODELS)
        comp (str): Component 'Z', '1', or '2'
        fullpaz (PAZ): Instance of PAZ that should be convolved with input signal

    Returns:
        samp_rate_lf (float or int): Sample rate of LF calibration data
        start_time_lf (obspy.UTCDateTime): Start time of LF input calibration stream
        lf_inp_wth_resp (numpy.ndarray): LF input timeseries convolved with fullpaz response
        lf_out (numpy-ndarray): LF measured output timeseries
        samp_rate_hf (): Sample rate of HF calibration data
        start_time_hf (): Start time of HF input calibration stream
        hf_inp_wth_resp (): HF input timeseries convolved with fullpaz response
        hf_out (): HF measured output timeseries
        freqs_lf (numpy.ndarray): LF Frequency array
        freqs_hf (numpy.ndarray): HF Frequency array
        lf_snr (float): LF Signal-to-noise of input+full and measured timeseries. Defined as 1/stdev(lfinput-meas)
        hf_snr (float): LF Signal-to-noise of input+full and measured timeseries. Defined as 1/stdev(lfinput-meas)


    """

    if lfpath:
        ms_fpath_lf = join(lfpath, lffile + '.ms')
        log_fpath_lf = join(lfpath, lffile + '.log')

        strm_lf, log_lf = ida.calibration.qcal_utils.read_qcal_files(ms_fpath_lf, log_fpath_lf)

        # trim settling and trailing times from traces
        ida.signals.utils.trim_stream(strm_lf, left=log_lf['settling_time'], right=log_lf['trailing_time'])

        # reverse polarity for those seismometers that are funky...
        ida.signals.utils.check_and_fix_polarities(strm_lf, sensor.upper())

        # split cal data into enzi, 21zi (one, two, vertical, input) tuple
        cal_lf = ida.calibration.qcal_utils.split_qcal_traces(strm_lf)

        # for triaxial sensors, rotate to UVW, then back the 12Z with
        # horizontal amplitudes combined using absolute magnitudes
        if sensor.upper() in TRIAXIAL_SEIS_MODELS:
            cal_lf_tpl = triaxial_horizontal_magnitudes(cal_lf, sensor.upper())
        else:
            cal_lf_tpl = cal_lf

        samp_rate_lf = cal_lf_tpl.input.stats.sampling_rate
        start_time_lf = cal_lf_tpl.input.stats.starttime
        npts_lf = cal_lf_tpl.input.stats.npts

        fraction = 0.1  # each side Creating tapers...
        taper_lf = tukey(npts_lf, alpha=fraction * 2, sym=True)
        taper_bin_cnt_lf = int(ceil(npts_lf * fraction))

        # make fitting paz response objs and reset gain to 1.0 for fitting
        lf_fit_paz = fullpaz.make_partial2(1.0, partial_mode=PAZ.PARTIAL_FITTING_LF)
        lf_fit_paz.h0 = 1.0

        freqs_lf = linspace(0, samp_rate_lf/2, npts_lf//2 + 1)  # count is to match behavior of np.fft.rfft below
        # compute frequency response in ACC units
        resp_tmp_lf = ida.signals.utils.compute_response(freqs_lf, lf_fit_paz, mode='acc')
        # normalize ferquency response
        resp_lf, _, _ = ida.signals.utils.normalize_response(resp_tmp_lf, freqs_lf, 0.05)

        # prep output channels
        if comp == 'Z':
            lf_out = cal_lf_tpl.vertical.data.copy()
        elif comp == '1':
            lf_out = cal_lf_tpl.one.data.copy()
        elif comp == '2':
            lf_out = cal_lf_tpl.two.data.copy()
        else:
            raise ValueError('Invalid component: ' + comp)

        # trim of amount used to input time-series
        lf_out = lf_out[taper_bin_cnt_lf:-taper_bin_cnt_lf]

        # normalize and de-mean
        lf_out.__itruediv__(lf_out.std())
        lf_out.__isub__(lf_out.mean())

        # Multiply input with taper and with starting frequency response...
        input_fft           = rfft(multiply(cal_lf_tpl.input.data[:npts_lf], taper_lf))
        inp_freqs_cnv_resp  = multiply(input_fft, resp_lf)
        lf_inp_wth_resp     = irfft(inp_freqs_cnv_resp, npts_lf)

        # remove taper region, normalize and de-mean
        lf_inp_wth_resp     = lf_inp_wth_resp[taper_bin_cnt_lf:-taper_bin_cnt_lf]
        lf_inp_wth_resp.__itruediv__(lf_inp_wth_resp.std())
        lf_inp_wth_resp.__isub__(lf_inp_wth_resp.mean())

        # compute an 'snr' to get sense of quality of 'starting' response.
        # imperically, values over ~75 will lead to a successful 'fitting' process
        lf_snr = 1/subtract(lf_inp_wth_resp, lf_out).std()
        print('lf snr:', lf_snr)

    else:
        samp_rate_lf, start_time_lf, lf_inp_wth_resp, lf_out, freqs_lf, lf_snr = None, None, None, None, None, None

    if hfpath:
        ms_fpath_hf = join(hfpath, hffile + '.ms')
        log_fpath_hf = join(hfpath, hffile + '.log')

        strm_hf, log_hf = ida.calibration.qcal_utils.read_qcal_files(ms_fpath_hf, log_fpath_hf)

        # trim settling and trailing times from traces
        ida.signals.utils.trim_stream(strm_hf, left=log_hf['settling_time'], right=log_hf['trailing_time'])
        # for those seismometers that are funky...
        ida.signals.utils.check_and_fix_polarities(strm_hf, sensor.upper())

        # split cal data into enzi, 21zi (one, two, vertical, input)
        cal_hf = ida.calibration.qcal_utils.split_qcal_traces(strm_hf)

        if sensor.upper() in TRIAXIAL_SEIS_MODELS:
            cal_hf_tpl = triaxial_horizontal_magnitudes(cal_hf, sensor.upper())
        else:
            cal_hf_tpl = cal_hf

        samp_rate_hf = cal_hf_tpl.input.stats.sampling_rate
        start_time_hf = cal_hf_tpl.input.stats.starttime
        npts_hf = cal_hf_tpl.input.stats.npts

        fraction = 0.1  # each side Creating tapers...
        taper_hf = tukey(npts_hf, alpha=fraction * 2, sym=True)
        taper_bin_cnt_hf = int(ceil(npts_hf * fraction))

        # make fitting paz response objs and reset gain to 1.0 for fitting
        hf_fit_paz = fullpaz.make_partial2(1.0, partial_mode=PAZ.PARTIAL_FITTING_HF)
        hf_fit_paz.h0 = 1.0

        freqs_hf = linspace(0, samp_rate_hf/2, npts_hf//2 + 1)  # count is to match behavior of np.fft.rfft below
        resp_tmp_hf = ida.signals.utils.compute_response(freqs_hf, hf_fit_paz, mode='acc')
        resp_hf, _, _ = ida.signals.utils.normalize_response(resp_tmp_hf, freqs_hf, 0.05)

        # prep output channels
        if comp == 'Z':
            hf_out = cal_hf_tpl.vertical.data.copy()
        elif comp == '1':
            hf_out = cal_hf_tpl.one.data.copy()
        elif comp == '2':
            hf_out = cal_hf_tpl.two.data.copy()
        else:
            raise ValueError('Invalid component: ' + comp)

        # new filtering
        hf_out.__isub__(hf_out.mean())
        hf_out = multiply(hf_out, taper_hf)
        hf_out = osf.bandpass(hf_out, 0.45, opsr/2, samp_rate_hf, zerophase=True)
        hf_out = hf_out[taper_bin_cnt_hf:-taper_bin_cnt_hf]

        hf_out.__itruediv__(hf_out.std())
        hf_out = ss.detrend(hf_out, type='constant')

        cal_hf_tpl.input.data = ss.detrend(cal_hf_tpl.input.data, type='constant')

        input_tapered       = multiply(cal_hf_tpl.input.data, taper_hf)
        input_fft           = rfft(input_tapered)
        inp_freqs_cnv_resp  = multiply(input_fft, resp_hf)
        hf_inp_wth_resp     = irfft(inp_freqs_cnv_resp, npts_hf)

        # apply bandpass filter
        hf_inp_wth_resp.__isub__(hf_inp_wth_resp.mean())
        hf_inp_wth_resp = osf.bandpass(hf_inp_wth_resp, 0.45, opsr/2, samp_rate_hf, zerophase=True)

        # remove taper region, normalize and de-mean
        hf_inp_wth_resp     = hf_inp_wth_resp[taper_bin_cnt_hf:-taper_bin_cnt_hf]
        hf_inp_wth_resp.__itruediv__(hf_inp_wth_resp.std())
        hf_inp_wth_resp = ss.detrend(hf_inp_wth_resp, type='constant')

        # plt.figure(111, figsize=(12, 8))
        # plt.plot(hf_inp_wth_resp[50000:52000], 'g')
        # plt.plot(hf_out[50000:52000], 'b')
        # plt.plot(hf_inp_wth_resp[50000:52000]-hf_out[50000:52000], 'r')
        # plt.show()

        hf_snr = 1/subtract(hf_inp_wth_resp, hf_out).std()
        print('hf snr:', hf_snr)

    else:
        samp_rate_hf, start_time_hf, hf_inp_wth_resp, hf_out, freqs_hf, hf_snr = None, None, None, None, None, None

    return samp_rate_lf, start_time_lf, lf_inp_wth_resp, lf_out, \
           samp_rate_hf, start_time_hf, hf_inp_wth_resp, hf_out, \
           freqs_lf, freqs_hf, lf_snr, hf_snr


def triaxial_horizontal_magnitudes(cal_tpl, seis_model):
    """Transform XYZ components to UVW, then to absolute values for ENZ.
    This obtains the true absolute value time series for each horizontal component

    :param cal_tpl: QCal Tuple with component output teimseries plus input
    :type cal_tpl: ida.calibration.qcal_utils.QCalData
    :param seis_model: Seismometer model key
    :type seis_model: str
    :return: New, transformed QCal tuple
    :rtype: ida.calibration.qcal_utils.QCalData
    """

    if seis_model in TRIAXIAL_SEIS_MODELS:
        uvw = ida.signals.utils.channel_xform((cal_tpl.two,
                                               cal_tpl.one,
                                               cal_tpl.vertical),
                                              TRIAXIAL_TRANSFORMS[seis_model][XFRM_TYPE_XYZ2UVW])
        enz = ida.signals.utils.channel_xform(uvw, TRIAXIAL_TRANSFORMS[seis_model][XFRM_TYPE_UVW2ENZ_ABS])
        new_cal_tpl = ida.calibration.qcal_utils.QCalData(two=enz[0],
                                                          one=enz[1],
                                                          vertical=enz[2],
                                                          input=cal_tpl.input)

    else:
        msg = 'The seismomemter model {} is unsupported or not a triaxial instrument.'.format(seis_model)
        logging.error(msg)
        raise Exception(msg)

    return new_cal_tpl


