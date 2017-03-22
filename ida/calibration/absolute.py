#######################################################################################
# Copyright (C) 2016  Regents of the University of California
#
# This is free software: you can redistribute it and/or modify it under the terms of
# the GNU General Public License (GNU GPL) as published by the Free Software Foundation,
# either version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# A copy of the GNU General Public License can be found in LICENSE.TXT in the root of
# the source code repository. It can also be found at http://www.gnu.org/licenses/.
#
# NOTES: Per GNU GPLv3 terms:
#   * This notice must be kept in this source file
#   * Changes to the source must be clearly noted with date & time of change
#
# If you use this software in a product, an explicit acknowledgment in the product
# documentation of the contribution by Project IDA, Institute of Geophysics and
# Planetary Physics, UCSD would be appreciated but is not required.
#######################################################################################
from datetime import datetime
import os.path
# from os import remove
from pathlib import Path # , PurePath
import yaml
from collections import namedtuple
import logging

import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
plt.ion()

#from fabulous.color import red, green, blue, bold
from obspy import read, UTCDateTime, Stream  # Trace
from obspy.signal.invsim import evalresp
import obspy.signal.filter as osf
from scipy.signal import tukey
# from numpy import array, cos, sin, subtract, angle, complex128, absolute, less
from numpy import array, ones, pi, arctan2, float64, multiply, divide, log10
from numpy import std, sqrt, dot, insert
from numpy.fft import rfft, irfft
import numpy.linalg as la
import scipy.signal as ss
from fabulous.color import red

from ida import IDA_PKG_VERSION_HASH_STR, IDA_PKG_VERSION_DATETIME
from ida.utils import i10get, pimseed
from ida.signals.utils import time_offset, taper_high_freq_resp, \
        dynlimit_resp_min
from ida.calibration.shaketable import rename_chan

class APSurveyComponentResult(object):

    def __init__(self, comp):
        self.component = comp
        self.seg_results = []
        self._amp_mean = None
        self._amp_std = None
        self._ang_mean = None
        self._ang_std = None
        self._ang_resid = None
        self._lrms_mean = None
        self._var_mean = None
        self._recalc_needed = False
        self._usable_count = 0

    def add_segment(self, apssegres):
        self.seg_results.append(apssegres)
        self._recalc_needed = True
        if apssegres.can_use:
            self._usable_count += 1

    def _recalc(self):
        if self.usable_count > 0:
            self._amp_mean = array([seg.amp for seg in self.seg_results if seg.can_use]).mean()
            self._amp_std = array([seg.amp for seg in self.seg_results if seg.can_use]).std()
            self._ang_mean = array([seg.ang for seg in self.seg_results if seg.can_use]).mean()
            self._ang_std = array([seg.ang for seg in self.seg_results if seg.can_use]).std()
            self._lrms_mean = array([seg.lrms for seg in self.seg_results if seg.can_use]).mean()
            self._var_mean = array([seg.var for seg in self.seg_results if seg.can_use]).mean()
        else:
            self._amp_mean = None
            self._amp_std = None
            self._ang_mean = None
            self._ang_std = None
            self._ang_resid = None
            self._lrms_mean = None
            self._var_mean = None
        self._recalc_needed = False

    @property
    def amp_mean(self):
        if self._recalc_needed:
            self._recalc()
        return self._amp_mean

    @property
    def amp_std(self):
        if self._recalc_needed:
            self._recalc()
        return self._amp_std

    @property
    def ang_mean(self):
        if self._recalc_needed:
            self._recalc()
        return self._ang_mean

    @property
    def ang_std(self):
        if self._recalc_needed:
            self._recalc()
        return self._ang_std

    @property
    def ang_resid(self):
        if self._recalc_needed:
            self._recalc()
        return self._ang_resid

    @property
    def lrms_mean(self):
        if self._recalc_needed:
            self._recalc()
        return self._lrms_mean

    @property
    def var_mean(self):
        if self._recalc_needed:
            self._recalc()
        return self._var_mean

    @property
    def usable_count(self):
        return self._usable_count

    @property
    def total_count(self):
        return len(self.seg_results)


class APSurveySegmentResult(object):

    def __init__(self, start_t_utc, ang, ang_resid, amp, lrms, var, coh, use_segment):
        self._start_t_utc = start_t_utc
        self._ang = ang
        self._ang_resid = ang_resid
        self._amp = amp
        self._lrms = lrms
        self._var = var
        self._coh = coh
        self._use_segment = use_segment

    @property
    def start_utc(self):
        return self._start_t_utc

    @property
    def start_epoch(self):
        return self._start_t_utc.timestamp

    @property
    def ang(self):
        return self._ang

    @property
    def ang_resid(self):
        return self._ang_resid

    @property
    def amp(self):
        return self._amp

    @property
    def lrms(self):
        return self._lrms

    @property
    def var(self):
        return self._var

    @property
    def coh(self):
        return self._coh

    @property
    def can_use(self):
        return self._use_segment

class APSurvey(object):

    ChanTpl = namedtuple('ChanTuple', 'z n e')

    def __init__(self, fn, debug=False):

        self.analdate = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        self.ida_cal_raw_dir = os.environ.get('IDA_CAL_RAW_DIR', '')
        self.debug = debug
        self.waveform_files = []

        self.ok = True

        self.logger = logging.getLogger(__name__)

        # in-memory streams and source miniseed data
        self.streams = {
            'abs': {
                'ref': Stream(),
                'pri': Stream(),
                'sec': Stream()
            },
            'azi': {
                'ref': Stream(),
                'pri': Stream(),
                'sec': Stream()
            }
        }
        self.msfiles = {
            'abs': {
                'ref': '',
                'pri': '',
                'sec': ''
            },
            'azi': {
                'ref': '',
                'pri': '',
                'sec': ''
            }
        }
        self.trtpls= {
            'abs': {
                'ref': None,
                'pri': None,
                'sec': None
            },
            'azi': {
                'ref': None,
                'pri': None,
                'sec': None
            }
        }
        # in-memory responses for 2 or 3 sensors
        # note that '_v_' responses are because responses have to account for
        # potentially different sample-rates  & nyquist freqs
        self.responses = {
            # responses for 'ref' sensor, adjusted for dif sampling rates of sta sensors
            'ref': {
                'pri': None,
                'sec': None,
            },
            'sec': {
                'pri': None, # resp of SEC sta sensor adjusted for dif sampling rate of PRI
                'sec': None  # regular response for SEC station sensor
            },
            # now do regular responses for PRI station sensors
            'pri': {
                'pri': None
            },
        }
        # time correction flags so only do it once
        self.ref_clock_adjustment = 0.0
        self.ref_clock_adjusted = False
        self.ref_clock_adjustment_sensor = 'n/a'

        self.ref_data_unusable = False  # assume will be ok, set True if probs on first read
        self.pri_data_unusable = False  # assume will be ok, set True if probs on first read
        self.sec_data_unusable = False  # assume will be ok, set True if probs on first read
        self.results = None

        with open(fn, 'rt') as cfl:
            config_txt = cfl.read()
        try:
            self._config = yaml.load(config_txt)
            self._config_dir = Path(fn).parent
        except:
            self._config = {}
            print(red('Error parsing YAML config file: ' + fn))
            self.ok = False
        else:
            self._process_config()

    def _process_config(self):

        if 'correlation_segment_size_secs' not in self._config:
            self.ok = False
            self.logmsg(logging.ERROR, 'Missing key in configuration file: ' + 'correlation_segment_size_secs')
        if 'segment_size_secs' not in self._config:
            self.ok = False
            self.logmsg(logging.ERROR, 'Missing key in configuration file: ' + 'segment_size_secs')
        if 'segment_size_trim_secs' not in self._config:
            self.ok = False
            self.logmsg(logging.ERROR, 'Missing key in configuration file: ' + 'segment_size_trim_secs')
        if 'coherence_cutoff' not in self._config:
            self.ok = False
            self.logmsg(logging.ERROR, 'Missing key in configuration file: ' + 'coherence_cutoff')
        if 'analysis_sample_rate_hz' not in self._config:
            self.ok = False
            self.logmsg(logging.ERROR, 'Missing key in configuration file: ' + 'analysis_sample_rate_hz')
        if 'analysis_bandpass_hz' not in self._config:
            self.ok = False
            self.logmsg(logging.ERROR, 'Missing key in configuration file: ' + 'analysis_bandpass_hz')

        if 'station' not in self._config:
            self.ok = False
            self.logmsg(logging.ERROR, 'Missing key in configuration file: ' + 'station')
        if 'pri_sensor_installed' not in self._config:
            self.ok = False
            self.logmsg(logging.ERROR, 'Missing key in configuration file: ' + 'pri_sensor_installed')
        elif self._config['pri_sensor_installed']:
            if 'pri_sensor_chans' not in self._config:
                self.ok = False
                self.logmsg(logging.ERROR, 'Missing key in configuration file: ' + 'pri_sensor_chans')
            if 'pri_sensor_loc' not in self._config:
                self.ok = False
                self.logmsg(logging.ERROR, 'Missing key in configuration file: ' + 'pri_sensor_loc')
        if 'sec_sensor_installed' not in self._config:
            self.ok = False
            self.logmsg(logging.ERROR, 'Missing key in configuration file: ' + 'sec_sensor_installed')
        elif self._config['sec_sensor_installed']:
            if 'sec_sensor_chans' not in self._config:
                self.ok = False
                self.logmsg(logging.ERROR, 'Missing key in configuration file: ' + 'sec_sensor_chans')
            if 'sec_sensor_loc' not in self._config:
                self.ok = False
                self.logmsg(logging.ERROR, 'Missing key in configuration file: ' + 'sec_sensor_loc')

        if 'ref_azimuth_data' not in self._config:
            self.ok = False
            self.logmsg(logging.ERROR, 'Missing key in configuration file: ' + 'ref_azimuth_data')

        if 'ref_absolute_data' not in self._config:
            self.ok = False
            self.logmsg(logging.ERROR, 'Missing key in configuration file: ' + 'ref_absolute_data')

        for sect in ['ref_azimuth_data', 'ref_absolute_data']:
            if sect in self._config:
                if 'process' not in self._config[sect]:
                    self.ok = False
                    self.logmsg(logging.ERROR,
                                'Missing key in configuration file: ' +
                                '{}/process'.format(sect))
                elif self._config[sect]['process']:
                    if 'ms_file' not in self._config[sect]:
                        self.ok = False
                        self.logmsg(logging.ERROR,
                                    'Missing key in configuration file: ' +
                                    '{}/ms_files'.format(sect))
                    if 'starttime_iso' not in self._config[sect]:
                        self.ok = False
                        self.logmsg(logging.ERROR,
                                    'Missing key in configuration file: ' +
                                    '{}/ms_files'.format(sect))
                    if 'starttime_iso' not in self._config[sect]:
                        self.ok = False
                        self.logmsg(logging.ERROR,
                                    'Missing key in configuration file: ' +
                                    '{}/endtime_iso'.format(sect))

        if 'arc_raw_dir' not in self._config:
            self.ok = False
            self.logmsg(logging.ERROR, 'Missing key in configuration file: ' + 'arc_raw_dir')
        else:
            fpath = self._config['arc_raw_dir']
            if not (os.path.exists(fpath) or not os.path.isdir(fpath)):
                self.logmsg(logging.ERROR, 'Raw Archive directory does not exist: {}'.format(fpath))

        if 'resp_file_dir' not in self._config:
            self.ok = False
            self.logmsg(logging.ERROR, 'Missing key in configuration file: ' + 'resp_file_dir')
        else:
            fpath = self._config['resp_file_dir']
            if not (os.path.exists(fpath) or not os.path.isdir(fpath)):
                self.logmsg(logging.ERROR, 'RESP file directory does not exist: {}'.format(fpath))

        if self.ok:
            # check azimuth settings
            if self.process_azimuth:
                if os.path.isabs(self._config['ref_azimuth_data']['ms_file']):
                    fpath = os.path.normpath(self._config['ref_azimuth_data']['ms_file'])
                else:
                    fpath = os.path.join(self.ida_cal_raw_dir, self._config['ref_azimuth_data']['ms_file'])

                if not (os.path.exists(fpath) and os.path.isfile(fpath)):
                    self.logmsg(logging.ERROR, 'Azimuth file for reference sensor not found: {}'.format(fpath))
                    self.ok = False
                else:
                    self.msfiles['azi']['ref'] = fpath

                try:
                    self._config['ref_azimuth_data']['starttime_iso'] = UTCDateTime(self.starttime('azi'))
                except:
                    self.logmsg(logging.ERROR, 'Error parsing starttime_iso for azimuth reference data: {}'.format(
                        self._config['ref_azimuth_data']['starttime_iso']
                    ))
                    self.ok = False
                try:
                    self._config['ref_azimuth_data']['endtime_iso'] = UTCDateTime(self.endtime('azi'))
                except:
                    self.logmsg(logging.ERROR, 'Error parsing endtime_iso for azimuth reference data: {}'.format(
                        self._config['ref_azimuth_data']['endtime_iso']
                    ))
                    self.ok = False

            # check absolute settings
            if self.process_absolute:
                if os.path.isabs(self._config['ref_absolute_data']['ms_file']):
                    fpath = os.path.normpath(self._config['ref_absolute_data']['ms_file'])
                else:
                    fpath = os.path.join(self.ida_cal_raw_dir, self._config['ref_absolute_data']['ms_file'])

                if not (os.path.exists(fpath) and os.path.isfile(fpath)):
                    self.logmsg(logging.ERROR, 'Absolute file for reference sensor not found: {}'.format(fpath))
                else:
                    self.msfiles['abs']['ref'] = fpath

                try:
                    self._config['ref_absolute_data']['starttime_iso'] = UTCDateTime(self.starttime('abs'))
                except:
                    self.logmsg(logging.ERROR, 'Error parsing starttime_iso for absolute reference data: {}'.format(
                        self._config['ref_absolute_data']['starttime_iso']
                    ))
                    self.ok = False
                try:
                    self._config['ref_absolute_data']['endtime_iso'] = UTCDateTime(self.endtime('abs'))
                except:
                    self.logmsg(logging.ERROR, 'Error parsing endtime_iso for absolute reference data: {}'.format(
                        self._config['ref_absolute_data']['endtime_iso']
                    ))
                    self.ok = False


    # if not using fabulous.color, don't really nee this abstraction
    def logmsg(self, loglevel, errmsg):
        if self.logger:
            if loglevel == logging.DEBUG:
                self.logger.debug(errmsg)
            elif loglevel == logging.INFO:
                self.logger.info(errmsg)
            elif loglevel == logging.WARN:
                self.logger.warning(errmsg)
            elif loglevel == logging.ERROR:
                self.logger.error(errmsg)
            elif loglevel == logging.CRITICAL:
                self.logger.critical(errmsg)

    @property
    def station(self):
        return self._config['station'].lower()

    @property
    def segment_size_secs(self):
        return self._config['segment_size_secs']

    @property
    def segment_size_trim(self):
        return self._config['segment_size_trim_secs']

    @property
    def correlation_segment_size(self):
        return self._config['correlation_segment_size_secs']

    @property
    def analysis_sample_rate(self):
        return self._config['analysis_sample_rate_hz']

    @property
    def process_azimuth(self):
        return self._config['ref_azimuth_data']['process']

    @property
    def process_absolute(self):
        return self._config['ref_absolute_data']['process']

    def dataset_enabled(self, dataset):
        if dataset.lower() == 'azi':
            return self.process_azimuth
        elif dataset.lower() == 'abs':
            return self.process_absolute
        else:
            return False

    @property
    def pri_sensor_installed(self):
        return self._config['pri_sensor_installed']

    @property
    def sec_sensor_installed(self):
        return self._config['sec_sensor_installed']

    @property
    def bp_start(self):
        return self._config['analysis_bandpass_hz'][0]

    @property
    def bp_stop(self):
        return self._config['analysis_bandpass_hz'][1]

    @property
    def coherence_cutoff(self):
        return self._config['coherence_cutoff']


    def starttime(self, datatype):
        if datatype == 'azi':
            return self._config['ref_azimuth_data']['starttime_iso']
        elif datatype == 'abs':
            return self._config['ref_absolute_data']['starttime_iso']

    def endtime(self, datatype):
        if datatype == 'azi':
            return self._config['ref_azimuth_data']['endtime_iso']
        elif datatype == 'abs':
            return self._config['ref_absolute_data']['endtime_iso']

    def starttime_datetime(self, datatype):
        return self.starttime(datatype).datetime

    def endtime_datetime(self, datatype):
        return self.endtime(datatype).datetime

    def respfilename(self, net, sta, chn, loc):
        resp_file = 'RESP.{}.{}.{}.{}'.format(net, sta, loc, chn)
        return os.path.join(self.resp_file_dir, resp_file)

    def station_sensor_loc(self, sensor):
        if sensor == 'pri':
            return self._config['pri_sensor_loc']
        elif sensor == 'sec':
            return self._config['sec_sensor_loc']
    @property
    def arc_raw_dir(self):
        return self._config['arc_raw_dir']

    @property
    def resp_file_dir(self):
        return self._config['resp_file_dir']

    def correct_ref_time(self, datatype, src1, src2):

        datatype = datatype.lower()
        src1 = src1.lower()
        src2 = src2.lower()

        if src2 not in ['pri', 'sec']:
            raise ValueError('read_sensor_data: sensor must be "pri" or "sec".')

        if datatype not in ['azi', 'abs']:
            raise ValueError('read_sensor_data: datatype must be "azi" or "abs".')

        if (src2 == 'pri') and not self.pri_sensor_installed:
            raise ValueError('Primary sensor processing not enabled. ' +
                             ' Can not use it to correct reference clock. ' +
                             ' Check configuration.')
        if (src2 == 'sec') and not self.sec_sensor_installed:
            raise ValueError('Secondary sensor processing not enabled. ' +
                             ' Can not use it to correct reference clock. ' +
                             ' Check configuration.')
        if (datatype == 'azi') and not self.process_azimuth:
            raise ValueError('Azimuth data processing not enabled. ' +
                             ' Can not use it to correct reference clock. ' +
                             ' Check configuration.')
        if (datatype == 'abs') and not self.process_absolute:
            raise ValueError('Absolute data processing not enabled. ' +
                             ' Can not use it to correct reference clock. ' +
                             ' Check configuration.')

        offset = 0
        cval = 0.0
        cfunc = []
        emsg = ''

        ref_z = self.trtpls[datatype][src1].z.copy()
        sensor_z = self.trtpls[datatype][src2].z.copy()

        # need to interpolate if ref data sampling rate > sensor sampling rate
        if ref_z.stats.sampling_rate > sensor_z.stats.sampling_rate:
            sensor_z.interpolate(ref_z.stats.sampling_rate, method='linear')

        if sensor_z and ref_z:
            ref_z.filter('bandpass', freqmin=self.bp_start, freqmax=2.0) # self.bp_stop)
            sensor_z.filter('bandpass', freqmin=self.bp_start, freqmax=2.0) # self.bp_stop)
            # take middle for time series correlation
            dur = sensor_z.stats.endtime - sensor_z.stats.starttime
            start_t = sensor_z.stats.starttime + dur/2 - self.correlation_segment_size/2
            end_t = start_t + self.correlation_segment_size
            ref_z.trim(start_t, end_t)
            sensor_z.trim(start_t, end_t)

            offset, cval, cfunc, emsg = time_offset(sensor_z, ref_z)
            if src1 == 'ref':
                for tr in self.streams[datatype][src1]:
                    tr.stats.starttime = tr.stats.starttime + offset
                    # print('Adjusting REF start time by: ', offset)
            else:
                print('{}/{} Offset: {}'.format(src1, src2, offset))

            if cval < 0.9:
                self.logmsg(logging.WARN, 'WARNING (time_offset): Correlation value < 0.9: ' + str(cval))

        return offset, cval, cfunc, emsg

    def read_ref_data(self, datatype):
        """
        Reads AZI or ABS miniseed data for ref sensor.
        """

        if (datatype == 'abs') and not self.process_absolute:
                self.logmsg(logging.WARN, 'Processing Absolute data is not enabled.')
                return False

        if (datatype == 'azi') and not self.process_azimuth:
                self.logmsg(logging.WARN, 'Processing Azimuth data is not enabled.')
                return False

        self.logmsg(logging.INFO, 'Reading REF sensor {} data...'.format(datatype.upper()))

        try:
            self.streams[datatype]['ref']= read(self.msfiles[datatype]['ref'], format='MSEED')
        except Exception as e:
            self.logmsg(logging.ERROR, 'Error reading reference data: ' + self.msfiles[datatype]['ref'])
            self.logmsg(logging.ERROR, 'Can not use REFERENCE data in sensor comparisons.')
            print()
            print(e)
            print()
            return False

        self.streams[datatype]['ref'].trim(starttime=self.starttime(datatype),
                                endtime=self.endtime(datatype))

        gaps = self.streams[datatype]['ref'].get_gaps()
        if gaps:
            self.logmsg(logging.ERROR, 'REFERENCE data has gaps within time span indicated.')
            self.streams[datatype]['ref'].print_gaps()
            self.logmsg(logging.ERROR, 'Can not use REFERENCE data in sensor comparisons.')
            self.streams[datatype]['ref'] = Stream()
            return False

        self.streams[datatype]['ref'].merge()

        for tr in self.streams[datatype]['ref']:
            tr.stats.channel = rename_chan(tr.stats.channel)
            if self._config['ref_kit_metadata']['network'].strip():
                tr.stats.network = self._config['ref_kit_metadata']['network']
            if self._config['ref_kit_metadata']['station'].strip():
                tr.stats.station = self._config['ref_kit_metadata']['station']
            if self._config['ref_kit_metadata']['location'].strip():
                tr.stats.location = self._config['ref_kit_metadata']['location']
        try:
            tr_z = self.streams[datatype]['ref'].select(component='Z')[0]
            tr_1 = self.streams[datatype]['ref'].select(component='1')[0]
            tr_2 = self.streams[datatype]['ref'].select(component='2')[0]
        except Exception as e:
            self.logmsg(logging.ERROR, 'Z12 components not found in reference data.')
            print(self.streams[datatype]['ref'])
            self.logmsg(logging.ERROR, 'Can not use REFERENCE data in sensor comparisons.')
            self.streams[datatype]['ref'] = Stream()
            return False

        self.trtpls[datatype]['ref'] = self.ChanTpl(z=tr_z, n=tr_1, e=tr_2)

        return True

    def read_sensor_data(self, datatype, sensor):
        """
        Retrieve IDA10 data from IDA Archive and comvert to miniseed for
        azi/abs time periods
        """

        datatype = datatype.lower()
        sensor = sensor.lower()

        if sensor not in ['pri', 'sec']:
            raise ValueError('read_sensor_data: sensor must be "pri" or "sec".')

        if datatype not in ['azi', 'abs']:
            raise ValueError('read_sensor_data: datatype must be "azi" or "abs".')

        # see if data previously read in...
        if self.streams[datatype][sensor]:
            self.logmsg(logging.INFO, 'Using {} sensor {} data already in memory'.format(
                sensor.upper(), datatype.upper()))
            return True

        if (sensor == 'pri') and not self.pri_sensor_installed:
            return False
        if (sensor == 'sec') and not self.sec_sensor_installed:
            return False
        if datatype == 'azi' and not self.process_azimuth:
            return False
        if datatype == 'abs' and not self.process_absolute:
            return False

        if os.path.isfile(self.arc_raw_dir):
            #assume miniseed file with ONLY needed channels and time period
            ms_name = self.arc_raw_dir
            self.waveform_files.append(ms_name)

        else:
            outname = './{}_{}_{}'.format(self.station,
                                          self.station_sensor_loc(sensor),
                                          datatype)
            ms_name = outname + '.ms'
            i10_name = outname + '.i10'
            if os.path.exists(i10_name): os.remove(i10_name)
            if os.path.exists(ms_name): os.remove(ms_name)

            self.waveform_files.append(i10_name)
            self.waveform_files.append(ms_name)

            self.logmsg(logging.INFO, 'Retrieving {} sensor {} data from archive' \
                        ' and saving in {}'.format(sensor.upper(), 
                                                   datatype.upper(),
                                                   ms_name))

            try:  # this is really too much in single try:
                i10get(self.arc_raw_dir,
                       self.station,
                       self.chanloc_codes(sensor),
                       self.starttime(datatype), self.endtime(datatype),
                       outfn=i10_name)
                pimseed(self.station, i10_name, ms_name)
            except:
                self.logmsg(logging.ERROR,
                            'Error reading and converting IDA10 data for sensor ({}/{})'.format(datatype, sensor))
                self.logmsg(logging.ERROR,
                            'Can not use {} sensor data in sensor comparisons.'.format(sensor))
                return False

        try:
            self.streams[datatype][sensor] = read(ms_name)
            self.streams[datatype][sensor].trim(starttime=self.starttime(datatype),
                                                endtime=self.endtime(datatype))
            self.streams[datatype][sensor].merge()
            gaps = self.streams[datatype][sensor].get_gaps()
            if gaps:
                self.logmsg(logging.ERROR,
                            '{} Sensor data has gaps within the time span indicated.'.format(sensor))
                self.streams[datatype][sensor].print_gaps()
                self.logmsg(logging.ERROR,
                            'Can not use {} sensor data in sensor comparisons.'.format(sensor))
                return False
            tr_z = self.streams[datatype][sensor].select(component='Z')[0] #  .copy()
            tr_1 = self.streams[datatype][sensor].select(component='1')[0] #  .copy()
            tr_2 = self.streams[datatype][sensor].select(component='2')[0] #  .copy()
            self.trtpls[datatype][sensor] = self.ChanTpl(z=tr_z, n=tr_1, e=tr_2)
            self.msfiles[datatype][sensor] = os.path.abspath(ms_name)
            self.streams[datatype][sensor].write(self.msfiles[datatype][sensor], format='MSEED')
        except:
            self.logmsg(logging.ERROR,
                        'Error reading processing miniseed data for sensor ({}/{})'.format(datatype, sensor))
            self.logmsg(logging.ERROR,
                        'Can not use {} sensor data in sensor comparisons.'.format(sensor))
            return False


        return True

    def read_responses(self, sens1, sens2):

        sens1 = sens1.lower()
        sens2 = sens2.lower()

        if sens1 not in ['ref', 'sec']:
            raise ValueError('read_responses: sens1 must be "ref" or "sec"')

        if sens2 not in ['pri', 'sec']:
            raise ValueError('read_responses: sens2 must be "pri" or "sec".')

        reftpl = self.trtpls['azi']['ref'] if self.trtpls['azi']['ref'] else self.trtpls['abs']['ref']
        pritpl = self.trtpls['azi']['pri'] if self.trtpls['azi']['pri'] else self.trtpls['abs']['pri']
        sectpl = self.trtpls['azi']['sec'] if self.trtpls['azi']['sec'] else self.trtpls['abs']['sec']

        if (sens1 == 'ref') and (not reftpl):
            self.logmsg(logging.ERROR, 'Can not retrieve REFERENCE sensor ' + 
                           'response before reading reference data.')
            return False, False
        if ((sens1 == 'pri') or (sens2 == 'pri')) and (not pritpl):
            self.logmsg(logging.ERROR, 'Can not retrieve PRIMARY sensor ' + 
                           'response before reading primary sensor data.')
            return False, False
        if ((sens1 == 'sec') or (sens2 == 'sec')) and (not sectpl):
            self.logmsg(logging.ERROR, 'Can not retrieve SECONDARY sensor ' + 
                           'response before reading secondary sensor data.')
            return False, False

        if sens1 == 'ref':
            pass
        elif sens1 == 'sec':
            reftpl = sectpl

        if sens2 == 'pri':
            sentpl = pritpl
        else:
            sentpl = sectpl

        if not self.responses[sens1][sens2]:
            self.logmsg(logging.INFO, 'Computing {} response for {}/{} comparison...'.format(sens1.upper(),
                                                                         sens1.upper(),
                                                                         sens2.upper()))
            ref_z, freqs, ref_ok = self.response_tr(reftpl.z.stats.starttime,
                                            sentpl.z.stats.sampling_rate*self.segment_size_secs,
                                            sentpl.z.stats.sampling_rate,
                                            reftpl.z.stats.network,
                                            reftpl.z.stats.station,
                                            reftpl.z.stats.channel,
                                            reftpl.z.stats.location,
                                            units='VEL')
            ref_1, freqs, ref_ok = self.response_tr(reftpl.n.stats.starttime,
                                            sentpl.n.stats.sampling_rate*self.segment_size_secs,
                                            sentpl.n.stats.sampling_rate,
                                            reftpl.n.stats.network,
                                            reftpl.n.stats.station,
                                            reftpl.n.stats.channel,
                                            reftpl.n.stats.location,
                                            units='VEL')
            ref_2, freqs, ref_ok = self.response_tr(reftpl.e.stats.starttime,
                                            sentpl.e.stats.sampling_rate*self.segment_size_secs,
                                            sentpl.e.stats.sampling_rate,
                                            reftpl.e.stats.network,
                                            reftpl.e.stats.station,
                                            reftpl.e.stats.channel,
                                            reftpl.e.stats.location,
                                            units='VEL')
            if ref_ok:
                dynlimit_resp_min(ref_z, 100.0)
                dynlimit_resp_min(ref_1, 100.0)
                dynlimit_resp_min(ref_2, 100.0)
                # this is the response of hte 'ref' sensors prepared for convolution
                # with sens2. Will NOT be regular 'ref' response if sample rates are different
                self.responses[sens1][sens2] = self.ChanTpl(z=ref_z, n=ref_1, e=ref_2)

        else:
            ref_ok = True
            self.logmsg(logging.INFO, 'Using {} response already in memory...'.format(sens1.upper()))

        if not self.responses[sens2][sens2]:
            self.logmsg(logging.INFO, 'Computing {} response for {}/{} comparison...'.format(sens2.upper(),
                                                                         sens1.upper(),
                                                                         sens2.upper()))
            sen_z, freqs, sen_ok = self.response_tr(sentpl.z.stats.starttime,
                                                sentpl.z.stats.sampling_rate*self.segment_size_secs,
                                                sentpl.z.stats.sampling_rate,
                                                sentpl.z.stats.network,
                                                sentpl.z.stats.station,
                                                sentpl.z.stats.channel,
                                                sentpl.z.stats.location,
                                                units='VEL')
            sen_1, freqs, sen_ok = self.response_tr(sentpl.n.stats.starttime,
                                                sentpl.n.stats.sampling_rate*self.segment_size_secs,
                                                sentpl.n.stats.sampling_rate,
                                                sentpl.n.stats.network,
                                                sentpl.n.stats.station,
                                                sentpl.n.stats.channel,
                                                sentpl.n.stats.location,
                                                units='VEL')
            sen_2, freqs, sen_ok = self.response_tr(sentpl.e.stats.starttime,
                                                sentpl.e.stats.sampling_rate*self.segment_size_secs,
                                                sentpl.e.stats.sampling_rate,
                                                sentpl.e.stats.network,
                                                sentpl.e.stats.station,
                                                sentpl.e.stats.channel,
                                                sentpl.e.stats.location,
                                                units='VEL')
            if sen_ok:
                # but lets exponentially taper off highest 5% of freqs before nyquist
                #  taper off high freq response before deconvolving strm2 resp
                # per idaresponse/resp.c
                taper_high_freq_resp(sen_z, 0.95)
                # clip minimum amp resp to 1/100.0 of max amp
                # per idaresponse/resp.c
                dynlimit_resp_min(sen_z, 100.0)
                taper_high_freq_resp(sen_1, 0.95)
                dynlimit_resp_min(sen_1, 100.0)
                taper_high_freq_resp(sen_2, 0.95)
                dynlimit_resp_min(sen_2, 100.0)
                # note this is the regular responses for either 'pri' or 'sec' sensors
                self.responses[sens2][sens2] = self.ChanTpl(z=sen_z, n=sen_1, e=sen_2)
        else:
            sen_ok = True
            self.logmsg(logging.INFO, 'Using {} response already in memory...'.format(sens2.upper()))

        return ref_ok, sen_ok

    def sensor_sample_rate_str(self, datatype, sensor):
        if self.streams[datatype][sensor]:
            return str(self.streams[datatype][sensor][0].stats.sampling_rate)
        else:
            return ''

    def chanloc_codes(self, sensor):

        if sensor == 'pri':
            chans = self._config.get('pri_sensor_chans', '').lower().split(',')
            loc = self._config.get('pri_sensor_loc', '').upper()
            if chans :
                return ','.join([chan+loc for chan in chans])
            else:
                return ''

        elif sensor == 'sec':
            chans = self._config.get('sec_sensor_chans', '').lower().split(',')
            loc = self._config.get('sec_sensor_loc', '').upper()
            if chans :
                return ','.join([chan+loc for chan in chans])
            else:
                return ''

    def response_tr(self, respdate, npts, sr, net, sta, chn, loc, units='VEL'):
        resp, f = None, None
        respfile = self.respfilename(net, sta, chn, loc)

        if os.path.exists(respfile):
            try:
                resp, f = evalresp(1/sr,
                                   npts,
                                   self.respfilename(net, sta, chn, loc),
                                   UTCDateTime(respdate),
                                   station=sta,
                                   channel=chn,
                                   network=net,
                                   locid=loc,
                                   units=units, freq=True)
                ok = True
            except:
                self.logmsg(logging.ERROR,
                            'Error running EVALRESP with RESP file {}.'.format(respfile))
                ok = False
        else:
            self.logmsg(logging.ERROR, 'RESP file not found: ' + respfile)
            ok = False

        return resp, f, ok

    def ms_filename(self, datatype, src):
        return Path(self.msfiles[datatype][src]).name


    def save_header(self, datatype, sumf, detf):
        header = '#'*144 + '\n'
        header += '# ANALYSIS PARAMETERS\n'
        header += '# =====================\n'
        header += '#   ida module version: {} ({})\n'.format(IDA_PKG_VERSION_HASH_STR,
                                                             IDA_PKG_VERSION_DATETIME)
        header += '#          analysis at: {}\n'.format(self.analdate)
        header += '#              dataset: {}\n'.format(datatype.upper())
        header += '#              station: {}\n'.format(self.station.upper())
        header += '#           pri sensor: {}\n'.format(self.pri_sensor_installed)
        header += '#     pri data sr (hz): {}\n'.format(self.sensor_sample_rate_str(datatype, 'pri'))
        header += '#            pri chans: {}\n'.format(self.chanloc_codes('pri').upper())
        header += '#           sec sensor: {}\n'.format(self.sec_sensor_installed)
        header += '#     sec data sr (hz): {}\n'.format(self.sensor_sample_rate_str(datatype, 'sec'))
        header += '#            sec chans: {}\n'.format(self.chanloc_codes('sec').upper())
        header += '#          ref ms file: {}\n'.format(self.msfiles[datatype]['ref'])
        header += '#     ref data sr (hz): {}\n'.format(self.sensor_sample_rate_str(datatype, 'ref'))
        header += '#   ref clock base sen: {}\n'.format(self.ref_clock_adjustment_sensor)
        header += '#   ref clock adjusted: {}\n'.format(self.ref_clock_adjusted)
        header += '# ref clock adj (secs): {}\n'.format(self.ref_clock_adjustment)
        header += '#       cal start time: {}\n'.format(self.starttime(datatype))
        header += '#         cal end time: {}\n'.format(self.endtime(datatype))
        header += '#      seg size (secs): {}\n'.format(self.segment_size_secs)
        header += '#      seg trim (secs): {}\n'.format(self.segment_size_trim)
        header += '#           coh cutoff: {}\n'.format(self.coherence_cutoff)
        header += '#     analysis sr (hz): {}\n'.format(self.analysis_sample_rate)
        header += '#     analysis bp (hz): {}\n'.format(self._config['analysis_bandpass_hz'])
        header += '#\n'
        sumhdr = '#H {:<4} {:4} {:4} {:4} '\
                 '{:>7} {:>7} {:>8} {:>7} {:>7} {:>7} {:>8} {:>6} {:>6} '\
                 '{:>10} {:<14} {:<14}\n'.format(
                      'sta', 'sen1', 'sen2', 'comp',
                      'ampmn', 'ampstd', 'angmn', 'angstd', 'lrmsmn', 'varmn',
                      'cohcut', 'segcnt', 'segtot', 'analyzedon',
                      'sen1_ms_file', 'sen2_ms_file')
        dethdr = '#H {:<4} {:4} {:4} {:4} {:27} {:17} ' \
                 '{:>7} {:>8} {:>7} {:>7} ' \
                '{:>7} {:>7} {:6} {:>10} '\
                 '{:<14} {:<14}\n'.format(
                     'sta','sen1', 'sen2', 'comp', 'start_time', 'start_epoch',
                     'amp', 'angle', 'lrms', 'var',
                     'coh', 'cohcut', 'status', 'analyzedon',
                     'sen1_ms_file', 'sen2_ms_file')
        sumf.write(header + sumhdr)
        detf.write(header + dethdr)

    def save_footer(self, datatype, sumf, detf):
        sumf.write('#'*144 + '\n')
        detf.write('#'*144 + '\n')

    def save_results(self, datatype, src1, src2, results, sumf, detf):
        analday = datetime.now().strftime('%Y-%m-%d')
        if results:
            for compres in results:
                if compres.usable_count == 0:
                    sumres = '{:<4} {}/{} : {} No usable segments. {} {}'.format(
                        self.station.upper(), src1.upper(), src2.upper(),
                        compres.component, self.msfiles[datatype][src1],
                        self.msfiles[datatype][src2])
                else:
                    sumres = '   {:<4} {:<4} {:<4} {:<4} '\
                             '{:7.4f} {:7.3f} '\
                             '{:8.3f} {:7.3f} {:7.3f} {:7.3f} '\
                             '{:8.3f} {:>6} {:>6} '\
                             '{:>10} {:<14} {:<14}'.format(
                                 self.station.upper(), src1.upper(), src2.upper(),
                                 compres.component, compres.amp_mean, compres.amp_std,
                                 (compres.ang_mean * 180./pi) % 360.,
                                 (compres.ang_std * 180./pi) % 360.,
                                 compres.lrms_mean, compres.var_mean,
                                 self.coherence_cutoff, compres.usable_count,
                                 compres.total_count, analday, self.ms_filename(datatype, src1),
                                 self.ms_filename(datatype, src2))
                self.logmsg(logging.INFO, sumres)
                sumf.write(sumres+'\n')

                for res in compres.seg_results:
                    excl = 'EXCL' if not res.can_use else 'Ok  '
                    detres = '   {:<4} {:<4} {:<4} {:<4} {} {:17.6f} {:7.4f} {:8.3f} '\
                            '{:7.3f} {:7.3f} {:7.4f} {:7.3f} {:<6} {:>10} {:<14} {:<14}'.format(
                        self.station.upper(), src1.upper(), src2.upper(), compres.component,
                        res.start_utc, res.start_epoch, res.amp, (res.ang * 180./pi + 360.) % 360.,
                        res.lrms, res.var, res.coh, self.coherence_cutoff, excl, analday,
                        self.ms_filename(datatype, src1), self.ms_filename(datatype, src2))
                    detf.write(detres+'\n')

    def analyze(self, datatype):

        self.waveform_files = []
        sensors_available = 0  # need 2+ to do any comparisons

        data_ok = self.read_ref_data(datatype)
        compare_ref = data_ok
        if not data_ok:
            self.logmsg(logging.WARN, 'Unable to process {} REFERENCE sensor data.'.format(
                datatype.upper()))
            self.logmsg(logging.WARN, 'Skipping comparisons with REFERENCE sensor')
        else:
            sensors_available += 1

        if self.sec_sensor_installed:
            data_ok = self.read_sensor_data(datatype, 'sec')
            compare_sec = data_ok
            if not data_ok:
                self.logmsg(logging.WARN, 'Unable to process {} SECONDARY sensor data.'.format(
                    datatype.upper()))
                self.logmsg(logging.WARN, 'Skipping comparisons with SECONDARY sensor')
            else:
                sensors_available += 1
                offset, _, _, _ = self.correct_ref_time(datatype, 'ref', 'sec')
                self.ref_clock_adjustment = offset
                self.ref_clock_adjusted = True
                self.ref_clock_adjustment_sensor = 'secondary'
        else:
            compare_sec = False  # not installed

        if self.pri_sensor_installed:
            data_ok = self.read_sensor_data(datatype, 'pri')
            compare_pri = data_ok
            if not data_ok:
                self.logmsg(logging.WARN, 'Unable to process {} PRIMARY sensor data.'.format(
                    datatype.upper()))
                self.logmsg(logging.WARN, 'Skipping comparisons with PRIMARY sensor')
            else:
                sensors_available += 1
                if not self.ref_clock_adjusted:
                    offset, _, _, _ = self.correct_ref_time(datatype, 'ref', 'pri')
                    self.ref_clock_adjustment = offset
                    self.ref_clock_adjusted = True
                    self.ref_clock_adjustment_sensor = 'primary'
        else:
            compare_pri = False  # not installed

        dataday = self.starttime_datetime(datatype).strftime('%Y-%m-%d')
        detfn = '{}_{}_{}_details.txt'.format(self.station.lower(), dataday, datatype)
        sumfn = '{}_{}_{}_summary.txt'.format(self.station.lower(), dataday, datatype)

        with open(sumfn, 'at') as sumf:
            with open(detfn, 'at') as detf:

                self.save_header(datatype, sumf, detf)

                if sensors_available >= 2:
                    if compare_ref and compare_sec:
                        src1, src2 = 'ref', 'sec'
                        results = self.compare_streams(datatype, src1, src2)
                        if results:
                            self.save_results(datatype, src1, src2, results, sumf, detf)
                        else:
                            sumf.write('# ERROR: No comparison results returned!\n')
                            detf.write('# ERROR: No comparison results returned!\n')

                    if compare_ref and compare_pri:
                        src1, src2 = 'ref', 'pri'
                        results = self.compare_streams(datatype, src1, src2)
                        if results:
                            self.save_results(datatype, src1, src2, results, sumf, detf)
                        else:
                            sumf.write('# ERROR: No comparison results returned!\n')
                            detf.write('# ERROR: No comparison results returned!\n')

                    if compare_pri and compare_sec:
                        src1, src2 = 'sec', 'pri'
                        results = self.compare_streams(datatype, src1, src2)
                        if results:
                            self.save_results(datatype, src1, src2, results, sumf, detf)
                        else:
                            sumf.write('# ERROR: No comparison results returned!\n')
                            detf.write('# ERROR: No comparison results returned!\n')
                else:
                    detf.write('# No comparisons possible for {} data sets\n'.format(datatype.upper()))
                    sumf.write('# No comparisons possible for {} data sets\n'.format(datatype.upper()))


                self.save_footer(datatype, detf, sumf)

        return sumfn, detfn, self.waveform_files

    def compare_streams(self, datatype, src1, src2):
        """
        Processes each pair of traces:
            trim sample count to multiple of analysis frequency rate
            decimate to analysis sample rate
            bandpass filter to analysis bandpass
            deconvolve strm2 sensor resp from strm2 and convolve strm1 response
            loop over segments performing analysis and compute
            angle, amp, and variance measures
        """

        src1_resp_ok, src2_resp_ok = self.read_responses(src1, src2)
        if not src1_resp_ok:
            self.logmsg(logging.WARN, 'Unable to read {} response.'.format(src1))
            self.logmsg(logging.WARN, 'Can not perform comparison: {}/{}'.format(
                src1.upper(), src2.upper()))
            return False
        elif not src2_resp_ok:
            self.logmsg(logging.WARN, 'Unable to read {} response.'.format(src2))
            self.logmsg(logging.WARN, 'Can not perform comparison: {}/{}'.format(
                src1.upper(), src2.upper()))
            return False

        if src1.lower() == 'ref':

            strm1_z = self.trtpls[datatype][src1].z
            strm1_1 = self.trtpls[datatype][src1].n
            strm1_2 = self.trtpls[datatype][src1].e
            strm2_z = self.trtpls[datatype][src2].z
            strm2_1 = self.trtpls[datatype][src2].n
            strm2_2 = self.trtpls[datatype][src2].e

        else:
            # if 'ref' not first sensor, assume pri vs sec
            # call as ('pri', 'sec'), but sec sensor should be the 'ref' sensor
            # since it is higher sample rate, and would not be able to convolve
            # response of lower rate sensor into high rate timeseries
            # and set clock_adj to 0.0. Assumed sensor locked together.
            strm1_z = self.trtpls[datatype]['sec'].z
            strm1_1 = self.trtpls[datatype]['sec'].n
            strm1_2 = self.trtpls[datatype]['sec'].e
            strm2_z = self.trtpls[datatype]['pri'].z
            strm2_1 = self.trtpls[datatype]['pri'].n
            strm2_2 = self.trtpls[datatype]['pri'].e

        strm1_z_resp = self.responses[src1][src2].z
        strm1_1_resp = self.responses[src1][src2].n
        strm1_2_resp = self.responses[src1][src2].e
        strm2_z_resp = self.responses[src2][src2].z
        strm2_1_resp = self.responses[src2][src2].n
        strm2_2_resp = self.responses[src2][src2].e

        self.results = self.ChanTpl(z=APSurveyComponentResult('Z'),
                                    n=APSurveyComponentResult('1'),
                                    e=APSurveyComponentResult('2'))
        self.compare_verticals(datatype, src1, src2, strm1_z, strm2_z,
                               strm1_z_resp, strm2_z_resp, self.results.z)
        self.compare_horizontals(datatype, src1, src2, strm1_1, strm1_2, strm2_1,
                                 strm1_1_resp, strm2_1_resp, self.results.n)
        self.compare_horizontals(datatype, src1, src2, strm1_1, strm1_2, strm2_2,
                                 strm1_2_resp, strm2_2_resp, self.results.e)

        return self.results

    def compare_horizontals(self, datatype, src1, src2, tr1_n, tr1_e, tr2, tr1_resp, tr2_resp, results):

        dbg_cnt = 0
        start_t = max(tr1_n.stats.starttime, tr2.stats.starttime)
        end_t = min(tr1_n.stats.endtime, tr2.stats.endtime)
        tr1_sr = tr1_n.stats.sampling_rate
        tr2_sr = tr2.stats.sampling_rate
        tr1_time_delta = 1.0/tr1_sr
        tr2_time_delta = 1.0/tr2_sr

        while start_t + self.segment_size_secs < end_t:
            dbg_cnt += 1

            # print('satrting segment analysis...')
            tr1_n_seg = tr1_n.slice(starttime=start_t,
                                    endtime=start_t + self.segment_size_secs - tr1_time_delta)
            tr1_e_seg = tr1_e.slice(starttime=start_t,
                                    endtime=start_t + self.segment_size_secs - tr1_time_delta)
            tr2_seg = tr2.slice(starttime=start_t,
                                endtime=start_t + self.segment_size_secs - tr2_time_delta)

            tr1_n_seg = tr1_n_seg.data.astype(float64)
            tr1_e_seg = tr1_e_seg.data.astype(float64)
            tr2_seg = tr2_seg.data.astype(float64)

            # need to demean and taper strm2, deconvolve it's resp, convolve strm1 resp
            tr1_n_seg = ss.detrend(tr1_n_seg, type='constant')
            tr1_e_seg = ss.detrend(tr1_e_seg, type='constant')
            tr2_seg = ss.detrend(tr2_seg, type='constant')

            # each side Creating tapers...
            fraction = (self.segment_size_trim/self.segment_size_secs)
            taper = tukey(len(tr2_seg), alpha=fraction * 2, sym=True)

            tr2_fft = rfft(multiply(tr2_seg, taper))
            tr2_fft_dcnv = divide(tr2_fft[1:], tr2_resp[1:])
            tr2_fft_dcnv = insert(tr2_fft_dcnv, 0, 0.0)
            tr2_dcnv = irfft(tr2_fft_dcnv)
            # demean again before convolving strm2 resp
            tr2_dcnv = ss.detrend(tr2_dcnv, type='constant')
            tr2_fft = rfft(multiply(tr2_dcnv, taper))
            tr2_fft_cnv = multiply(tr2_fft, tr1_resp)
            tr2_cnv = irfft(tr2_fft_cnv)
            del tr2_fft, tr2_fft_dcnv, tr2_dcnv, tr2_fft_cnv

            # trim segment_size_trim in secs
            trim_cnt = round(tr1_sr * self.segment_size_trim)
            tr1_n_seg = tr1_n_seg[trim_cnt:-trim_cnt]
            tr1_e_seg = tr1_e_seg[trim_cnt:-trim_cnt]
            trim_cnt = round(tr2_sr * self.segment_size_trim)
            tr2_cnv = tr2_cnv[trim_cnt:-trim_cnt]

            tr1_n_seg = ss.resample(tr1_n_seg, round(len(tr1_n_seg) / (tr1_sr / self.analysis_sample_rate)))
            tr1_e_seg = ss.resample(tr1_e_seg, round(len(tr1_e_seg) / (tr1_sr / self.analysis_sample_rate)))
            tr2_cnv = ss.resample(tr2_cnv, round(len(tr2_cnv) / (tr2_sr / self.analysis_sample_rate)))

            tr1_n_seg = ss.detrend(tr1_n_seg, type='linear')
            tr1_e_seg = ss.detrend(tr1_e_seg, type='linear')
            tr2_cnv = ss.detrend(tr2_cnv, type='linear')

            fraction = 1/12
            taper = tukey(len(tr1_n_seg), fraction * 2, sym=True)
            tr1_n_seg *= taper
            tr1_e_seg *= taper
            tr2_cnv *= taper

            tr1_n_seg = osf.bandpass(tr1_n_seg,
                                     self.bp_start, self.bp_stop,
                                     self.analysis_sample_rate, zerophase=True)
            tr1_e_seg = osf.bandpass(tr1_e_seg,
                                     self.bp_start, self.bp_stop,
                                     self.analysis_sample_rate, zerophase=True)
            tr2_cnv = osf.bandpass(tr2_cnv,
                                   self.bp_start, self.bp_stop,
                                   self.analysis_sample_rate, zerophase=True)

            if self.debug and (dbg_cnt == 5):
                fig = plt.figure(figsize=(8.5, 11))

                plt.plot(tr2_cnv)
                plt.plot(tr1_n_seg)
                plt.plot(tr1_e_seg)
                plt.show()
                input('hit any ley')


            # set up matrix to solve for cos(wn) & cos(we) where:
            #   wn is angle of tr2_cnv from ref North
            #   we is angle of tr2_cnv from ref East
            # Treat cos(we) as sin(pi/2 - we) where pi/2 - we is angle with North,
            # then take arctan cos(wn)/cos(we)
            dc = ones(len(tr1_n_seg), dtype=float64)  # to take care of any non-zero means
            mat = array([dc, tr1_n_seg, tr1_e_seg])
            matinv = mat.transpose()
            solution, resid, _, _ = la.lstsq(matinv, tr2_cnv)
            w = arctan2(solution[2], solution[1])

            amp_ratio =  sqrt(solution[1]*solution[1] + solution[2]*solution[2])

            lrms = log10(sqrt(multiply(tr2_cnv, tr2_cnv).sum() / len(tr2_cnv)))
            syn  = dot(matinv, solution)
            res = tr2_cnv - syn
            myvar = std(res) / std(tr2_cnv)
            coh = dot(tr2_cnv, syn) / sqrt(dot(tr2_cnv, tr2_cnv) * dot(syn, syn))

            results.add_segment(APSurveySegmentResult(start_t, w, resid, amp_ratio,
                                                      lrms, myvar, coh,
                                                      coh >= self.coherence_cutoff))

            # start time for next segment
            start_t += self.segment_size_secs


    def compare_verticals(self, datatype, src1, src2, tr1, tr2, tr1_resp, tr2_resp, results):

        dbg_cnt = 0
        start_t = max(tr1.stats.starttime, tr2.stats.starttime)
        tr1_sr = tr1.stats.sampling_rate
        tr2_sr = tr2.stats.sampling_rate
        tr1_time_delta = 1.0/tr1_sr
        tr2_time_delta = 1.0/tr2_sr
        while start_t + self.segment_size_secs < tr1.stats.endtime:
            dbg_cnt += 1

            # print('satrting segment analysis...')
            tr1_seg = tr1.slice(starttime=start_t,
                                        endtime=start_t + self.segment_size_secs - tr1_time_delta)
            tr2_seg = tr2.slice(starttime=start_t,
                                        endtime=start_t + self.segment_size_secs - tr2_time_delta)

            tr1_seg = tr1_seg.data.astype(float64)
            tr2_seg = tr2_seg.data.astype(float64)

            # need to demean and taper strm2, deconvolve it's resp, convolve strm1 resp
            tr1_seg = ss.detrend(tr1_seg, type='constant')
            tr2_seg = ss.detrend(tr2_seg, type='constant')

            # each side Creating tapers...
            fraction = (self.segment_size_trim/self.segment_size_secs)
            taper = tukey(len(tr2_seg), alpha=fraction * 2, sym=True)

            tr2_fft = rfft(multiply(tr2_seg, taper))
            tr2_fft_dcnv = divide(tr2_fft[1:], tr2_resp[1:])
            tr2_fft_dcnv = insert(tr2_fft_dcnv, 0, 0.0)
            tr2_dcnv = irfft(tr2_fft_dcnv)
            # demean again before convolving strm2 resp
            tr2_dcnv = ss.detrend(tr2_dcnv, type='constant')
            tr2_fft = rfft(multiply(tr2_dcnv, taper))
            tr2_fft_cnv = multiply(tr2_fft, tr1_resp)
            tr2_cnv = irfft(tr2_fft_cnv)
            del tr2_fft, tr2_fft_dcnv, tr2_dcnv, tr2_fft_cnv

            # trim segment_size_trim in secs
            trim_cnt = round(tr1_sr * self.segment_size_trim)
            tr1_seg = tr1_seg[trim_cnt:-trim_cnt]
            trim_cnt = round(tr2_sr * self.segment_size_trim)
            tr2_cnv = tr2_cnv[trim_cnt:-trim_cnt]

            tr1_seg = ss.resample(tr1_seg, round(len(tr1_seg) / (tr1_sr / self.analysis_sample_rate)))
            tr2_cnv = ss.resample(tr2_cnv, round(len(tr2_cnv) / (tr2_sr / self.analysis_sample_rate)))

            tr1_seg = ss.detrend(tr1_seg, type='linear')
            tr2_cnv = ss.detrend(tr2_cnv, type='linear')

            fraction = 1/12
            taper = tukey(len(tr1_seg), fraction * 2, sym=True)
            tr1_seg *= taper
            tr2_cnv *= taper

            tr1_seg = osf.bandpass(tr1_seg, self.bp_start, self.bp_stop, self.analysis_sample_rate, zerophase=True)
            tr2_cnv = osf.bandpass(tr2_cnv, self.bp_start, self.bp_stop, self.analysis_sample_rate, zerophase=True)

            amp_ratio = tr2_cnv.std() / tr1_seg.std()

            lrms = log10(sqrt(multiply(tr2_cnv, tr2_cnv).sum() / len(tr2_cnv)))

            syn  = tr1_seg * amp_ratio
            coh = dot(tr2_cnv, syn) / sqrt(dot(tr2_cnv, tr2_cnv) * dot(syn, syn))

            res = tr2_cnv - syn
            myvar = res.std() / tr2_cnv.std()

            results.add_segment(APSurveySegmentResult(start_t, 0.0, 0.0, amp_ratio,
                                                      lrms, myvar,
                                                      coh, coh >= self.coherence_cutoff))

            # start time for next segment
            start_t += self.segment_size_secs


