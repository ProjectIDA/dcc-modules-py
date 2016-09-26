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
# from datetime import datetime
import os.path
# from os import remove
from pathlib import Path # , PurePath
import yaml
from collections import namedtuple
# import matplotlib.pyplot as plt

from fabulous.color import red, bold
from obspy import read, UTCDateTime, Stream  # Trace
from obspy.signal.invsim import evalresp
# import obspy.signal.filter as ops
from scipy.signal import tukey
# from numpy import array, cos, sin, subtract, angle, complex128, absolute, less
from numpy import array, ones, pi, arctan2, float64, multiply, divide, log10
from numpy import mean, std, sqrt, dot, insert
from numpy.fft import rfft, irfft
import numpy.linalg as la
import scipy.signal as ss

from ida.utils import i10get, pimseed
from ida.signals.utils import time_offset, decimate_factors_425, taper_high_freq_resp, \
        dynlimit_resp_min
from ida.calibration.shaketable import rename_chan


class APSurvey(object):

    ChanTpl = namedtuple('ChanTuple', 'z n e')

    def __init__(self, fn):

        self.errs = []

        # in-memory streams of source miniseed data
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
        # note that '_v_' responses are because responses have to accound for
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
        # decimation factors for each sensor
        self.decifactors = {
            'ref': [1],
            'pri': [1],
            'sec': [1]
        }
        # time correction flags so only do it once 
        self._ref_time_corrected = {
            'abs': False,
            'azi': False
        }

        with open(fn, 'rt') as cfl:
            config_txt = cfl.read()
        try:
            self._config = yaml.load(config_txt)
            self._config_dir = Path(fn).parent
        except:
            self._config = {}
            print(red(bold('Error parsing YAML config file: ' + fn)))
            self.errs.append('Error parsing YAML config file: ' + fn)
        else:
            self._process_config()

        self._config['resp_dir'] = os.environ.get('SEEDRESP') or ''

    def _process_config(self):

        # check ENV
        if not os.environ.get('IDA_CAL_ANALYSIS_DIR'):
            self.errs.append('The env var IDA_CAL_ANALYSIS_DIR must be set.')
        if not os.environ.get('IDA_CAL_RAW_DIR'):
            self.errs.append('The env var IDA_CAL_RAW_DIR must be set.')
        if not os.environ.get('IDA_ARCHIVE_RAW_DIR'):
            self.errs.append('The env var IDA_ARCHIVE_RAW_DIR must be set.')
        if not os.environ.get('SEEDRESP'):
            self.errs.append('The env var SEEDRESP must be set.')

        raw_path = os.path.join(os.environ['IDA_CAL_RAW_DIR'], self.station, 'APSurvey')
        self.analysis_output_path = os.path.join(os.environ['IDA_CAL_ANALYSIS_DIR'],
                                                 self.station, 'APSurvey')
        # check azimuth settings
        if self.process_azimuth:

            fpath = os.path.join(raw_path, self.ref_seed_file('azi'))
            self._config['ref_azimuth_data']['ms_file'] = fpath

            if not (os.path.exists(fpath) and os.path.isfile(fpath)):
                self.errs.append('Azimuth file for reference sensor not found: {}'.format(fpath))

            try:
                self._config['ref_azimuth_data']['starttime_iso'] = UTCDateTime(self.starttime('azi'))
            except:
                self.errs.append('Error parsing starttime_iso for azimuth reference data: {}'.format(
                    self._config['ref_azimuth_data']['starttime_iso']
                ))
            try:
                self._config['ref_azimuth_data']['endtime_iso'] = UTCDateTime(self.endtime('azi'))
            except:
                self.errs.append('Error parsing endtime_iso for azimuth reference data: {}'.format(
                    self._config['ref_azimuth_data']['endtime_iso']
                ))

        # check absolute settings
        if self.process_absolute:

            fpath = os.path.join(raw_path, self.ref_seed_file('abs'))
            self._config['ref_absolute_data']['ms_file'] = fpath

            if not (os.path.exists(fpath) and os.path.isfile(fpath)):
                self.errs.append('Absolute file for reference sensor not found: {}'.format(fpath))

            try:
                self._config['ref_absolute_data']['starttime_iso'] = UTCDateTime(self.starttime('abs'))
            except:
                self.errs.append('Error parsing starttime_iso for absolute reference data: {}'.format(
                    self._config['ref_absolute_data']['starttime_iso']
                ))
            try:
                self._config['ref_absolute_data']['endtime_iso'] = UTCDateTime(self.endtime('abs'))
            except:
                self.errs.append('Error parsing endtime_iso for absolute reference data: {}'.format(
                    self._config['ref_absolute_data']['endtime_iso']
                ))

    @property
    def station(self):
        return self._config['station_info']['station'].lower()

    def ref_seed_file(self, datatype):
        if datatype == 'azi':
            return self._config['ref_azimuth_data']['ms_file']
        elif datatype == 'abs':
            return self._config['ref_absolute_data']['ms_file']

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

    @property
    def segment_size_secs(self):
        return self._config['analysis']['segment_size_secs']

    @property
    def segment_size_trim(self):
        return self._config['analysis']['segment_size_trim']

    @property
    def minimum_segment_cnt(self):
        return self._config['analysis']['minimum_segment_cnt']

    @property
    def correlation_segment_size(self):
        return self._config['analysis']['correlation_segment_size_secs']

    @property
    def analysis_sample_rate(self):
        return self._config['analysis']['analysis_sample_rate']

    @property
    def process_azimuth(self):
        return self._config['ref_azimuth_data']['process_azimuth']

    @property
    def process_absolute(self):
        return self._config['ref_absolute_data']['process_absolute']

    @property
    def pri_sensor_installed(self):
        return self._config['station_info']['pri_sensor_installed']

    @property
    def sec_sensor_installed(self):
        return self._config['station_info']['sec_sensor_installed']

    @property
    def bp_start(self):
        return self._config['analysis']['analysis_bandpass'][0]

    @property
    def bp_stop(self):
        return self._config['analysis']['analysis_bandpass'][1]

    @property
    def coherence_cutoff(self):
        return self._config['analysis']['coherence_cutoff']

    def respfilename(self, net, sta, chn, loc):
        resp_file = 'RESP.{}.{}.{}.{}'.format(net, sta, loc, chn)
        return os.path.join(self._config['resp_dir'], resp_file)

    def station_sensor_loc(self, sensor):
        if sensor == 'pri':
            return self._config['station_info']['pri_sensor_loc']
        elif sensor == 'sec':
            return self._config['station_info']['sec_sensor_loc']

    def correct_ref_time(self, datatype, sensor):

        datatype = datatype.lower()
        sensor = sensor.lower()

        if sensor not in ['pri', 'sec']:
            raise ValueError('read_sensor_data: sensor must be "pri" or "sec".')

        if datatype not in ['azi', 'abs']:
            raise ValueError('read_sensor_data: datatype must be "azi" or "abs".')

        if (sensor == 'pri') and not self.pri_sensor_installed:
            raise ValueError('Primary sensor processing not enabled. ' \
                             ' Can not use it to correct reference clock. ' \
                             ' Check configuration.')
        if (sensor == 'sec') and not self.sec_sensor_installed:
            raise ValueError('Secondary sensor processing not enabled. ' \
                             ' Can not use it to correct reference clock. ' \
                             ' Check configuration.')
        if (datatype == 'azi') and not self.process_azimuth:
            raise ValueError('Azimuth data processing not enabled. ' \
                             ' Can not use it to correct reference clock. ' \
                             ' Check configuration.')
        if (datatype == 'abs') and not self.process_absolute:
            raise ValueError('Absolute data processing not enabled. ' \
                             ' Can not use it to correct reference clock. ' \
                             ' Check configuration.')

        offset = 0
        cval = 0.0
        cfunc = []
        emsg = ''

        ref_z = self.trtpls[datatype]['ref'].z.copy()
        sensor_z = self.trtpls[datatype][sensor].z.copy()

        # need to interpolate if ref data sampling rate > sensor sampling rate
        if ref_z.stats.sampling_rate > sensor_z.stats.sampling_rate:
            sensor_z.interpolate(ref_z.stats.sampling_rate, method='linear')

        if sensor_z and ref_z:
            ref_z.filter('bandpass', freqmin=0.1, freqmax=1.0)
            sensor_z.filter('bandpass', freqmin=0.1, freqmax=1.0)
            # take middle for time series correlation
            dur = sensor_z.stats.endtime - sensor_z.stats.starttime
            start_t = sensor_z.stats.starttime + dur/2 - self.correlation_segment_size/2
            end_t = start_t + self.correlation_segment_size
            ref_z.trim(start_t, end_t)
            sensor_z.trim(start_t, end_t)

            offset, cval, cfunc, emsg = time_offset(sensor_z, ref_z)
            for tr in self.streams[datatype]['ref']:
                tr.stats.starttime = tr.stats.starttime + offset

            if cval < 0.9:
                print(red('WARNING (time_offset): Correlation value < 0.9'))

        return offset, cval, cfunc, emsg

    def read_ref_data(self, datatype):
        """
        Reads AZI or ABS miniseed data for ref sensor.
        """

        if (datatype == 'abs') and not self.process_absolute:
                print(red(bold('Processing Reference Absolute data not enabled.')))
                return False

        if (datatype == 'azi') and not self.process_azimuth:
                print(red(bold('Processing Reference Azimuth data not enabled.')))
                return False

        if self.streams[datatype]['ref']:
            print('Using REF sensor {} data already in memory.'.format(datatype.upper()))
            return True

        print('Retrieving REF sensor {} data...'.format(datatype.upper()))
        try:
            self.streams[datatype]['ref']= read(self.ref_seed_file(datatype), format='MSEED')
        except:
            print(red(bold('Error reading reference data: ' + self.ref_seed_file(datatype))))
            return False

        self.streams[datatype]['ref'].trim(starttime=self.starttime(datatype),
                                endtime=self.endtime(datatype))
        self.streams[datatype]['ref'].merge()

        for tr in self.streams[datatype]['ref']:
            tr.stats.channel = rename_chan(tr.stats.channel)
            tr.stats.network = self._config['ref_kit_metadata']['network']
            tr.stats.station = self._config['ref_kit_metadata']['station']
            tr.stats.location = self._config['ref_kit_metadata']['location']
        tr_z = self.streams[datatype]['ref'].select(component='Z')[0]#.copy()
        tr_1 = self.streams[datatype]['ref'].select(component='1')[0]#.copy()
        tr_2 = self.streams[datatype]['ref'].select(component='2')[0]#.copy()
        self.trtpls[datatype]['ref'] = self.ChanTpl(z=tr_z, n=tr_1, e=tr_2)

        self.decifactors['ref'] = decimate_factors_425(tr_z.stats.sampling_rate,
                                                       self.analysis_sample_rate)
        return True

    def read_sensor_data(self, datatype, sensor):
        """
        Retrieve IDA10 data from IDA Archive and comvert to miniseed for
        azi/abs time periods
        """

        # see if data previously read in...
        if self.streams[datatype][sensor]:
            print('Using {} sensor {} data already in memory'.format(sensor.upper(),
                                                                     datatype.upper()))
            return True

        datatype = datatype.lower()
        sensor = sensor.lower()

        if sensor not in ['pri', 'sec']:
            raise ValueError('read_sensor_data: sensor must be "pri" or "sec".')

        if datatype not in ['azi', 'abs']:
            raise ValueError('read_sensor_data: datatype must be "azi" or "abs".')

        if (sensor == 'pri') and not self.pri_sensor_installed:
            print(red(bold('Primary sensor processing not enabled. Check configuration.')))
            return False
        if (sensor == 'sec') and not self.sec_sensor_installed:
            print(red(bold('Secondary sensor processing not enabled. Check configuration.')))
            return False
        if datatype == 'azi' and not self.process_azimuth:
            print(red(bold('Azimuth data processing not enabled. Check configuration.')))
            return False
        if datatype == 'abs' and not self.process_absolute:
            print(red(bold('Absolute data processing not enabled. Check configuration.')))
            return False

        print('Retrieving {} sensor {} data...'.format(sensor.upper(),
                                                       datatype.upper()))
        outname = './{}_{}_{}'.format(self.station,
                                      self.station_sensor_loc(sensor),
                                      datatype)
        if os.path.exists(outname+'.i10'): os.remove(outname+'.i10')
        if os.path.exists(outname+'.ms'): os.remove(outname+'.ms')
        i10get(self.station,
               self._chanloc_codes(sensor),
               self.starttime(datatype), self.endtime(datatype),
               outfn=outname+'.i10')
        pimseed(self.station, outname+'.i10', outname+'.ms')
        self.streams[datatype][sensor] = read(outname + '.ms')
        self.streams[datatype][sensor].trim(starttime=self.starttime(datatype),
                                            endtime=self.endtime(datatype))
        self.streams[datatype][sensor].merge()
        tr_z = self.streams[datatype][sensor].select(component='Z')[0].copy()
        tr_1 = self.streams[datatype][sensor].select(component='1')[0].copy()
        tr_2 = self.streams[datatype][sensor].select(component='2')[0].copy()
        self.trtpls[datatype][sensor] = self.ChanTpl(z=tr_z, n=tr_1, e=tr_2)
        self.decifactors[sensor] = decimate_factors_425(tr_z.stats.sampling_rate,
                                                     self.analysis_sample_rate)

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
            print(red(bold('Can not retrieve reference sensor ' \
                           'response before reading reference data.')))
            return False
        if ((sens1 == 'pri') or (sens2 == 'pri')) and (not pritpl):
            print(red(bold('Can not retrieve primary sensor ' \
                           'response before reading primary sensor data.')))
            return False
        if ((sens1 == 'sec') or (sens2 == 'sec')) and (not sectpl):
            print(red(bold('Can not retrieve secondary sensor ' \
                           'response before reading secondary sensor data.')))
            return False


        if sens1 == 'ref':
            pass
        elif sens1 == 'sec':
            reftpl = sectpl

        if sens2 == 'pri':
            sentpl = pritpl
        else:
            sentpl = sectpl

        if not self.responses[sens1][sens2]:
            print('Computing {} response for {}/{} comparison...'.format(sens1.upper(),
                                                                         sens1.upper(),
                                                                         sens2.upper()))
            ref_z, freqs = self.response_tr(reftpl.z.stats.starttime,
                                            sentpl.z.stats.sampling_rate*self.segment_size_secs,
                                            sentpl.z.stats.sampling_rate,
                                            reftpl.z.stats.network,
                                            reftpl.z.stats.station,
                                            reftpl.z.stats.channel,
                                            reftpl.z.stats.location,
                                            units='VEL')
            dynlimit_resp_min(ref_z, 100.0)
            ref_1, freqs = self.response_tr(reftpl.n.stats.starttime,
                                            sentpl.n.stats.sampling_rate*self.segment_size_secs,
                                            sentpl.n.stats.sampling_rate,
                                            reftpl.n.stats.network,
                                            reftpl.n.stats.station,
                                            reftpl.n.stats.channel,
                                            reftpl.n.stats.location,
                                            units='VEL')
            dynlimit_resp_min(ref_1, 100.0)
            ref_2, freqs = self.response_tr(reftpl.e.stats.starttime,
                                            sentpl.e.stats.sampling_rate*self.segment_size_secs,
                                            sentpl.e.stats.sampling_rate,
                                            reftpl.e.stats.network,
                                            reftpl.e.stats.station,
                                            reftpl.e.stats.channel,
                                            reftpl.e.stats.location,
                                            units='VEL')
            dynlimit_resp_min(ref_2, 100.0)
            # this is the response of hte 'ref' sensors prepared for convolution
            # with sens2. Will NOT be regular 'ref' response if sample rates are different
            self.responses[sens1][sens2] = self.ChanTpl(z=ref_z, n=ref_1, e=ref_2)
        else:
            print('Using {} response already in memory...'.format(sens1.upper()))

        if not self.responses[sens2][sens2]:
            print('Computing {} response for {}/{} comparison...'.format(sens2.upper(),
                                                                         sens1.upper(),
                                                                         sens2.upper()))
            sen_z, freqs = self.response_tr(sentpl.z.stats.starttime,
                                            sentpl.z.stats.sampling_rate*self.segment_size_secs,
                                            sentpl.z.stats.sampling_rate,
                                            sentpl.z.stats.network,
                                            sentpl.z.stats.station,
                                            sentpl.z.stats.channel,
                                            sentpl.z.stats.location,
                                            units='VEL')
            # but lets exponentially taper off highest 5% of freqs before nyquist
            #  taper off high freq response before deconvolving strm2 resp
            # per idaresponse/resp.c
            taper_high_freq_resp(sen_z, 0.95)
            # clip minimum amp resp to 1/100.0 of max amp
            # per idaresponse/resp.c
            dynlimit_resp_min(sen_z, 100.0)
            sen_1, freqs = self.response_tr(sentpl.n.stats.starttime,
                                            sentpl.n.stats.sampling_rate*self.segment_size_secs,
                                            sentpl.n.stats.sampling_rate,
                                            sentpl.n.stats.network,
                                            sentpl.n.stats.station,
                                            sentpl.n.stats.channel,
                                            sentpl.n.stats.location,
                                            units='VEL')
            taper_high_freq_resp(sen_1, 0.95)
            dynlimit_resp_min(sen_1, 100.0)
            sen_2, freqs = self.response_tr(sentpl.e.stats.starttime,
                                            sentpl.e.stats.sampling_rate*self.segment_size_secs,
                                            sentpl.e.stats.sampling_rate,
                                            sentpl.e.stats.network,
                                            sentpl.e.stats.station,
                                            sentpl.e.stats.channel,
                                            sentpl.e.stats.location,
                                            units='VEL')
            taper_high_freq_resp(sen_2, 0.95)
            dynlimit_resp_min(sen_2, 100.0)
            # note this is the regular responses for either 'pri' or 'sec' sensors
            self.responses[sens2][sens2] = self.ChanTpl(z=sen_z, n=sen_1, e=sen_2)
        else:
            print('Using {} response already in memory...'.format(sens2.upper()))

        return True

    def _chanloc_codes(self, sensor):

        if sensor == 'pri':
            chans = self._config['station_info']['pri_sensor_chans'].split(',')
            loc = self._config['station_info']['pri_sensor_loc']
            if chans :
                return ','.join([chan+loc for chan in chans])
            else:
                return ''

        elif sensor == 'sec':
            chans = self._config['station_info']['sec_sensor_chans'].split(',')
            loc = self._config['station_info']['sec_sensor_loc']
            if chans :
                return ','.join([chan+loc for chan in chans])
            else:
                return ''

    def response_tr(self, respdate, npts, sr, net, sta, chn, loc, units='VEL'):
        resp, f = evalresp(1/sr,
                           npts,
                           self.respfilename(net, sta, chn, loc),
                           UTCDateTime(respdate),
                           station=sta,
                           channel=chn,
                           network=net,
                           locid=loc,
                           units=units, freq=True)
        return resp, f

    # def compare_streams(self, strm1, strm2):
    def compare_streams(self, datatype='abs', src1='ref', src2='sec'):
        """
        Loops through component traces and runs azi and abs analysis on segments

        Processes each pair of traces:
            trim sample count to multiple of analysis frequency rate
            decimate to analysis sample rate
            bandpass filter to analysis bandpass
            deconvolve strm2 sensor resp from strm2 and convolve strm1 response
            loop over segments performing analysis and compute
            angle, amp, and variance measures


        """
        self.read_ref_data(datatype)

        if 'pri' in [src1, src2]:
            self.read_sensor_data(datatype, 'pri')
        if 'sec' in [src1, src2]:
            self.read_sensor_data(datatype, 'sec')

        self.read_responses(src1, src2)

        if src1.lower() == 'ref':

            self.correct_ref_time(datatype, src2)

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

        self.compare_verticals(src1, src2, strm1_z, strm2_z, strm1_z_resp, strm2_z_resp)
        self.compare_horizontals(src1, src2, strm1_1, strm1_2, strm2_1, strm1_1_resp, strm2_1_resp)
        self.compare_horizontals(src1, src2, strm1_1, strm1_2, strm2_2, strm1_2_resp, strm2_2_resp)

    def compare_horizontals(self, src1, src2, tr1_n, tr1_e, tr2, tr1_resp, tr2_resp):

        # # lets find decimation factors
        strm1_factors = self.decifactors[src1]
        strm2_factors = self.decifactors[src2]

        dbg_cnt = 0
        amp_list = []  # list of amp values for each segment beating coh_cutoff
        ang_list = []  # list of amp values for each segment beating coh_cutoff
        dc_list = []  # list of amp values for each segment beating coh_cutoff
        start_t = max(tr1_n.stats.starttime, tr2.stats.starttime)
        end_t = min(tr1_n.stats.endtime, tr2.stats.endtime)
        tr1_sr = tr1_n.stats.sampling_rate
        tr2_sr = tr2.stats.sampling_rate
        tr1_time_delta = 1.0/tr1_sr
        tr2_time_delta = 1.0/tr2_sr
        comp = tr2.stats.channel[-1].upper()

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

            # decimate, detrend and bandpass filter
            for factor in strm1_factors:
                tr1_n_seg = ss.decimate(tr1_n_seg, factor, ftype='fir', zero_phase=True)
                tr1_e_seg = ss.decimate(tr1_e_seg, factor, ftype='fir', zero_phase=True)
            for factor in strm2_factors:
                tr2_cnv = ss.decimate(tr2_cnv, factor, ftype='fir', zero_phase=True)

            tr1_n_seg = ss.detrend(tr1_n_seg, type='linear')
            tr1_e_seg = ss.detrend(tr1_e_seg, type='linear')
            tr2_cnv = ss.detrend(tr2_cnv, type='linear')

            fraction = 1/12
            taper = tukey(len(tr1_n_seg), fraction * 2, sym=True)
            tr1_n_seg *= taper
            tr1_e_seg *= taper
            tr2_cnv *= taper

            # print('bandpass..')
            b, a = ss.iirdesign(wp=[.1, .3], ws=[0.05, .4], gstop=80, gpass=1, ftype='cheby2')
            tr1_n_seg = ss.lfilter(b, a, tr1_n_seg)
            tr1_e_seg = ss.lfilter(b, a, tr1_e_seg)
            tr2_cnv = ss.lfilter(b, a, tr2_cnv)

            # set up matrix to solve for cos(w) & sin(w) where w is angle of tr2_cnv from North
            # set up with [N, E] so solution comes out [cos, sin]
            dc = ones(len(tr1_n_seg), dtype=float64)  # to take care of any non-zero means
            mat = array([dc, tr1_n_seg, tr1_e_seg])
            matinv = mat.transpose()
            solution, _, _, _ = la.lstsq(matinv, tr2_cnv)
            # print('solution:', solution)
            w = arctan2(solution[2], solution[1])
            if w < 0:
                w += 2*pi
            w_deg = w*180/pi
            amp_ratio =  sqrt(solution[1]*solution[1] + solution[2]*solution[2])

#       import numpy.linalg as la
#       # Defining matrixes
#       n = np.array([1.,3,4])
#       e = np.array([2.,6,8])
#       b = np.array([2.2, 6.7, 8.9])
#
#       a1= np.array([n,e])
#       a1t = a1.transpose()
#       print(a1)
#       print(a1t)
#
#       x, _, _, _ = la.lstsq(a1t, b)
#       print(x)
            lrms = log10(sqrt(multiply(tr2_cnv, tr2_cnv).sum() / len(tr2_cnv)))
            syn  = dot(matinv, solution)
            res = tr2_cnv - syn
            myvar = std(res) / std(tr2_cnv)
            coh = dot(tr2_cnv, syn) / sqrt(dot(tr2_cnv, tr2_cnv) * dot(syn, syn))

            res = '{} ({})  Comp: {}  ang: {:6.3f}; amp: {:5.3f}; lrms: {:5.3f}; ' \
                    'var: {:5.3f}; coh: {:5.3f}'.format(
                    start_t, start_t.timestamp, comp, w_deg, amp_ratio, lrms,  myvar, coh)

            if True: # coh >= self.coherence_cutoff:
                amp_list.append(amp_ratio)
                ang_list.append(w_deg)
                dc_list.append(solution[0])
                print(res)
            else:
                print(red(bold(res)))

            # start time for next segment
            start_t += self.segment_size_secs


        print('Comp:', comp, 'seg cnt >= ', self.coherence_cutoff, ':', len(amp_list), 'of', dbg_cnt)
        if len(amp_list) > 0:
            print(comp, 'AMP RATIO AVERAGE:', mean(amp_list))
            print(comp, 'AMP RATIO  STDDEV:', std(amp_list))
            print(comp, 'ANG AVERAGE:', mean(ang_list))
            print(comp, 'ANG STDDEV:', std(ang_list))
        else:
            print(red(bold('No segments with coh >= cutoff ')))

    def compare_verticals(self, src1, src2, tr1, tr2, tr1_resp, tr2_resp):

        # # lets find decimation factors
        strm1_factors = self.decifactors[src1]
        strm2_factors = self.decifactors[src2]

        dbg_cnt = 0
        amp_list = []  # list of amp values for each segment beating coh_cutoff
        start_t = max(tr1.stats.starttime, tr2.stats.starttime)
        while start_t + self.segment_size_secs < tr1.stats.endtime:
            dbg_cnt += 1

            # print('satrting segment analysis...')
            tr1_seg = tr1.slice(starttime=start_t,
                                        endtime=start_t + self.segment_size_secs - 1/tr1.stats.sampling_rate)
            tr2_seg = tr2.slice(starttime=start_t,
                                        endtime=start_t + self.segment_size_secs - 1/tr2.stats.sampling_rate)

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
            trim_cnt = round(tr1.stats.sampling_rate * self.segment_size_trim)
            tr1_seg = tr1_seg[trim_cnt:-trim_cnt]
            trim_cnt = round(tr2.stats.sampling_rate * self.segment_size_trim)
            tr2_cnv = tr2_cnv[trim_cnt:-trim_cnt]
#            if dbg_cnt == 9: #  in range(10, 10):
#                plt.plot(tr1_seg)
#                plt.plot(tr2_dcnv, 'r')
#                plt.plot(tr2_cnv, 'm')
#                plt.show()

            # decimate, detrend and bandpass filter
            for factor in strm1_factors:
                tr1_seg = ss.decimate(tr1_seg, factor, ftype='fir', zero_phase=True)
            for factor in strm2_factors:
                tr2_cnv = ss.decimate(tr2_cnv, factor, ftype='fir', zero_phase=True)
#           if dbg_cnt in range(10,10):
#               plt.plot(tr1_seg)
#               plt.plot(tr2_cnv)
#               plt.show()

            # print('strm1 mean:{}; strm2 mean:{}'.format(tr1_seg.mean(), tr2_cnv.mean()))
            tr1_seg = ss.detrend(tr1_seg, type='linear')
            tr2_cnv = ss.detrend(tr2_cnv, type='linear')
            # print('strm1 mean:{}; strm2 mean:{}'.format(tr1_seg.mean(), tr2_cnv.mean()))

            fraction = 1/12
            taper = tukey(len(tr1_seg), fraction * 2, sym=True)
            tr1_seg *= taper
            tr2_cnv *= taper

            # print('bandpass..')
            b, a = ss.iirdesign(wp=[.1, .3], ws=[0.05, .4], gstop=80, gpass=1, ftype='cheby2')
            tr1_seg = ss.lfilter(b, a, tr1_seg)
            tr2_cnv = ss.lfilter(b, a, tr2_cnv)
#           if dbg_cnt == 10:
#               plt.plot(tr1_seg)
#               plt.plot(tr2_cnv)
#               plt.show()

            amp_ratio = tr2_cnv.std() / tr1_seg.std()
            lrms = log10(sqrt(multiply(tr2_cnv, tr2_cnv).sum() / len(tr2_cnv)))
            syn  = tr1_seg * amp_ratio
            res = tr2_cnv - syn
            myvar = res.std() / tr2_cnv.std()
            coh = dot(tr2_cnv, syn) / sqrt(dot(tr2_cnv, tr2_cnv) * dot(syn, syn))

            res = '{}:  ang: 0.0; amp: {}; lrms: {}; var: {}; coh: {}'.format(start_t, amp_ratio,
                                                                   lrms, myvar, coh)

            if coh >= self.coherence_cutoff:
                amp_list.append(amp_ratio)
                print(res)
            else:
                print(red(bold(res)))

            # start time for next segment
            start_t += self.segment_size_secs


        print('SEG CNT >= ', self.coherence_cutoff, ':', len(amp_list), 'of', dbg_cnt)
        if len(amp_list) > 0:
            print('AMP RATIO AVERAGE:', mean(amp_list))
            print('AMP RATIO  STDDEV:', std(amp_list))
        else:
            print(red(bold('No segments with coh >= cutoff ')))


    def adjust_ref_clock(self, ref_st, sen_st):
        """
        Adjust timing of ref_st based on correlation with sen_tr.

        It is assumed that sen_tr has GPS lock and ref_tr does not.
        Will typically use Z component.

        Args:
            ref_st (Stream): reference stream,
            sen_st (Stream): deployed sensor stream with locked clock.

        Returns:
            new_ref_st (Stream): copy of ref_st with starttime adjusted

        """

class APSurveryComponentResult(object):

    def __init__(self, ):
        self.seg_results = []
        self._amp_mean = None
        self._amp_std = None
        self._ang_mean = None
        self._ang_std = None
        self._ang_resid = None
        self._lrms_mean = None
        self._lrms_std = None
        self._var_mean = None
        self._var_std = None

    def _recalc(self):
        pass

class APSurveySegmentResult(object):

    def __init__(self, refstrachnloc, targ_stachnloc, start_t, ang, ang_resid, amp, lrms, var, coh):

        self._ang = ang
        self._ang_resid = ang_resid
        self._amp = amp
        self._lrms = lrms
        self._var = var
        self._coh = coh

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

