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
#from datetime import datetime
import os.path
#from os import remove
from pathlib import Path, PurePath
import yaml
#import collections

from fabulous.color import red, bold
from obspy import read, UTCDateTime, Stream#, Trace
from obspy.signal.invsim import evalresp
import obspy.signal.filter as ops
from numpy import float32, std, mean
from numpy.fft import rfft, irfft
import scipy.signal as ss

from ida.utils import i10get, pimseed
from ida.signals.utils import time_offset
from ida.calibration.shaketable import rename_chan


class AbsOnsiteConfig(object):

    def __init__(self, fn):

        self.errs = []

        # in-memory streams of source miniseed data
        self._ref_azi_strm = Stream()
        self._ref_abs_strm = Stream()
        self._pri_azi_strm = Stream()
        self._pri_abs_strm = Stream()
        self._sec_azi_strm = Stream()
        self._sec_abs_strm = Stream()

        # in-memory streams of miniseed data prepped for analysis


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
            self.process_config()

        self._config['respfile_dir'] = os.environ.get('SEEDRESP') or ''

        if self.errs:
            print(red(bold('\nThe following problems were encountered while parsing the config file: \n')))
            print(red('\n'.join(self.errs)))


    def process_config(self):

        # check ENV
        if not os.environ.get('IDA_CAL_ABS_SITEDATA_DIR'):
            self.errs.append('The env var IDA_CAL_ABS_SITEDATA_DIR must be set to the root directory of the onsite reference data.')
        if not os.environ.get('IDA_ARCHIVE_RAW_DIR'):
            self.errs.append('The env var IDA_ARCHIVE_RAW_DIR must be set to the root directory IDA10 waveform data.')

        # check azimuth settings
        if self._config['process_azimuth']:

            if not PurePath(self._config['azimuth_ref_data']['ms_file']).is_absolute():
                fpath = os.path.join(os.environ.get('IDA_CAL_ABS_SITEDATA_DIR'),
                                     self._config['azimuth_ref_data']['ms_file'])
                self._config['azimuth_ref_data']['ms_file'] = fpath
            else:
                fpath = self._config['azimuth_ref_data']['ms_file']

            if not (os.path.exists(fpath) and os.path.isfile(fpath)):
                self.errs.append('Azimuth file for reference sensor not found: {}'.format(fpath))

            try:
                self._config['azimuth_ref_data']['starttime_iso'] = \
                    UTCDateTime(self._config['azimuth_ref_data']['starttime_iso'])
            except:
                self.errs.append('Error parsing starttime_iso for azimuth reference data: {}'.format(
                    self._config['azimuth_ref_data']['starttime_iso']
                ))
            try:
                self._config['azimuth_ref_data']['endtime_iso'] = \
                    UTCDateTime(self._config['azimuth_ref_data']['endtime_iso'])
            except:
                self.errs.append('Error parsing endtime_iso for azimuth reference data: {}'.format(
                    self._config['azimuth_ref_data']['endtime_iso']
                ))

        # check absolute settings
        if self._config['process_absolute']:

            if not PurePath(self._config['absolute_ref_data']['ms_file']).is_absolute():
                fpath = os.path.join(os.environ.get('IDA_CAL_ABS_SITEDATA_DIR'),
                                     self._config['absolute_ref_data']['ms_file'])
                self._config['absolute_ref_data']['ms_file'] = fpath
            else:
                fpath = self._config['absolute_ref_data']['ms_file']

            if not (os.path.exists(fpath) and os.path.isfile(fpath)):
                self.errs.append('Absolute file for reference sensor not found: {}'.format(fpath))

            try:
                self._config['absolute_ref_data']['starttime_iso'] = \
                    UTCDateTime(self._config['absolute_ref_data']['starttime_iso'])
            except:
                self.errs.append('Error parsing starttime_iso for absolute reference data: {}'.format(
                    self._config['absolute_ref_data']['starttime_iso']
                ))
            try:
                self._config['absolute_ref_data']['endtime_iso'] = \
                    UTCDateTime(self._config['absolute_ref_data']['endtime_iso'])
            except:
                self.errs.append('Error parsing endtime_iso for absolute reference data: {}'.format(
                    self._config['absolute_ref_data']['endtime_iso']
                ))


    @property
    def ref_azi_seed_file(self):
        return  self._config['azimuth_ref_data']['ms_file']

    @property
    def ref_abs_seed_file(self):
        return  self._config['absolute_ref_data']['ms_file']

    @property
    def azi_starttime(self):
        return self._config['azimuth_ref_data']['starttime_iso']

    @property
    def azi_endtime(self):
        return self._config['azimuth_ref_data']['endtime_iso']

    @property
    def abs_starttime(self):
        return self._config['absolute_ref_data']['starttime_iso']

    @property
    def abs_endtime(self):
        return self._config['absolute_ref_data']['endtime_iso']

    @property
    def segment_size_secs(self):
        return self._config['segment_size_secs']

    @property
    def minimum_segment_cnt(self):
        return self._config['minimum_segment_cnt']

    @property
    def correlation_segment_size(self):
        return self._config['correlation_segment_size_secs']

    def respfilename(self, tr):
        resp_file = 'RESP.{}.{}.{}.{}'.format(
            tr.stats.network, tr.stats.station, tr.stats.loc, tr.stats.channel
        )
        return os.path.join(self._config['resp_filedir'], resp_file)

    def response_tr(self, tr):
        resp = evalresp(1/tr.stats.sampling_rate,
                        tr.stats.npts,
                        self.respfilename(tr),
                        tr.stats.starttime,
                        station=tr.stats.station,
                        channel=tr.stats.channel,
                        network=tr.stats.network,
                        locid=tr.stats.loc,
                        units='VEL')
        return resp


    def correct_ref_abs_time(self):

        offset = 0
        cval = 0.0
        cfunc = []
        emsg = ''

        if self._ref_abs_strm.select(channel='BHZ'):
           ref_z = self._ref_abs_strm.select(channel='BHZ').copy()
        else:
           emsg = 'Error: Can not correlate reference and site sensor time series. ' + \
                  'No vertical channel found in reference sensor data'
           return offset, cval, cfunc, emsg

        if self._sec_abs_strm.select(channel='BHZ'):
           sensor_z = self._sec_abs_strm.select(channel='BHZ').copy()
        elif self._pri_abs_strm.select(channel='BHZ'):
           sensor_z = self._pri_abs_strm.select(channel='BHZ').copy()
        else:
           emsg = 'Error: Can not correlate reference and site sensor time series. ' + \
                          'No vertical channel found in primary or secondary sensor data'
           return offset, cval, cfunc, emsg


       # may need to decimate ref data if using 20Hz Primary
        ref_z.decimate(round(ref_z[0].stats.sampling_rate / sensor_z[0].stats.sampling_rate))

        if sensor_z and ref_z:
           self.ref_z = ref_z
           self.sensor_z = sensor_z
           ref_z.filter('bandpass',
                        freqmin=self._config['analysis_bandpass'][0],
                        freqmax=self._config['analysis_bandpass'][1])
           ref_z.normalize()
           sensor_z.filter('bandpass',
                        freqmin=self._config['analysis_bandpass'][0],
                        freqmax=self._config['analysis_bandpass'][1])
           sensor_z.normalize()
           # take middle two hours for time series correlation
           dur = sensor_z[0].stats.endtime - sensor_z[0].stats.starttime
           start_t = sensor_z[0].stats.starttime + dur/2 - self.correlation_segment_size/2
           #start_t = sensor_z[0].stats.starttime + 60*60*4
           end_t = start_t + self.correlation_segment_size
           ref_z.trim(start_t, end_t)
           sensor_z.trim(start_t, end_t)

           offset, cval, cfunc, emsg = time_offset(sensor_z[0], ref_z[0])
           for tr in self._ref_abs_strm:
               tr.stats.starttime = tr.stats.starttime + offset

        return offset, cval, cfunc, emsg

    def read_ref_data(self):
        """
        Reads AZI and ABS miniseed data for ref sensor.
        """
        if self._config['process_azimuth'] and self.ref_azi_seed_file:
            try:
                self._ref_azi_strm = read(self.ref_azi_seed_file, format='MSEED')
            except:
                print(red(bold('Error reading reference azimuth data: ' + self.ref_azi_seed_file)))
            else:
                self._ref_azi_strm.trim(starttime=self.azi_starttime, endtime=self.azi_endtime)
                self._ref_azi_strm.merge()


        if self._config['process_absolute'] and self.ref_abs_seed_file:
            try:
                self._ref_abs_strm = read(self.ref_abs_seed_file, format='MSEED')
            except:
                print(red(bold('Error reading reference absolute data: ' + self.ref_abs_seed_file)))
            else:
                self._ref_abs_strm.trim(starttime=self.abs_starttime, endtime=self.abs_endtime)
                self._ref_abs_strm.merge()

        self.clean_ref_channels()

    def clean_ref_channels(self):
        """
        Cleanup ref data channel codes not configured properly on digitizer
        """

        if self._ref_azi_strm:
            for tr in self._ref_azi_strm:
                tr.stats.channel = rename_chan(tr.stats.channel)
                tr.stats.network = self._config['field_kit_metadata']['network']
                tr.stats.station = self._config['field_kit_metadata']['station']
                tr.stats.location = self._config['field_kit_metadata']['location']
        if self._ref_abs_strm:
            for tr in self._ref_abs_strm:
                tr.stats.channel = rename_chan(tr.stats.channel)
                tr.stats.network = self._config['field_kit_metadata']['network']
                tr.stats.station = self._config['field_kit_metadata']['station']
                tr.stats.location = self._config['field_kit_metadata']['location']

    def read_sensor_data(self):
        """
        Retrieve IDA10 data from IDA Archive and comvert to miniseed for azi/abs time periods
        """

        if self._config['process_azimuth']:

            print('Retrieving primary sensor data for azimuth period...')
            outname = './{}_{}_{}'.format(self._config['site_info']['station'],
                                           self._config['site_info']['pri_sensor_loc'],
                                           'azi')
            if os.path.exists(outname+'.i10'): os.remove(outname+'.i10')
            if os.path.exists(outname+'.ms'): os.remove(outname+'.ms')
            i10get(self._config['site_info']['station'],
                   self._pri_chanloc_codes(),
                   self.azi_starttime, self.azi_endtime,
                   outfn=outname+'.i10')
            pimseed(self._config['site_info']['station'], outname+'.i10', outname+'.ms')
            self._pri_azi_strm = read(outname + '.ms')
            self._pri_azi_strm.trim(starttime=self.azi_starttime, endtime=self.azi_endtime)
            self._pri_azi_strm.merge()

            print('Retrieving secondary sensor data for azimuth period...')
            outname = './{}_{}_{}'.format(self._config['site_info']['station'],
                                             self._config['site_info']['sec_sensor_loc'],
                                             'azi')
            if os.path.exists(outname+'.i10'): os.remove(outname+'.i10')
            if os.path.exists(outname+'.ms'): os.remove(outname+'.ms')
            i10get(self._sec_azi_strm[0].stats.station,
                   self._sec_chanloc_codes(),
                   self.azi_starttime, self.azi_endtime,
                   outfn=outname + '.i10')
            pimseed(self._config['site_info']['station'], outname + '.i10', outname + '.ms')
            self._sec_azi_strm = read(outname + '.ms')
            self._sec_azi_strm.trim(starttime=self.azi_starttime, endtime=self.azi_endtime)
            self._sec_azi_strm.merge()

        if self._config['process_absolute'] and self._ref_abs_strm:
            print('Retrieving primary sensor data for absolute period...')
            outname = './{}_{}_{}'.format(self._config['site_info']['station'],
                                           self._config['site_info']['pri_sensor_loc'],
                                           'abs')
            if os.path.exists(outname+'.i10'): os.remove(outname+'.i10')
            if os.path.exists(outname+'.ms'): os.remove(outname+'.ms')
            i10get(self._config['site_info']['station'],
                   self._pri_chanloc_codes(),
                   self.abs_starttime, self.abs_endtime,
                   outfn=outname+'.i10')
            pimseed(self._config['site_info']['station'], outname+'.i10', outname+'.ms')
            self._pri_abs_strm = read(outname + '.ms')
            self._pri_abs_strm.trim(starttime=self.abs_starttime, endtime=self.abs_endtime)
            self._pri_abs_strm.merge()

            print('Retrieving secondary sensor data for absolute period...')
            outname = './{}_{}_{}'.format(self._config['site_info']['station'],
                                           self._config['site_info']['sec_sensor_loc'],
                                           'abs')
            if os.path.exists(outname+'.i10'): os.remove(outname+'.i10')
            if os.path.exists(outname+'.ms'): os.remove(outname+'.ms')
            i10get(self._config['site_info']['station'],
                   self._sec_chanloc_codes(),
                   self.abs_starttime, self.abs_endtime,
                   outfn=outname+'.i10')
            pimseed(self._config['site_info']['station'], outname+'.i10', outname+'.ms')
            self._sec_abs_strm = read(outname + '.ms')
            self._sec_abs_strm.trim(starttime=self.abs_starttime, endtime=self.abs_endtime)
            self._sec_abs_strm.merge()


    def _pri_chanloc_codes(self):
        chans = self._config['site_info']['pri_sensor_chans'].split(',')
        loc = self._config['site_info']['pri_sensor_loc']
        if chans :
            return ','.join([chan+loc for chan in chans])
        else:
            return ''

    def _sec_chanloc_codes(self):
        chans = self._config['site_info']['sec_sensor_chans'].split(',')
        loc = self._config['site_info']['sec_sensor_loc']
        if chans :
            return ','.join([chan+loc for chan in chans])
        else:
            return ''

    def compare_streams(self, strm1, strm2):
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
        strm1_z = strm1.select(component='Z')[0].copy()
        strm1_1 = strm1.select(component='1')[0].copy()
        strm1_2 = strm1.select(component='2')[0].copy()
        strm2_z = strm2.select(component='Z')[0].copy()
        strm2_1 = strm2.select(component='1')[0].copy()
        strm2_2 = strm2.select(component='2')[0].copy()

        strm1_z_resp = self.response_tr(strm1_z)
        strm1_1_resp = self.response_tr(strm1_1)
        strm1_2_resp = self.response_tr(strm1_2)
        strm2_z_resp = self.response_tr(strm2_z)
        strm2_1_resp = self.response_tr(strm2_1)
        strm2_2_resp = self.response_tr(strm2_2)

        # lets do vertical first...
        dbg_cnt = 0
        start_t = strm1_z.stats.starttime
        while start_t + self.segment_size_secs < strm1_z.stats.endtime:

            strm1_z_seg = strm1_z.slice(starttime=start_t, endtime=start_t + self.segment_size_secs).data.astype(float32)
            strm2_z_seg = strm2_z.slice(starttime=start_t, endtime=start_t + self.segment_size_secs).data.astype(float32)

            mn = mean(strm1_z_seg)
            strm1_z_seg -= mn
            mn = mean(strm2_z_seg)
            strm2_z_seg -= mn

            # move strm2 into strm1 response space
            strm2_z_fft = rfft(strm2_z_seg)
            strm2_z_fft /= strm2_z_resp
            strm2_z_fft *= strm1_z_resp
            strm2_z_cnv = irfft(strm2_z_seg, len(strm2_z_seg))

            # trim ends
            trim_cnt = strm1_z.stats.sampling_rate * self._config['segment_size_trim']
            strm1_z_seg = strm1_z_seg[trim_cnt:-trim_cnt]
            trim_cnt = strm2_z.stats.sampling_rate * self._config['segment_size_trim']
            strm2_z_seg = strm2_z_seg[trim_cnt:-trim_cnt]

            # decimate, detrend and bandpass filter
            strm1_z_seg = ss.decimate(strm1_z_seg,
                                      int(round(strm1_z[0].stats.sampling_rate / self._config['analysis_sample_rate'])),
                                      ftype='fir')
            strm2_z_cnv = ss.decimate(strm2_z_cnv,
                                      int(round(strm2_z[0].stats.sampling_rate / self._config['analysis_sample_rate'])),
                                      ftype='fir')

            strm1_z_seg = ss.detrend(strm1_z_seg,
                                     type='linear')
            strm2_z_cnv = ss.detrend(strm2_z_cnv,
                                     type='linear')

            strm1_z_seg = ops.bandpass(strm1_z_seg,
                                       self._config['analysis_bandpass'][0],
                                       self._config['analysis_bandpass'][1],
                                       self._config['analysis_sample_rate'])
            strm2_z_cnv = ops.bandpass(strm2_z_cnv,
                                       self._config['analysis_bandpass'][0],
                                       self._config['analysis_bandpass'][1],
                                       self._config['analysis_sample_rate'])

            amp_ratio = std(strm2_z.cnv) / std(strm1_z_seg)

            print('{}: rel amp: {}'.format(start_t, amp_ratio))

            start_t += self.segment_size_secs
            dbg_cnt += 1

            if dbg_cnt > 5:
                break

            # GET DATA ARRAYS, conv/deconv strm1 resp from strm2,
            # rfft/irfft


            # CHANGE BELOW TO USE SCIPY FUNCTIONS
            # strm1_z_seg.decimate(round(strm1_z_seg.stats.sampling_rate / self._config['analysis_sample_rate']), no_filter=True)
            # strm2_z_seg.decimate(round(strm2_z_seg.stats.sampling_rate / self._config['analysis_sample_rate']), no_filter=True)
            # strm1_1_seg = strm1_1.slice(starttime=start_t, endtime=start_t + self.segment_size_secs).data.astype(float32)
            # strm1_2_seg = strm1_2.slice(starttime=start_t, endtime=start_t + self.segment_size_secs).data.astype(float32)
            # strm2_1_seg = strm2_1.slice(starttime=start_t, endtime=start_t + self.segment_size_secs).data.astype(float32)
            # strm2_2_seg = strm2_2.slice(starttime=start_t, endtime=start_t + self.segment_size_secs).data.astype(float32)


    def compare_traces(self, ref_tr, sen_tr, segment_size_samples=1024*40):
        """These """

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

