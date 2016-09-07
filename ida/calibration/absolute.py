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
from pathlib import Path, PurePath
import yaml
import functools
#import collections

from fabulous.color import red, bold
from obspy import read, Stream, Trace, UTCDateTime

from ida.utils import i10get
from ida.calibration.shaketable import rename_chan


class AbsOnsiteConfig(object):

    def __init__(self, fn):

        self.errs = []

        # in-memory streams of source miniseed data
        self._ref_azi_strm = None
        self._ref_abs_strm = None
        self._pri_azi_strm = None
        self._pri_abs_strm = None
        self._sec_azi_strm = None
        self._sec_abs_strm = None
        self._ref_azi_strm_proc = None
        self._ref_abs_strm_proc = None
        self._pri_azi_strm_proc = None
        self._pri_abs_strm_proc = None
        self._sec_azi_strm_proc = None
        self._sec_abs_strm_proc = None

        # in-memory streams of miniseed data prepped for analysis


        with open(fn, 'rt') as cfl:
            config_txt = cfl.read()

        try:
            self._config = yaml.load(config_txt)
            self._config_dir = Path(fn).parent
        except:
            print(red(bold('Error parsing YAML config file: ' + fn)))
            self.errs.append('Error parsing YAML config file: ' + fn)
        else:
            self.process_config()


    def process_config(self):

        # check ENV
        if not os.environ.get('IDA_CAL_SITEDATA_DIR'):
            self.errs.append('The env var IDA_CAL_SITEDATA_DIR must be set to the root directory of the onsite reference data.')
        if not os.environ.get('IDA_ARCHIVE_RAW_DIR'):
            self.errs.append('The env var IDA_ARCHIVE_RAW_DIR must be set to the root directory IDA10 waveform data.')

        # check azimuth settings
        if self._config['process_azimuth'] == 1:

            if not PurePath(self._config['azimuth_ref_data']['ms_file']).is_absolute():
                fpath = os.path.join(os.environ.get('IDA_CAL_SITEDATA_DIR'),
                                     self._config['azimuth_ref_data']['ms_file'])
                self._config['azimuth_ref_data']['ms_file'] = fpath
            else:
                fpath = self._config['azimuth_ref_data']['ms_file']

            if not (os.path.exists(fpath) and os.path.isfile(fpath)):
                self.errs.append('Reference azimuth file not found: {}'.format(fpath))

            try:
                self._config['azimuth_ref_data']['starttime_iso'] = \
                    UTCDateTime(self._config['azimuth_ref_data']['starttime_iso'])
            except:
                self.errs.append('Error parsing azimuth_ref_data starttime_iso: {}'.format(
                    self._config['azimuth_ref_data']['starttime_iso']
                ))
            try:
                self._config['azimuth_ref_data']['endtime_iso'] = \
                    UTCDateTime(self._config['azimuth_ref_data']['endtime_iso'])
            except:
                self.errs.append('Error parsing azimuth_ref_data endtime_iso: {}'.format(
                    self._config['azimuth_ref_data']['endtime_iso']
                ))

        # check absolute settings
        if self._config['process_absolute'] == 1:

            if not PurePath(self._config['absolute_ref_data']['ms_file']).is_absolute():
                fpath = os.path.join(os.environ.get('IDA_CAL_SITEDATA_DIR'),
                                     self._config['absolute_ref_data']['ms_file'])
                self._config['absolute_ref_data']['ms_file'] = fpath
            else:
                fpath = self._config['absolute_ref_data']['ms_file']

            if not (os.path.exists(fpath) and os.path.isfile(fpath)):
                self.errs.append('Reference absolute file not found: {}'.format(fpath))

            try:
                self._config['absolute_ref_data']['starttime_iso'] = \
                    UTCDateTime(self._config['absolute_ref_data']['starttime_iso'])
            except:
                self.errs.append('Error parsing azimuth_ref_data starttime_iso: {}'.format(
                    self._config['absolute_ref_data']['starttime_iso']
                ))
            try:
                self._config['absolute_ref_data']['endtime_iso'] = \
                    UTCDateTime(self._config['absolute_ref_data']['endtime_iso'])
            except:
                self.errs.append('Error parsing azimuth_ref_data endtime_iso: {}'.format(
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
    def segment_size_samples(self):
        return self._config['segment_size_samples']

    @property
    def minimum_segment_cnt(self):
        return self._config['minimum_segment_cnt']

    def read_ref_data(self):
        """
        Reads AZI and ABS miniseed data for ref sensor.

        Performs trim (in-place) based on start/end times in config and saves pointer to streams

        """
        if self.ref_azi_seed_file:
            try:
                self._ref_azi_strm = read(self.ref_azi_seed_file, format='MSEED')
            except:
                print(red(bold('Error reading reference azimuth data: ' + self.ref_azi_seed_file)))

        if self.ref_abs_seed_file:
            try:
                self._ref_abs_strm = read(self.ref_abs_seed_file, format='MSEED')
            except:
                print(red(bold('Error reading reference absolute data: ' + self.ref_abs_seed_file)))

    def clean_ref_channels(self):
        """
        Cleanup ref data channel codes not configured properly on digitizer
        """

        if self._ref_azi_strm:
            for tr in self._ref_azi_strm:
                tr.stats.channel = rename_chan(tr.stats.channel)
        if self._ref_abs_strm:
            for tr in self._ref_abs_strm:
                tr.stats.channel = rename_chan(tr.stats.channel)


    def trim_data(self):
        """
        Trim all existing streams to azi/abs start/end times
        """

        if self._ref_azi_strm: self._ref_azi_strm.trim(self.azi_starttime, self.azi_endtime)
        if self._pri_azi_strm: self._pri_azi_strm.trim(self.azi_starttime, self.azi_endtime)
        if self._sec_azi_strm: self._sec_azi_strm.trim(self.azi_starttime, self.azi_endtime)
        if self._ref_abs_strm: self._ref_abs_strm.trim(self.abs_starttime, self.abs_endtime)
        if self._pri_abs_strm: self._pri_abs_strm.trim(self.abs_starttime, self.abs_endtime)
        if self._sec_abs_strm: self._sec_abs_strm.trim(self.abs_starttime, self.abs_endtime)


    def get_sensor_data(self):
        """
        Retrieve IDA10 data from IDA Archive and comvert to miniseed for azi/abs time periods
        """

        if self._pri_azi_strm:
            print('Retrieving primary sensor data for azimuth period...')
            i10get(self._pri_azi_strm[0].stats.station,
                   ','.join(self._pri_chanloc_codes()),
                   self.azi_starttime, self.azi_endtime,
                   outfn='./{}_{}_{}.ms'.format(self._pri_azi_strm[0].stats.station,
                                                self._config['site_settings']['pri_sensor_loc'],
                                                'azi.i10'))


    def _pri_chanloc_codes(self):
        if self._pri_azi_strm:
            return [tr.stats.channel+self._config['site_settings']['pri_sensor_loc'] for tr in self._pri_azi_strm]
        else:
            return ''


    def compare_streams(self, ref_st, sen_st):
        """
        Loops through component traces and runs azi and abs analysis on segments

        Processes each pair of traces:
            trim sample count to multiple of analysis frequency rate
            decimate to analysis sample rate
            bandpass filter to analysis bandpass
            deconvolve inst resp and convolve ref response to sens_st
            loop over segments performing analysis and compute
            angle, amp, and variance measures


        """
        pass

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

