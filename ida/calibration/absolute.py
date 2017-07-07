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
from pathlib import Path
import yaml
from collections import namedtuple
import logging

# import matplotlib
# matplotlib.use('Qt4Agg')
# import matplotlib.pyplot as plt
# plt.ion()

from obspy import read, UTCDateTime, Stream
from obspy.signal.invsim import evalresp
import obspy.signal.filter as osf
from scipy.signal import tukey
from numpy import array, ones, pi, arctan2, float64, multiply, divide, log10
from numpy import std, sqrt, dot, insert, inf
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
    """
    Class representing the results for a single component (Z12).
    It accumulates individual APSurveySegmentResult objects that each hold 
    for an individual segment
    """

    def __init__(self, comp):
        """
        Constructor for APSurveyComponentResult
        
        Args:
            comp (str): 'Z', '1', or '2' 
        """
        self.component = comp
        self.seg_results = []
        self._amp_mean = None  # mean of segments ampl ratios
        self._amp_std = None  # stdev of segments ampl ratios
        self._ang_mean = None  # mean of segments relative angles
        self._ang_std = None  # stdev of segments relative angles
        self._lrms_mean = None
        self._var_mean = None
        self._recalc_needed = False
        self._usable_count = 0

    def add_segment(self, apssegres):
        """
        Adds an individual segment result object to accumulating component result
        and keeps track of number of 'usable' segments
        
        Args:
            apssegres (APSurveySegmentResult): Segment result object to add  

        Returns: None

        """
        self.seg_results.append(apssegres)
        self._recalc_needed = True
        if apssegres.can_use:
            self._usable_count += 1

    def _recalc(self):
        """
        Recalculate overall compionent results from individual segment results.
        
        Returns: None

        """
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
            self._lrms_mean = None
            self._var_mean = None
        self._recalc_needed = False

    @property
    def amp_mean(self):
        """
        Property returning the overall mean relative amplitude of the individual segments for this component
        
        Returns: (float) Mean relative amplitude 

        """
        if self._recalc_needed:
            self._recalc()
        return self._amp_mean

    @property
    def amp_std(self):
        """
        Property returning the standard deviation of the relative amplitude of the 
        individual segments for this component
        
        Returns: (float) STD deviations of the segment relative amplitude values 

        """
        if self._recalc_needed:
            self._recalc()
        return self._amp_std

    @property
    def ang_mean(self):
        """
        Property returning the overall mean relative angle of the individual segments for this component
        WITH RESPECT to the NORTH component of the baseline sensor

        Returns: (float) Mean relative angle 

        """
        if self._recalc_needed:
            self._recalc()
        return self._ang_mean

    @property
    def ang_std(self):
        """
        Property returning the standard deviation of the relative angle of the 
        individual segments for this component WITH RESPECT to the NORTH component of the baseline sensor

        Returns: (float) STD deviations of the segment relative angle values 

        """
        if self._recalc_needed:
            self._recalc()
        return self._ang_std

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
        self._ang_resid = ang_resid  # currently not used for overall component results
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
    """
    Performs relative azimuth and sensitivity calculations for two or more sensors in pairs. Structured on 
    IDA historical approach that can analysze both 2 sensors for a station at one time wrt to one or two sets 
    of reference sensor data obtained with the Azimuth Pointing System Kit, Compact Trillium portable sensor and 
    Taurus and Centaur portable digitizers.
    
    In addition to measuring azimuth and sensitivity of each station sensor wrt to the reference sensor data,
    the second (loc=10) sensor is also compared against the first (loc=00) sensor. The relative reslts can then
    be examined to see that the results of the three pairs fo sensors are internally consistent.
    
    A typical IDA APS calibration in the field will opbtain data for two lengthy time periods (15+ hours). 
        1) 'Azimuth' data, with the reference sensor in a known orientation (from APS GPS infomation. This generally can only 
    be done in open air, and does not permit the most accurate sensitivity comparisons between the reference sensor 
    and the station sensors. 
        2) 'Absolute' data is acquired with the reference sensor very near the station sensors,
    on the same pier, for example. In this instance, the azimuth and clock on the reference sensor will not be 
    very well known, but the location allows for a better comparison of sensitivities. 
    
    When both sets of data can be obtained, the first is used to compute the relative azimuth of the 
    station sensors with the reference sensor, adn teh second is ude to compute the relative sensitivities.
    
    The analysis process is driven by a YAML configuration file that is supplied to the APSurvey constructor. 
    The configuration file contains the following information (represented with Python data structures):
    
         {
         # segment size. The entire timeseries is split into segments of this length. Analysis is performed on
         each one and the results aggregated for an overall result,
         'correlation_segment_size_secs': 7200,

         # filtering bandpass used in analysis of each segment
         'analysis_bandpass_hz': [0.1, 0.3],

         # sample rate that both segments are decimated to before analysis
         'analysis_sample_rate_hz': 5,

         # size of individual segments 
         'segment_size_secs': 1024,

         # trim/taper size used at each end when deconvolving/convolving responses.
         'segment_size_trim_secs': 128,

         # coherence cutoff. If the pari of segments don't have a coherence of this value or greater 
         # they are not used in the overall calculation for htat component. Start at 0.99 and go as low as 0.95
         # for noisy data.
         'coherence_cutoff': 0.99,

          # station sensor meta data used to retrieve data from IDA archives
          # and to find the correct sensor RESP file for convolution/deconvolution
         'station': 'DGAR'
         'pri_sensor_installed': True,
         # channels should be in Z, 1, 2 order
         'pri_sensor_chans': 'BHZ,BH1,BH2',
         'pri_sensor_loc': '00',
         'pri_sensor_installed': True,
         'sec_sensor_chans': 'BHZ,BH1,BH2',
         'sec_sensor_loc': '10',

         # location and timeseries start/end info. Time ar used when retrieving station sensor data from IDA archives.
         # the 'process' flag indicates whether to include the dataset in the analysis. One of the two must always 
         # be set to True
         'ref_absolute_data': {'endtime_iso': '',
                               'ms_file': '',
                               'process': False,
                               'starttime_iso': ''},
         'ref_azimuth_data': {'endtime_iso': '2016-12-03 02:00:00',
                              'ms_file': '/ida/cal/raw/dgar/apsurvey/2016-12-02-azi2/STN01_centaur-3_0831_20161202_030400.seed',
                              'process': True,
                              'starttime_iso': '2016-12-02 05:15:00'},

         # to override Trillium reference sensor timeseries metadata if metadata does not confirm to IDA conventions.
         'ref_kit_metadata': {'location': '01', 'network': 'II', 'station': 'TRI'},

         # location of root of IDA IDA10 archives
         'arc_raw_dir': '/ida/archive/raw',
         # NOTE: This may also be a MINISEED file, if a non-IDA sensor is being analyzed. 
         # in this case, the metadata from the 'pri_sensor' above is used in output of results'

         # location of RESP files used for comuting sensors responses
         'resp_file_dir': '/ida/dcc/response/RESP',
         
         }
                
    """

    ChanTpl = namedtuple('ChanTuple', 'z n e')

    def __init__(self, fn, debug=False):

        self.analdate = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        self.ida_cal_raw_dir = os.environ.get('IDA_CAL_RAW_DIR', '')
        self.debug = debug

        # keep a list of waveform files. Include i10 and ms data, just to show user at end
        self.waveform_files = []

        self.ok = True

        self.logger = logging.getLogger(__name__)

        # in-memory streams
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
        # source miniseed data
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
        # 3-component ChanTpl's. 'abs' or 'azi' may be left empty if not running both analyses
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
            # responses for 'ref' sensor adjusted for dif sampling rates of sta sensors
            # these are used to convolve reference response into station sensor timeseries
            'ref': {
                'pri': None,
                'sec': None,
            },
            'sec': {
                # response for 'sec' sensor (if it exists) adjusted for dif sampling rates of pri sensor
                # used to convolve sec response into pri timeseries when comparing pri sensors against sec.
                'pri': None,

                # regular response for sec station sensor used to deconvolve sta sensor response
                'sec': None
            },
            'pri': {
                # regular response for sec station sensor used to deconvolve sta sensor response
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

        # will be a ChanTpl containing a APSurveyComponentResults for each Z12 component
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
        """
        Checks for existance of keys and some existance of some file/directory values 
        in parsed yaml config file.
        
        Sets self.ok flag to False if it runs into any problems. 
        Flag needs to be check from calling code
        
        Returns: None

        """

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
            if not (os.path.exists(fpath)) and (',' not in fpath):
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

    # TODO: if not using fabulous.color, don't really need this abstraction
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

    def starttime_datetime(self, dataset):
        return self.starttime(dataset).datetime

    def endtime_datetime(self, dataset):
        return self.endtime(dataset).datetime

    def _respfilename(self, net, sta, chn, loc):
        """
        Construct the absolute or relative path of RESP file for given net, sta, chn, loc 
        using RESP directory supplied in config variable resp_file_dir
        
        Args:
            net (str): Network code 
            sta (str): Station code
            chn (str): Channel code
            loc (str): Location code

        Returns:
            (str): Absolute path of RESP file

        """
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

    def dataset_enabled(self, dataset):
        """
        Called externally to make sure dataset is enabled in config file before try calling self.analyze(dataset)

        Args:
            dataset (str): 'azi' or 'abs' 

        Returns:
            (bool): True if data set 'dataset' is enabled in the configuration. False otherwise.

        """
        if dataset.lower() == 'azi':
            return self.process_azimuth
        elif dataset.lower() == 'abs':
            return self.process_absolute
        else:
            return False

    def _correct_ref_time(self, dataset, sens1, sens2):
        """
        Cross correlate two timeseries to indentify any clock offset using the 'Z' copmonent data 
        from each sensor. It takes a segment of the timeseries that is correlation_segment_size long 
        (specified in config) and in the center of the timeseries.
        
        Teh two timeseries are bandpass filter from config bp_start (default is 0.1hz) to 2hz.
        
        The correlation calculation is performed by the obspy.signal.xcorr function underneath.
        
        Adjusts the starttime of the traces from sens1 to remove clock discrepancy. Positive offset 
        adjustment value indicates that the sens1 clock is slow (behind) sens2. Negative indicates the opposite. 
        
        We adjust sens1 because when comparing with a true reference sensor timeseries the clock 
        may have been set by hand, not GPS. So we assume the sta sensor clocks have better time 
        using the Q330 and GPS.
        
        Args:
            dataset (str): 
            sens1 (str): should be 'ref' 
            sens2 (str): 'pri' or 'sec' 

        Returns:
            (float, float, np.array, str): (Time offset in seconds, 
                                            max correlation value,
                                            correlation function values
                                            error msg returned from underlying routines
                                            )
        """

        sens1 = sens1.lower()
        sens2 = sens2.lower()

        if sens2 not in ['pri', 'sec']:
            raise ValueError('_correct_ref_time: sens2 must be "pri" or "sec".')

        if sens1 not in ['ref']:
            raise ValueError('_correct_ref_time: sens1 must be "ref".')

        if (sens2 == 'pri') and not self.pri_sensor_installed:
            raise ValueError('Primary sensor processing not enabled. ' +
                             ' Can not use it to correct reference clock. ' +
                             ' Check configuration.')
        if (sens2 == 'sec') and not self.sec_sensor_installed:
            raise ValueError('Secondary sensor processing not enabled. ' +
                             ' Can not use it to correct reference clock. ' +
                             ' Check configuration.')
        if (dataset == 'azi') and not self.process_azimuth:
            raise ValueError('Azimuth data processing not enabled. ' +
                             ' Can not use it to correct reference clock. ' +
                             ' Check configuration.')
        if (dataset == 'abs') and not self.process_absolute:
            raise ValueError('Absolute data processing not enabled. ' +
                             ' Can not use it to correct reference clock. ' +
                             ' Check configuration.')

        offset = 0
        cval = 0.0
        cfunc = []
        emsg = ''

        ref_z = self.trtpls[dataset][sens1].z.copy()
        sensor_z = self.trtpls[dataset][sens2].z.copy()

        # need to interpolate if ref data sampling rate > sensor sampling rate
        if ref_z.stats.sampling_rate > sensor_z.stats.sampling_rate:
            sensor_z.interpolate(ref_z.stats.sampling_rate, method='linear')

        if sensor_z and ref_z:
            ref_z.filter('bandpass', freqmin=self.bp_start, freqmax=2.0)
            sensor_z.filter('bandpass', freqmin=self.bp_start, freqmax=2.0)

            # take middle for time series correlation
            dur = sensor_z.stats.endtime - sensor_z.stats.starttime
            start_t = sensor_z.stats.starttime + dur/2 - self.correlation_segment_size/2
            end_t = start_t + self.correlation_segment_size
            ref_z.trim(start_t, end_t)
            sensor_z.trim(start_t, end_t)

            # compute time offset and adjust starttime of sens1 traces
            # use window isze of 1 minute at 40hz
            offset, cval, cfunc, emsg = time_offset(sensor_z, ref_z, 
                                                    winsize=2400)
            for tr in self.streams[dataset][sens1]:
                tr.stats.starttime = tr.stats.starttime + offset

            # log a messages if max correlation function value is < 0.9
            if cval < 0.9:
                self.logmsg(logging.WARN, 'WARNING (time_offset): Correlation value < 0.9: ' + str(cval))

        return offset, cval, cfunc, emsg

    def _read_ref_data(self, dataset):
        """
        Read reference sensor miniseed data for given dataset, 
        and store ChanTpl of 3 component traces.
        
        Channels with old-style IDA reference sensor codes are reanmed to Z12.
        Channels codes can be overriden with reference sensors metadata supplied 
        in the config file, if it is set.
        
        Ultimately, the refernce sensor channel codes MUST end in Z, 1 and 2 to succeed.
        
        Args:
            dataset (str): 'azi' or 'abs'

        Returns:
            (bool): Success flag. True if Reference data read successfully, False if not.

        """

        if (dataset == 'abs') and not self.process_absolute:
                self.logmsg(logging.WARN, 'Processing Absolute data is not enabled.')
                return False

        if (dataset == 'azi') and not self.process_azimuth:
                self.logmsg(logging.WARN, 'Processing Azimuth data is not enabled.')
                return False

        self.logmsg(logging.INFO, 'Reading REF sensor {} data...'.format(dataset.upper()))

        try:
            self.streams[dataset]['ref']= read(self.msfiles[dataset]['ref'], format='MSEED')
        except Exception as e:
            self.logmsg(logging.ERROR, 'Error reading reference data: ' + self.msfiles[dataset]['ref'])
            self.logmsg(logging.ERROR, 'Can not use REFERENCE data in sensor comparisons.')
            print()
            print(e)
            print()
            return False

        self.streams[dataset]['ref'].trim(starttime=self.starttime(dataset),
                                          endtime=self.endtime(dataset))

        gaps = self.streams[dataset]['ref'].get_gaps()
        if gaps:
            self.logmsg(logging.ERROR, 'REFERENCE data has gaps within time span indicated.')
            self.streams[dataset]['ref'].print_gaps()
            self.logmsg(logging.ERROR, 'Can not use REFERENCE data in sensor comparisons.')
            self.streams[dataset]['ref'] = Stream()
            return False

        # merge multiple traces together in one contiguous trace
        self.streams[dataset]['ref'].merge()

        for tr in self.streams[dataset]['ref']:
            # rename refernece sensors traces with old-style IDA channel codes
            tr.stats.channel = rename_chan(tr.stats.channel)

            # override channel codes with those in config file, if they are set.
            if self._config['ref_kit_metadata']['network'].strip():
                tr.stats.network = self._config['ref_kit_metadata']['network']
            if self._config['ref_kit_metadata']['station'].strip():
                tr.stats.station = self._config['ref_kit_metadata']['station']
            if self._config['ref_kit_metadata']['location'].strip():
                tr.stats.location = self._config['ref_kit_metadata']['location']
        try:
            # split out traces based on channel Z12 codes
            tr_z = self.streams[dataset]['ref'].select(component='Z')[0]
            tr_1 = self.streams[dataset]['ref'].select(component='1')[0]
            tr_2 = self.streams[dataset]['ref'].select(component='2')[0]
        except Exception as e:
            self.logmsg(logging.ERROR, 'Z12 components not found in reference data.')
            print(self.streams[dataset]['ref'])
            self.logmsg(logging.ERROR, 'Can not use REFERENCE data in sensor comparisons.')
            self.streams[dataset]['ref'] = Stream()
            return False

        # finally, store ref sensor traces for this dataset
        self.trtpls[dataset]['ref'] = self.ChanTpl(z=tr_z, n=tr_1, e=tr_2)

        return True

    def _read_sensor_data(self, dataset, sensor):
        """
        Retrieve IDA10 data from IDA Archive and convert to miniseed for indicated
        azi or abs time period. 
        
        ALTERNATELY, if the ARC_RAW_DIR config setting is a file, it is ASSUMED to be a miniseed file.
        This is ONLY to run single analysis (not both 'azi' and 'abs') on an existing miniseed file
        instead of pulling data from IDA10 archive. IN this miniseed file there MUST be only 3 traces
        with Z12 channel codes.
        
        Args:
            dataset (str): 'azi' or 'abs' to indicate which time period to read from archive 
            sensor (str): 'pri' or 'sec' to indicate which sensor's data to read form archive

        Returns: 
            (bool): Success status: True sensor data read, or False if not.
        
        """

        dataset = dataset.lower()
        sensor = sensor.lower()

        if dataset not in ['azi', 'abs']:
            raise ValueError('_read_sensor_data: dataset must be "azi" or "abs".')

        if sensor not in ['pri', 'sec']:
            raise ValueError('_read_sensor_data: sensor must be "pri" or "sec".')

        # see if data previously read in.
        # For cases when running both 'azi' and 'abs' in one pass
        if self.streams[dataset][sensor]:
            self.logmsg(logging.INFO, 'Using {} sensor {} data already in memory'.format(
                sensor.upper(), dataset.upper()))
            return True

        # user error here... command line conflicts with config file
        # TODO: better logging/feedback
        if (sensor == 'pri') and not self.pri_sensor_installed:
            return False
        if (sensor == 'sec') and not self.sec_sensor_installed:
            return False
        if dataset == 'azi' and not self.process_azimuth:
            return False
        if dataset == 'abs' and not self.process_absolute:
            return False

        # check for non IDA sensor data supplied as miniseed file.
        # TODO: this is a kludge and this situation needs to be reworked
        # if os.path.isfile(self.arc_raw_dir):
        if not os.path.isdir(self.arc_raw_dir):
            # assume pri, and possibly sec,  miniseed file names with ONLY needed 
            # channels (Z12/ZNE) and time period

            fns = self.arc_raw_dir.split(',')
            fns = [fn.strip() for fn in fns]
            if (sensor == 'pri') and (len(fns) > 0) and os.path.isfile(fns[0]):
                ms_name = fns[0]
                self.waveform_files.append(fns[0])
            elif (sensor == 'sec') and (len(fns) > 1) and os.path.isfile(fns[1]):
                ms_name = fns[1]
                self.waveform_files.append(fns[1])
            else:
                self.logmsg(logging.ERROR, 
                            'Error finding ms file [{}] for {} sensor.'.format(
                                self.arc_raw_dir, sensor
                            ))

                return False
        else:
            # retrieve data from IDA IDA10 archive

            # construct root of output file names
            outname = './{}_{}_{}'.format(self.station,
                                          self.station_sensor_loc(sensor),
                                          dataset)
            ms_name = outname + '.ms'
            i10_name = outname + '.i10'
            if os.path.exists(i10_name):
                os.remove(i10_name)
            if os.path.exists(ms_name):
                os.remove(ms_name)

            self.waveform_files.append(i10_name)
            self.waveform_files.append(ms_name)

            self.logmsg(logging.INFO, 'Retrieving {} sensor {} data from archive' \
                        ' and saving in {}'.format(sensor.upper(),
                                                   dataset.upper(),
                                                   ms_name))

            try:
                # get IDA10 data and convert to miniseed
                i10get(self.arc_raw_dir,
                       self.station,
                       self._chanloc_codes(sensor),
                       self.starttime(dataset), self.endtime(dataset),
                       outfn=i10_name)
                pimseed(self.station, i10_name, ms_name)
            except:
                self.logmsg(logging.ERROR,
                            'Error reading and converting IDA10 data for sensor ({}/{})'.format(dataset, sensor))
                self.logmsg(logging.ERROR,
                            'Can not use {} sensor data in sensor comparisons.'.format(sensor))
                return False

        try:
            # TODO: Split up, too much in this block
            # read ms file, merge traces, check for gaps and split trces into ChanTpl for processing
            self.streams[dataset][sensor] = read(ms_name)
            self.streams[dataset][sensor].trim(starttime=self.starttime(dataset),
                                               endtime=self.endtime(dataset))
            self.streams[dataset][sensor].merge()
            gaps = self.streams[dataset][sensor].get_gaps()
            if gaps:
                self.logmsg(logging.ERROR,
                            '{} Sensor data has gaps'.format(sensor))
                self.streams[dataset][sensor].print_gaps()
                self.logmsg(logging.ERROR,
                            'Can not use {} sensor data in sensor comparisons.'.format(sensor))
                return False
            tr_z = self.streams[dataset][sensor].select(component='Z')[0]
            tr_1 = self.streams[dataset][sensor].select(component='1')
            if not tr_1:
                tr_1 = self.streams[dataset][sensor].select(component='N')
            if tr_1:
                tr_1 = tr_1[0]
            else:
                self.logmsg(logging.ERROR, 'Could not find NORTH component trace in miniseed data for sensor {}'.format(sensor))

            tr_2 = self.streams[dataset][sensor].select(component='2')
            if not tr_2:
                tr_2 = self.streams[dataset][sensor].select(component='E')
            if tr_2:
                tr_2 = tr_2[0]
            else:
                self.logmsg(logging.ERROR, 'Could not find EAST component trace in miniseed data for sensor {}'.format(sensor))

            self.trtpls[dataset][sensor] = self.ChanTpl(z=tr_z, n=tr_1, e=tr_2)
            self.msfiles[dataset][sensor] = os.path.abspath(ms_name)
            self.streams[dataset][sensor].write(self.msfiles[dataset][sensor], format='MSEED')
        except:
            self.logmsg(logging.ERROR,
                        'Error reading processing miniseed data for sensor ({}/{})'.format(dataset, sensor))
            self.logmsg(logging.ERROR,
                        'Can not use {} sensor data in sensor comparisons.'.format(sensor))
            return False

        return True

    def _read_responses(self, sens1, sens2):
        """
        Computes a pair of responses, for each component, for comparing a specific pair of sensors.

        One is the response for sens1, but computed with the sample rate of sens2. This is used to convolve 
        the sen1 response into the sens2 timeseries.

        The other is the response of sens, at it's normal operating sample rate. This is used to deconvolving
        sens2 response from it's timeseries


        sens1 is considered the 'reference' sensor. and sens2 timeseries is

        For each component:
            The


        azimuth and abs sensitivity of sens2 wrt sens1.

        Args:
            sens1 (str): 'ref' or 'sec'. These are the two sensors that can be used as baseline 
            sens2 (str):  'pri' or 'sec'. These are the two sensors that are compared against a baseline sensor 

        Returns:
            (bool, bool): Success flags. (Ref response OK, Sensor response OK)

        """

        sens1 = sens1.lower()
        sens2 = sens2.lower()

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

        # compute responses for sens2. For deconcolving response later...
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
                # exponentially taper off highest 5% of freqs before nyquist
                # taper off high freq response before deconvolving sens2 resp
                # per idaresponse/resp.c
                taper_high_freq_resp(sen_z, 0.95)
                taper_high_freq_resp(sen_1, 0.95)
                taper_high_freq_resp(sen_2, 0.95)
                # clip minimum amp resp to 1/100.0 of max amp
                # per idaresponse/resp.c
                dynlimit_resp_min(sen_z, 100.0)
                dynlimit_resp_min(sen_1, 100.0)
                dynlimit_resp_min(sen_2, 100.0)
                # note this is the regular responses for either 'pri' or 'sec' sensors
                self.responses[sens2][sens2] = self.ChanTpl(z=sen_z, n=sen_1, e=sen_2)
        else:
            sen_ok = True
            self.logmsg(logging.INFO, 'Using {} response already in memory...'.format(sens2.upper()))

        # compute responses for sens1 with sens2 sample rate, for convolving sens1 response into sens2, later...
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
                # clip minimum amp resp to 1/100.0 of max amp
                # per idaresponse/resp.c
                # TODO: why idaresponse/resp.c has no high freq taper here...?
                dynlimit_resp_min(ref_z, 100.0)
                dynlimit_resp_min(ref_1, 100.0)
                dynlimit_resp_min(ref_2, 100.0)
                self.responses[sens1][sens2] = self.ChanTpl(z=ref_z, n=ref_1, e=ref_2)

        else:
            ref_ok = True
            self.logmsg(logging.INFO, 'Using {} response already in memory...'.format(sens1.upper()))

        return ref_ok, sen_ok

    def _sensor_sample_rate_str(self, dataset, sensor):
        """
        Reads sample rate from the first Trace in dataset/sensor Stream and returns it formatted as a string
        
        Args:
            dataset (str): 'azi' or 'abs' 
            sensor (str): 'ref' or 'pri' or 'sec'

        Returns:
            (str): sample rate value represented as a string. Empty string if dataset/sensor does not exist

        """
        if self.streams[dataset][sensor]:
            return str(self.streams[dataset][sensor][0].stats.sampling_rate)
        else:
            return ''

    def _chanloc_codes(self, sensor):
        """
        
        Args:
            sensor (str): 'pri' or 'sec' 

        Returns:
            (str): Comma separated list of CHAN+LOC codes as set in config for 'sensor'. 
                    Returns empty string if not set in config or if invalid sensor requested 

        """

        if sensor == 'pri':
            chans = self._config.get('pri_sensor_chans', '').lower().split(',')
            loc = self._config.get('pri_sensor_loc', '').upper()
            if chans:
                return ','.join([chan+loc for chan in chans])
            else:
                return ''

        elif sensor == 'sec':
            chans = self._config.get('sec_sensor_chans', '').lower().split(',')
            loc = self._config.get('sec_sensor_loc', '').upper()
            if chans:
                return ','.join([chan+loc for chan in chans])
            else:
                return ''
        else:
            return ''

    def response_tr(self, respdate, npts, sr, net, sta, chn, loc, units='VEL'):
        """
        Obtain channel/loc frequency response using evalresp (from obspy) and RESP file on disk.
         
        Args:
            respdate (datetime or UTCDateTime): Time at for response should be returned 
            npts (int): NUmber of points in the response 
            sr (int or float): sample rate at which response should be computed 
            net (str): Network code
            sta (str): Station code
            chn (str): Channel code
            loc (str): Location code
            units (str): 'DIS', 'VEL', or 'ACC'

        Returns:
            (np.array, np.array, bool): 
                (Complex array with response, array with frequencies, True/False Success flag)

        """
        resp, f = None, None
        respfn = self._respfilename(net, sta, chn, loc)

        if os.path.exists(respfn):
            try:
                resp, f = evalresp(1 / sr,
                                   npts,
                                   respfn,
                                   UTCDateTime(respdate),
                                   station=sta,
                                   channel=chn,
                                   network=net,
                                   locid=loc,
                                   units=units, freq=True)
                ok = True
            except:
                self.logmsg(logging.ERROR,
                            'Error running EVALRESP with RESP file {}.'.format(respfn))
                ok = False
        else:
            self.logmsg(logging.ERROR, 'RESP file not found: ' + respfn)
            ok = False

        return resp, f, ok

    def ms_filename(self, dataset, sensor):
        """
        Args:
            dataset (str): 'azi' or 'abs' 
            sensor (str): 'pri' or 'sec'

        Returns:
            (str): miniseed filename for dataset and sensor

        """
        return Path(self.msfiles[dataset][sensor]).name


    def _get_result_headers(self, datatype):
        """
        Construct headers with analysis parameters for summary and detailed results files.
        
        Args:
            datatype (str): 'azi' or 'abs' 

        Returns:
            (str, str): (Summary Hdr, Detailed Hdr) 

        """

        header = '#'*144 + '\n'
        header += '# ANALYSIS PARAMETERS\n'
        header += '# =====================\n'
        header += '#   ida module version: {} ({})\n'.format(IDA_PKG_VERSION_HASH_STR,
                                                             IDA_PKG_VERSION_DATETIME)
        header += '#          analysis at: {}\n'.format(self.analdate)
        header += '#              dataset: {}\n'.format(datatype.upper())
        header += '#              station: {}\n'.format(self.station.upper())
        header += '#           pri sensor: {}\n'.format(self.pri_sensor_installed)
        header += '#     pri data sr (hz): {}\n'.format(self._sensor_sample_rate_str(datatype, 'pri'))
        header += '#            pri chans: {}\n'.format(self._chanloc_codes('pri').upper())
        header += '#           sec sensor: {}\n'.format(self.sec_sensor_installed)
        header += '#     sec data sr (hz): {}\n'.format(self._sensor_sample_rate_str(datatype, 'sec'))
        header += '#            sec chans: {}\n'.format(self._chanloc_codes('sec').upper())
        header += '#          ref ms file: {}\n'.format(self.msfiles[datatype]['ref'])
        header += '#     ref data sr (hz): {}\n'.format(self._sensor_sample_rate_str(datatype, 'ref'))
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

        return header + sumhdr, header + dethdr

    def _get_result_text(self, dataset, sens1, sens2, results):
        """
        Construct text for results summary and detail result files form 'results'
        Args:
            dataset (str): 'azi' or 'abs' 
            sens1 (str): 'ref' or 'sec'
            sens2 (str): 'pri' or 'sec'
            results (ChanTpl): ChanTpl with APSurveyComponentResults for each component

        Returns:
            (str, str): (Summary result text, Detailed result text). 
                        Returns empty strings if no results (should never happen)

        """

        sumres = ''
        detres = ''

        analday = datetime.now().strftime('%Y-%m-%d')
        if results:
            # loop through each component in results ChanTpl
            for compres in results:
                if compres.usable_count == 0:
                    # so sad...
                    sumres += '{:<4} {}/{} : {} No usable segments. {} {}'.format(
                        self.station.upper(), sens1.upper(), sens2.upper(),
                        compres.component, self.msfiles[dataset][sens1],
                        self.msfiles[dataset][sens2])
                else:
                    # construct summary result text record from component form aggregate calculations
                    sumres += '   {:<4} {:<4} {:<4} {:<4} '\
                             '{:7.4f} {:7.3f} '\
                             '{:8.3f} {:7.3f} {:7.3f} {:7.3f} '\
                             '{:8.3f} {:>6} {:>6} '\
                             '{:>10} {:<14} {:<14}'.format(
                                 self.station.upper(), sens1.upper(), sens2.upper(),
                                 compres.component, compres.amp_mean, compres.amp_std,
                                 (compres.ang_mean * 180./pi) % 360.,
                                 (compres.ang_std * 180./pi) % 360.,
                                 compres.lrms_mean, compres.var_mean,
                                 self.coherence_cutoff, compres.usable_count,
                                 compres.total_count, analday, self.ms_filename(dataset, sens1),
                                 self.ms_filename(dataset, sens2))
                self.logmsg(logging.INFO, sumres)
                sumres += '\n'

                # build detailed result text, one text line for the result for each segment
                for res in compres.seg_results:
                    excl = 'EXCL' if not res.can_use else 'Ok  '
                    detres += '   {:<4} {:<4} {:<4} {:<4} {} {:17.6f} {:7.4f} {:8.3f} '\
                            '{:7.3f} {:7.3f} {:7.4f} {:7.3f} {:<6} {:>10} {:<14} {:<14}'.format(
                        self.station.upper(), sens1.upper(), sens2.upper(), compres.component,
                        res.start_utc, res.start_epoch, res.amp, (res.ang * 180./pi + 360.) % 360.,
                        res.lrms, res.var, res.coh, self.coherence_cutoff, excl, analday,
                        self.ms_filename(dataset, sens1), self.ms_filename(dataset, sens2))
                    detres += '\n'

        return sumres, detres

    def analyze(self, dataset):
        """
        Called externally to perform the analysis and write out results for a given 'dataset'.
        
        The refernce data is read first. Then the secondary station sensor is read. The secondary sensor is read
        first because in most cases it will be sampling at the same rate (40hz) as the reference sensor and is
        therefor better to use to run a correlation to identify and time offset between the reference and 
        station sensors' timeseries
        
        Once any time offset is correction
        Args:
            dataset (str): 'azi' or 'abs 

        Returns:
            (str, str, list): Summary, detail and list of waveform files names, respectively.

        """

        dataset = dataset.lower()
        if dataset not in ['azi', 'abs']:
            raise ValueError('analyze: dataset must be "azi" or "abs".')

        self.waveform_files = []
        sensors_available = 0  # need 2+ to do any comparisons

        # read reference sensor data
        data_ok = self._read_ref_data(dataset)
        compare_ref = data_ok
        if not data_ok:
            self.logmsg(logging.WARN, 'Unable to process {} REFERENCE sensor data.'.format(
                dataset.upper()))
            self.logmsg(logging.WARN, 'Skipping comparisons with REFERENCE sensor')
        else:
            sensors_available += 1

        # read 'sec' sta sensor if enabled/active in config
        if self.sec_sensor_installed:
            data_ok = self._read_sensor_data(dataset, 'sec')
            compare_sec = data_ok
            if not data_ok:
                self.logmsg(logging.WARN, 'Unable to process {} SECONDARY sensor data.'.format(
                    dataset.upper()))
                self.logmsg(logging.WARN, 'Skipping comparisons with SECONDARY sensor')
            else:
                sensors_available += 1
                # compute and correct ref time offset
                offset, _, _, _ = self._correct_ref_time(dataset, 'ref', 'sec')
                # save offset info just to report out in results
                self.ref_clock_adjustment = offset
                self.ref_clock_adjusted = True
                self.ref_clock_adjustment_sensor = 'secondary'
        else:
            compare_sec = False  # not installed

        # read 'pri' sta sensor if enabled/active in config
        if self.pri_sensor_installed:
            data_ok = self._read_sensor_data(dataset, 'pri')
            compare_pri = data_ok
            if not data_ok:
                self.logmsg(logging.WARN, 'Unable to process {} PRIMARY sensor data.'.format(
                    dataset.upper()))
                self.logmsg(logging.WARN, 'Skipping comparisons with PRIMARY sensor')
            else:
                # compute time offset IFF not already done with 'sec' sensor
                sensors_available += 1
                if not self.ref_clock_adjusted:
                    # compute and correct ref time offset
                    offset, _, _, _ = self._correct_ref_time(dataset, 'ref', 'pri')
                    # save offset info just to report out in results
                    self.ref_clock_adjustment = offset
                    self.ref_clock_adjusted = True
                    self.ref_clock_adjustment_sensor = 'primary'
        else:
            compare_pri = False  # not installed

        # build result file names from 'dataset' and date at start of timeseries
        dataday = self.starttime_datetime(dataset).strftime('%Y-%m-%d')
        detfn = '{}_{}_{}_details.txt'.format(self.station.lower(), dataday, dataset)
        sumfn = '{}_{}_{}_summary.txt'.format(self.station.lower(), dataday, dataset)

        with open(sumfn, 'at') as sumf:
            with open(detfn, 'at') as detf:

                sumhdr, dethdr = self._get_result_headers(dataset)
                sumf.write(sumhdr)
                detf.write(dethdr)

                comparisons = []  # list of sensor paairs that can be compared
                if sensors_available >= 2:
                    if compare_ref and compare_sec:
                        comparisons.append(('ref', 'sec'))

                    if compare_ref and compare_pri:
                        comparisons.append(('ref', 'pri'))

                    if compare_pri and compare_sec:
                        comparisons.append(('sec', 'pri'))

                    # loop through pairs of sensors, comparing data streams and writing out results...
                    for sens1, sens2 in comparisons:
                        results = self._compare_streams(dataset, sens1, sens2)
                        if results:
                            sumres, detres = self._get_result_text(dataset, sens1, sens2, results)
                            sumf.write(sumres)
                            detf.write(detres)
                        else:
                            sumf.write('# ERROR: No comparison results returned!\n')
                            detf.write('# ERROR: No comparison results returned!\n')
                else:
                    detf.write('# No comparisons possible for {} data sets\n'.format(dataset.upper()))
                    sumf.write('# No comparisons possible for {} data sets\n'.format(dataset.upper()))

                # write footers to results files
                sumf.write('#'*144 + '\n')
                detf.write('#'*144 + '\n')

        return sumfn, detfn, self.waveform_files

    def _compare_streams(self, datatype, sens1, sens2):
        """
        Prepares for analysis of wach of the three components.
        Identifies proper respoonses for de/convolution

        Args:
            datatype (str): 'azi' or 'abs 
            sens1 (str): 'ref' or 'sec' 
            sens2 (str): 'pri' or 'sec'

        Returns: 
            (ChanTpl): Contains tuple of APSurveyComponentResults containing the results of 
            this sensor pair comparisons.

        """

        if sens1 not in ['ref', 'sec']:
            raise ValueError('_compare_streams: sens1 must be "ref" or "sec"')
        if sens2 not in ['pri', 'sec']:
            raise ValueError('_compare_streams: sens2 must be "pri" or "sec".')
        if sens1 == sens2:
            raise ValueError('_compare_streams: sens1 and sens2 must different sensors.')

        # maek sure responses for the the pari of sensors being compared have been loaded successfully
        sens1_resp_ok, sens2_resp_ok = self._read_responses(sens1, sens2)
        if not sens1_resp_ok:
            self.logmsg(logging.WARN, 'Unable to read {} response.'.format(sens1))
            self.logmsg(logging.WARN, 'Can not perform comparison: {}/{}'.format(
                sens1.upper(), sens2.upper()))
            return False
        elif not sens2_resp_ok:
            self.logmsg(logging.WARN, 'Unable to read {} response.'.format(sens2))
            self.logmsg(logging.WARN, 'Can not perform comparison: {}/{}'.format(
                sens1.upper(), sens2.upper()))
            return False

        # grab Traces for both sensors
        strm1_z = self.trtpls[datatype][sens1].z
        strm1_1 = self.trtpls[datatype][sens1].n
        strm1_2 = self.trtpls[datatype][sens1].e
        strm2_z = self.trtpls[datatype][sens2].z
        strm2_1 = self.trtpls[datatype][sens2].n
        strm2_2 = self.trtpls[datatype][sens2].e

        # get responnse of sens1 computed for SR of sens2
        # to convolve into timeseries of sens2
        strm1_z_resp = self.responses[sens1][sens2].z
        strm1_1_resp = self.responses[sens1][sens2].n
        strm1_2_resp = self.responses[sens1][sens2].e

        # get response of sens2 for deconvolving response
        strm2_z_resp = self.responses[sens2][sens2].z
        strm2_1_resp = self.responses[sens2][sens2].n
        strm2_2_resp = self.responses[sens2][sens2].e

        # create a set of empty component results
        self.results = self.ChanTpl(z=APSurveyComponentResult('Z'),
                                    n=APSurveyComponentResult('1'),
                                    e=APSurveyComponentResult('2'))

        # analyze components
        self._compare_verticals(strm1_z, strm2_z, strm1_z_resp, strm2_z_resp, self.results.z)
        self._compare_horizontals(strm1_1, strm1_2, strm2_1, strm1_1_resp, strm2_1_resp, self.results.n)
        self._compare_horizontals(strm1_1, strm1_2, strm2_2, strm1_2_resp, strm2_2_resp, self.results.e)

        return self.results

    def _compare_horizontals(self, tr1_n, tr1_e, tr2, tr1_resp, tr2_resp, results):
        """
        Compares tr2 from one sensor with the north and east (tr1_n, tr1_e) traces from the another sensor 
        to determine the relative orientation of tr2 WRT tr1_n and relative sensitivities.
        The traces is are analyzed in segments of length specified in the config.
        The current tr2 response is removed and the response from tr1 is applied

        Result is expressed as tr2 angle WRT to tr1_n. 

        Positive relative angles are clock-wise adn saved in radians 

        Args:
            tr1_n (Trace): sensor 1 ('baseline') north trace
            tr1_e (Trace): sensor 1 ('baseline') east trace
            tr2 (Trace): sensor 2 trace
            tr1_resp (np.array): sensor 1 response to apply to tr2
            tr2_resp (np.array): tr2 response to deconvolve
            results (APSurveyComponentResult): COmponent results to assumulating individual segments results

        Returns: None

        """

        start_t = max(tr1_n.stats.starttime, tr2.stats.starttime)
        end_t = min(tr1_n.stats.endtime, tr2.stats.endtime)
        tr1_sr = tr1_n.stats.sampling_rate
        tr2_sr = tr2.stats.sampling_rate
        tr1_time_delta = 1.0/tr1_sr
        tr2_time_delta = 1.0/tr2_sr

        # loop through traces segment by segment
        while start_t + self.segment_size_secs < end_t:

            tr1_n_seg = tr1_n.slice(starttime=start_t,
                                    endtime=start_t + self.segment_size_secs - tr1_time_delta)
            tr1_e_seg = tr1_e.slice(starttime=start_t,
                                    endtime=start_t + self.segment_size_secs - tr1_time_delta)
            tr2_seg = tr2.slice(starttime=start_t,
                                endtime=start_t + self.segment_size_secs - tr2_time_delta)

            if (round(tr1_n_seg.std()) > 0) and (round(tr1_e_seg.std()) > 0) and (round(tr2_seg.std()) > 0):

                # get raw float data
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

                # remove tr2 response from timeseries
                tr2_fft = rfft(multiply(tr2_seg, taper))
                tr2_fft_dcnv = divide(tr2_fft[1:], tr2_resp[1:])
                tr2_fft_dcnv = insert(tr2_fft_dcnv, 0, [0.0])
                tr2_dcnv = irfft(tr2_fft_dcnv)

                # demean again before convolving strm2 resp
                tr2_dcnv = ss.detrend(tr2_dcnv, type='constant')

                # convolve tr2 with tr1 response
                tr2_fft = rfft(multiply(tr2_dcnv, taper))
                tr2_fft_cnv = multiply(tr2_fft, tr1_resp)
                tr2_cnv = irfft(tr2_fft_cnv)
                del tr2_fft, tr2_fft_dcnv, tr2_dcnv, tr2_fft_cnv

                # trim off taper segment_size_trim in secs
                trim_cnt = round(tr1_sr * self.segment_size_trim)
                tr1_n_seg = tr1_n_seg[trim_cnt:-trim_cnt]
                tr1_e_seg = tr1_e_seg[trim_cnt:-trim_cnt]
                trim_cnt = round(tr2_sr * self.segment_size_trim)
                tr2_cnv = tr2_cnv[trim_cnt:-trim_cnt]

                # resample all 3 timeseries to analsysi sample rate
                tr1_n_seg = ss.resample(tr1_n_seg, round(len(tr1_n_seg) / (tr1_sr / self.analysis_sample_rate)))
                tr1_e_seg = ss.resample(tr1_e_seg, round(len(tr1_e_seg) / (tr1_sr / self.analysis_sample_rate)))
                tr2_cnv = ss.resample(tr2_cnv, round(len(tr2_cnv) / (tr2_sr / self.analysis_sample_rate)))

                # demean and detrend
                tr1_n_seg = ss.detrend(tr1_n_seg, type='linear')
                tr1_e_seg = ss.detrend(tr1_e_seg, type='linear')
                tr2_cnv = ss.detrend(tr2_cnv, type='linear')

                # apply taper before filtering
                fraction = 1/12
                taper = tukey(len(tr1_n_seg), fraction * 2, sym=True)
                tr1_n_seg *= taper
                tr1_e_seg *= taper
                tr2_cnv *= taper

                # apply band pass filter with stops from config
                tr1_n_seg = osf.bandpass(tr1_n_seg,
                                         self.bp_start, self.bp_stop,
                                         self.analysis_sample_rate, zerophase=True)
                tr1_e_seg = osf.bandpass(tr1_e_seg,
                                         self.bp_start, self.bp_stop,
                                         self.analysis_sample_rate, zerophase=True)
                tr2_cnv = osf.bandpass(tr2_cnv,
                                       self.bp_start, self.bp_stop,
                                       self.analysis_sample_rate, zerophase=True)


                dc = ones(len(tr1_n_seg), dtype=float64)  # to take care of any non-zero means

                # set up matrix of constant 1, and tr1 north and east traces
                mat = array([dc, tr1_n_seg, tr1_e_seg])
                matinv = mat.transpose()

                # find least squares solution for:
                #     tr2 = solution[2] * tr1_e + solution[1] * tr1_n + solution[0]
                # ration of solution[2] and solution[1] are effectively the tan of angle omega
                # between tr2 and tr1_n.
                solution, resid, _, _ = la.lstsq(matinv, tr2_cnv)
                w = arctan2(solution[2], solution[1])
                # the lsq solution coeeficients are used to determine whether the tr2 values are scaled
                # up or down from tr1. SInce:
                #
                #   tan(w) = solution[1] / solution [2]
                #
                # we can define a scaling factor f such that:
                #
                #   solution[1] = f * cos(w), and
                #   solution[2] = f * sin(w)
                #
                # then
                #
                #   solution[1]**2 + solution[2]**2 = f**2 * ( cos**2(w) + sin**2(w) )
                #
                # simplifying & solving for f
                #
                #   f = sqrt(solution[1]**2 + solution[2]**2)
                #
                #
                # This factor represents the factor by which the current response for tr2 (sensor 2) is off.
                # When factor is > 1.0, the current response is underestimating actual sensitivity by this factor.
                # When factor is < 1.0, the current response is overestimating actual sensitivity by this factor.
                #
                # whew.
                #
                amp_ratio =  sqrt(solution[1]*solution[1] + solution[2]*solution[2])

                # log root mean square of the tr2 timseries amplitude
                lrms = log10(sqrt(multiply(tr2_cnv, tr2_cnv).sum() / len(tr2_cnv)))

                # 'synthetic' timeseries based on lsq solution, and calcualte the residauls
                syn  = dot(matinv, solution)
                res = tr2_cnv - syn

                # variance indicator based on stdev of residuals and tr2 timeseries.
                myvar = std(res) / std(tr2_cnv)

                # coherence of waveforms tr2 and synthetic
                coh = dot(tr2_cnv, syn) / sqrt(dot(tr2_cnv, tr2_cnv) * dot(syn, syn))

            else:

                amp_ratio = 0.0
                resid = array([])

                # log root mean square of the tr2 timseries amplitude
                lrms = 0.0

                # variance indicator based on stdev of residuals and tr2 timeseries.
                myvar = inf

                # coherence of waveforms tr2 and synthetic
                coh = 0.0


            # Add this segment's results to component results
            results.add_segment(APSurveySegmentResult(start_t, w, resid, amp_ratio,
                                                      lrms, myvar, coh,
                                                      coh >= self.coherence_cutoff))

            # start time for next segment
            start_t += self.segment_size_secs


    def _compare_verticals(self, tr1, tr2, tr1_resp, tr2_resp, results):
        """
        Compares the vertical component amplitudes of two traces from different sensors to determine the 
        relative sensitivities of hte two components. 
        The traces is are analyzed in segments of length specified in the config.
        The current tr2 response is removed and the response from tr1 is applied

        Args:
            tr1 (Trace): sensor 1 ('baseline') trace
            tr2 (Trace): sensor 2 trace
            tr1_resp (np.array): sensor 1 response to apply to tr2
            tr2_resp (np.array): tr2 response to deconvolve
            results (APSurveyComponentResult): COmponent results to assumulating individual segments results

        Returns: None

        """

        # find latest astarttime of the two traces.
        # The trace starttimes will be different when due to non-zero clock offset adjustment.
        start_t = max(tr1.stats.starttime, tr2.stats.starttime)
        tr1_sr = tr1.stats.sampling_rate
        tr2_sr = tr2.stats.sampling_rate
        tr1_time_delta = 1.0/tr1_sr
        tr2_time_delta = 1.0/tr2_sr

        # loop through traces segment by segment
        while start_t + self.segment_size_secs < tr1.stats.endtime:

            tr1_seg = tr1.slice(starttime=start_t,
                                        endtime=start_t + self.segment_size_secs - tr1_time_delta)
            tr2_seg = tr2.slice(starttime=start_t,
                                        endtime=start_t + self.segment_size_secs - tr2_time_delta)

            # get raw float data
            tr1_seg = tr1_seg.data.astype(float64)
            tr2_seg = tr2_seg.data.astype(float64)

            # need to demean and taper strm2, deconvolve it's resp, convolve strm1 resp
            tr1_seg = ss.detrend(tr1_seg, type='constant')
            tr2_seg = ss.detrend(tr2_seg, type='constant')

            # each side Creating tapers...
            fraction = (self.segment_size_trim/self.segment_size_secs)
            taper = tukey(len(tr2_seg), alpha=fraction * 2, sym=True)

            # remove tr2 response from timeseries
            tr2_fft = rfft(multiply(tr2_seg, taper))
            tr2_fft_dcnv = divide(tr2_fft[1:], tr2_resp[1:])
            tr2_fft_dcnv = insert(tr2_fft_dcnv, 0, [0.0])
            tr2_dcnv = irfft(tr2_fft_dcnv)

            # demean again before convolving strm2 resp
            tr2_dcnv = ss.detrend(tr2_dcnv, type='constant')
            tr2_fft = rfft(multiply(tr2_dcnv, taper))
            tr2_fft_cnv = multiply(tr2_fft, tr1_resp)
            tr2_cnv = irfft(tr2_fft_cnv)
            del tr2_fft, tr2_fft_dcnv, tr2_dcnv, tr2_fft_cnv

            # trim off taper segment_size_trim in secs
            trim_cnt = round(tr1_sr * self.segment_size_trim)
            tr1_seg = tr1_seg[trim_cnt:-trim_cnt]
            trim_cnt = round(tr2_sr * self.segment_size_trim)
            tr2_cnv = tr2_cnv[trim_cnt:-trim_cnt]

            # resample all 3 timeseries to analsysi sample rate
            tr1_seg = ss.resample(tr1_seg, round(len(tr1_seg) / (tr1_sr / self.analysis_sample_rate)))
            tr2_cnv = ss.resample(tr2_cnv, round(len(tr2_cnv) / (tr2_sr / self.analysis_sample_rate)))

            # demean and detrend
            tr1_seg = ss.detrend(tr1_seg, type='linear')
            tr2_cnv = ss.detrend(tr2_cnv, type='linear')

            # apply taper before filtering
            fraction = 1/12
            taper = tukey(len(tr1_seg), fraction * 2, sym=True)
            tr1_seg *= taper
            tr2_cnv *= taper

            # apply band pass filter with stops from config
            tr1_seg = osf.bandpass(tr1_seg, self.bp_start, self.bp_stop, self.analysis_sample_rate, zerophase=True)
            tr2_cnv = osf.bandpass(tr2_cnv, self.bp_start, self.bp_stop, self.analysis_sample_rate, zerophase=True)

            tr1_stdev = tr1_seg.std()
            tr2_stdev = tr2_cnv.std()

            if (round(tr1_stdev) == 0) or (round(tr2_stdev) == 0):
                # compute amplitude ratio
                amp_ratio = 0.0

                # log root mean square of the tr2 timseries amplitude
                lrms = 0.0

                # variance indicator based on stdev of residuals and tr2 timeseries
                myvar = inf

                # coherence of waveforms tr2 and synthetic
                coh = 0.0
            else:
                # compute amplitude ratio
                amp_ratio = tr2_stdev / tr1_stdev

                # log root mean square of the tr2 timseries amplitude
                lrms = log10(sqrt(multiply(tr2_cnv, tr2_cnv).sum() / len(tr2_cnv)))

                # 'synthetic' timeseries based on lsq solution, and calcualte the residauls
                syn  = tr1_seg * amp_ratio
                res = tr2_cnv - syn

                # variance indicator based on stdev of residuals and tr2 timeseries
                myvar = res.std() / tr2_cnv.std()

                # coherence of waveforms tr2 and synthetic
                coh = dot(tr2_cnv, syn) / sqrt(dot(tr2_cnv, tr2_cnv) * dot(syn, syn))

            # Add this segment's results to component results
            results.add_segment(APSurveySegmentResult(start_t, 0.0, 0.0, amp_ratio,
                                                      lrms, myvar,
                                                      coh, coh >= self.coherence_cutoff))

            # start time for next segment
            start_t += self.segment_size_secs

