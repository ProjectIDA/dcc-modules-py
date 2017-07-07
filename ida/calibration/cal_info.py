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
import glob
from os import getcwd, environ
from os.path import exists, abspath, join, dirname
from pathlib import Path
import datetime
import yaml
#import pprint

from fabulous.color import bold, blue, red, green
from pandas.core.frame import DataFrame

from ida.tui import select, SelectResult
from ida.instruments import CALTYPE_RBLF, CALTYPE_RBHF, SEISMOMETER_MODELS
from ida.calibration import nom_resp_for_model, cur_resp_for_model_station_comp, local_resp_files
from ida.signals.paz import PAZ
from ida.db.io import read
from ida.db.query import get_stages

VERSION = '0.1'


class CalInfo():

    STARTING_GROUP = 1
    LF_FITTING_GROUP = 2
    HF_FITTING_GROUP = 3

    CAL_TYPE_TITLE = {
        CALTYPE_RBLF: 'LOW Freq',
        CALTYPE_RBHF: 'HIGH Freq'
    }

    tui_indent = 4
    tui_indent_str = ' '*tui_indent

    def __init__(self, sta=None, loc=None,
                 config_file=None, comp=None, chancode=None,
                 cal_raw_dir=None, cur_paz_dir=None,
                 nom_paz_dir=None, db_dir=None):

        self.config_file = config_file
        self.cal_raw_dir = cal_raw_dir
        self.nom_paz_dir = nom_paz_dir
        self.cur_paz_dir = cur_paz_dir
        self.db_dir = db_dir

        self._info = {
            'sta': None,
            'loc': None,
            'sensor': None,
            'comp': None,
            'ctbto': None,
            'ctbto_calper': None,
            'chan': None,
            'opsr': None,
            'lfdatedir': None,
            'hfdatedir': None,
            'lfdatestr': None,
            'hfdatestr': None,
            'lffile': None,
            'hffile': None,
            'lfpath': None,
            'hfpath': None,
            'respfn': None,
            'fullpaz': None,
            'newpaz': None,
            'lfpert': None,
            'hfpert': None,
            # below 4 keys only used for CTBTO calculations
            'q330_nom_sens': None,
            'q330_calib_factor': None,
            'sensor_nom_sens': None,
            'sensor_calib_z': None,
            'sensor_calib_1': None,
            'sensor_calib_2': None,
        }

        self.omit_lf = False
        self.omit_hf = False

#        _, self.stages_df = read('datascope', self.db_dir, 'stage')
        self.stages_df = None
        self._lf_chn_stages = None
        self._hf_chn_stages = None

        self.ok = True

        if not config_file:
            self.mode = 'interactive'
            self.sta = sta
            self.loc = loc
            self.cal_raw_dir = cal_raw_dir
            self.nom_paz_dir = nom_paz_dir
            self.cur_paz_dir = cur_paz_dir
            self.db_dir = db_dir

            # if not using configuration file, must access datascopedb
            # from db_dir and then to read stages into pandas dataframe
            self.db_dir = db_dir
            if db_dir:
                if exists(abspath(db_dir)):
                    self.db_dir = db_dir
                    _, self.stages_df = read('datascope', self.db_dir, 'stage')
                else:
                    self.db_dir = None
                    raise ValueError(self.tui_indent_str +
                                     bold(red(
                                         'Error: Database directory [' + db_dir + '] does not exist.'
                                     )))
        else:
            self.mode = 'config-file'
            with open(self.config_file, 'rt') as cfl:
                config_txt = cfl.read()
            try:
                self._config = yaml.load(config_txt)
                self._config_dir = Path(self.config_file).parent
            except:
                self._config = {}
                print(red('Error parsing YAML config file: ' + self.config_file))
                self.ok = False
            else:
                # check for any fatal error conditions in successfully parsed config file
                self.ok = self.config_file_isvalid()
                self.sta = self._config['station'].lower()
                self.loc = self._config['loc']
                self.run_comp = comp.lower()
                self.run_chan = chancode
                self.nom_paz_dir = self._config['nom_paz_dir']
                self.cur_paz_dir = self._config['cur_paz_dir']

    def config_file_isvalid(self):

        if 'station' not in self._config:
            print(red('Error in configuration file: "station" key is missing.'))
            valid = False

        if 'loc' not in self._config:
            print(red('Error in configuration file: "loc" key is missing.'))
            valid = False

        valid = True
        if 'process_lf' not in self._config:
            print(red('Error in configuration file: "process_lf" key missing.'))
            valid = False
        else:
            process_lf = self._config['process_lf']
            if process_lf:
                if 'lf_ms_file' not in self._config:
                    print(red('Error in configuration file: "lf_ms_file" key is missing.'))
                if 'lf_qcal_file' not in self._config:
                    print(red('Error in configuration file: "lf_qcal_file" key is missing.'))
                if 'lf_cal_date' not in self._config:
                    print(red('Error in configuration file: "lf_cal_date" key is missing.'))
        if 'process_hf' not in self._config:
            print(red('Error in configuration file: "process_hf" key missing.'))
            valid = False
        else:
            process_hf = self._config['process_hf']
            if process_hf:
                if 'hf_ms_file' not in self._config:
                    print(red('Error in configuration file: "hf_ms_file" key is missing.'))
                if 'hf_qcal_file' not in self._config:
                    print(red('Error in configuration file: "hf_qcal_file" key is missing.'))
                if 'hf_cal_date' not in self._config:
                    print(red('Error in configuration file: "hf_cal_date" key is missing.'))
            else:
                if 'process_lf' in self._config and not process_lf:
                    print(red('Error in configuration file: neither LF or HF procesing is enabled.'))
                    valid = False

        if 'sensor_code' not in self._config:
            print(red('Error in configuration file: "sensor_code" key is missing.'))
            valid = False

        if 'sampling_rate' not in self._config:
            print(red('Error in configuration file: "sampling_rate" key is missing.'))
            valid = False

        if 'ctbto_analysis' not in self._config:
            print(red('Error in configuration file: "ctbto_analysis" key is missing.'))
            valid = False
        elif self._config['ctbto_analysis']:

#           if 'q330_sensor_port' not in self._config:
#               print(red('Error in configuration file: "q330_sensor_port" key is missing.'))
#               valid = False

            if 'sensor_nom_sens' not in self._config:
                print(red('Error in configuration file: "sensor_nom_sens" key is missing.'))
                valid = False

            if 'q330_nom_sens' not in self._config:
                print(red('Error in configuration file: "q330_nom_sens" key is missing.'))
                valid = False

            if 'q330_calib_factor' not in self._config:
                print(red('Error in configuration file: "q330_calib_factor" key is missing.'))
                valid = False

            if 'sensor_calib_factor_z' not in self._config:
                print(red('Error in configuration file: "sensor_calib_factor_z" key is missing.'))
                valid = False

            if 'sensor_calib_factor_1' not in self._config:
                print(red('Error in configuration file: "sensor_calib_factor_1" key is missing.'))
                valid = False

            if 'sensor_calib_factor_2' not in self._config:
                print(red('Error in configuration file: "sensor_calib_factor_2" key is missing.'))
                valid = False

            if 'ctbto_calper' not in self._config:
                print(red('Error in configuration file: "ctbto_calper" key is missing.'))
                valid = False

        if 'nom_paz_dir' not in self._config:
            print(red('Error in configuration file: "nom_paz_dir" key is missing.'))
            valid = False
        else:
            if not exists(self._config['nom_paz_dir']):
                print(red('Error in configuration: nom_paz_dir does not exist.' + \
                          self._config['nom_paz_dir']))
                valid = False

        if 'cur_paz_dir' not in self._config:
            print(red('Error in configuration file: "cur_paz_dir" key is missing.'))
            valid = False
        else:
            if not exists(self._config['cur_paz_dir']):
                print(red('Error in configuration: cur_paz_dir does not exist: ' + \
                          self._config['cur_paz_dir']))
                valid = False

        return valid

    def reset(self, key_list):
        for key in key_list:
            self._info[key] = None

    def reset_all(self):
        for key in self._info:
            self._info[key] = None

    def reset_all_except(self, key_list):
        for key in self._info:
            if key not in key_list:
                self._info[key] = None

    def check_type(self, val, valtype, valname):
        if not isinstance(val, valtype):
            raise TypeError('Invalid type for ' + valname + ':' + val + \
                            'Should be ' + str(valtype))

    def sensor_gnom(self):
        if self.mode == 'interactive':
            return self.chn_stages.iloc[0].gnom
        elif self.mode == 'config-file':
            return self._info['sensor_nom_sens']

    def sensor_gcalib(self):
        if self.mode == 'interactive':
            return self.chn_stages.iloc[0].gcalib
        elif self.mode == 'config-file':
            return self._info['sensor_calib_'+ self.comp.lower()]

    def q330_gnom(self):
        if self.mode == 'interactive':
            return self.chn_stages.iloc[2].gnom
        elif self.mode == 'config-file':
            return self._info['q330_nom_sens']

    def q330_gcalib(self):
        if self.mode == 'interactive':
            return self.chn_stages.iloc[2].gcalib
        elif self.mode == 'config-file':
            return self._info['q330_calib_factor']

    @property
    def sta(self):
        return self._info['sta'] if self._info['sta'] else ''

    @sta.setter
    def sta(self, value):
        if value:
            self.check_type(value, str, 'Station code')
            self._info['sta'] = value
        else:
            self._info['sta'] = None

        self.reset_all_except(['sta'])

    @property
    def loc(self):
        return self._info['loc'] if self._info['loc'] else ''

    @loc.setter
    def loc(self, value):
        if value:
            self.check_type(value, str, 'Location code')
            self._info['loc'] = value
        else:
            self._info['loc'] = None

        self.reset_all_except(['sta', 'loc'])

    @property
    def sensor(self):
        return self._info['sensor'] if self._info['sensor'] else ''

    @sensor.setter
    def sensor(self, value):

        if value:
            self.check_type(value, str, 'sensor')
            if self.mode == 'config-file':
                self._info['sensor'] = value
            else:
                self._info['sensor'] = value if exists(join(self.cal_raw_dir,
                                                            self.sta,
                                                            self.loc,
                                                            value)) else None
        else:
            self._info['sensor'] = None

        self.reset_all_except(['sta', 'loc', 'sensor', 'comp',' chan'])

    @property
    def comp(self):
        return self._info['comp'] if self._info['comp'] else ''

    @comp.setter
    def comp(self, value):

        if value:
            self.check_type(value, str, 'component')
            value = value.upper()
            if value in ['Z', '1', '2']:
                self._info['comp']= value
            else:
                raise ValueError('Invalid value for Component: ' + value + ". Should be one of ['Z', '1', '2']")
        else:
            self._info['comp'] = None

    @property
    def ctbto(self):
        return self._info['ctbto'] if self._info['ctbto'] else ''

    @ctbto.setter
    def ctbto(self, value):

        if value:
            self.check_type(value, str, 'sensor ctbto Y/N flag')
            value = value.upper()
            if value in ['Y', 'N']:
                self._info['ctbto']= value
            else:
                raise ValueError('Invalid value for CTBTO Y/N flag: ' + value + ". Should be one of ['Y', 'N']")
        else:
            self._info['ctbto'] = None

    @property
    def ctbto_calper(self):
        return self._info['ctbto_calper'] if self._info['ctbto_calper'] else ''

    @ctbto_calper.setter
    def ctbto_calper(self, value):

        if value:
            if isinstance(value, (float, int)) and value > 0.0:
                self._info['ctbto_calper']= float(value)
            else:
                raise ValueError('Invalid value for CTBTO CALPER: ' + value + ". Must be a numeric value > 0")
        else:
            self._info['ctbto_calper'] = None

    def is_ctbto(self):
        return self.ctbto == 'Y'

    @property
    def chan(self):
        return self._info['chan'] if self._info['chan'] else ''

    @chan.setter
    def chan(self, value):

        if value and (len(value) == 3):
            self.check_type(value, str, 'channel code')
            value = value.upper()
            self._info['chan'] = value
        else:
            self._info['chan'] = None

    @property
    def opsr(self):
        return self._info['opsr'] if self._info['opsr'] else ''

    @opsr.setter
    def opsr(self, value):

        if value:
            self.check_type(value, (int), 'Operating sample rate must be an integer')
            self._info['opsr'] = value
        else:
            self._info['opsr'] = None

    @property
    def lfdatestr(self):
        return self.lfdatedir[:10]

    @property
    def lfdatedir(self):
        return self._info['lfdatedir'] if self._info['lfdatedir'] else ''

    @lfdatedir.setter
    def lfdatedir(self, value):

        if self.sta and self.loc and self.sensor:
            if value:
                self.check_type(value, str, 'low freq cal date')
                if self.mode == 'config-file':
                    # really only the ISO date from config file, not a dir segment
                    self._info['lfdatedir'] = value
                else:
                    rbpath = join(self.cal_raw_dir,
                                  self.sta,
                                  self.loc,
                                  self.sensor,
                                  CALTYPE_RBLF,
                                  value)
                    if exists(rbpath):
                        if self.lfdatedir != value:
                            self.reset(['lffile'])

                        self._info['lfdatedir'] = value
                        self.reset(['lffile', 'respfn', 'fullpaz', 'lfpert', 'hfpert'])

                        ms_files = glob.glob(join(rbpath, '*.ms'))
                        log_files = glob.glob(join(rbpath, '*.log'))
                        ms_stems = [Path(msfn).stem for msfn in ms_files]
                        log_stems = [Path(logfn).stem for logfn in log_files]
                        paired_files = [stem for stem in ms_stems if stem in log_stems]
                        if len(paired_files) == 1:
                            self.lffile = paired_files[0]

            else:
                self._info['lfdatedir'] = None
                self.reset(['lffile', 'respfn', 'fullpaz', 'lfpert', 'hfpert'])


    @property
    def hfdatestr(self):
        return self.hfdatedir[:10]

    @property
    def hfdatedir(self):
        return self._info['hfdatedir'] if self._info['hfdatedir'] else ''

    @hfdatedir.setter
    def hfdatedir(self, value):


        if self.sta and self.loc and self.sensor:
            if value:
                self.check_type(value, str, 'high freq cal date')
                if self.mode == 'config-file':
                    # really only the ISO date from config file, not a dir segment
                    self._info['hfdatedir'] = value
                else:
                    rbpath = join(self.cal_raw_dir,
                                  self.sta,
                                  self.loc,
                                  self.sensor,
                                  CALTYPE_RBHF,
                                  value)
                    if exists(rbpath):
                        if self.hfdatedir != value:
                            self.reset(['hffile'])

                        self._info['hfdatedir'] = value
                        self.reset(['hffile'])

                        ms_files = glob.glob(join(rbpath, '*.ms'))
                        log_files = glob.glob(join(rbpath, '*.log'))
                        ms_stems = [Path(msfn).stem for msfn in ms_files]
                        log_stems = [Path(logfn).stem for logfn in log_files]
                        paired_files = [stem for stem in ms_stems if stem in log_stems]
                        if len(paired_files) == 1:
                            self.hffile = paired_files[0]

            else:
                self._info['hfdatedir'] = None
                self.reset(['hffile'])

    @property
    def lffile(self):
        return self._info['lffile'] if self._info['lffile'] else ''

    @lffile.setter
    def lffile(self, value):

        if self.sta and self.loc and self.sensor and self.lfdatedir:
            if value:
                self.check_type(value, str, 'low freq cal data filename')
                if self.mode == 'config-file':
                    # value is full path if coming from config file
                    value = abspath(value)
                    if exists(value):
                        self._info['lffile'] = Path(value).stem
                        self._info['lfpath'] = dirname(value)
                    else:
                        print(red('Error: LF file does not exist: ' + value))
                        print(red('Disabling LF analysis.'))
                        self.omit_lf = True
                        self._info['lffile'] = None
                        self._info['lfpath'] = None

                else:
                    fpath = join(self.cal_raw_dir, self.sta, self.loc, self.sensor,
                                 CALTYPE_RBLF, self.lfdatedir)
                    if exists(join(fpath, value + '.ms')):
                        self._info['lffile'] = value
                        self._info['lfpath'] = fpath
                    else:
                        self._info['lffile'] = None
                        self._info['lfpath'] = None
            else:
                self._info['lffile'] = None
                self._info['lfpath'] = None

    @property
    def lfpath(self):
        return self._info['lfpath'] if self._info['lfpath'] else ''

    @property
    def hffile(self):
        return self._info['hffile'] if self._info['hffile'] else ''

    @hffile.setter
    def hffile(self, value):

        if self.sta and self.loc and self.sensor and self.hfdatedir:
            if value:
                self.check_type(value, str, 'high freq cal data filename')
                if self.mode == 'config-file':
                    # value is full path if coming from config file
                    value = abspath(value)
                    if exists(value):
                        self._info['hffile'] = Path(value).stem
                        self._info['hfpath'] = dirname(value)
                    else:
                        print(red('Error: HF file does not exist: ' + value))
                        print(red('Disabling HF analysis.'))
                        self.omit_hf = True
                        self._info['hffile'] = None
                        self._info['hfpath'] = None
                else:
                    fpath = join(self.cal_raw_dir, self.sta, self.loc, self.sensor,
                                 CALTYPE_RBHF, self.hfdatedir)
                    if exists(join(fpath, value + '.ms')):
                        self._info['hffile'] = value
                        self._info['hfpath'] = fpath
                    else:
                        self._info['hffile'] = None
                        self._info['hfpath'] = None
            else:
                self._info['hffile'] = None
                self._info['hfpath'] = None

    @property
    def hfpath(self):
        return self._info['hfpath'] if self._info['hfpath'] else ''

    @property
    def respfn(self):
        return self._info['respfn'] if self._info['respfn'] else ''

    @respfn.setter
    def respfn(self, value):

        if value:
            self.check_type(value, str, 'Response filename')
            self._info['respfn'] = value if exists(value) else None
        else:
            self._info['respfn'] = None

    @property
    def fullpaz(self):
        return self._info['fullpaz']

    @fullpaz.setter
    def fullpaz(self, value):

        if value:
            self.check_type(value, PAZ, 'Poles and Zeros response object')
            if self.respfn:
                self._info['fullpaz'] = value
            else:
                self._info['fullpaz'] = None
        else:
            self._info['fullpaz'] = None

    @property
    def newpaz(self):
        return self._info['newpaz']

    @newpaz.setter
    def newpaz(self, value):

        if value:
            self.check_type(value, PAZ, 'Poles and Zeros response object')
            self._info['newpaz'] = value
        else:
            self._info['newpaz'] = None

    @property
    def lfpert(self):
        return self._info['lfpert']

    @lfpert.setter
    def lfpert(self, value):

        if value:
            self.check_type(value, tuple, 'Response fitting perturb LF p/z')
            if self.fullpaz:
                self._info['lfpert'] = value
            else:
                self._info['lfpert'] = None
        else:
            self._info['lfpert'] = None

    @property
    def hfpert(self):
        return self._info['hfpert']

    @hfpert.setter
    def hfpert(self, value):

        if value:
            self.check_type(value, tuple, 'Response fitting perturb HF p/z')
            if self.fullpaz:
                self._info['hfpert'] = value
            else:
                self._info['hfpert'] = None
        else:
            self._info['hfpert'] = None

    def __repr__(self):
        return str(self._info)

    def select_raw_cal_sensordir(self):

        staloc_path = join(self.cal_raw_dir, self.sta, self.loc)
        sensordirlist = glob.glob(join(staloc_path, '*'))
        sensordirlist = sorted([Path(item).stem for item in sensordirlist if Path(item).stem.upper() in SEISMOMETER_MODELS])
        list_len = len(sensordirlist)

        errmsg = ''
        if list_len > 0:
            result, _, user_choice_groups, _ = select([sensordirlist],
                                                      title='Select Sensor Directory',
                                                      prompt='Select # of desired Sensor directory ("q" => quit): ',
                                                      implicit_quit_q=True,
                                                      menu_on_error=True, err_message='Invalid selection. Please try again.',
                                                      indent_width=self.tui_indent)
            if result == SelectResult.ok:
                ndx = user_choice_groups[0][0]
                self.sensor = sensordirlist[ndx]

        else:
            result = SelectResult.error
            errmsg = 'No directories for supported sensor models found in: ' + abspath(staloc_path)
            errmsg = errmsg + '\n' + self.tui_indent_str + 'Supported sensor models: ' + str(SEISMOMETER_MODELS)

        return result, errmsg

    def select_component(self):

        if self.mode == 'interactive':
            prompt = 'Enter selection ("q" => quit, "b" => back): '
            allow_back = True
        elif self.mode == 'config-file':
            prompt = 'Enter selection ("q" => quit): '
            allow_back = False

        complist = [('Z', 'Vertical'),
                    ('1', '1'),
                    ('2', '2')]
        result, _, user_choice_groups, _ = select([complist], title='Select Component',
                                                  prompt=prompt, implicit_quit_q=True,
                                                  implicit_back_b=allow_back,
                                                  menu_on_error=True,
                                                  err_message='Invalid choice. Please try again',
                                                  indent_width=self.tui_indent)

        if result == SelectResult.ok:
            ndx = user_choice_groups[0][0]
            self.comp = complist[ndx][0]

        return result

    def select_ctbto_flag(self):

        flaglist = [('Y', 'Yes: Sensor requires a CTBTO Calibration Report'),
                    ('N', 'No: Not a CTBTO sensor')]

        # flaglist = [item for item in ctbto_dict.items()]
        result, _, user_choice_groups, _ = select([flaglist], 'CTBTO Status',
                                                  prompt='Is this a CTBTO sensor ("q" => quit, "b" => back)? ',
                                                  implicit_quit_q=True, implicit_back_b=True,
                                                  menu_on_error=True, err_message='Invalid choice. Please try again',
                                                  indent_width=self.tui_indent)

        if result == SelectResult.ok:
            ndx = user_choice_groups[0][0]
            self.ctbto = flaglist[ndx][0]

        return result

    def enter_ctbto_calper(self):

        result = SelectResult.noop
        while result not in [SelectResult.quit,
                             SelectResult.goback,
                             SelectResult.ok]:
            print()
            calper_str = input(bold(blue(self.tui_indent_str + ' Enter CTBTO Calibration Period (CALPER) value [1.0]: ')))

            try:
                calper = float(calper_str)
                if calper <= 0.0:
                    print()
                    print(self.tui_indent_str + bold(red(
                        'ERROR: CALPER must be > 0.0.'
                    )))
                    print()
                else:
                    self.ctbto_calper = calper
                    result = SelectResult.ok
            except:
                if calper_str.upper() == 'Q':
                    result = SelectResult.quit
                elif calper_str.upper() == 'B':
                    result = SelectResult.goback
                else:
                    print()
                    print(self.tui_indent_str +
                          bold(red('ERROR: CTBTO Calibration Period (CALPER) must be a positive numeric value.')))
                    print()

        return result

    def enter_chan(self):

        result = SelectResult.noop
        while result not in [SelectResult.quit,
                             SelectResult.goback,
                             SelectResult.ok]:
            print()
            chan = input(bold(blue(self.tui_indent_str + ' Enter 3 character channel code: ')))
            chan = chan.upper()
            if (len(chan) == 3):
                if chan[2] in ['Z', '1', '2']:
                    self.chan = chan
                    result = SelectResult.ok
                else:
                    result = SelectResult.error

            elif len(chan) == 1:
                if chan == 'Q':
                    result = SelectResult.quit
                elif chan == 'B':
                    result = SelectResult.goback
                else:
                    result = SelectResult.error

            if result == SelectResult.error:
                print()
                print(self.tui_indent_str + bold(red('Invalid channel code.')))
                print(self.tui_indent_str + bold(red('Make sure code is 3 characters and ends in Z, 1 or 2.')))
                print()

        return result

    def enter_opsr(self):

        result = SelectResult.noop
        while result not in [SelectResult.quit,
                             SelectResult.goback,
                             SelectResult.ok]:
            print()
            opsr_str = input(bold(blue(self.tui_indent_str + ' Enter sensor operating frequency [20, 40]: ')))

            try:
                opsr = round(float(opsr_str))
                if opsr not in [20, 40]:
                    print()
                    print(self.tui_indent_str + bold(red(
                        'ERROR: Unsupported frequency. 20 and 40 are the supported frrequencies.'
                    )))
                    print()
                else:
                    self.opsr = opsr
                    result = SelectResult.ok
            except:
                if opsr_str.upper() == 'Q':
                    result = SelectResult.quit
                elif opsr_str.upper() == 'B':
                    result = SelectResult.goback
                else:
                    print()
                    print(self.tui_indent_str + bold(red('ERROR: Invalid frequency. Please enter frequency in Hz.')))
                    print()

        return result

    def select_raw_cal_date(self, cal_type):

        if cal_type not in [CALTYPE_RBLF, CALTYPE_RBHF]:
            raise ValueError('Invalid CAL TYPE supplied: ' + cal_type)

        datedirlist = glob.glob(join(self.cal_raw_dir, self.sta, self.loc, self.sensor, cal_type, '*'))
        datedirlist = sorted([Path(item).stem for item in datedirlist])

        omit_option = [('S', 'Skip {} processing'.format(CalInfo.CAL_TYPE_TITLE[cal_type]))]

        errmsg = ''
        result, choices, user_choice_groups, _ = select([datedirlist, omit_option],
                                                        title=CalInfo.CAL_TYPE_TITLE[cal_type] + ' dates for {} {}: '.format(self.sta,
                                                                                                          self.loc,
                                                                                                          self.sensor),
                                                        prompt='Select # of desired DATE ("q" => quit, "b" => back): ',
                                                        implicit_quit_q=True, implicit_back_b=True,
                                                        menu_on_error=True, err_message='Invalid selection. Please try again.',
                                                        indent_width=self.tui_indent)
        if result == SelectResult.ok:
            choice = choices[0].upper()
            if cal_type == CALTYPE_RBLF:
                if choice == 'S':
                    self.omit_lf = True
                else:
                    self.lfdatedir = datedirlist[user_choice_groups[0][0]]
            else:
                if choice == 'S':
                    self.omit_hf = True
                else:
                    self.hfdatedir = datedirlist[user_choice_groups[0][0]]

        return result, errmsg

    def select_starting_response_file(self):

        # run for IDA deployment ONLY

        pick_groups = []
        group_titles = []

        sensor_fn_list = []
        if self.mode == 'interactive':
            if isinstance(self.chn_stages, DataFrame):
                resp_file = self.chn_stages.iloc[0].dfile
                if resp_file:
                    sensor_fn_list = [resp_file]
            pick_groups.append(sensor_fn_list)
            group_titles.append('Response on Date of Calibration')
        elif self.mode == 'config-file':
            sensor_fn_list = cur_resp_for_model_station_comp(self.cur_paz_dir,
                                                             self.sensor.lower(),
                                                            self.sta,
                                                            self.comp)
            pick_groups.append(sensor_fn_list)
            group_titles.append('Existing Responses for Station, Sensor, Component')

        nom_resps = nom_resp_for_model(self.nom_paz_dir, self.sensor.lower())
        nom_resps_names = sorted([Path(respfn).name for respfn in nom_resps])
        pick_groups.append(nom_resps_names)
        group_titles.append('Nominal Responses for sensor: ' + self.sensor.upper())

        local_resp = local_resp_files()
        local_resp_names = sorted([Path(respfn).name for respfn in local_resp])
        pick_groups.append(local_resp_names)
        group_titles.append('Responses Found in CWD')

        errmsg = ''
        if len(pick_groups) > 0:
            result, user_choices, _, choice_tpls = select(pick_groups, 'Select Response file to use as starting model',
                                                          prompt='Select # of desired response file ("q" => quit, "b" => back): ',
                                                          group_titles=group_titles,
                                                          implicit_quit_q=True, implicit_back_b=True,
                                                          menu_on_error=True, err_message='Invalid choice. Please try again',
                                                          indent_width=self.tui_indent)
            if result == SelectResult.ok:
                fn = pick_groups[choice_tpls[0][0]][choice_tpls[0][1]]
                if choice_tpls[0][0] == 0:     #  first group is single current resp file as in DB
                    self.respfn = join(self.cur_paz_dir, fn)
                elif choice_tpls[0][0] == 1:   #  second group are files from nominal response dir
                    self.respfn = join(self.nom_paz_dir, fn)
                elif choice_tpls[0][0] == 2:   #  third group are respo files found in cwd
                    self.respfn = join(getcwd(), fn)

                # read in resp data and make sure file is valid (will raise an eception if not)
                self.fullpaz = PAZ(pzfilename=self.respfn, fileformat='ida', mode='vel', units='hz')

        else:
            result = SelectResult.error
            errmsg = 'No response files found in DB, Resp Dir or cwd.'

        return result, errmsg

    def select_perturb_map(self, cal_type, paz=None):

        if paz:
            starting_paz = paz
        else:
            starting_paz = self.fullpaz

        prompt = 'Enter comma separated list ("q" => quit): '
        grp_titles = []
        pert_choices = []

        if not isinstance(starting_paz, PAZ):
            raise TypeError('PAZ must be a populated PAZ object')

        poles_pert_def, zeros_pert_def = self.fullpaz.perturb_defaults()
        if cal_type == CALTYPE_RBLF:
            paz_fit = starting_paz.make_partial2(norm_freq=1.0, partial_mode=self.fullpaz.PARTIAL_FITTING_LF)
            ppert_def = poles_pert_def[0]
            zpert_def = zeros_pert_def[0]
        else:
            paz_fit = starting_paz.make_partial2(norm_freq=1.0, partial_mode=self.fullpaz.PARTIAL_FITTING_HF)
            ppert_def = poles_pert_def[1]
            zpert_def = zeros_pert_def[1]

        if ppert_def or zpert_def:
            pert_choices.append([('D', 'Use Defaults (indicated by "<==")')])
            grp_titles.append('')

        # make list of LF p/z to user perturbing choices
        # Do LOW FREQ FIRST
        zero_pert_choices = []
        pole_pert_choices = []
        for ndx, val in enumerate(paz_fit.zeros()):
            if ndx in zpert_def:
                zero_pert_choices.append(str(val) + ' <==')
            else:
                zero_pert_choices.append(str(val))
        for ndx, val in enumerate(paz_fit.poles()):
            if ndx in ppert_def:
                pole_pert_choices.append(str(val) + ' <==')
            else:
                pole_pert_choices.append(str(val))

        errmsg = ''
        pert_choices.extend([zero_pert_choices, pole_pert_choices])
        grp_titles.extend(['Zeros', 'Poles'])
        result, choices, pert_choice_groups, _ = select(pert_choices,
                                                        'Select {} zeros & poles to perturb'.format(CalInfo.CAL_TYPE_TITLE[cal_type]),
                                                        prompt=prompt, group_titles=grp_titles,
                                                        multiple_choice=True,
                                                        implicit_quit_q=True,
                                                        menu_on_error=True, err_message='Invalid entry. Please try again',
                                                        indent_width=self.tui_indent)

        if result == SelectResult.ok:
            if choices[0].upper() == 'D':  # using defaults
                pert = (ppert_def, zpert_def)  # beware, put poles then zeros in this map tuple
            else:
                pert = (pert_choice_groups[-1], pert_choice_groups[-2])

            if cal_type == CALTYPE_RBLF:
                self.lfpert = pert
            else:
                self.hfpert = pert

        return result, errmsg

    @property
    def chn_stages(self):

        if isinstance(self._lf_chn_stages, DataFrame) and not self._lf_chn_stages.empty:
            stages_df = self._lf_chn_stages
        elif isinstance(self._hf_chn_stages, DataFrame) and not self._hf_chn_stages.empty:
            stages_df = self._hf_chn_stages
        else:
            stages_df = None

        return stages_df

    def sensor_cnt(self):

        if self.sta and self.loc:
            return len(glob.glob(join(self.cal_raw_dir, self.sta, self.loc, '*')))
        else:
            return 0

    def date_cnt(self, cal_type):
        if cal_type not in [CALTYPE_RBLF, CALTYPE_RBHF]:
            raise ValueError('Invalid CAL TYPE supplied: ' + cal_type)

        if self.sta and self.loc and self.sensor:
            return len(glob.glob(join(self.cal_raw_dir, self.sta, self.loc, self.sensor, cal_type, '*')))
        else:
            return 0

    def find_qcal_files(self, cal_type):

        if cal_type not in [CALTYPE_RBLF, CALTYPE_RBHF]:
            raise ValueError('Invalid CAL TYPE supplied: ' + cal_type)

        errmsg = ''
        if cal_type == CALTYPE_RBLF:
            adate = self.lfdatedir
        else:
            adate = self.hfdatedir

        result = SelectResult.noop
        if self.sta and self.loc and self.sensor and adate:
            dpath = abspath(join(self.cal_raw_dir, self.sta, self.loc, self.sensor, cal_type, adate))

            if exists(dpath):
                msfiles = glob.glob(join(dpath, '*.ms'))
                logfiles = glob.glob(join(dpath, '*.log'))

                if (len(msfiles) == 1) and (len(logfiles) == 1):
                    if cal_type == CALTYPE_RBLF:
                        self.lffile = Path(msfiles[0]).stem
                        self._info['lfpath'] = dpath
                    else:
                        self.hffile = Path(msfiles[0]).stem
                        self._info['hfpath'] = dpath
                    result = SelectResult.ok
                else:
                    result = SelectResult.error
                    errmsg = 'Error: Missing or multiple MS and LOG file(s) in directory: ' + dpath
            else:
                result = SelectResult.error
                errmsg = 'Directory does not exist: ' + dpath

        return result, errmsg

    def new_filename_stem(self):
        if self.sta and self.loc and self.sensor and self.comp and (self.lfdatedir or self.hfdatedir):
            adate = self.lfdatedir if self.lfdatedir else self.hfdatedir
            stem = self.sensor.lower() + '_' + self.sta.lower() + '_' + self.loc + '_' + \
                   adate.replace('-', '') + "_" + self.comp.upper()
        else:
            stem = 'INCOMPLETE_CAL_INFO'

        return stem

    def is_complete(self):  # , group):

        done = self.sta and self.loc and self.sensor and \
               self.opsr and self.comp and self.chan and \
               (self.lfdatedir or self.omit_lf) and (self.hfdatedir or self.omit_hf) and \
               self.respfn and self.fullpaz and self.ctbto

        return done

    def collect_backup(self):

        # in reverse order of querying user
        done = False

        if not done and self.ctbto:
            self.ctbto = None
            self._info['q330_nom_sens'] = None
            self._info['q330_calib_factor'] = None
            self._info['sensor_nom_sens'] = None
            self._info['sensor_calib_z'] = None
            self._info['sensor_calib_1'] = None
            self._info['sensor_calib_2'] = None
            done = True

        # if not done and self.hfpert:
        #     self.hfpert = None
        #     done = True
        #
        # if not done and self.lfpert:
        #     self.lfpert = None
        #     done = True
        #
        if not done and self.respfn:
            self.respfn = None
            self.fullpaz = None
            done = True

        if not done and self.opsr:
            self.opsr = None
            done = True
            if isinstance(self.chn_stages, DataFrame): # assume opsr pulled from DB, so back to chan
                self.chan = None

        if not done and self.chan:
            self.chan = None
            self._lf_chn_stages = None  # clear stages since channel changing
            self._hf_chn_stages = None  # clear stages since channel changing
            done = True

        if not done and self.comp:
            self.comp = None
            done = True

        if not done and (self.hfdatedir or self.omit_hf):
            self.hfdatedir = None
            self.hffile = None
            self.omit_hf = False
            self._hf_chn_stages = None  # clear stages since lfdate (epoch) changing
            self._info['hfpath'] = None
            done = True

        if not done and (self.lfdatedir or self.omit_lf):
            self.lfdatedir = None
            self.lffile = None
            self.omit_lf = False
            self._lf_chn_stages = None  # clear stages since lfdate (epoch) changing
            self._info['lfpath'] = None
            done = True

        if not done and self.sensor:
            self.sensor = None
            done = self.sensor_cnt() > 1

    def collect_next(self):

        col_res = SelectResult.noop
        col_msg = ''

        if not self.sensor:
            if self.mode == 'config-file':
                self.sensor = self._config['sensor_code'].lower()
                col_res = SelectResult.ok
            else:   # 'interactive':
                col_res, col_msg = self.select_raw_cal_sensordir()

        elif not self.lfdatedir and not self.omit_lf:
            if self.mode == 'config-file':
                self.omit_lf = not self._config['process_lf']
                if not self.omit_lf:
                    self.lfdatedir = self._config['lf_cal_date']
                    self.lffile = abspath(self._config['lf_ms_file'])
                col_res = SelectResult.ok
            else:
                col_res, col_msg = self.select_raw_cal_date(CALTYPE_RBLF)
                if col_res == SelectResult.goback:
                    self.collect_backup()
                    col_res = SelectResult.ok
                elif col_res == SelectResult.ok:
                    if not self.omit_lf:
                        col_res, col_msg = self.find_qcal_files(CALTYPE_RBLF)

        elif not self.hfdatedir and not self.omit_hf:
            if self.mode == 'config-file':
                self.omit_hf = not self._config['process_hf']
                if not self.omit_hf:
                    self.hfdatedir = self._config['hf_cal_date']
                    self.hffile = abspath(self._config['hf_ms_file'])
                col_res = SelectResult.ok
            else:
                col_res, col_msg = self.select_raw_cal_date(CALTYPE_RBHF)
                if col_res == SelectResult.goback:
                    self.collect_backup()
                    col_res = SelectResult.ok
                elif col_res == SelectResult.ok:
                    if not self.omit_hf:
                        col_res, col_msg = self.find_qcal_files(CALTYPE_RBHF)


        elif self.omit_lf and self.omit_hf:
            print(bold(red(CalInfo.tui_indent_str + ' ERROR: You can not skip both LF and HF processing.\n')))
            self.omit_lf = False # reset so both dates prompted for again
            self.omit_hf = False
            col_res = SelectResult.ok

        elif not self.comp:
            if self.mode == 'config-file':
                self.comp = self.run_comp
                col_res = SelectResult.ok
            else:
                col_res = self.select_component()
                if col_res == SelectResult.goback:
                    self.collect_backup()
                    col_res = SelectResult.ok

        elif not self.chan:
            if self.mode == 'config-file':
                self.chan = self.run_chan
                col_res = SelectResult.ok
            else:
                col_res = self.enter_chan()
                if col_res == SelectResult.goback:
                    self.collect_backup()
                    col_res = SelectResult.ok
                elif col_res == SelectResult.ok:
                    # good channel, see if we can pull sample rate out of DB, if available
                    if isinstance(self.stages_df, DataFrame):
                        if not self.omit_lf:
                            self._lf_chn_stages = get_stages(self.stages_df, self.sta, self.loc, self.chan, self.lfdatedir)
                        elif not self.omit_hf:
                            self._hf_chn_stages = get_stages(self.stages_df, self.sta, self.loc, self.chan, self.hfdatedir)

                        if isinstance(self.chn_stages, DataFrame) and not (len(self.chn_stages) < 3):
                            # get channel freq from stage 3 srate column
                            # ToDo: bury the srate ref in public acces method to make db independent
                            try:
                                self.opsr = round(float(self.chn_stages.iloc[2].srate))
                            except:
                                pass

        elif not self.opsr:
            if self.mode == 'config-file':
                self.opsr = self._config['sampling_rate']
                col_res = SelectResult.ok
            else:
                col_res = self.enter_opsr()
                if col_res == SelectResult.goback:
                    self.collect_backup()
                    col_res = SelectResult.ok

        elif not self.respfn:
            if self.mode == 'config-file':
                # for 'config-file' mode,print info after incorporating config file
                self.print_info()

            col_res, col_msg = self.select_starting_response_file()
            if col_res == SelectResult.goback:
                self.collect_backup()
                col_res = SelectResult.ok

        elif not self.ctbto:
            if self.mode == 'config-file':
                self.ctbto = 'Y' if self._config['ctbto_analysis'] else 'N'
                col_res = SelectResult.ok
                if self.is_ctbto():
#                    self._info['q330_sensor_port'] = self._config['q330_sensor_port']
# TODO: add KEY checks here
                    self._info['q330_nom_sens'] = self._config['q330_nom_sens']
                    self._info['q330_calib_factor'] = self._config['q330_calib_factor']
                    self._info['sensor_nom_sens'] = self._config['sensor_nom_sens']
                    self._info['sensor_calib_z'] = self._config['sensor_calib_factor_z']
                    self._info['sensor_calib_1'] = self._config['sensor_calib_factor_1']
                    self._info['sensor_calib_2'] = self._config['sensor_calib_factor_2']
            else:
                col_res= self.select_ctbto_flag()
                if col_res == SelectResult.goback:
                    self.collect_backup()
                    col_res = SelectResult.ok

        elif self.ctbto and not self.ctbto_calper:
            if self.mode == 'config-file':
                self.ctbto_calper = self._config['ctbto_calper']
                col_res = SelectResult.ok
            else:
                col_res = self.enter_opsr()
                if col_res == SelectResult.goback:
                    self.collect_backup()
                    col_res = SelectResult.ok

        return col_res, col_msg

    def collect_info(self):

        col_res = SelectResult.noop
        while col_res not in [SelectResult.quit, SelectResult.error] and not self.is_complete():
            col_res, col_msg = self.collect_next()
            if col_res == SelectResult.error:
                print(bold(red(self.tui_indent_str + col_msg)))
                self.collect_backup()

        return col_res == SelectResult.ok

    def print_info(self):

        staloc = self.sta.upper() + '-' + self.loc
        staloc = bold('{:<8}'.format(staloc))
        opsr =   bold('{:<4}'.format(self.opsr))
        comp =   bold('{:<3}'.format(self.comp))
        chan =   bold('{:<6}'.format(self.chan))
        ctbto =  bold('{:<3}/{}'.format(self.ctbto, self.ctbto_calper))
        sen =    bold('{:<10}'.format(self.sensor))
        lffile = bold('{}'.format(self.lffile))
        hffile = bold('{}'.format(self.hffile))
        lfpath = bold('{}'.format(self.lfpath))
        hfpath = bold('{}'.format(self.hfpath))
        respfn = bold('{}'.format(self.respfn))
        lfpert = bold('{} / {}'.format(self.fullpaz.zeros()[self.lfpert[1]] if self.fullpaz and self.lfpert else '',
                                  self.fullpaz.poles()[self.lfpert[0]] if self.fullpaz and self.lfpert else '',))
        hfpert = bold('{} / {}'.format(self.fullpaz.zeros()[self.hfpert[1]] if self.fullpaz and self.hfpert else '',
                                  self.fullpaz.poles()[self.hfpert[0]] if self.fullpaz and self.hfpert else '',))

        # have to do the format first because fabulous.color strings not always handled properly by format().
        stalbl =    '{:>8}'.format('Sta-Loc:')
        opsrlbl =   '{:>3}'.format('SR:')
        complbl =   '{:>4}'.format('Comp:')
        chanlbl =   '{:>4}'.format('Chan:')
        senlbl =    '{:>4}'.format('Sen:')
        lffilelbl = '{:>8}'.format('LF File:')
        hffilelbl = '{:>8}'.format('HF File:')
        lfpathlbl = '{:>8}'.format('LF Path:')
        hfpathlbl = '{:>8}'.format('HF Path:')
        ctbtolbl =  '{:>5}'.format('CTBTO:')
        respfnlbl = '{:>16}'.format('Starting Resp:')
        lfzpertlbl ='{:>16}'.format('LF Pert [z]/[p]:')
        hfzpertlbl ='{:>16}'.format('HF Pert [z]/[p]:')

        stalbl =    str(bold(green(stalbl))) if self.sta and self.loc else str(bold(red(stalbl)))
        opsrlbl =    str(bold(green(opsrlbl))) if self.opsr else str(bold(red(opsrlbl)))
        senlbl =    str(bold(green(senlbl))) if self.sensor else str(bold(red(senlbl)))
        complbl =   str(bold(green(complbl))) if self.comp else str(bold(red(complbl)))
        chanlbl =   str(bold(green(chanlbl))) if (self.chan or self.ctbto == 'N') else str(bold(red(chanlbl)))
        ctbtolbl =  str(bold(green(ctbtolbl))) if self.ctbto else str(bold(red(ctbtolbl)))
        lffilelbl = str(bold(green(lffilelbl))) if self.lffile else str(bold(red(lffilelbl)))
        hffilelbl = str(bold(green(hffilelbl))) if self.hffile else str(bold(red(hffilelbl)))
        lfpathlbl = str(bold(green(lfpathlbl))) if self.lfpath else str(bold(red(lfpathlbl)))
        hfpathlbl = str(bold(green(hfpathlbl))) if self.hfpath else str(bold(red(hfpathlbl)))
        respfnlbl = str(bold(green(respfnlbl))) if self.respfn else str(bold(red(respfnlbl)))
        lfzpertlbl = str(bold(green(lfzpertlbl))) if self.lfpert else str(bold(red(lfzpertlbl)))
        hfzpertlbl = str(bold(green(hfzpertlbl))) if self.hfpert else str(bold(red(hfzpertlbl)))

        banner = self.tui_indent_str + bold(blue('='*70))
        print(banner)
        print(self.tui_indent_str, stalbl, staloc,
              opsrlbl, opsr,
              complbl, comp,
              chanlbl, chan,
              ctbtolbl, ctbto,
              senlbl, sen)
        print(self.tui_indent_str, lfpathlbl, lfpath)
        print(self.tui_indent_str, lffilelbl, lffile)
        print(self.tui_indent_str, hfpathlbl, hfpath)
        print(self.tui_indent_str, hffilelbl, hffile)
        print(self.tui_indent_str, respfnlbl, respfn)
        print(self.tui_indent_str, lfzpertlbl, lfpert)
        print(self.tui_indent_str, hfzpertlbl, hfpert)
        print(banner)


    def __str__(self):

        staloc = '{:>8}'.format('Sta-Loc:') + '{:<8}'.format(self.sta.upper() + '-' + self.loc)
        opsr =   '{:>3}'.format('SR:') + '{:<4}'.format(self.opsr)
        comp =   '{:>4}'.format('Comp:') + '{:<3}'.format(self.comp)
        chan =   '{:>4}'.format('Chan:') + '{:<6}'.format(self.chan)
        ctbto =  '{:>5}'.format('CTBTO:') + '{:<8}'.format(self.ctbto)
        sen =    '{:>4}'.format('Sen:') + '{:<10}'.format(self.sensor.upper())
        lffile = '{:>8}'.format('LF File:') + '{}'.format(self.lffile)
        hffile = '{:>8}'.format('HF File:') + '{}'.format(self.hffile)
        lfpath = '{:>8}'.format('LF Path:') + '{}'.format(self.lfpath)
        hfpath = '{:>8}'.format('HF Path:') + '{}'.format(self.hfpath)
        respfn = '{:>16}'.format('Starting Resp:') + '{}'.format(self.respfn)
        lfpert = '{:>16}'.format('LF Pert [z]/[p]:') + \
                  '{} / {}'.format(self.fullpaz.zeros()[self.lfpert[1]] if self.fullpaz and self.lfpert else '',
                                  self.fullpaz.poles()[self.lfpert[0]] if self.fullpaz and self.lfpert else '',)
        hfpert = '{:>16}'.format('HF Pert [z]/[p]:') + \
                 '{} / {}'.format(self.fullpaz.zeros()[self.hfpert[1]] if self.fullpaz and self.hfpert else '',
                                  self.fullpaz.poles()[self.hfpert[0]] if self.fullpaz and self.hfpert else '',)

        now = datetime.datetime.now(datetime.timezone(datetime.timedelta(0)))
        banner = '=' * 70 + '\n'
        txt = banner
        # noinspection PyArgumentList
        txt = txt + 'Processed on: ' + now.astimezone().strftime('%Y-%m-%d %H:%M:%S %Z') + \
              'by user: ' + environ['USER'] + '\n'
        txt = txt + staloc + opsr + comp + chan + ctbto + sen + '\n'
        txt = txt + lfpath + '\n'
        txt = txt + lffile + '\n'
        txt = txt + hfpath + '\n'
        txt = txt + hffile + '\n'
        txt = txt + respfn + '\n'
        txt = txt + lfpert + '\n'
        txt = txt + hfpert + '\n'
        txt = txt + banner + '\n'

        return txt
