#!/usr/bin/env python

# from time import strptime
import argparse
import glob
from os import getcwd
from os.path import exists, abspath, join
from pathlib import Path

from fabulous.color import bold, blue, red, green
from pandas.core.frame import DataFrame

from ida.tui import pick2, PickResult
from ida.instruments import CALTYPE_RBLF, CALTYPE_RBHF
from ida.calibration import nom_resp_for_model, local_resp_files
from ida.signals.paz import PAZ
from ida.db.io import read
from ida.db.query import find_sensor_file, get_stages

VERSION = '0.1'


class CalInfo():

    STARTING_GROUP = 1
    FITTING_GROUP = 2

    tui_indent = 4
    tui_indent_str = ' '*tui_indent

    def __init__(self, cal_raw_dir, resp_cur_dir, resp_nom_dir, cal_analysis_dir):
        if exists(abspath(cal_raw_dir)):
            self.cal_raw_dir = cal_raw_dir
        else:
            self.cal_raw_dir = ''
            raise ValueError(self.tui_indent_str +
                             bold(red('Error: IDA Cal Raw directory [' + cal_raw_dir + '] does not exist.')))

        if exists(abspath(resp_cur_dir)):
            self.resp_cur_dir = resp_cur_dir
        else:
            raise ValueError(self.tui_indent_str +
                             bold(red('Error: Current IDA Response directory [' + resp_cur_dir + '] does not exist.')))

        if exists(abspath(resp_nom_dir)):
            self.resp_nom_dir = resp_nom_dir
        else:
            raise ValueError(self.tui_indent_str +
                             bold(red('Error: IDA NOMINAL Response irectory [' + resp_nom_dir + '] does not exist.')))

        if exists(abspath(cal_analysis_dir)):
            self.cal_analysis_dir = cal_analysis_dir
        else:
            self.cal_analysis_dir = ''
            raise ValueError(self.tui_indent_str +
                             bold(red('Error: IDA Cal Analysis directory [' + cal_analysis_dir + '] does not exist.')))

        self._info = {
            'staloc': None,
            'sensor': None,
            'comp': None,
            'ctbto': None,
            'chan': None,
            'opsr': None, #todo: query user or DB for this
            'lfdate': None,
            'hfdate': None,
            'lffile': None,
            'hffile': None,
            'lfpath': None,
            'hfpath': None,
            'respfn': None,
            'fullpaz': None,
            'newpaz': None,
            'lfpert': None,
            'hfpert': None,
        }

        _, self.stages_df = read('datascope', 'stage')
        self.chn_stages = []

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
            raise TypeError('Invalid type for ' + valname + '. Should be ' + str(valtype))

    @property
    def staloc(self):
        return self._info['staloc'] if self._info['staloc'] else ''

    @staloc.setter
    def staloc(self, value):

        self.reset_all()
        if value:
            self.check_type(value, str, 'staloc')
            self._info['staloc'] = value if exists(join(self.cal_raw_dir, value)) else None
        else:
            self._info['staloc'] = None

    @property
    def sta(self):
        if self.staloc and len(self.staloc) > 2:
            return self.staloc[:-2]
        else:
            return ''

    @property
    def loc(self):
        if self.staloc and len(self.staloc) > 2:
            return self.staloc[-2:]
        else:
            return ''

    @property
    def sensor(self):
        return self._info['sensor'] if self._info['sensor'] else ''

    @sensor.setter
    def sensor(self, value):

        if self._info['staloc']:
            if value:
                self.check_type(value, str, 'sensor')
                self._info['sensor'] = value if exists(join(self.cal_raw_dir,
                                                            self.staloc,
                                                            value)) else None
            else:
                self._info['sensor'] = None

            self.reset_all_except(['staloc', 'sensor', 'comp',' chan'])

    @property
    def comp(self):
        return self._info['comp'] if self._info['comp'] else ''

    @comp.setter
    def comp(self, value):

        if value:
            self.check_type(value, str, 'component')
            value = value.upper()
            if value in ['Z', 'N', 'E']:
                self._info['comp']= value
                self.reset(['chan'])
            else:
                raise ValueError('Invalid value for Component: ' + value + ". Should be one of ['Z', 'N', 'E']")
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

        self.check_type(value, (int, float), 'operating sample rate')
        self._info['opsr'] = value

    @property
    def lfdate(self):
        return self._info['lfdate'] if self._info['lfdate'] else ''

    @lfdate.setter
    def lfdate(self, value):

        if self.staloc and self.sensor:
            if value:
                self.check_type(value, str, 'low freq cal date')
                rbpath = join(self.cal_raw_dir,
                              self.staloc,
                              self.sensor,
                              CALTYPE_RBLF,
                              value)
                if exists(rbpath):
                    if self.lfdate != value:
                        self.reset(['lffile'])

                    self._info['lfdate'] = value
                    self.reset(['lffile', 'respfn', 'fullpaz', 'lfpert', 'hfpert'])

                    ms_files = glob.glob(join(rbpath, '*.ms'))
                    log_files = glob.glob(join(rbpath, '*.log'))
                    ms_stems = [Path(msfn).stem for msfn in ms_files]
                    log_stems = [Path(logfn).stem for logfn in log_files]
                    paired_files = [stem for stem in ms_stems if stem in log_stems]
                    if len(paired_files) == 1:
                        self.lffile = paired_files[0]

            else:
                self._info['lfdate'] = None
                self.reset(['lffile', 'respfn', 'fullpaz', 'lfpert', 'hfpert'])


    @property
    def hfdate(self):
        return self._info['hfdate'] if self._info['hfdate'] else ''

    @hfdate.setter
    def hfdate(self, value):

        self.check_type(value, str, 'high freq cal date')

        if self.staloc and self.sensor:
            if value:
                rbpath = join(self.cal_raw_dir,
                                      self.staloc,
                                      self.sensor,
                                      CALTYPE_RBHF,
                                      value)
                if exists(rbpath):
                    if self.hfdate != value:
                        self.reset(['hffile'])

                    self._info['hfdate'] = value
                    self.reset(['hffile'])

                    ms_files = glob.glob(join(rbpath, '*.ms'))
                    log_files = glob.glob(join(rbpath, '*.log'))
                    ms_stems = [Path(msfn).stem for msfn in ms_files]
                    log_stems = [Path(logfn).stem for logfn in log_files]
                    paired_files = [stem for stem in ms_stems if stem in log_stems]
                    if len(paired_files) == 1:
                        self.hffile = paired_files[0]

            else:
                self._info['hfdate'] = None
                self.reset(['hffile'])

    @property
    def lffile(self):
        return self._info['lffile']if self._info['lffile'] else ''

    @lffile.setter
    def lffile(self, value):

        if self.staloc and self.sensor and self.lfdate:

            if value:
                self.check_type(value, str, 'low freq cal data filename')
                fpath = join(self.cal_raw_dir, self.staloc, self.sensor,
                             CALTYPE_RBLF, self.lfdate)
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

        if self.staloc and self.sensor and self.hfdate:

            if value:
                self.check_type(value, str, 'high freq cal data filename')
                fpath = join(self.cal_raw_dir, self.staloc, self.sensor,
                             CALTYPE_RBHF, self.hfdate)
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

        if self.staloc and self.comp and self.sensor and self.lfdate:

            if value:
                self.check_type(value, str, 'response filename')
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

    def select_raw_cal_staloc(self, staloc_stub):

        staloc_stub = staloc_stub.lower()
        stalocdirlist = glob.glob(join(self.cal_raw_dir, staloc_stub + '*[01][01]'))
        stalocdirlist = sorted([Path(item).stem for item in stalocdirlist])

        errmsg = ''
        if stalocdirlist:
            self.print_info()
            result, _, user_choice_groups, _ = pick2([stalocdirlist],
                                                              title='Select Raw Cal Sta+Loc',
                                                              prompt='Select # of desired StaLoc (or "q" to quit): ',
                                                              implicit_quit_q=True, menu_on_error=True,
                                                              err_message='Invalid selection. Please try again.',
                                                              indent_width=self.tui_indent)
            if result == PickResult.collect_ok:
                ndx = user_choice_groups[0][0]
                self.staloc = stalocdirlist[ndx]

        else:
            result = PickResult.collect_error
            errmsg = 'No station directories found in: ' + abspath(self.cal_raw_dir)

        return result, errmsg

    def select_raw_cal_sensor(self):

        staloc = self._info['staloc']
        sensordirlist = glob.glob(join(self.cal_raw_dir, staloc, '*'))
        sensordirlist = sorted([Path(item).stem for item in sensordirlist])
        list_len = len(sensordirlist)

        errmsg = ''
        if list_len > 1:
            self.print_info()
            result, user_choices, user_choice_groups, _ = pick2([sensordirlist],
                                                              title='Select Sensor',
                                                              prompt='Select # of desired SENSOR ("q" => quit, "b" => back): ',
                                                              implicit_quit_q=True, implicit_back_b=True,
                                                              menu_on_error=True, err_message='Invalid selection. Please try again.',
                                                              indent_width=self.tui_indent)
            if result == PickResult.collect_ok:
                ndx = user_choice_groups[0][0]
                self.sensor = sensordirlist[ndx]

        elif list_len == 1:
            self.sensor = sensordirlist[0]
            result = PickResult.collect_ok
        else:
            result = PickResult.collect_error
            errmsg = 'No sensor directories found in: ' + abspath(join(self.cal_raw_dir, staloc))

        return result, errmsg

    def select_component(self):

        self.print_info()

        complist = [('Z', 'Vertical'),
                    ('N', 'North/South'),
                    ('E', 'East/West')]
        result, _, user_choice_groups, _ = pick2([complist],
                                               title='Select Component',
                                               prompt='Enter selection ("q" => quit, "b" => back): ',
                                               implicit_quit_q=True, implicit_back_b=True,
                                               menu_on_error=True, err_message='Invalid choice. Please try again',
                                               indent_width=self.tui_indent)

        if result == PickResult.collect_ok:
            ndx = user_choice_groups[0][0]
            self.comp = complist[ndx][0]

        return result

    def select_ctbto_flag(self):

        self.print_info()

        ctbto_dict = {
            'Y': 'Yes: Sensor requires a CTBTO Calibration Report',
            'N': 'No: Not a CTBTO sensor'
        }

        flaglist = [item for item in ctbto_dict.items()]
        result, _, user_choice_groups, _ = pick2([flaglist], 'CTBTO Status',
                                               prompt='Is this a CTBTO sensor ("q" => quit, "b" => back)? ',
                                               implicit_quit_q=True, implicit_back_b=True,
                                               menu_on_error=True, err_message='Invalid choice. Please try again',
                                               indent_width=self.tui_indent)

        if result == PickResult.collect_ok:
            ndx = user_choice_groups[0][0]
            self.ctbto = flaglist[ndx][0]

        return result

    def enter_chan(self):

        self.print_info()

        result = PickResult.collect_noop
        while result not in [PickResult.collect_quit,
                             PickResult.collect_back,
                             PickResult.collect_ok]:
            print()
            chan = input(bold(blue(self.tui_indent_str + 'Enter 3 character channel code: ')))
            chan = chan.upper()
            if (len(chan) == 3):
                is_valid = (chan[2] == self.comp) or \
                           (chan[2] in ['1', 'N'] and self.comp == 'N') or \
                           (chan[2] in ['2', 'E'] and self.comp == 'E')
                if is_valid:
                    self.chan = chan
                    result = PickResult.collect_ok
                else:
                    result = PickResult.collect_error

            elif len(chan) == 1:
                if chan == 'Q':
                    result = PickResult.collect_quit
                elif chan == 'B':
                    result = PickResult.collect_back
                else:
                    result = PickResult.collect_error

            if result == PickResult.collect_error:
                print()
                print(self.tui_indent_str + bold(red('Invalid channel code.')))
                print(self.tui_indent_str + bold(red('Make sure 3 characters, matches component "{}" and try again'.format(self.comp))))
                print()

        return result

    def select_raw_cal_date(self, cal_type):

        title_dict = {
            CALTYPE_RBLF: 'Low Freq',
            CALTYPE_RBHF: 'High Freq'
        }

        if cal_type not in [CALTYPE_RBLF, CALTYPE_RBHF]:
            raise ValueError('Invalid CAL TYPE supplied: ' + cal_type)

        datedirlist = glob.glob(join(self.cal_raw_dir, self.staloc, self.sensor, cal_type, '*'))
        datedirlist = sorted([Path(item).stem for item in datedirlist])

        errmsg = ''
        if len(datedirlist) > 1:
            self.print_info()
            result, user_choices, user_choice_groups, _ = pick2([datedirlist],
                                                    title=title_dict[cal_type] + ' dates for {} {}: '.format(self.staloc,
                                                                                                            self.sensor),
                                                    prompt='Select # of desired DATE ("q" => quit, "b" => back): ',
                                                    implicit_quit_q=True,
                                                    menu_on_error=True, err_message='Invalid selection. Please try again.',
                                                    indent_width=self.tui_indent)
            if result == PickResult.collect_ok:
                    ndx = user_choice_groups[0][0]
                    if cal_type == CALTYPE_RBLF:
                        self.lfdate = datedirlist[ndx]
                    else:
                        self.hfdate = datedirlist[ndx]

        elif len(datedirlist) == 1:
            result = PickResult.collect_ok
            if cal_type == CALTYPE_RBLF:
                self.lfdate = datedirlist[0]
            else:
                self.hfdate = datedirlist[0]

        else:
            result = PickResult.collect_error
            errmsg = 'No raw {} calibration dates found for sensor {}.'.format(title_dict[cal_type].upper(),
                                                                               self.sensor.upper())

        return result, errmsg

    def select_starting_response_file(self):

        # run for IDA deployment ONLY

        pick_groups = []
        group_titles = []

        if isinstance(self.chn_stages, DataFrame) and not self.chn_stages.empty:
            resp_file = self.chn_stages.iloc[0].dfile
            if resp_file:
                sensor_fn_list = [resp_file]
            else:
                sensor_fn_list = []
        else:
            sensor_fn_list = []
        pick_groups.append(sensor_fn_list)
        group_titles.append('Response on Date of Calibration')

        nom_resps = nom_resp_for_model(self.sensor.lower())
        nom_resps_names = sorted([Path(respfn).name for respfn in nom_resps])
        pick_groups.append(nom_resps_names)
        group_titles.append('Nominal Responses for sensor: ' + self.sensor.upper())

        local_resp = local_resp_files()
        local_resp_names = sorted([Path(respfn).name for respfn in local_resp])
        pick_groups.append(local_resp_names)
        group_titles.append('Responses Found in CWD')

        errmsg = ''
        if len(pick_groups) > 0:
            self.print_info()
            result, user_choices, _, choice_tpls = pick2(pick_groups, 'Select Response file to use as starting model',
                                                         prompt='Select # of desired DATE ("q" => quit, "b" => back): ',
                                                         group_titles=group_titles,
                                                         implicit_quit_q=True, implicit_back_b=True,
                                                         menu_on_error=True, err_message='Invalid choice. Please try again',
                                                         indent_width=self.tui_indent)

            if result == PickResult.collect_ok:
                fn = pick_groups[choice_tpls[0][0]][choice_tpls[0][1]]
                if choice_tpls[0][0] == 0:     #  first group is single current resp file as in DB
                    self.respfn = join(self.resp_cur_dir, fn)
                elif choice_tpls[0][0] == 1:   #  second group are files from nominal response dir
                    self.respfn = join(self.resp_nom_dir, fn)
                elif choice_tpls[0][0] == 2:   #  third group are respo files found in cwd
                    self.respfn = join(getcwd(), fn)

                # read in resp data and make sure file is valid (will raise an eception if not)
                self.fullpaz = PAZ(pzfilename=self.respfn, fileformat='ida', mode='vel', units='hz')

        else:
            result = PickResult.collect_error
            errmsg = 'No response files found in DB, Resp Dir or cwd.'

        return result, errmsg

    def select_perturb_map(self):

        grp_titles = ['', 'Zeros', 'Poles']
        prompt = 'Enter comma separated list ("q" => quit, "b" => back): '

        if not isinstance(self.fullpaz, PAZ):
            raise TypeError('fullpaz must be a populated PAZ object')

        paz_fit_lf = self.fullpaz.make_partial2(norm_freq=1.0, partial_mode=self.fullpaz.PARTIAL_FITTING_LF)

        poles_pert_def, zeros_pert_def = self.fullpaz.perturb_defaults()
        defchoice = [('D', 'Use Defaults (indicated by "<==")')]

        # make list of LF p/z to user perturbing choices
        # Do LOW FREQ FIRST
        zero_pert_choices = []
        pole_pert_choices = []
        for ndx, val in enumerate(paz_fit_lf.zeros()):
            if ndx in zeros_pert_def[0]:
                zero_pert_choices.append(str(val) + ' <==')
            else:
                zero_pert_choices.append(str(val))
        for ndx, val in enumerate(paz_fit_lf.poles()):
            if ndx in poles_pert_def[0]:
                pole_pert_choices.append(str(val) + ' <==')
            else:
                pole_pert_choices.append(str(val))

        errmsg = ''
        pert_choices = [defchoice, zero_pert_choices, pole_pert_choices]
        result, choices, pert_choice_groups, _ = pick2(pert_choices,
                                                        'Select LOW FREQ zeros & poles to perturb',
                                                        prompt=prompt, group_titles=grp_titles,
                                                        multiple_choice=True,
                                                        implicit_quit_q=True, implicit_back_b=True,
                                                        menu_on_error=True, err_message='Invalid entry. Please try again',
                                                        indent_width=self.tui_indent)

        if result == PickResult.collect_ok:
            if choices[0].upper() == 'D':  # using defaults
                self.lfpert = (poles_pert_def[0], zeros_pert_def[0])  # beware, put poles then zeros in this map tuple
            else:
                self.lfpert = (pert_choice_groups[2], pert_choice_groups[1])

            # NOW HIGH FREQ
            paz_fit_hf = self.fullpaz.make_partial2(norm_freq=1.0, partial_mode=self.fullpaz.PARTIAL_FITTING_HF)

            zero_pert_choices = []
            pole_pert_choices = []
            for ndx, val in enumerate(paz_fit_hf.zeros()):
                if ndx in zeros_pert_def[1]:
                    zero_pert_choices.append(str(val) + ' <==')
                else:
                    zero_pert_choices.append(str(val))
            for ndx, val in enumerate(paz_fit_hf.poles()):
                if ndx in poles_pert_def[1]:
                    pole_pert_choices.append(str(val) + ' <==')
                else:
                    pole_pert_choices.append(str(val))

            pert_choices = [defchoice, zero_pert_choices, pole_pert_choices]
            success, choices, pert_choice_groups, _ = pick2(pert_choices,
                                                            'Select HIGH FREQ zeros & poles to perturb',
                                                            prompt=prompt, group_titles=grp_titles,
                                                            multiple_choice=True,
                                                            implicit_quit_q=True, implicit_back_b=True,
                                                            menu_on_error=True,
                                                            err_message='Invalid entry. Please try again',
                                                            indent_width=self.tui_indent)
            # print(success, pert_choice_groups)
            if result == PickResult.collect_ok:
                if choices[0].upper() == 'D':  # using defaults
                    self.hfpert = (poles_pert_def[1], zeros_pert_def[1])  # beware, put poles then zeros in this map tuple
                else:
                    self.hfpert = (pert_choice_groups[2], pert_choice_groups[1])

        return result, errmsg

    def sensor_cnt(self):

        if self.staloc:
            return len(glob.glob(join(self.cal_raw_dir, self.staloc, '*')))
        else:
            return 0

    def date_cnt(self, cal_type):
        if cal_type not in [CALTYPE_RBLF, CALTYPE_RBHF]:
            raise ValueError('Invalid CAL TYPE supplied: ' + cal_type)

        if self.staloc and self.sensor:
            return len(glob.glob(join(self.cal_raw_dir, self.staloc, self.sensor, cal_type, '*')))
        else:
            return 0

    def find_qcal_files(self, cal_type):

        if cal_type not in [CALTYPE_RBLF, CALTYPE_RBHF]:
            raise ValueError('Invalid CAL TYPE supplied: ' + cal_type)

        errmsg = ''
        if cal_type == CALTYPE_RBLF:
            adate = self.lfdate
        else:
            adate = self.hfdate

        result = PickResult.collect_noop
        if self.staloc and self.sensor and adate:
            dpath = abspath(join(self.cal_raw_dir, self.staloc, self.sensor, cal_type, adate))

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
                    result = PickResult.collect_ok
                else:
                    result = PickResult.collect_error
                    errmsg = 'Error: Missing or multiple MS and LOG file(s) in directory: ' + dpath
            else:
                result = PickResult.collect_error
                errmsg = 'Directory does not exist: ' + dpath

        return result, errmsg

    def new_resp_filename_stem(self):
        if self.is_complete(self.FITTING_GROUP):
            stem = self.sensor.lower() + '_' + self.staloc.lower() + '_' + self.lfdate + "_" + self.comp.upper()
        else:
            stem = 'INCOMPLETE_CAL_INFO'

        return stem

    def is_complete(self, group):

        if group == self.STARTING_GROUP:
            done = self.staloc and self.sensor and self.opsr and \
                   self.comp and self.chan and \
                   self.lfdate and self.hfdate and \
                   self.lffile and self.hffile and \
                   self.respfn and self.fullpaz
        elif group == self.FITTING_GROUP:
            done = self.staloc and self.sensor and self.opsr and \
                   self.comp and self.chan and \
                   self.lfdate and self.hfdate and \
                   self.lffile and self.hffile and \
                   self.respfn and self.fullpaz and \
                   self.lfpert and self.hfpert and self.ctbto

        return done

    def collect_backup(self):

        # in reverse order of querying user
        done = False

        if not done and self.ctbto:
            self.ctbto = None
            done = True

        if self.lfpert or self.hfpert:
            self.lfpert = None
            self.hfpert = None
            done = True

        if self.respfn:
            self.respfn = None
            self.fullpaz = None
            done = True

        if not done and self.chan:
            self.chan = None
            self.chn_stages = None  # clear stages since channel changing
            done = True

        if not done and self.comp:
            self.comp = None
            done = True

        if not done and self.hfdate:
            self.hfdate = None
            self.hffile = None
            self._info['hfpath'] = None
            done = self.date_cnt(CALTYPE_RBHF) > 1

        if not done and self.lfdate:
            self.lfdate = None
            self.lffile = None
            self._info['lfpath'] = None
            done = self.date_cnt(CALTYPE_RBLF) > 1

        if not done and self.sensor:
            self.sensor = None
            done = self.sensor_cnt() > 1

        if not done:
            self.staloc = None
            done = True


    def collect_next(self):

        col_res = PickResult.collect_noop
        col_msg = ''

        if not self.staloc:
            col_res, col_msg = self.select_raw_cal_staloc('')

        elif not self.sensor:
            col_res, col_msg = self.select_raw_cal_sensor()
            if col_res == PickResult.collect_back:
                self.staloc = ''
                col_res = PickResult.collect_ok

        elif not self.lfdate:
            col_res, col_msg = self.select_raw_cal_date(CALTYPE_RBLF)
            if col_res == PickResult.collect_back:
                self.collect_backup()
                col_res = PickResult.collect_ok

        elif not self.hfdate:
            col_res, col_msg = self.select_raw_cal_date(CALTYPE_RBHF)
            if col_res == PickResult.collect_back:
                self.collect_backup()
                col_res = PickResult.collect_ok

        elif not self.lffile:
            col_res, col_msg = self.find_qcal_files(CALTYPE_RBLF)

        elif not self.hffile:
            col_res, col_msg = self.find_qcal_files(CALTYPE_RBHF)

        elif not self.comp:
            col_res = self.select_component()
            if col_res == PickResult.collect_back:
                self.collect_backup()
                col_res = PickResult.collect_ok

        elif not self.chan:
            col_res = self.enter_chan()
            if col_res == PickResult.collect_back:
                self.collect_backup()
                col_res = PickResult.collect_ok
            else:
                # good channel, now lets get stages info form db for later use
                self.chn_stages = get_stages(self.stages_df, self.sta, self.loc, self.chan, self.lfdate)
                if isinstance(self.chn_stages, DataFrame) and not self.chn_stages.empty:
                    # get channel freq from stage 3 srate column
                    # ToDo: bury the srate ref in public acces method to make db independent
                    self.opsr = self.chn_stages.iloc[2].srate

        elif not self.respfn:
            col_res, col_msg = self.select_starting_response_file()
            if col_res == PickResult.collect_back:
                self.collect_backup()
                col_res = PickResult.collect_ok

        elif not (self.lfpert and self.hfpert):
            col_res, col_msg = self.select_perturb_map()
            if col_res == PickResult.collect_back:
                self.collect_backup()
                col_res = PickResult.collect_ok

        elif not self.ctbto:
            col_res= self.select_ctbto_flag()
            if col_res == PickResult.collect_back:
                self.collect_backup()
                col_res = PickResult.collect_ok


        return col_res, col_msg

    def collect(self, group):

        if not self.cal_raw_dir:
            return PickResult.collect_error, \
                   bold(red(self.tui_indent_str + 'Calibration raw data directory root not set.'))
        elif not exists(self.cal_raw_dir):
            return PickResult.collect_error, \
                   bold(red(self.tui_indent_str + 'Calibration raw data directory does not exist: '+self.cal_raw_dir))

        col_res = PickResult.collect_noop
        while col_res != PickResult.collect_quit and not self.is_complete(group):
            col_res, col_msg = self.collect_next()
            if col_res == PickResult.collect_error:
                print(bold(red(self.tui_indent_str + col_msg)))
                self.collect_backup()

        if col_res != PickResult.collect_quit:
            self.print_info()
            return True
        else:
            print('Looser User quit.')

    def log_msg(self):
        pass

    def print_info(self):

        if self.staloc:
            staloc = self.staloc[:-2].upper() + '-' + self.staloc[-2:]
        else:
            sta, loc = '', ''
            staloc = ''

        staloc = bold('{:<8}'.format(staloc))
        opsr =   bold('{:<4}'.format(self.opsr))
        comp =   bold('{:<3}'.format(self.comp))
        chan =   bold('{:<6}'.format(self.chan))
        ctbto =  bold('{:<3}'.format(self.ctbto))
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

        stalbl =    str(bold(green(stalbl))) if self.staloc else str(bold(red(stalbl)))
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
        # print(self.tui_indent_str, stalbl, sta, loclbl, loc, senlbl, sen)
        # print(self.tui_indent_str, complbl, comp, chanlbl, chan, ctbtolbl, ctbto)
        # print(self.tui_indent_str, lfdatelbl, lfdate, hfdatelbl, hfdate)
        print(self.tui_indent_str, lfpathlbl, lfpath)
        print(self.tui_indent_str, lffilelbl, lffile)
        print(self.tui_indent_str, hfpathlbl, hfpath)
        print(self.tui_indent_str, hffilelbl, hffile)
        print(self.tui_indent_str, respfnlbl, respfn)
        print(self.tui_indent_str, lfzpertlbl, lfpert)
        print(self.tui_indent_str, hfzpertlbl, hfpert)
        print(banner)


    def save_log(self, log_msg):
        pass


