import glob
import os
from pathlib import Path

from ida.tui import pick
from ida import IDA_CAL_RAW_DIR, \
    IDA_RESPONSES_CUR_DIR, \
    IDA_RESPONSES_NOM_DIR, \
    IDA_DATASCOPEDB_DIR


def select_raw_cal_sensor(station):

    station = station.lower()
    loc, sensor_model, sensor_root_dir = None, None, None

    sensordirlist = glob.glob(os.path.join(IDA_CAL_RAW_DIR, station + '??', '*'))
    sensordirlist = sorted(['/'.join(item.split(os.sep)[-2:]) for item in sensordirlist])
    # rawdirlist.append('Quit')

    if sensordirlist > 0:
        success, ndx = pick(sensordirlist,
                            title='Sensor List for ' + station.upper(),
                            prompt='Select # of desired SENSOR (or "q" to quit): ',
                            allow_quit_q=True,
                            menu_on_error=True,
                            err_message='Invalid selection. Please try again.')
    else:
        success = False

    if success:
        sensor_dir = sensordirlist[ndx]
        loc = sensor_dir.split('/')[0][-2:]
        sensor_model = sensor_dir.split(os.sep)[1]
        sensor_root_dir = os.path.join(IDA_CAL_RAW_DIR, sensor_dir)

    return success, station, loc, sensor_model, sensor_root_dir


def select_raw_cal_date(sensor_root_dir, cal_type):

    if not os.path.exists(sensor_root_dir):
        raise ValueError('Directory does not exist: '+ sensor_root_dir)

    date_dir, date_str = None, None

    # get date dir list
    datedirlist = glob.glob(os.path.join(sensor_root_dir, cal_type, '*'))
    datedirlist = sorted([os.sep.join(item.split(os.sep)[-2:]) for item in datedirlist])

    if len(datedirlist) == 0:
        suc = False
        date_str = 'No raw calibration dates found for this sensor.'

    elif len(datedirlist) > 1:
        didquit, ndx = pick(datedirlist, 
                        title='RBLF dates for sensor: ' + os.sep.join(sensor_root_dir.split(os.sep)[-2:]),
                        prompt='Select # of desired DATE (or "q" to quit): ',
                        allow_quit_q=True, 
                        menu_on_error=True,
                        err_message='Invalid selection. Please try again.')
        suc = not didquit
        if suc:
            date_dir = datedirlist[ndx]
        else:
            date_str = 'User canceled.'

    elif len(datedirlist) == 1:
        suc = True
        date_dir = datedirlist[0]

    # make full path
    if suc: 
        date_str = date_dir.split(os.sep)[1]
        date_dir = os.path.join(sensor_root_dir, date_dir)

    return suc, date_str, date_dir


def select_raw_cal_files(rb_dir):

    if not os.path.exists(rb_dir):
        raise ValueError('Directory does not exist: '+ rb_dir)

    # get ms and log files and compare lists
    ms_files = sorted(glob.glob(os.path.join(rb_dir, '*.ms')))
    log_files = sorted(glob.glob(os.path.join(rb_dir, '*.log')))
    ms_stems = [Path(msfn).stem for msfn in ms_files]
    log_stems = [Path(logfn).stem for logfn in log_files]

    paired_cals = [stem for stem in ms_stems if stem in log_stems]

    rb_files = None

    if len(ms_files) == 0:
        rb_files = (None, None)
        suc = True

    elif len(ms_files) == 1:
        suc = True
        rb_files = (os.path.join(rb_dir, ms_files[0]), os.path.join(rb_dir, log_files[0]))

    else:  # need to pick
        didquit, ndx = pick(paired_cals, 
                            title='Miniseed files found in: ' + os.sep.join(rb_dir.split(os.sep)[-2:]),
                            prompt='Select # of the file to process (or "q" to quit): ',
                            allow_quit_q=True, 
                            menu_on_error=True,
                            err_message='Invalid selection. Please try again.')
        suc = not didquit
        if suc:
            rb_files = (os.path.join(rb_dir, paired_cals[ndx]+'.ms'), os.path.join(rb_dir, paired_cals[ndx]+'.log'))


    return suc, rb_files


def select_perturb_map(paz):

    if not isinstance(paz, PAZ):
        raise TypeError('paz must be a populated PAZ object')

    paz_fit_lf = paz.make_partial2(norm_freq=1.0, partial_mode=paz.PARTIAL_FITTING_LF)
    paz_fit_hf = paz.make_partial2(norm_freq=1.0, partial_mode=paz.PARTIAL_FITTING_HF)

    poles_pert_def, zeros_pert_def = paz.perturb_defaults()
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

    pert_choices = [defchoice, zero_pert_choices, pole_pert_choices]
    success, choices, pert_choice_groups = pick2(pert_choices, 'Select LOW Freq zeros & poles to perturb',
                                                 prompt='Enter selection (or "q" to quit): ',
                                                 group_titles=['',
                                                               'LOW Freq Zeros',
                                                               'LOW Freq Poles'],
                                                 multiple_choice=True,
                                                 implicit_quit_q=True, menu_on_error=True)

    # print(success, pert_choice_groups)
    if not success:
        return False, None, None

    if choices[0].upper() == 'D':  # using defaults
        lf_map = (poles_pert_def[0], zeros_pert_def[0])  # beware, put poles then zeros in this map tuple
    else:
        lf_map = (pert_choice_groups[2], pert_choice_groups[1])


    # NOW HIGH FREQ
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
    success, choices, pert_choice_groups = pick2(pert_choices, 'Select HIGH Freq zeros & poles to perturb',
                                                 prompt='Enter selection (or "q" to quit): ',
                                                 group_titles=['',
                                                               'HIGH Freq Zeros',
                                                               'HIGH Freq Poles'],
                                                 multiple_choice=True,
                                                 implicit_quit_q=True, menu_on_error=True)

    # print(success, pert_choice_groups)

    if not success:
        return False, None, None

    if choices[0].upper() == 'D':  # using defaults
        hf_map = (poles_pert_def[1], zeros_pert_def[1])  # beware, put poles then zeros in this map tuple
    else:
        hf_map = (pert_choice_groups[2], pert_choice_groups[1])


    return success, lf_map, hf_map  # each in (poles, zeros) order
