import glob
import os
import re
from sys import exit
from pathlib import Path

from ida.utils import pick
from ida.calibration import IDA_CAL_RAW_DIR, IDA_RESPONSES_CUR_DIR, IDA_RESPONSES_NOM_DIR

def get_rb_cal_raw_dirs(station):

    station = station.lower()

    sensordirlist = glob.glob(os.path.join(IDA_CAL_RAW_DIR, station + '??', '*'))
    sensordirlist = sorted(['/'.join(item.split(os.sep)[-2:]) for item in sensordirlist])
    # rawdirlist.append('Quit')

    didquit, ndx = pick(sensordirlist, 
                        title='Sensor List for ' + station.upper(), 
                        prompt='Select # of desired SENSOR (or "q" to quit): ',
                        allow_quit_q=True, 
                        menu_on_error=True,
                        err_message='Invalid selection. Please try again.',
                        indent_width=8)

    if didquit: return None

    sensor_dir = sensordirlist[ndx]
    sensor_root_dir = os.path.join(IDA_CAL_RAW_DIR, sensor_dir)
    # select rblf date
    datedirlist = glob.glob(os.path.join(sensor_root_dir, 'rblf', '*'))
    datedirlist = sorted(['/'.join(item.split(os.sep)[-2:]) for item in datedirlist])

    didquit, ndx = pick(datedirlist, 
                        title='RBLF dates for sensor: ' + sensor_dir, 
                        prompt='Select # of desired DATE (or "q" to quit): ',
                        allow_quit_q=True, 
                        menu_on_error=True,
                        err_message='Invalid selection. Please try again.',
                        indent_width=8)

    if didquit: 
        return None
    else:
        rblf_dir = os.path.join(sensor_root_dir, datedirlist[ndx])

    # select rbhf date
    datedirlist = glob.glob(os.path.join(sensor_root_dir, 'rbhf', '*'))
    datedirlist = sorted(['/'.join(item.split(os.sep)[-2:]) for item in datedirlist])

    didquit, ndx = pick(datedirlist, 
                        title='RBHF dates for sensor: ' + sensor_dir, 
                        prompt='Select # of desired DATE (or "q" to quit): ',
                        allow_quit_q=True, 
                        menu_on_error=True,
                        err_message='Invalid selection. Please try again.',
                        indent_width=8)

    if didquit: 
        return None
    else:
        rbhf_dir = os.path.join(sensor_root_dir, datedirlist[ndx])

    return (rblf_dir, rbhf_dir)


def get_cal_raw_files(rb_dir):

    if not os.path.exists(rb_dir):
        return 'Directory not found: ' + rb_dir, True

    # get ms and log files and compare lists
    ms_files = sorted(glob.glob(os.path.join(rb_dir, '*.ms')))
    log_files = sorted(glob.glob(os.path.join(rb_dir, '*.log')))

    paired_cals = []

    for ndx, msfn in enumerate(ms_files):
        p = Path(msfn)
        stem = p.stem
        if log_files[ndx].startswith(stem):
            paired_cals.append(stem)

    if len(ms_files) > 1:
        didquit, ndx = pick(paired_cals, 
                            title='Cal data found in: ' + rb_dir, 
                            prompt='Select # of cal data to process (or "q" to quit): ',
                            allow_quit_q=True, 
                            menu_on_error=True,
                            err_message='Invalid selection. Please try again.',
                            indent_width=8)

        if didquit: 
            return 'User canceled.', True
        else:
            rb_files = (os.path.join(rb_dir, paired_cals[ndx]+'.ms'), os.path.join(rb_dir, paired_cals[ndx]+'.log'))

    else:
            rb_files = (os.path.join(rb_dir, ms_files[0]), os.path.join(rb_dir, log_files[0]))

    return rb_files, False


def get_sensor_response_file(station, sensor_model=None, component=None):

    station = station.lower()

    resplist = glob.glob(os.path.join(IDA_RESPONSES_CUR_DIR, '*' + station + '*'))
    resplist = sorted([resp.split(os.sep)[-1] for resp in resplist], reverse=True)

    if sensor_model and isinstance(str, sensor_model):
        # include resp only if file starts with sensor_model
        resplist = [resp for resp in resplist if re.search('^' + sensor_model + '_', resp, flags=re.IGNORECASE)]

    if component and isinstance(str, component):
        # include only if matches compoentn or NO component indicator on file
       include = (re.search('_' + component + '$', resp,  flags=re.IGNORECASE) != None)
        # or 
        #                     (re.search('_[ZEN]$', resp, flags=re.IGNORECASE) not None)
        resplist = [resp for resp in resplist if include]


    didquit, ndx = pick(resplist, 
                        title='Sensor responses for: ' + station.upper(), 
                        prompt='Select # of desired response (or "q" to quit): ',
                        allow_quit_q=True, 
                        menu_on_error=True,
                        err_message='Invalid selection. Please try again.',
                        indent_width=8)

    if didquit: 
        resp_file = None
    else:
        resp_file = os.path.join(IDA_RESPONSES_CUR_DIR, resplist[ndx])

    return resp_file


