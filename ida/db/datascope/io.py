import os
import pandas as pd
from ida.db.datascope import STAGE_COLS, CHAN_COLS, DATE_CNVTRS


def read(table_name):

    db_dir = os.environ.get('')

    if table_name.upper() == 'STAGE':
        return _read_stage()
    elif table_name.upper() == 'CHAN':
        return _read_chan
    else:
        raise ValueError('Unsupported datascope table: ' + table_name)


def _read_stage():

    stage_colspecs = [(stage_col[1], stage_col[1] + stage_col[2]) for stage_col in STAGE_COLS]
    stage_names  = [stage_col[0] for stage_col in STAGE_COLS]

    stage_df = pd.read_fwf('IDA.stage', 
        names=stage_names, 
        colspecs=stage_colspecs, 
        header=None, 
        converters=DATE_CNVTRS)

    return stage_df


def _read_chan():

    chan_colspecs = [(chan_col[1], chan_col[1] + chan_col[2]) for chan_col in CHAN_COLS]
    chan_names  = [chan_col[0] for chan_col in CHAN_COLS]

    chan_df = pd.read_fwf('IDA.chan', 
        names=chan_names, 
        colspecs=chan_colspecs, 
        header=None, 
        converters=DATE_CNVTRS)

    return chan_df

