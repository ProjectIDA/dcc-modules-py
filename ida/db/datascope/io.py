import os
import pandas as pd
from ida import IDA_DATASCOPEDB_DIR
from ida.db.datascope import STAGE_COLS, CHAN_COLS, DATE_CNVTRS, DB_TABLES


def read(table_name):

    if table_name not in DB_TABLES:
        raise ValueError('Invalid DB_TABLE: '+ table_name)

    if table_name.upper() == 'STAGE':
        return _read_stage()
    elif table_name.upper() == 'CHAN':
        return _read_chan


def _read_stage():

    stage_colspecs = [(stage_col[1], stage_col[1] + stage_col[2]) for stage_col in STAGE_COLS]
    stage_names  = [stage_col[0] for stage_col in STAGE_COLS]

    table_path = os.path.join(IDA_DATASCOPEDB_DIR, 'IDA.stage')

    stage_df = pd.read_fwf(table_path, 
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

