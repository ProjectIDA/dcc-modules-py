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
import os
import pandas as pd
from ida.db.datascope import STAGE_COLS, CHAN_COLS, DATE_CNVTRS, DB_TABLES

def read(db_dir, table_name):

    if table_name not in DB_TABLES:
        raise ValueError('Invalid DB_TABLE: '+ table_name)

    if table_name.upper() == 'STAGE':
        return _read_stage(db_dir)
    elif table_name.upper() == 'CHAN':
        return _read_chan(db_dir)


def _read_stage(db_dir):

    stage_colspecs = [(stage_col[1], stage_col[1] + stage_col[2]) for stage_col in STAGE_COLS]
    stage_names  = [stage_col[0] for stage_col in STAGE_COLS]
    table_path = os.path.join(db_dir, 'IDA.stage')

    if os.path.exists(table_path):

        stage_df = pd.read_fwf(table_path,
            names=stage_names,
            colspecs=stage_colspecs,
            header=None,
            converters=DATE_CNVTRS)
        results = True

    else:
        results = False
        stage_df = None

    return results, stage_df


def _read_chan(db_dir):

    chan_colspecs = [(chan_col[1], chan_col[1] + chan_col[2]) for chan_col in CHAN_COLS]
    chan_names  = [chan_col[0] for chan_col in CHAN_COLS]

    table_path = os.path.join(db_dir, 'IDA.chan')

    if os.path.exists(table_path):

        chan_df = pd.read_fwf(table_path,
            names=chan_names,
            colspecs=chan_colspecs,
            header=None,
            converters=DATE_CNVTRS)

        results = True

    else:
        results = False
        chan_df = None

    return results, chan_df


