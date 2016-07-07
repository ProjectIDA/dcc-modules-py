import datetime
from pandas.core.frame import DataFrame

import ida.db.datascope.query

def find_sensor_file(df, station, loca, comp, date_8601_str):

    if isinstance(df, DataFrame):
        dfile_list = ida.db.datascope.query.find_sensor_file(df, station, loca, comp, date_8601_str)

    return dfile_list

