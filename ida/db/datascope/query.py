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
import datetime

from pandas.core.frame import DataFrame

def find_sensor_file(df, station, loca, comp, date_8601_str):

    if not isinstance(df, DataFrame):
        raise ValueError('df is not Pandas Data Frame')
    if not isinstance(station, str) or len(station) < 1:
        raise ValueError('Invalid station value: {}\nShould be the alphanumeric STATION code.'.format(station))
    if not isinstance(loca, str) or len(loca) != 2:
        raise ValueError('Invalid location value: {}\nShould be the 2 character location code.'.format(loca))
    if not isinstance(comp, str) or len(comp) != 1:
        raise ValueError('Invalid component value: {}\nShould be one of {}.'.format(station, ['Z', 'N', 'E']))
    if not isinstance(date_8601_str, str) or len(date_8601_str) != 10:
        raise ValueError('Invalid date value: {}\nShould be a date in format YYYY-MM-DD.'.format(date_8601_str))

    dt = datetime.datetime.strptime(date_8601_str, '%Y-%m-%d')

    comp = comp.lower()
    comp_list = ['vh'+comp.lower()]
    if comp == 'n':
        comp_list.append('vh1')
    elif comp == 'e':
        comp_list.append('vh2')

    condition = (df.stageid == 1) & (df.sta == station.upper()) & (df.loca == loca) & \
                (df.chn.isin(comp_list)) & (df.begt <= dt) & (df.endt >= dt)
    file_df = df[condition].dfile

    return list(file_df)

def get_stages(df, station, loc, chn, date_8601_str):
    '''Query db and return 0 or more stage records for given channel & date'''


    if not isinstance(df, DataFrame):
        raise ValueError('df is not Pandas Data Frame')

    try:
        dt = datetime.datetime.strptime(date_8601_str, '%Y-%m-%dT%H:%M:%S')
    except:
        try:
            dt = datetime.datetime.strptime(date_8601_str, '%Y-%m-%dT%H:%M')
        except:
            # last chance
            dt = datetime.datetime.strptime(date_8601_str, '%Y-%m-%d')


    condition = (df.sta == station.upper()) & (df.loca == loc) & \
                (df.chn == chn.lower()) & (df.begt <= dt) & (df.endt >= dt)

    stage_list = df[condition].sort_values('stageid', ascending=True)

    return stage_list
