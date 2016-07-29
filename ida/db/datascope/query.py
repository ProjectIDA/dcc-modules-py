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
