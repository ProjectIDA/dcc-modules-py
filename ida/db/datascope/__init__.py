import datetime

STAGE_COLS = [
    
    ('sta',       0,   6),
    ('chn',       7,   8),
    ('loc',      16,   2),
    ('begt',     19,  17),
    ('endt',     37,  17),
    ('stageid',  55,   8),
    ('ssident',  64,  16),
    ('gnom',     81,  11),
    ('gcalib',   93,  10),
    ('iunits',  104,  16),
    ('ounits',  121,  16),
    ('izero',   138,  8),
    ('decifac', 147,  8),
    ('srate',   156,  11),
    ('leadfac', 168,  11),
    ('dir',     180,  64),
    ('dfile',   245,  32),
    ('lddate',  278,  17),
]



CHAN_COLS = [
    # name     sta wid  
    ('sta',      0,  6),
    ('chn',      7,  8),
    ('loc',     16,  2),
    ('begt',    19, 17),
    ('endt',    37, 17),
    ('edepth',  55,  9),
    ('hang',    65,  6),
    ('vang',    72,  6),
    ('flag',    79,  2),
    ('instype', 82,  6),
    ('nomfreq', 89, 16),
]

DATE_CNVTRS = {
    'begt': parse_dt,
    'endt': parse_dt,
    'lddate': parse_dt
}

COL_DELIM = ' '

HEADER_ROW = None

def parse_dt(dt_str):
    return datetime.datetime.fromtimestamp(float(dt_str))
