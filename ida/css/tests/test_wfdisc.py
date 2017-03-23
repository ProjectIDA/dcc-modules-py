import os.path
from pathlib import Path
import pytest
import tempfile
import ida.css.wfdisc

@pytest.fixture
def wfdisc_rec():
    return 'TKL BHZ 1488931200.01900 -1 -1 2017066 1488945599.99400 576000 40.0000000 0.063238 1.000000 - - s3 - ./ TKL_ALL.20170308.0000.mwf 3400 -1 1489009080.00000'

@pytest.fixture
def wfdisc_filename_18recs():
    return 'TKL_ALL.20170308.0000.wfdisc'
    
@pytest.fixture
def wfdisc_filepath():
    return os.path.join(os.path.dirname(__file__), 'data')
    
@pytest.fixture
def wfdisc_fullpath():
    return os.path.join(os.path.dirname(__file__), 'data', 'TKL_ALL.20170308.0000.wfdisc')
    
@pytest.fixture
def segment_first_values():
    return [5704,6054,5677,5364,6076,5314,1820,1929,1304,1978,1453,1463,3243,3724,3683,3779,3945,3511]

@pytest.fixture
def segment_last_values():
    return [6059,5671,5382,6034,5289,5097,1918,1326,1973,1464,1442,1939,3723,3664,3774,3888,3528,3745]

def test_wfdisc_segment_create(wfdisc_rec, wfdisc_fullpath):

    asegment = ida.css.wfdisc.WfdiscSegment(wfdisc_fullpath, wfdisc_rec)
    assert isinstance(asegment, ida.css.wfdisc.WfdiscSegment)


def test_wfdisc_create_bad_file():

    with pytest.raises(ida.css.wfdisc.WfdiscFileNotFoundException):
        asegment = ida.css.wfdisc.WfdiscFile('missingfile.txt')


@pytest.mark.parametrize("wfdiscrec, errcls", [
    ('TKL BHZ X488931200.01900 -1 -1 2017066 1488945599.99400 576000 40.0000000 0.063238 1.000000 - - s3 - ./ TKL_ALL.20170308.0000.mwf 3400 -1 1489009080.00000', ida.css.wfdisc.WfdiscTimeFormatError),
    ('TKL BHZ 1488931200.01900 -1.9 -1 2017066 1488945599.99400 576000 40.0000000 0.063238 1.000000 - - s3 - ./ TKL_ALL.20170308.0000.mwf 3400 -1 1489009080.00000', ida.css.wfdisc.WfdiscWfidFormatError),
    ('TKL BHZ 1488931200.01900 -1 -X 2017066 1488945599.99400 576000 40.0000000 0.063238 1.000000 - - s3 - ./ TKL_ALL.20170308.0000.mwf 3400 -1 1489009080.00000', ida.css.wfdisc.WfdiscChanidFormatError),
    ('TKL BHZ 1488931200.01900 -1 -1 20s17066 1488945599.99400 576000 40.0000000 0.063238 1.000000 - - s3 - ./ TKL_ALL.20170308.0000.mwf 3400 -1 1489009080.00000', ida.css.wfdisc.WfdiscJdateFormatError),
    ('TKL BHZ 1488931200.01900 -1 -1 2017066 X488945599.99400 576000 40.0000000 0.063238 1.000000 - - s3 - ./ TKL_ALL.20170308.0000.mwf 3400 -1 1489009080.00000', ida.css.wfdisc.WfdiscEndtimeFormatError),
    ('TKL BHZ 1488931200.01900 -1 -1 2017066 1488945599.99400 576000.5 40.0000000 0.063238 1.000000 - - s3 - ./ TKL_ALL.20170308.0000.mwf 3400 -1 1489009080.00000', ida.css.wfdisc.WfdiscNsampFormatError),
    ('TKL BHZ 1488931200.01900 -1 -1 2017066 1488945599.99400 576000 40s.0000000 0.063238 1.000000 - - s3 - ./ TKL_ALL.20170308.0000.mwf 3400 -1 1489009080.00000', ida.css.wfdisc.WfdiscSamprateFormatError),
    ('TKL BHZ 1488931200.01900 -1 -1 2017066 1488945599.99400 576000 40.0000000 0.x063238 1.000000 - - s3 - ./ TKL_ALL.20170308.0000.mwf 3400 -1 1489009080.00000', ida.css.wfdisc.WfdiscCalibFormatError),
    ('TKL BHZ 1488931200.01900 -1 -1 2017066 1488945599.99400 576000 40.0000000 0.063238 1x.000000 - - s3 - ./ TKL_ALL.20170308.0000.mwf 3400 -1 1489009080.00000', ida.css.wfdisc.WfdiscCalperFormatError),
    ('TKL BHZ 1488931200.01900 -1 -1 2017066 1488945599.99400 576000 40.0000000 0.063238 1.000000 - - s3 - ./ TKL_ALL.20170308.0000.mwf x3400 -1 1489009080.00000', ida.css.wfdisc.WfdiscFoffFormatError),
    ('TKL BHZ 1488931200.01900 -1 -1 2017066 1488945599.99400 576000 40.0000000 0.063238 1.000000 - - s3 - ./ TKL_ALL.20170308.0000.mwf 3400 -x1 1489009080.00000', ida.css.wfdisc.WfdiscCommidFormatError),
    ('TKL BHZ 1488931200.01900 -1 -1 2017066 1488945599.99400 576000 40.0000000 0.063238 1.000000 - - s3 - ./ TKL_ALL.20170308.0000.mwf 3400 -1 148900908s0', ida.css.wfdisc.WfdiscLddateFormatError),
    ('TKL BHZ 1488931200.01900 -1 -1 2017066 1488945599.99400 576000 40.0000000 0.063238 1.000000 - - X3 - ./ TKL_ALL.20170308.0000.mwf 3400 -1 1489009080', ida.css.wfdisc.WfdiscDatatypeValueError),
])
def test_wfdisc_segment_bad_input(wfdiscrec, errcls, wfdisc_fullpath):

    with pytest.raises(errcls):
        _ = ida.css.wfdisc.WfdiscSegment(wfdisc_fullpath, wfdiscrec)


def test_wfdisc_read_recs(wfdisc_fullpath):
    
    wfdics_data = ida.css.wfdisc.WfdiscFile(wfdisc_fullpath)
    assert wfdics_data.segment_cnt == 18
  

def test_wfdisc_first_last_sample_values(wfdisc_fullpath, segment_first_values, segment_last_values):
    
    wfdisc = ida.css.wfdisc.WfdiscFile(wfdisc_fullpath)
    for ndx, seg in enumerate(wfdisc.segments):
        assert seg.samples[0] == segment_first_values[ndx]
        assert seg.samples[-1] == segment_last_values[ndx]
        
        
def test_wfdisc_convert_raw_wf_s3():
    
    data = b'\x00\x17\xa6\x00\x17\xa2\x00\x17\x9a\x00\x17\x8b\xFF\x17\x82'
    newdata = ida.css.wfdisc.WfdiscSegment.convert_s3(data)
    assert list(newdata) == [6054, 6050, 6042, 6027, 16717698]

