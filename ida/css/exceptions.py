"""
ida.css.exceptions
~~~~~~~~~~~~~~~~~~~

This module contains the CSS exceptions.
"""

class WfdiscError(Exception):
    """Baae class for WFDISC file dan data errors"""


class WfdiscFileNotFoundException(WfdiscError):
    def __init__(self, filename, msg=None):
        if msg is None:
            msg = "File not found %s" % filename
        super(WfdiscError, self).__init__(msg)


class WfdiscSegmentRecordInvalidSize(WfdiscError):
    def __init__(self, msg):
        if msg is None:
            msg = "Invalid WFDISC record"
        super(WfdiscError, self).__init__(msg)


class WfdiscSegmentTimeFormatError(WfdiscError):
    """"wfdisc Rec TIME format conversion error"""


class WfdiscSegmentWfidFormatError(WfdiscError):
    """wfdisc Rec WFID format conversion error"""


class WfdiscSegmentChanidFormatError(WfdiscError):
    "wfdisc Rec CHANID format conversion error"


class WfdiscSegmentJdateFormatError(WfdiscError):
    "wfdisc Rec JDATE format conversion error"


class WfdiscSegmentEndtimeFormatError(WfdiscError):
    "wfdisc Rec ENDTIME format conversion error"


class WfdiscSegmentNsampFormatError(WfdiscError):
    "wfdisc Rec NSAMP format conversion error"


class WfdiscSegmentSamprateFormatError(WfdiscError):
    "wfdisc Rec SAMPRATE format conversion error"


class WfdiscSegmentCalibFormatError(WfdiscError):
    "wfdisc Rec CALIB format conversion error"


class WfdiscSegmentCalperFormatError(WfdiscError):
    "wfdisc Rec CALPER format conversion error"


class WfdiscSegmentFileNotFoundError(WfdiscError):
    "wfdisc DFILE File Not Found error"


class WfdiscSegmentDatatypeValueError(WfdiscError):
    "wfdisc DFILE File Not Found error"


class WfdiscSegmentFoffFormatError(WfdiscError):
    "wfdisc Rec FOFF format conversion error"


class WfdiscSegmentCommidFormatError(WfdiscError):
    "wfdisc Rec COMMID format conversion error"


class WfdiscSegmentLddateFormatError(WfdiscError):
    "wfdisc Rec LDDATE format conversion error"

