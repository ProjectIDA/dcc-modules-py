"""
IDA CSS Library
~~~~~~~~~~~~~~~~~~~~~

The IDA CSS library is a basic pibrary for reading CSS 3.0 Wfdisc and waveform data.

Basic usage for converting wfdisc to miniseed:
   >>> import wfdisc
   >>> wfd = wfdisc.WfdiscFile('stationdata.wfdisc')
   >>> wfd.save_miniseed(filename='stationdata.ms', network='II')

"""
from .wfdisc import WfdiscFile, WfdiscSegment
from .exceptions import (
    WfdiscError, WfdiscFileNotFoundException, 
    WfdiscSegmentRecordInvalidSize, WfdiscSegmentTimeFormatError, WfdiscSegmentWfidFormatError, 
    WfdiscSegmentChanidFormatError, WfdiscSegmentJdateFormatError, WfdiscSegmentEndtimeFormatError, 
    WfdiscSegmentNsampFormatError, WfdiscSegmentSamprateFormatError, WfdiscSegmentCalibFormatError, 
    WfdiscSegmentCalperFormatError, WfdiscSegmentFileNotFoundError, WfdiscSegmentDatatypeValueError, 
    WfdiscSegmentFoffFormatError, WfdiscSegmentCommidFormatError, WfdiscSegmentLddateFormatError
)