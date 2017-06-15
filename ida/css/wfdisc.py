import os.path
from inspect import currentframe
import sys
import numpy as np
from obspy import Trace, Stream, UTCDateTime
from obspy.core.trace import Stats
from .exceptions import (
    WfdiscFileNotFoundException,
    WfdiscSegmentRecordInvalidSize, WfdiscSegmentTimeFormatError, WfdiscSegmentWfidFormatError,
    WfdiscSegmentChanidFormatError, WfdiscSegmentJdateFormatError, WfdiscSegmentEndtimeFormatError,
    WfdiscSegmentNsampFormatError, WfdiscSegmentSamprateFormatError, WfdiscSegmentCalibFormatError,
    WfdiscSegmentCalperFormatError, WfdiscSegmentFileNotFoundError, WfdiscSegmentDatatypeValueError,
    WfdiscSegmentFoffFormatError, WfdiscSegmentCommidFormatError, WfdiscSegmentLddateFormatError
)

class WfdiscFile(object):

    def __init__(self, filename, *args, **kwargs):

        if not os.path.exists(filename):
            raise WfdiscFileNotFoundException(filename or '<no filename supplied>')

        self._fpath = os.path.split(filename)[0]
        self._segments = []
        # self._samples = []
        with open(filename, 'rt') as wffil:
            for segrec in wffil:
                newwfseg = WfdiscSegment(filename, segrec, *args, **kwargs)
                ok, err = newwfseg.load_wf_data()
                if ok:
                    self._segments.append(newwfseg)
                else:
                    print(err, file=sys.stderr)


    @property
    def segment_cnt(self):
        return len(self._segments)

    @property
    def segments(self):
        return self._segments

    def __len__(self):
        return len(self._segments)

    def __getitem__(self, item):
        return self._segments[item]

    def write_miniseed(self, filename: str, network: str) -> (bool, str):

        st = Stream()
        #create Trace for each segment
        for seg in self.segments:
            header = Stats()
            header['network'] = network
            header['station'] = seg.seginfo[seg.KEY_STA]
            if len(seg.seginfo[seg.KEY_CHAN]) == 5:
                header['channel'] = seg.seginfo[seg.KEY_CHAN][:3]
                header['location'] = seg.seginfo[seg.KEY_STA][3:5]
            else:
                header['channel'] = seg.seginfo[seg.KEY_CHAN]
                header['location'] = ''
            header['delta'] = round(1.0 / seg.seginfo[seg.KEY_SAMPRATE], 3)
            header['calib'] = seg.seginfo[seg.KEY_CALIB]
            header['npts'] = seg.samples.size
            header['starttime'] = UTCDateTime(seg.seginfo[seg.KEY_TIME])
            header['dataquality'] = 'R'
            header['mseed'] = {}
            if sys.byteorder == 'big':
                header['mseed']['byteorder'] = '>'
            elif sys.byteorder == 'little':
                header['mseed']['byteorder'] = '<'

            tr = Trace(data=seg.samples, header=header)
            st.append(tr)

        if not hasattr(filename, 'write'):
            fd = open(filename, 'wb')
        else:
            fd = filename

        st.write(fd, format='MSEED')

        return True, ''


class WfdiscSegment(object):
    """Object representing the data and metadata of a single record in a Wfdisc file"""

    SEGMENT_REC_FIELDCNT = 20
    KEY_STA = 'sta'
    KEY_CHAN = 'chan'
    KEY_TIME = 'time'
    KEY_WFID = 'wfid'
    KEY_CHANID = 'chanid'
    KEY_JDATE = 'jdate'
    KEY_ENDTIME = 'endtime'
    KEY_NSAMP = 'nsamp'
    KEY_SAMPRATE = 'samprate'
    KEY_CALIB = 'calib'
    KEY_CALPER = 'calper'
    KEY_INSTYPE = 'instype'
    KEY_SEGTYPE = 'segtype'
    KEY_DATATYPE = 'datatype'
    KEY_CLIP = 'clip'
    KEY_DIR = 'dir'
    KEY_DFILE = 'dfile'
    KEY_FOFF = 'foff'
    KEY_COMMID = 'commid'
    KEY_LDDATE = 'lddate'

    KEYS_DATATYPES = {
        's3': {'size': 3}
    }

    def __init__(self, wfdisc_fn: str, segrec: str, *args, **kwargs):
        """Constructor for WfdiscSegment object

        Args:
            wfdisc_fn: Source wfdisc file path. Needed to obtain absolute location of WF data file
            segrec: Record from Wfdisc fiel with which to create this segment
        """

        self._skip_bindata = kwargs.get('skip_bindata', False)

        self.seginfo = {}
        self._samples = None
        self._fpath = wfdisc_fn

        flds = segrec.split()
        if len(flds) == self.SEGMENT_REC_FIELDCNT:
            self.seginfo[self.KEY_STA] = flds[0]
            self.seginfo[self.KEY_CHAN] = flds[1]
            try:
                self.seginfo[self.KEY_TIME] = float(flds[2])
            except ValueError as e:
                raise WfdiscSegmentTimeFormatError("Error converting TIME to float: {}".format(flds[2]))
            try:
                self.seginfo[self.KEY_WFID] = int(flds[3])
            except ValueError as e:
                raise WfdiscSegmentWfidFormatError("Error converting WFID to integer: {}".format(flds[3]))
            try:
                self.seginfo[self.KEY_CHANID] = int(flds[4])
            except ValueError as e:
                raise WfdiscSegmentChanidFormatError("Error converting CHANID to integer: {}".format(flds[4]))
            try:
                self.seginfo[self.KEY_JDATE] = int(flds[5])
            except ValueError as e:
                raise WfdiscSegmentJdateFormatError("Error converting JDATE to integer: {}".format(flds[5]))
            try:
                self.seginfo[self.KEY_ENDTIME] = float(flds[6])
            except ValueError as e:
                raise WfdiscSegmentEndtimeFormatError("Error converting ENDTIME to float: {}".format(flds[6]))
            try:
                self.seginfo[self.KEY_NSAMP] = int(flds[7])
            except ValueError as e:
                raise WfdiscSegmentNsampFormatError("Error converting NSAMP to integer: {}".format(flds[7]))
            try:
                self.seginfo[self.KEY_SAMPRATE] = float(flds[8])
            except ValueError as e:
                raise WfdiscSegmentSamprateFormatError("Error converting SAMPRATE to float: {}".format(flds[8]))
            try:
                self.seginfo[self.KEY_CALIB] = float(flds[9])
            except ValueError as e:
                raise WfdiscSegmentCalibFormatError("Error converting CALIB to float: {}".format(flds[9]))
            try:
                self.seginfo[self.KEY_CALPER] = float(flds[10])
            except ValueError as e:
                raise WfdiscSegmentCalperFormatError("Error converting CALPER to float: {}".format(flds[10]))

            self.seginfo[self.KEY_INSTYPE] = flds[11]
            self.seginfo[self.KEY_SEGTYPE] = flds[12]

            if self._skip_bindata or (flds[13] in self.KEYS_DATATYPES):
                self.seginfo[self.KEY_DATATYPE] = flds[13]
            else:
                raise WfdiscSegmentDatatypeValueError('Invalid data type: {}'.format(flds[13]))

            self.seginfo[self.KEY_CLIP] = flds[14]
            self.seginfo[self.KEY_DIR] = flds[15]
            self.seginfo[self.KEY_DFILE] = flds[16]

            try:
                self.seginfo[self.KEY_FOFF] = int(flds[17])
            except ValueError as e:
                raise WfdiscSegmentFoffFormatError("Error converting FOFF to integer: {}".format(flds[17]))
            try:
                self.seginfo[self.KEY_COMMID] = int(flds[18])
            except ValueError as e:
                raise WfdiscSegmentCommidFormatError("Error converting COMMID to integer: {}".format(flds[18]))
            try:
                self.seginfo[self.KEY_LDDATE] = float(flds[19])
            except ValueError as e:
                raise WfdiscSegmentLddateFormatError("Error converting LDDATE to float: {}".format(flds[19]))

        else:
            raise WfdiscSegmentRecordInvalidSize(len(flds))


    @property
    def samples(self):
        return self._samples

    @samples.setter
    def samples(self, value):
        self._samples = value

    def load_wf_data(self) -> (bool, str):
        """Read WF datam from WF binary file into self.samples

        Returns:
            Success: flag to indicate load success (True) or failure (False)
            Errmsg: to indicate nature of error. '' if successful.
        """

        if not self._skip_bindata:
            start_byte = self.seginfo[self.KEY_FOFF]
            sample_size = self.KEYS_DATATYPES[self.seginfo[self.KEY_DATATYPE]]['size']
            bytes_to_read = sample_size * self.seginfo[self.KEY_NSAMP]

            wffn = os.path.join(os.path.split(self._fpath)[0], self.seginfo[self.KEY_DIR], self.seginfo[self.KEY_DFILE])
            if os.path.exists(wffn) and os.path.isfile(wffn):
                with open(wffn, 'rb') as wffl:
                    wffl.seek(start_byte)
                    rawdata = wffl.read(bytes_to_read)
            else:
                return False, '{}: File not found: {}'.format(currentframe().f_code.co_name, wffn)

            # convert data based on datatype
            if self.seginfo[self.KEY_DATATYPE] == 's3':
                wfdata = WfdiscSegment.convert_s3(rawdata)
            else:
                return False, '{}: Unsupported wfdisc datatype: {}'.format(currentframe().f_code.co_name,
                                                                           self.seginfo[self.KEY_DATATYPE])

            self.samples = wfdata

        return True, ''

    @classmethod
    def convert_s3(self, rawdata: bytes) -> np.ndarray:
        """Function to convert Wfdisc WF data in 's3' format.

        Args:
            rawdata: Sequence of raw bytes to convert.

        Returns:
            Converted bytes in numpy array of int32

        """
        wfdata = []

        for ndx in range(0, len(rawdata), 3):
            val = (rawdata[ndx] << 16) + (rawdata[ndx + 1] << 8) + (rawdata[ndx + 2])
            if val > pow(2, 23):
                val -= pow(2, 24)
            wfdata.append(val)

        return np.array(wfdata, np.int32)

