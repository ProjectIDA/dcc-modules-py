import os.path
import sys
import numpy as np
from obspy import Trace, Stream, UTCDateTime
from obspy.core.trace import Stats

class WfdiscError(Exception):
    """Baae class for WFDISC file dan data errors"""

class WfdiscFileNotFoundException(WfdiscError):
    pass

class WfdiscSegmentRecordInvalidSize(WfdiscError):
    pass

class WfdiscTimeFormatError(WfdiscError):
    "wfdisc Rec TIME format conversion error"

class WfdiscWfidFormatError(WfdiscError):
    "wfdisc Rec WFIDformat conversion error"

class WfdiscChanidFormatError(WfdiscError):
    "wfdisc Rec CHANID format conversion error"

class WfdiscJdateFormatError(WfdiscError):
    "wfdisc Rec JDATE format conversion error"

class WfdiscEndtimeFormatError(WfdiscError):
    "wfdisc Rec ENDTIME format conversion error"

class WfdiscNsampFormatError(WfdiscError):
    "wfdisc Rec NSAMP format conversion error"

class WfdiscSamprateFormatError(WfdiscError):
    "wfdisc Rec SAMPRATE format conversion error"

class WfdiscCalibFormatError(WfdiscError):
    "wfdisc Rec CALIB format conversion error"

class WfdiscCalperFormatError(WfdiscError):
    "wfdisc Rec CALPER format conversion error"

class WfdiscDfileNotFoundError(WfdiscError):
    "wfdisc DFILE File Not Found error"
    
class WfdiscDatatypeValueError(WfdiscError):
    "wfdisc DFILE File Not Found error"
    
class WfdiscFoffFormatError(WfdiscError):
    "wfdisc Rec FOFF format conversion error"

class WfdiscCommidFormatError(WfdiscError):
    "wfdisc Rec COMMID format conversion error"

class WfdiscLddateFormatError(WfdiscError):
    "wfdisc Rec LDDATE format conversion error"

class WfdiscFile(object):
    
    def __init__(self, filename=None):
        
        if filename and not os.path.exists(filename):
            raise WfdiscFileNotFoundException('WFDISC file not fond: {}'. format(filename))
        
        self._fpath = os.path.split(filename)[0] 
        self._segments = []
        # self._samples = []
        with open(filename, 'rt') as wffil:
            for segrec in wffil:
                newwfseg = WfdiscSegment(filename, segrec)
                # newwfseg.samples = self.load_wf_data(newwfseg)
                self._segments.append(newwfseg)
                # self._samples.append(self.load_wf_data(newwfseg))


    @property
    def segment_cnt(self):
        return len(self._segments)
    
    @property
    def segments(self):
        return self._segments
    
    # @property
    # def samples(self):
    #     return self._samples
    
   
    def write_miniseed(self, filename, network):
        
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
            # header['endtime'] = UTCDateTime(seg.seginfo[seg.KEY_ENDTIME])

            header['mseed'] = {}            
            if sys.byteorder == 'big':
                header['mseed']['byteorder'] = '>'
            elif sys.byteorder == 'little':
                header['mseed']['byteorder'] = '<'

            header['dataquality'] = 'R'
            # header['encoding'] = 'STEIM2'
            # header['filesize'] = 8192
            # header['number_of_records'] = 2
            # header['record_length'] = 4096
            
            tr = Trace(data=seg.samples, header=header)
            st.append(tr)

        if not hasattr(filename, 'write'):
            fd = open(filename, 'wb')
        else:
            fd = filename
        
        st.write(fd, format='MSEED')

    
    
class WfdiscSegment(object):
    
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
        's3': {
            'endian': 'big',
            'size': 3
        }
    }
        
    def __init__(self, wfdisc_fn, segrec):
        
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
                raise WfdiscTimeFormatError("Error converting TIME to float: {}".format(flds[2]))
            try:
                self.seginfo[self.KEY_WFID] = int(flds[3])
            except ValueError as e:
                raise WfdiscWfidFormatError("Error converting WFID to integer: {}".format(flds[3]))
            try:
                self.seginfo[self.KEY_CHANID] = int(flds[4])
            except ValueError as e:
                raise WfdiscChanidFormatError("Error converting CHANID to integer: {}".format(flds[4]))
            try:
                self.seginfo[self.KEY_JDATE] = int(flds[5])
            except ValueError as e:
                raise WfdiscJdateFormatError("Error converting JDATE to integer: {}".format(flds[5]))
            try:
                self.seginfo[self.KEY_ENDTIME] = float(flds[6])
            except ValueError as e:
                raise WfdiscEndtimeFormatError("Error converting ENDTIME to float: {}".format(flds[6]))
            try:
                self.seginfo[self.KEY_NSAMP] = int(flds[7])
            except ValueError as e:
                raise WfdiscNsampFormatError("Error converting NSAMP to integer: {}".format(flds[7]))
            try:
                self.seginfo[self.KEY_SAMPRATE] = float(flds[8])
            except ValueError as e:
                raise WfdiscSamprateFormatError("Error converting SAMPRATE to float: {}".format(flds[8]))
            try:
                self.seginfo[self.KEY_CALIB] = float(flds[9])
            except ValueError as e:
                raise WfdiscCalibFormatError("Error converting CALIB to float: {}".format(flds[9]))
            try:
                self.seginfo[self.KEY_CALPER] = float(flds[10])
            except ValueError as e:
                raise WfdiscCalperFormatError("Error converting CALPER to float: {}".format(flds[10]))
                
            self.seginfo[self.KEY_INSTYPE] = flds[11]
            self.seginfo[self.KEY_SEGTYPE] = flds[12]
            
            if flds[13] in self.KEYS_DATATYPES:
                self.seginfo[self.KEY_DATATYPE] = flds[13]
            else:
                raise WfdiscDatatypeValueError('Invalid data type: {}'.format(flds[13]))
            
            self.seginfo[self.KEY_CLIP] = flds[14]
            self.seginfo[self.KEY_DIR] = flds[15]
            self.seginfo[self.KEY_DFILE] = flds[16]

            try:
                self.seginfo[self.KEY_FOFF] = int(flds[17])
            except ValueError as e:
                raise WfdiscFoffFormatError("Error converting FOFF to integer: {}".format(flds[17]))
            try:
                self.seginfo[self.KEY_COMMID] = int(flds[18])
            except ValueError as e:
                raise WfdiscCommidFormatError("Error converting COMMID to integer: {}".format(flds[18]))
            try:
                self.seginfo[self.KEY_LDDATE] = float(flds[19])
            except ValueError as e:
                raise WfdiscLddateFormatError("Error converting LDDATE to float: {}".format(flds[19]))
            
            self.load_wf_data()

        else:
            raise WfdiscSegmentRecordInvalidSize(len(flds))
               
        
    @property
    def samples(self):
        return self._samples
    
    @samples.setter
    def samples(self, value):
        self._samples = value

    def load_wf_data(self):

        start_byte = self.seginfo[self.KEY_FOFF]
        sample_size = self.KEYS_DATATYPES[self.seginfo[self.KEY_DATATYPE]]['size']
        bytes_to_read = sample_size * self.seginfo[self.KEY_NSAMP]

        wffn = os.path.join(os.path.split(self._fpath)[0], self.seginfo[self.KEY_DIR], self.seginfo[self.KEY_DFILE])
        with open(wffn, 'rb') as wffl:
            wffl.seek(start_byte)
            rawdata = wffl.read(bytes_to_read)

        # convert data based on datatype
        wfdata = np.array([])
        if self.seginfo[self.KEY_DATATYPE] == 's3':
            wfdata = WfdiscSegment.convert_s3(rawdata)

        self.samples = wfdata

    @classmethod
    def convert_s3(self, rawdata):

        wfdata = []

        for ndx in range(0, len(rawdata), 3):
            val = (rawdata[ndx] << 16) + (rawdata[ndx + 1] << 8) + (rawdata[ndx + 2])
            if val > pow(2, 23):
                val -= pow(2, 24) 
            wfdata.append(val)

        return np.array(wfdata, dtype=np.int32)

