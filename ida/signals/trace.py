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

import numpy as np
import copy

class IDATrace(object):

    def __init__(self, header, data=None):
        self._header = copy.deepcopy(header)
        if isinstance(data, np.ndarray):
            self._data = data
        elif isinstance(data, list):
            self._data = np.array(data)
        elif data:
            raise TypeError('IDATrace data must be None or of type list or numpy.ndarray')


    @property
    def data(self):
        return self._data


    @data.setter
    def data(self, data):
        self._data = data
        self._header['npts'] = len(data)


    @property
    def header(self):
            return self._header


    @property
    def sampling_rate(self):
        return self.header['sampling_rate']

    @property
    def station(self):
        return self._header['station']


    @property
    def starttime(self):
        return self._header['starttime']


    @property
    def network(self):
        return self._header['network']


    @property
    def location(self):
        return self._header['location']



    @property
    def npts(self):
        return self._header['npts']


    @property
    def channel(self):
        return self._header['channel']


    @channel.setter
    def channel(self, chan):
        self.header['channel'] = chan

    def mseed(self):
        return self._header['mseed']


    @property
    def byteorder(self):
        return self._header['mseed']['byteorder']


    @property
    def record_length(self):
        return self._header['mseed']['record_length']


    @property
    def encoding(self):
        return self._header['mseed']['encoding']


    @property
    def dataquality(self):
        return self._header['mseed']['dataquality']


    @property
    def filesize(self):
        return self._header['mseed']['filesize']


    @property
    def number_of_records(self):
        return self._header['mseed']['number_of_records']

    @property
    def endtime(self):
        return self.starttime + self.npts / self.sampling_rate


    def trim(self, starttime, endtime):

        if starttime < self.starttime:
            raise ValueError('Trimmed starttime must be >= current starttime.')

        if endtime > self.endtime:
            raise ValueError('Trimmed endtime must be <= current endtime.')

        if endtime < starttime:
            raise ValueError('Endtime can not be before starttime.')

        startdelta = starttime - self.starttime
        enddelta = endtime - self.endtime  # will be 0 or negative

        self._header['starttime'] = starttime
        sr = self.sampling_rate
        self.data = self.data[int(sr * startdelta):int(sr * enddelta)]


    def get_id(self):
        """
        Return a SEED compatible identifier of the trace.

        :rtype: str
        :return: SEED identifier

        The SEED identifier contains the network, station, location and channel
        code for the current Trace object.

        .. rubric:: Example

        >>> meta = {'station': 'MANZ', 'network': 'BW', 'channel': 'EHZ'}
        >>> tr = IDATrace(header=meta)
        >>> print(tr.get_id())
        BW.MANZ..EHZ
        >>> print(tr.id)
        BW.MANZ..EHZ
        """
        out = "{network}.{station}.{location}.{channel}"
        return out.format(**self.header)

    id = property(get_id)


    def __len__(self):
        """
        Return number of data samples of the current trace.

        :rtype: int
        :return: Number of data samples.

        .. rubric:: Example

        >>> trace = IDATrace(data=np.array([1, 2, 3, 4]))
        >>> trace.count()
        4
        >>> len(trace)
        4
        """
        return len(self.data)

    count = __len__


    def __str__(self, id_length=None):
        """
        Return short summary string of the current trace.

        :rtype: str
        :return: Short summary string of the current trace containing the SEED
            identifier, start time, end time, sampling rate and number of
            points of the current trace.

        .. rubric:: Example

        >>> tr = IDATrace(header={'station':'FUR', 'network':'GR'})
        >>> str(tr)  # doctest: +ELLIPSIS
        'GR.FUR.. | 1970-01-01T00:00:00.000000Z - ... | 1.0 Hz, 0 samples'
        """
        # set fixed id width
        if id_length:
            id_fmt = '{:' + str(id_length) + '}'
        else:
            id_fmt = '{:<}'

        trace_id = id_fmt.format(self.get_id())

        out = " | {} - {} | {:.1f} Hz, {:<} samples".format(str(self.starttime),
                                                            str(self.endtime),
                                                            self.sampling_rate,
                                                            self.npts
                                                            )

        # check for masked array
        if np.ma.count_masked(self.data):
            out += ' (masked)'
        return trace_id + out

