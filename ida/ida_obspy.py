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

from numpy import require
from ida.obspy.io.mseed.core import _read_mseed, _is_mseed
from ida.signals.trace import IDATrace
from ida.signals.stream import IDAStream

"""IDA Wrapper around subset of Obspy functionality"""

def read_mseed(file_name_or_object, dtype=None):
    """Read miniseed data from specified file and optionally cast to type dtype.

    :param file_name_or_object: filename or stream to read from
    :type file_name_or_object: str, file
    :param dtype: numpy datatype to cast raw sample values to
    :type dtype: np.dtype
    :return: IDAStream instance consting of IDATraces with miniseed raw data
    :rtype: IDAStream
    """

    if _is_mseed(file_name_or_object):

        tracelist = _read_mseed(file_name_or_object)

        idatracelist = []
        for trace in tracelist:
            header = trace[0]
            data = trace[1]
            if dtype:
                data = require(data, dtype)

            idatracelist.append(IDATrace(header, data=data))
    else:
        raise TypeError('File is not a valid miniseed object.')

    return IDAStream(traces=idatracelist)

