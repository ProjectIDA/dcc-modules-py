#######################################################################################################################
# Copyright (C) 2018  Regents of the University of California
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

import pytest
import argparse
import datetime
import sys
from utils import imageOutputFile, mseedInputFilename

# test image output file name generation

def test_image_output_filename():
    starttime = datetime.datetime.strptime("2018-05-13T00-00-00", "%Y-%m-%dT%H-%M-%S")
    endtime = datetime.datetime.strptime("2018-05-13T23-59-59", "%Y-%m-%dT%H-%M-%S")

    outputFilename = imageOutputFile("aak", "00", "bhz", starttime, endtime, '/ida/web/ppsd/page_images/2018/05')
    expectedFilename = '/ida/web/ppsd/page_images/2018/05/aak.00.bhz_20180513T000000-20180513T235959.png'

    assert outputFilename == expectedFilename


# test mseed input file name generation

def test_mseed_input_filename():
    inputFilename = mseedInputFilename('aak/2018/133/II.AAK.10.BH1.2018.133')
    expectedFilename = '/ida/archive/ms/aak/2018/133/II.AAK.10.BH1.2018.133'

    assert inputFilename == expectedFilename


