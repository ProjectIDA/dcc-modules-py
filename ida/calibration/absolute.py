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
#from datetime import datetime
#import os.path
from pathlib import Path
import yaml
#import collections

from fabulous.color import red, bold

class AbsOnsiteConfig(object):

    def __init__(self, fn):
        with open(fn, 'rt') as cfl:
            config_txt = cfl.read()

        try:
            self._config = yaml.load(config_txt)
            self._config_dir = Path(fn).parent
        except:
            print(red(bold('Error parsing config file: ' + fn)))
            print(config_txt)

    @property
    def ref_azi_seed_file(self):
        return  self._config['azimuth_processing']['ms_file']

    @property
    def ref_abs_seed_file(self):
        return  self._config['abssens_processing']['ms_file']

    @property
    def segment_size_samples(self):
        return self._config['segment_size_samples']

    @property   
    def minimum_segment_cnt(self):
        return self._config['minimum_segment_cnt']


